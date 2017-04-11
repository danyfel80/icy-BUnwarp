/*
 * Copyright 2010-2016 Institut Pasteur.
 * 
 * This file is part of Icy.
 * 
 * Icy is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Icy is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Icy. If not, see <http://www.gnu.org/licenses/>.
 */
package algorithms.danyfel80.bigimage;

import java.awt.Dimension;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;

import algorithms.danyfel80.icyBufferedImage.util.IcyBufferedImageCursor;
import icy.common.exception.UnsupportedFormatException;
import icy.common.listener.RichProgressListener;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import loci.common.services.ServiceException;
import loci.formats.FormatException;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageToTiffConverter implements Callable<Void> {

	List<RichProgressListener>	listeners;
	boolean											isInterrupted;

	File			inputFile;
	File			outputFile;
	boolean[]	channels;

	/**
	 * @param file
	 * @param convertChannels
	 */
	public BigImageToTiffConverter(File inputFile, File outputFile, boolean[] convertChannels) {
		this.listeners = new LinkedList<>();

		this.inputFile = inputFile;
		this.outputFile = outputFile;
		this.channels = convertChannels;
	}

	/**
	 * Adds a listener to the progress listeners subscribed to this converter.
	 * 
	 * @param listener
	 */
	public void addProgressListener(RichProgressListener listener) {
		this.listeners.add(listener);
	}

	/**
	 * Interrupts the current conversion.
	 */
	public void interrupt() {
		this.isInterrupted = true;
	}

	/*
	 * (non-Javadoc)
	 * @see java.util.concurrent.Callable#call()
	 */
	@Override
	public Void call() throws Exception {
		this.isInterrupted = false;
		String convertingMessage = "Converting Tiles";

		// Recover input image information
		int imgSizeC = 0;
		try {
			imgSizeC = BigImageUtil.getSequenceChannelCount(inputFile.getAbsolutePath());
		} catch (IOException | UnsupportedFormatException e) {
			throw e;
		}
		Dimension imgSize = BigImageUtil.getSequenceSize(inputFile.getAbsolutePath());
		DataType imgType = BigImageUtil.getSequenceDataType(inputFile.getAbsolutePath());

		// Compute reading tile size
		int numProc = Runtime.getRuntime().availableProcessors();
		long ram = Runtime.getRuntime().freeMemory();

		long maxTileSideSize = (ram / imgSizeC);
		maxTileSideSize /= numProc;
		maxTileSideSize /= imgType.getSize();
		maxTileSideSize = (long) Math.ceil(Math.sqrt(maxTileSideSize));

		int tileSideSize = 16;
		while (maxTileSideSize > tileSideSize * 16) {
			tileSideSize *= 16;
		}
		Dimension tileDimension = new Dimension(tileSideSize, tileSideSize);
		System.out.println(
				String.format("Reading image with tiles of %d by %d pixels.", tileDimension.width, tileDimension.height));

		Dimension imageTileGridDimension = new Dimension(
				(int) Math.ceil((double) imgSize.width / (double) tileDimension.width),
				(int) Math.ceil((double) imgSize.height / (double) tileDimension.height));

		// Open saver
		BigImageSaver saver = null;
		int numComponents = 0;
		for (boolean b: channels) {
			numComponents += b ? 1 : 0;
		}
		try {
			saver = new BigImageSaver(outputFile, imgSize, numComponents, DataType.UBYTE, tileDimension);
		} catch (ServiceException | FormatException | IOException e) {
			throw e;
		}

		ThreadPoolExecutor conversionTP = (ThreadPoolExecutor) Executors.newFixedThreadPool(numProc);
		ExecutorCompletionService<Rectangle> conversionCS = new ExecutorCompletionService<>(conversionTP);

		for (int y = 0; y < imageTileGridDimension.height; y++) {
			int tileHeight = tileDimension.height;
			tileHeight = (y * tileHeight + tileHeight < imgSize.height) ? tileHeight : imgSize.height - y * tileHeight;

			for (int x = 0; x < imageTileGridDimension.width; x++) {
				int tileWidth = tileDimension.width;
				tileWidth = (x * tileWidth + tileWidth < imgSize.width) ? tileWidth : imgSize.width - x * tileWidth;

				Rectangle tileRectangle = new Rectangle(x * tileDimension.width, y * tileDimension.height, tileWidth,
						tileHeight);

				conversionCS.submit(new TileConversionTask(saver, tileRectangle));
			}
		}
		conversionTP.shutdown();

		int totalProgress = imageTileGridDimension.width * imageTileGridDimension.height;
		int currentProgress = 0;

		try {
			while (currentProgress < totalProgress) {
				conversionCS.take().get();
				currentProgress++;
				for (RichProgressListener listener: listeners) {
					listener.notifyProgress(currentProgress, totalProgress,
							convertingMessage + String.format("(%d/%d)", currentProgress, totalProgress), null);
				}
			}
		} catch (ExecutionException e) {
			conversionTP.shutdownNow();
			throw e;
		} finally {
			saver.closeWriter();
		}

		return null;
	}

	/**
	 * @author Daniel Felipe Gonzalez Obando
	 */
	public class TileConversionTask implements Callable<Rectangle> {

		BigImageSaver	saver;
		Rectangle			tileRectangle;

		/**
		 * @param saver
		 * @param tileRectangle
		 */
		public TileConversionTask(BigImageSaver saver, Rectangle tileRectangle) {
			this.saver = saver;
			this.tileRectangle = tileRectangle;
		}

		/*
		 * (non-Javadoc)
		 * @see java.util.concurrent.Callable#call()
		 */
		@Override
		public Rectangle call() throws Exception {
			// Check interruption
			if (isInterrupted) return null;

			// load tile
			Sequence tile = null;
			BigImageLoader loader = new BigImageLoader();
			tile = loader.loadDownsampledImage(inputFile.getAbsolutePath(), tileRectangle, tileRectangle.width,
					tileRectangle.height, false);
			IcyBufferedImage tileImage = tile.getFirstImage();
			// Check interruption
			if (isInterrupted) return null;

			// Convert tile
			IcyBufferedImage resultTileImage;
			int numComponents = 0;
			for (boolean b: channels) {
				numComponents += b ? 1 : 0;
			}
			resultTileImage = new IcyBufferedImage(tileImage.getSizeX(), tileImage.getSizeY(), numComponents,
					tileImage.getDataType_());
			IcyBufferedImageCursor tileCursor = new IcyBufferedImageCursor(tile, 0, 0);
			IcyBufferedImageCursor resultTileCursor = new IcyBufferedImageCursor(resultTileImage);

			for (int c = 0, rc = 0; c < channels.length; c++) {
				if (channels[c]) {
					for (int x = 0; x < tileRectangle.width; x++) {
						for (int y = 0; y < tileRectangle.height; y++) {
							resultTileCursor.set(x, y, rc, tileCursor.get(x, y, c));
						}
					}
					rc++;
				}
			}
			resultTileCursor.commitChanges();
			// Check interruption
			if (isInterrupted) return null;

			// Save result
			saver.saveTile(resultTileImage, tileRectangle.getLocation());

			return tileRectangle;
		}

	}
}
