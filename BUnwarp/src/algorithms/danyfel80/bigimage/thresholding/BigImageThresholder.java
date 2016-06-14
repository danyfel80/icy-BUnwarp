package algorithms.danyfel80.bigimage.thresholding;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import algorithms.danyfel80.bigimage.BigImageLoader;
import algorithms.danyfel80.bigimage.BigImageSaver;
import algorithms.danyfel80.bigimage.BigImageUtil;
import icy.common.exception.UnsupportedFormatException;
import icy.image.lut.LUT;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;
import loci.common.services.ServiceException;
import loci.formats.FormatException;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageThresholder {

	boolean isInterrupted;

	final File inputFile;
	final File outputFile;
	final double minValue;
	final double maxValue;

	/**
	 * Constructs the thresholder. to execute the procedure use the method
	 * execute().
	 * 
	 * @param inputFile
	 *          Input image file.
	 * @param outputFile
	 *          Output image file.
	 * @param minValue
	 *          Threshold minimum value. Any value under minValue is set to 0.
	 * @param maxValue
	 *          Threshold maximum value. Any value over maxValue is set to the
	 *          maximum data type value.
	 */
	public BigImageThresholder(File inputFile, File outputFile, double minValue, double maxValue) {
		this.inputFile = inputFile;
		this.outputFile = outputFile;
		this.minValue = Math.min(minValue, maxValue);
		this.maxValue = Math.max(minValue, maxValue);
	}

	/**
	 * Executes the multithreaded thresholding
	 */
	public void execute() {
		this.isInterrupted = false;

		// Recover input image information
		Dimension imgSize = BigImageUtil.getSequenceSize(inputFile.getAbsolutePath());
		int imgSizeC = BigImageUtil.getSequenceChannelCount(inputFile.getAbsolutePath());
		//DataType imgDataType = BigImageUtil.getSequenceDataType(inputFile.getAbsolutePath());

		// Compute tile information
		int numProc = Runtime.getRuntime().availableProcessors();
		long ram = Runtime.getRuntime().freeMemory();

		long tileSideSize = ram / (imgSizeC + 1);
		tileSideSize /= numProc;
		tileSideSize = (long) Math.ceil(Math.sqrt(tileSideSize));

		Dimension tileCount = new Dimension(imgSize.width / (int) tileSideSize, imgSize.height / (int) tileSideSize);
		if (tileCount.width * tileSideSize < imgSize.width) {
			tileCount.width++;
		}
		if (tileCount.height * tileSideSize < imgSize.height) {
			tileCount.height++;
		}
		Dimension tileSize = new Dimension((int) tileSideSize, (int) tileSideSize);
		Dimension residualTileSize = new Dimension(imgSize.width % (int) tileSideSize, imgSize.height % (int) tileSideSize);
		
		// Open saver
		BigImageSaver saver = null;
		try {
			saver = new BigImageSaver(outputFile, imgSize, 1, DataType.UBYTE, tileSize);
		} catch (ServiceException | FormatException | IOException e) {
			e.printStackTrace();
			return;
		}

		ExecutorService threadPool = Executors.newFixedThreadPool(numProc);

		// Compute each tile
		for (int i = 0; i < tileCount.width; i++) {
			int tileSizeX = ((i + 1) * tileSize.width > imgSize.width) ? residualTileSize.width : tileSize.width;

			for (int j = 0; j < tileCount.height; j++) {
				int tileSizeY = ((j + 1) * tileSize.height > imgSize.height) ? residualTileSize.height : tileSize.height;

				Point currTilePos = new Point(i * tileSize.width, j * tileSize.height);
				Dimension currTileSize = new Dimension(tileSizeX, tileSizeY);
				Rectangle currTileRect = new Rectangle(currTilePos, currTileSize);

				threadPool.submit(new TileThresholder(inputFile.getAbsolutePath(), currTileRect, saver, minValue, maxValue));
			}
		}
		threadPool.shutdown();
		try {
			threadPool.awaitTermination(1, TimeUnit.HOURS);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		try {
			saver.closeWriter();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * This method sends a stop signal to the execution state and eventually stops
	 * the execution.
	 */
	public void stopExecution() {
		isInterrupted = true;
	}

	/**
	 * @author Daniel Felipe Gonzalez Obando
	 *
	 */
	private class TileThresholder implements Runnable {

		private final Rectangle currTileRect;
		private final BigImageSaver saver;
		private final double minValue;
		private final double maxValue;

		/**
		 * @param absolutePath
		 * @param currTileRect
		 * @param saver
		 * @param minValue
		 * @param maxValue
		 */
		public TileThresholder(String absolutePath, Rectangle currTileRect, BigImageSaver saver, double minValue,
		    double maxValue) {
			this.currTileRect = currTileRect;
			this.saver = saver;
			this.minValue = minValue;
			this.maxValue = maxValue;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run() {
			if (isInterrupted) {
				return;
			}
			Sequence tile = null;
			try {
				tile = BigImageLoader.loadDownsampledImage(inputFile.getAbsolutePath(), currTileRect, currTileRect.width,
				    currTileRect.height, false);
			} catch (UnsupportedFormatException | IOException e) {
				System.err.println("Error loading tile: " + currTileRect);
				e.printStackTrace();
				return;
			}
			
			if (isInterrupted) {
				return;
			}
			LUT lut = tile.createCompatibleLUT();
			for (int i = 0; i < lut.getNumChannel(); i++) {
				lut.getLutChannel(i).setMin(minValue);
				lut.getLutChannel(i).setMax(maxValue);
				lut.getLutChannel(i).setMinBound(minValue);
				lut.getLutChannel(i).setMaxBound(maxValue);
			}
			Sequence bi = SequenceUtil.convertColor(tile, BufferedImage.TYPE_BYTE_GRAY, lut);
			try {
				saver.saveTile(bi.getFirstImage(), currTileRect.getLocation());
			} catch (ServiceException | IOException | FormatException e) {
				e.printStackTrace();
				return;
			}

		}

	}
}
