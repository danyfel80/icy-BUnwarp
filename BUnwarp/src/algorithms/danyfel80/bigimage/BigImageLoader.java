package algorithms.danyfel80.bigimage;

import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Dimension;
import java.io.IOException;

import org.apache.commons.io.FilenameUtils;

import icy.common.exception.UnsupportedFormatException;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.sequence.Sequence;
import icy.type.DataType;
import loci.formats.ome.OMEXMLMetadataImpl;
import plugins.adufour.ezplug.EzGUI;
import plugins.adufour.ezplug.EzPlug;
import plugins.kernel.importer.LociImporterPlugin;

/**
 * This class allows easy big image loading through downsampling and tile
 * loading.
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageLoader {

	private static EzGUI pluginGUI;
	
	public static void setPluginGUI(EzGUI pluginGUI) {
		BigImageLoader.pluginGUI = pluginGUI;
	}
	
	private static void setProgress(double progress) {
		if (pluginGUI != null) {
			pluginGUI.setProgressBarValue(progress);
		}
	}
	
	private static void setStatusMessage(String message) {
		if (pluginGUI != null) {
			pluginGUI.setProgressBarMessage(message);
		}
	}
	
	/**
	 * Loads an image from a file and downsamples it according to the desired
	 * resulting size. This method is performed using tiles to manage memory
	 * usage. If resultWidth and resultHeight are 0 then the image is loaded in
	 * full size.
	 * 
	 * @param path
	 *          Image file path
	 * @param resultWidth
	 *          max width of the resulting image.
	 * @param resultHeight
	 *          max height of the resulting image.
	 * @return downsampled image.
	 * @throws IOException
	 * @throws UnsupportedFormatException
	 */
	public static Sequence loadDownsampledImage(String path, Rectangle tile, int resultMaxWidth, int resultMaxHeight)
	    throws UnsupportedFormatException, IOException {
		
		double progress = 0;
		setProgress(progress);
		setStatusMessage(String.format("Loading image: %d", (int)(progress*100)));
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadataImpl imgProps = importer.getMetaData();
			String imgName = FilenameUtils.getBaseName(path);
			int imgSizeX = imgProps.getPixelsSizeX(0).getValue();
			int imgSizeY = imgProps.getPixelsSizeY(0).getValue();
			int imgSizeC = imgProps.getPixelsSizeC(0).getValue();
			DataType imgDataType = DataType.getDataTypeFromPixelType(imgProps.getPixelsType(0));
			if (tile == null) {
				tile = new Rectangle(0, 0, imgSizeX, imgSizeY);
			} else {
				if (tile.x > imgSizeX || tile.y > imgSizeY)
					return null;
				if (tile.x + tile.width > imgSizeX)
					tile.width = imgSizeX - tile.x;
				if (tile.y + tile.height > imgSizeY)
					tile.height = imgSizeY - tile.y;
			}
			System.out.println("tile to extract: " + tile);

			if (resultMaxWidth == 0) {
				resultMaxWidth = tile.width;
			}
			if (resultMaxHeight == 0) {
				resultMaxHeight = tile.height;
			}

			Runtime.getRuntime().gc();
			long ram = Runtime.getRuntime().freeMemory();
			int nProc = Runtime.getRuntime().availableProcessors();
			System.out.println("Available memory: " + ram + " bytes, Available processors: " + nProc);
			ram /= imgSizeC;
			double szMax = Math.sqrt(ram);

			int tileSize = (int) Math.ceil(szMax / nProc);
			
			while (tile.width / resultMaxWidth > 2 * tileSize)
				tileSize *= 2;
			while (tile.height / resultMaxHeight > 2 * tileSize)
				tileSize *= 2;
			System.out.println("Image size: " + imgSizeX + "px*" + imgSizeY + "px. Tile size: " + tileSize);

			double imgScaledSizeX = imgSizeX;
			double imgScaledSizeY = imgSizeY;
			double tmpSizeX = tile.getWidth();
			double tmpSizeY = tile.getHeight();
			double tileTmpSizeX = tileSize;
			double tileTmpSizeY = tileSize;
			int resolution = 1;
			double scaleFactor = 1d;
			while (tmpSizeX > resultMaxWidth || tmpSizeY > resultMaxHeight) {
				imgScaledSizeX /= 2d;
				imgScaledSizeY /= 2d;
				tmpSizeX /= 2d;
				tmpSizeY /= 2d;
				tileTmpSizeX /= 2d;
				tileTmpSizeY /= 2d;
				resolution *= 2;
				scaleFactor /= 2d;
			}

			System.out.println("resolution scale: " + resolution);

			int outTileSizeX = (int) Math.round(tileTmpSizeX);
			int outTileSizeY = (int) Math.round(tileTmpSizeY);
			int outSizeX = ((int) (tile.width / tileSize)) * outTileSizeX;
			int outSizeY = ((int) (tile.height / tileSize)) * outTileSizeY;
			if (tile.width % tileSize > 0)
				outSizeX += (int) Math.round((double) (tile.width % tileSize) * scaleFactor);
			if (tile.height % tileSize > 0)
				outSizeY += (int) Math.round((double) (tile.height % tileSize) * scaleFactor);

			System.out.println("Result image size: " + outSizeX + "px*" + outSizeY + "px. Tile size: " + outTileSizeX + "px*"
			    + outTileSizeY + "px");

			Sequence result = new Sequence(new IcyBufferedImage(outSizeX, outSizeY, imgSizeC, imgDataType));
			result.beginUpdate();
			IcyBufferedImage resultImg = result.getFirstImage();

			TileLoaderThread[] threads = new TileLoaderThread[nProc];
			Point[] points = new Point[nProc];
			int currProc = 0;
			
			int totalTilesX = (tile.width / tileSize) + (tile.width % tileSize > 0? 1: 0);
			int totalTilesY = (tile.height / tileSize) + (tile.height % tileSize > 0? 1: 0);
			int totalTiles = totalTilesX * totalTilesY;
			int treatedTiles = 0;
			
			for (int posX = tile.x, posTileX = 0; posX < tile.x + tile.width && posX < imgSizeX;) {
				int tileSizeX = (posX + tileSize) <= imgSizeX ? tileSize : imgSizeX - posX;
				tileSizeX = (posX + tileSizeX) <= tile.x + tile.width ? tileSizeX : tile.x + tile.width - posX;
				int outCurrTileSizeX = (posTileX + outTileSizeX) <= outSizeX ? outTileSizeX : outSizeX - posTileX;

				for (int posY = tile.y, posTileY = 0; posY < tile.y + tile.height && posY < imgSizeY;) {
					// System.out.println("Processing: (" + posX + ", " + posY + ")");

					int tileSizeY = (posY + tileSize) <= imgSizeY ? tileSize : imgSizeY - posY;
					tileSizeY = (posY + tileSizeY) <= tile.y + tile.height ? tileSizeY : tile.y + tile.height - posY;
					int outCurrTileSizeY = (posTileY + outTileSizeY) <= outSizeY ? outTileSizeY : outSizeY - posTileY;

					threads[currProc] = new TileLoaderThread(path, new Rectangle(posX, posY, tileSizeX, tileSizeY),
					    new Dimension(outCurrTileSizeX, outCurrTileSizeY));
					points[currProc] = new Point(posTileX, posTileY);
					threads[currProc++].start();

					if (currProc >= nProc || ((posX + tileSizeX) >= imgSizeX && (posY + tileSizeY) >= imgSizeY)
					    || ((posX + tileSizeX) >= tile.x + tile.width && (posY + tileSizeY) >= tile.y + tile.height)) {
						for (int p = 0; p < currProc; p++) {
							try {
								threads[p].join();
								resultImg.copyData(threads[p].resultImage, null, new Point(points[p].x, points[p].y));
								threads[p] = null;
								points[p] = null;
								treatedTiles++;
								progress = (double)treatedTiles/(double)totalTiles;
								setProgress(progress);
								setStatusMessage(String.format("Loading image: %d, tile: %d / %d", (int)(progress*100), treatedTiles, totalTiles));
							} catch (InterruptedException e) {
								e.printStackTrace();
							}
						}
						currProc = 0;
					}

					posY += tileSizeY;
					posTileY += outCurrTileSizeY;
				}
				posX += tileSizeX;
				posTileX += outCurrTileSizeX;
			}
			System.out.println(treatedTiles);
			result.setImage(0, 0, resultImg);
			result.dataChanged();
			result.endUpdate();
			result.setName(imgName);
			Runtime.getRuntime().gc();
			return result;
		} catch (UnsupportedFormatException | IOException e) {
			throw e;
		} finally {
			importer.close();
		}
	}

	/**
	 * Class to load tiles with multi-threading
	 * 
	 * @author Daniel Felipe Gonzalez Obando
	 */
	private static class TileLoaderThread extends Thread {

		private final String path;
		private final Rectangle rect;
		private final Dimension dimension;
		private IcyBufferedImage resultImage;

		/**
		 * Constructor
		 * 
		 * @param path
		 *          Path of the image to be loaded
		 * @param rect
		 *          Tile to be loaded
		 * @param resultImage
		 *          final scaled dimensions of the result image. Used for
		 *          down/up-sampling.
		 */
		public TileLoaderThread(String path, Rectangle rect, Dimension dimension) {
			super();
			this.path = path;
			this.rect = rect;
			this.dimension = dimension;
			this.resultImage = null;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Thread#run()
		 */
		@Override
		public void run() {
			LociImporterPlugin importer = new LociImporterPlugin();
			try {
				importer.open(path, 0);
				resultImage = importer.getImage(0, 1, rect, 0, 0);
				resultImage = IcyBufferedImageUtil.scale(resultImage, (int) dimension.getWidth(), (int) dimension.getHeight());
				importer.close();
			} catch (UnsupportedFormatException | IOException e) {
				e.printStackTrace();
			}
		}
	}
}
