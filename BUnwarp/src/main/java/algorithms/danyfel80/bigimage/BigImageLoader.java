package algorithms.danyfel80.bigimage;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FilenameUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Node;

import icy.common.exception.UnsupportedFormatException;
import icy.gui.frame.progress.ProgressFrame;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import icy.type.point.Point5D;
import icy.util.XMLUtil;
import ome.xml.meta.OMEXMLMetadata;
import plugins.adufour.ezplug.EzGUI;
import plugins.kernel.importer.LociImporterPlugin;
import plugins.kernel.roi.roi2d.ROI2DArea;

/**
 * This class allows easy big image loading through downsampling and tile
 * loading.
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageLoader {

	private EzGUI pluginGUI = null;
	private ProgressFrame progressFrame = null;
	private boolean isInterrupted;

	public void setPluginGUI(EzGUI pluginGUI) {
		this.pluginGUI = pluginGUI;
	}

	private void setProgress(double progress) {
		if (pluginGUI != null) {
			this.pluginGUI.setProgressBarValue(progress);
		}
	}

	private void setStatusMessage(String message) {
		if (pluginGUI != null) {
			this.pluginGUI.setProgressBarMessage(message);
		}
	}

	public void interrupt() {
		this.isInterrupted = true;
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
	public Sequence loadDownsampledImage(String path, Rectangle tile, int resultMaxWidth, int resultMaxHeight,
			boolean showProgressBar) throws UnsupportedFormatException, IOException {
		this.isInterrupted = false;
		double progress = 0;
		setProgress(progress);
		if (showProgressBar) {
			this.progressFrame = new ProgressFrame("Loading image...");
		}

		String xmlPath = FilenameUtils.getFullPath(path) + FilenameUtils.getBaseName(path) + ".xml";
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadata imgProps = importer.getOMEXMLMetaData();
			String imgName = FilenameUtils.getBaseName(path);
			int imgSizeX = imgProps.getPixelsSizeX(0).getValue();
			int imgSizeY = imgProps.getPixelsSizeY(0).getValue();
			Dimension imgSize = new Dimension(imgSizeX, imgSizeY);
			int imgSizeC = imgProps.getPixelsSizeC(0).getValue();
			DataType imgDataType = DataType.getDataTypeFromPixelType(imgProps.getPixelsType(0));
			if (tile == null) {
				tile = new Rectangle(0, 0, imgSizeX, imgSizeY);
			} else {
				if (tile.x < 0) {
					tile.width += tile.x;
					tile.x = 0;
				}
				if (tile.y < 0) {
					tile.height += tile.y;
					tile.y = 0;
				}
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

			System.gc();
			long ram = Runtime.getRuntime().freeMemory();
			int nProc = Runtime.getRuntime().availableProcessors();
			if (showProgressBar) {
				System.out.println("Available memory: " + ram + " bytes, Available processors: " + nProc);
			}
			ram /= imgSizeC;
			double szMax = Math.sqrt(ram);
			// szMax /= 2;

			int tileSize = (int) Math.ceil(szMax / nProc);

			while (tile.width / resultMaxWidth > 2 * tileSize)
				tileSize *= 2;
			while (tile.height / resultMaxHeight > 2 * tileSize)
				tileSize *= 2;
			System.out.println("Tile size: " + tileSize);

			// double imgScaledSizeX = imgSizeX;
			// double imgScaledSizeY = imgSizeY;
			double tmpSizeX = tile.getWidth();
			double tmpSizeY = tile.getHeight();
			double tileTmpSizeX = tileSize;
			double tileTmpSizeY = tileSize;
			double scaleFactor = 1d;
			while (tmpSizeX > resultMaxWidth || tmpSizeY > resultMaxHeight) {
				// imgScaledSizeX /= 2d;
				// imgScaledSizeY /= 2d;
				tmpSizeX /= 2d;
				tmpSizeY /= 2d;
				tileTmpSizeX /= 2d;
				tileTmpSizeY /= 2d;
				scaleFactor /= 2d;
			}
			if (showProgressBar) {
				System.out.println("output resolution: " + scaleFactor);
			}
			int outTileSizeX = (int) Math.round(tileTmpSizeX);
			int outTileSizeY = (int) Math.round(tileTmpSizeY);
			int outSizeX = ((int) (tile.width / tileSize)) * outTileSizeX;
			int outSizeY = ((int) (tile.height / tileSize)) * outTileSizeY;
			if (tile.width % tileSize > 0)
				outSizeX += (int) Math.round((double) (tile.width % tileSize) * scaleFactor);
			if (tile.height % tileSize > 0)
				outSizeY += (int) Math.round((double) (tile.height % tileSize) * scaleFactor);
			if (showProgressBar) {
				System.out.println("Result image size: " + outSizeX + "px*" + outSizeY + "px. Tile size: " + outTileSizeX
						+ "px*" + outTileSizeY + "px");
			}
			Sequence result = new Sequence(new IcyBufferedImage(outSizeX, outSizeY, imgSizeC, imgDataType));
			result.beginUpdate();
			IcyBufferedImage resultImg = result.getFirstImage();

			TileLoaderThread[] threads = new TileLoaderThread[nProc];
			Point[] points = new Point[nProc];
			int currProc = 0;

			int totalTilesX = (tile.width / tileSize) + (tile.width % tileSize > 0 ? 1 : 0);
			int totalTilesY = (tile.height / tileSize) + (tile.height % tileSize > 0 ? 1 : 0);
			int totalTiles = totalTilesX * totalTilesY;
			int treatedTiles = 0;

			for (int posX = tile.x, posTileX = 0; posX < tile.x + tile.width && posX < imgSizeX;) {
				int tileSizeX = (posX + tileSize) <= imgSizeX ? tileSize : imgSizeX - posX;
				tileSizeX = (posX + tileSizeX) <= tile.x + tile.width ? tileSizeX : tile.x + tile.width - posX;
				int outCurrTileSizeX = (posTileX + outTileSizeX) <= outSizeX ? outTileSizeX : outSizeX - posTileX;

				for (int posY = tile.y, posTileY = 0; posY < tile.y + tile.height && posY < imgSizeY;) {
					// System.out.println("Processing: (" + posX + ", " + posY + ")");
					if (isInterrupted) {
						break;
					}

					int tileSizeY = (posY + tileSize) <= imgSizeY ? tileSize : imgSizeY - posY;
					tileSizeY = (posY + tileSizeY) <= tile.y + tile.height ? tileSizeY : tile.y + tile.height - posY;
					int outCurrTileSizeY = (posTileY + outTileSizeY) <= outSizeY ? outTileSizeY : outSizeY - posTileY;

					threads[currProc] = new TileLoaderThread(path, new Rectangle(posX, posY, tileSizeX, tileSizeY),
							new Dimension(outCurrTileSizeX, outCurrTileSizeY), imgSize);
					points[currProc] = new Point(posTileX, posTileY);
					threads[currProc++].start();

					if (currProc >= nProc || ((posX + tileSizeX) >= imgSizeX && (posY + tileSizeY) >= imgSizeY)
							|| ((posX + tileSizeX) >= tile.x + tile.width && (posY + tileSizeY) >= tile.y + tile.height)) {
						for (int p = 0; p < currProc; p++) {
							try {
								progress = (double) treatedTiles / (double) totalTiles;
								setProgress(progress);
								if (showProgressBar) {
									this.progressFrame.setPosition(progress * 100);
									this.progressFrame.setMessage(String.format("Loading image: %d%%, tile: %d/%d",
											(int) (progress * 100), treatedTiles, totalTiles));
								}
								setStatusMessage(String.format("Loading image: %d%%, tile: %d/%d", (int) (progress * 100), treatedTiles,
										totalTiles));
								threads[p].join();
								resultImg.copyData(threads[p].resultImage, null, new Point(points[p].x, points[p].y));
								threads[p] = null;
								points[p] = null;
								treatedTiles++;

							} catch (InterruptedException e) {
								e.printStackTrace();
							}
						}
						currProc = 0;
					}
					posY += tileSizeY;
					posTileY += outCurrTileSizeY;
				}
				if (isInterrupted) {
					break;
				}
				posX += tileSizeX;
				posTileX += outCurrTileSizeX;
			}
			if (showProgressBar) {
				System.out.println("Treated Tiles" + treatedTiles);
			}
			result.setImage(0, 0, resultImg);
			// Load ROIs
			result.addROIs(loadROIs(xmlPath, tile, scaleFactor), false);
			result.dataChanged();
			result.endUpdate();
			result.setName(imgName);
			Runtime.getRuntime().gc();
			return result;
		} catch (UnsupportedFormatException | IOException e) {
			throw e;
		} finally {
			importer.close();
			setStatusMessage("Done loading");
			setProgress(1);
			if (showProgressBar) {
				this.progressFrame.dispose();
			}
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
		private final Dimension outDimension;
		private final Dimension fullDimension;

		private IcyBufferedImage resultImage;
		private final double scaleX;
		private final double scaleY;

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
		public TileLoaderThread(String path, Rectangle rect, Dimension outDimension, Dimension fullDimension) {
			super();
			this.path = path;
			this.rect = rect;
			this.outDimension = outDimension;
			this.fullDimension = fullDimension;
			this.resultImage = null;
			this.scaleX = rect.width / outDimension.width;
			this.scaleY = rect.height / outDimension.height;
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

				// get bigger tile to avoid interpolation issues
				Dimension scaledExtractedDim = new Dimension(outDimension);
				Point scaledPosition = new Point(0, 0);
				Rectangle bigRect = new Rectangle(rect);
				if (rect.x - 3 >= 0) {
					bigRect.x -= 3;
					scaledExtractedDim.width += (int) (3 / scaleX);
					scaledPosition.x += (int) (3 / scaleX);
				}
				if (rect.y - 3 >= 0) {
					bigRect.y -= 3;
					scaledExtractedDim.height += (int) (3 / scaleY);
					scaledPosition.y += (int) (3 / scaleY);
				}
				if (rect.x + rect.width + 3 < fullDimension.width) {
					bigRect.width += 3;
					scaledExtractedDim.width += (int) (3 / scaleX);
				}
				if (rect.y + rect.height + 3 < fullDimension.height) {
					bigRect.height += 3;
					scaledExtractedDim.height += (int) (3 / scaleY);
				}

				resultImage = importer.getImage(0, 0, bigRect, 0, 0);
				if (resultImage.getWidth() != scaledExtractedDim.getWidth()
						|| resultImage.getHeight() != scaledExtractedDim.height) {
					resultImage = IcyBufferedImageUtil.scale(resultImage, scaledExtractedDim.width, scaledExtractedDim.height);
				}
				IcyBufferedImageUtil.getSubImage(resultImage, scaledPosition.x, scaledPosition.y, outDimension.width,
						outDimension.height);
				importer.close();
			} catch (UnsupportedFormatException | IOException e) {
				e.printStackTrace();
			}
		}
	}

	public ROI2D loadDownsampledMask(Sequence srcSeq, String srcPath, Rectangle tile, int resultMaxWidth,
			int resultMaxHeight, boolean showProgressBar) throws UnsupportedFormatException, IOException {

		String maskPath = FilenameUtils.getFullPath(srcPath);
		maskPath += FilenameUtils.getBaseName(srcPath) + "_mask.";
		maskPath += FilenameUtils.getExtension(srcPath);

		BooleanMask2D boolMask = new BooleanMask2D();

		if (new File(maskPath).exists()) {
			Sequence maskSeq = loadDownsampledImage(maskPath, tile, resultMaxWidth, resultMaxHeight, showProgressBar);
			double[] maskData = Array1DUtil.arrayToDoubleArray(maskSeq.getDataXY(0, 0, 0), maskSeq.isSignedDataType());
			boolMask.mask = new boolean[maskData.length];

			for (int i = 0; i < maskData.length; i++) {
				boolMask.mask[i] = maskData[i] > 0 ? true : false;
			}

			boolMask.bounds = new Rectangle(maskSeq.getBounds2D());
		} else {
			boolMask.mask = new boolean[srcSeq.getSizeX() * srcSeq.getSizeY()];

			for (int i = 0; i < boolMask.mask.length; i++) {
				boolMask.mask[i] = true;
			}

			boolMask.bounds = new Rectangle(srcSeq.getBounds2D());
		}

		return new ROI2DArea(boolMask);
	}

	public List<ROI> loadROIs(String xmlPath, Rectangle rect, double scale) {
		File xmlFile = new File(xmlPath);
		Document doc = XMLUtil.loadDocument(xmlFile);
		Node root = XMLUtil.getRootElement(doc);
		if (root == null)
			return new ArrayList<ROI>();
		Node rois = XMLUtil.getChild(root, "rois");
		List<ROI> roiList = (List<ROI>) ((rois != null) ? ROI.loadROIsFromXML(rois) : new ArrayList<ROI>());
		roiList = BigImageUtil.getROIsInTile(roiList, rect);
		for (ROI roi : roiList) {
			Point5D pos = roi.getPosition5D();
			pos.setX(pos.getX() - rect.x);
			pos.setY(pos.getY() - rect.y);
			roi.setPosition5D(pos);
			BigImageUtil.scaleROI(roi, scale);
		}

		return roiList;
	}
}
