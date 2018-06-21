/**
 * 
 */
package algorithms.danyfel80.registration.bunwarp.big;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.io.FilenameUtils;

import algorithms.danyfel80.bigimage.BigImageSaver;
import algorithms.danyfel80.registration.bunwarp.BSplineModel;
import icy.common.exception.UnsupportedFormatException;
import icy.common.listener.DetailedProgressListener;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;
import icy.type.collection.array.Array2DUtil;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import ome.xml.meta.OMEXMLMetadata;
import plugins.kernel.importer.LociImporterPlugin;

/**
 * Toolbox to make registration on big images. (working by tiles and downsampled
 * images)
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageTools {

	/**
	 * Gets the size of the image specified by the given path.
	 * 
	 * @param path
	 *          Path of the image
	 * @return Image dimension
	 */
	public static Dimension getSequenceSize(String path) {
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadata imgProps = importer.getOMEXMLMetaData();
			int imgSizeX = imgProps.getPixelsSizeX(0).getValue();
			int imgSizeY = imgProps.getPixelsSizeY(0).getValue();

			return new Dimension(imgSizeX, imgSizeY);

		} catch (UnsupportedFormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				importer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return null;
	}

	/**
	 * Gets the channel count of the image specified by the given path.
	 * 
	 * @param path
	 *          Path of the image
	 * @return Image dimension. 0 if image is not correctly read.
	 */
	public static int getSequenceChannelCount(String path) {
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadata imgProps = importer.getOMEXMLMetaData();
			return imgProps.getPixelsSizeC(0).getValue();

		} catch (UnsupportedFormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				importer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return 0;
	}

	/**
	 * Gets the data type of the image specified by the given path.
	 * 
	 * @param path
	 *          Path of the image
	 * @return Data type of the image.
	 */
	public static DataType getSequenceDataType(String path) {
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadata imgProps = importer.getOMEXMLMetaData();
			return DataType.getDataTypeFromPixelType(imgProps.getPixelsType(0));

		} catch (UnsupportedFormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				importer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return null;
	}

	/**
	 * Loads the sequence specified by path subsampled to fit into a 500px by
	 * 500px dimension.
	 * 
	 * @param path
	 *          Image file path
	 * @param name
	 *          Name of the resulting sequence
	 * @return subsampled sequence. null is returned if the sequence is not
	 *         correctly loaded.
	 */
	public static Sequence loadSubsampledSequence(String path, String name) {
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadata imgProps = importer.getOMEXMLMetaData();
			int imgSizeX = imgProps.getPixelsSizeX(0).getValue();
			int imgSizeY = imgProps.getPixelsSizeY(0).getValue();

			int maxResolution = Math.max(imgSizeX, imgSizeY);
			int resolution = 0;
			while (maxResolution > 2000) {
				maxResolution /= 2;
				resolution++;
			}

			IcyBufferedImage imageDS = importer.getImage(0, resolution, 0, 0);
			Sequence sequenceDS = new Sequence(imageDS);
			if (maxResolution > 1000) {
				sequenceDS = SequenceUtil.scale(sequenceDS, (int) Math.round(imageDS.getSizeX() / 2.0),
						(int) Math.round(imageDS.getSizeY() / 2.0));
			}
			sequenceDS.setName(name);
			return sequenceDS;

		} catch (UnsupportedFormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				importer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return null;
	}

	/**
	 * Loads a tile from the image stored at the specified path.
	 * 
	 * @param path
	 *          Image file path
	 * @param name
	 *          Result sequence name
	 * @param tile
	 *          Area to load
	 * @return Image tile or null if the image is not correctly loaded.
	 */
	public static Sequence loadImageTile(String path, String name, Rectangle tile) {
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			IcyBufferedImage imageDS = importer.getImage(0, 0, tile, 0, 0);
			Sequence sequenceDS = new Sequence(imageDS);
			sequenceDS.setName(name);
			return sequenceDS;
		} catch (UnsupportedFormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				importer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return null;
	}

	/**
	 * Computes the area (bounding box) used in the source image for the
	 * transformation given the B-spline model.
	 * 
	 * @param intervals
	 *          Number of intervals in the B-spline model.
	 * @param cx
	 *          X coordinate B-spline model coefficients.
	 * @param cy
	 *          Y coordinate B-spline model coefficients.
	 * @param targetFullSize
	 *          target image full size.
	 * @return Bounding box in the source image used to transform into the target
	 *         image.
	 */
	public static Rectangle computeTransformationUsedArea(int intervals, double[][] cx, double[][] cy,
			Dimension targetFullSize, DetailedProgressListener progressListener) {
		// Load B-Spline transformation model
		BSplineModel swx = new BSplineModel(cx, progressListener);
		BSplineModel swy = new BSplineModel(cy, progressListener);

		// Check the number of processors in the computer
		int nproc = Runtime.getRuntime().availableProcessors();

		// Threads (tiles) setup
		int blockHeight = targetFullSize.height / nproc;
		if (targetFullSize.height % 2 != 0)
			blockHeight++;

		int nThreads = nproc;

		Thread[] threads = new Thread[nThreads];
		Rectangle[] rects = new Rectangle[nThreads];
		Rectangle[] usedBoxes = new Rectangle[nThreads];

		// Threads execution
		for (int i = 0; i < nThreads; i++) {
			// last block size is the rest of the window
			int y0 = i * blockHeight;

			if (nThreads - 1 == i)
				blockHeight = targetFullSize.height - i * blockHeight;

			rects[i] = new Rectangle(0, y0, targetFullSize.width, blockHeight);
			usedBoxes[i] = new Rectangle();

			threads[i] = new Thread(
					new TransformBoundingBoxComputationOnTile(swx, swy, intervals, targetFullSize, rects[i], usedBoxes[i]));
			threads[i].start();
		}

		// Wait until all threads have finished.
		for (int i = 0; i < nThreads; i++) {
			try {
				threads[i].join();
				threads[i] = null;
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		// Compute full bounding box
		Rectangle usedBox = new Rectangle(usedBoxes[0]);

		for (int i = 1; i < usedBoxes.length; i++) {
			usedBox.x = Math.min(usedBox.x, usedBoxes[i].x);
			usedBox.y = Math.min(usedBox.y, usedBoxes[i].y);
			usedBox.width = Math.max(usedBox.width, usedBoxes[i].width);
			usedBox.height = Math.max(usedBox.height, usedBoxes[i].height);

			rects[i] = null;
			usedBoxes[i] = null;
		}

		return usedBox;
	}

	/**
	 * Class to compute the bounding box of the transformation on a tile. Used as
	 * a thread.
	 * 
	 * @author Daniel Felipe Gonzalez Obando
	 */
	private static class TransformBoundingBoxComputationOnTile implements Runnable {
		/** B-spline deformation in x */
		final BSplineModel swx;
		/** B-spline deformation in y */
		final BSplineModel swy;
		/** number of intervals between B-spline coefficients */
		final int intervals;
		/** target current size */
		final Dimension targetFullSize;
		/** area of the image to be transformed */
		final Rectangle tileBox;
		/** bounding box used for the transformation of this tile */
		final Rectangle tileBoundingBox;

		/**
		 * Constructor
		 * 
		 * @param swx
		 *          B-spline deformation in x
		 * @param swy
		 *          B-spline deformation in y
		 * @param intervals
		 *          Number of intervals between B-spline coefficients
		 * @param targetFullSize
		 *          Target current size
		 * @param tileBox
		 *          Area of the image to be transformed
		 * @param tileBoundingBox
		 *          Bounding box used for the transformation of this tile
		 */
		public TransformBoundingBoxComputationOnTile(BSplineModel swx, BSplineModel swy, int intervals,
				Dimension targetFullSize, Rectangle tileBox, Rectangle tileBoundingBox) {
			this.swx = swx;
			this.swy = swy;
			this.intervals = intervals;
			this.targetFullSize = targetFullSize;
			this.tileBox = tileBox;
			this.tileBoundingBox = tileBoundingBox;
		}

		/**
		 * Execution of the thread to compute warped coordinates.
		 * 
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run() {
			// Compute warped coordinates
			int tileEndHeight = tileBox.y + tileBox.height;
			int tileEndWidth = tileBox.x + tileBox.width;

			for (int v_rect = 0, v = tileBox.y; v < tileEndHeight; v++, v_rect++) {
				final double tv = (double) (v * intervals) / (double) (targetFullSize.height - 1) + 1.0F;

				for (int u_rect = 0, u = tileBox.x; u < tileEndWidth; u++, u_rect++) {
					final double tu = (double) (u * intervals) / (double) (targetFullSize.width - 1) + 1.0F;

					final double x = swx.prepareForInterpolationAndInterpolateI(tu, tv, false, false);
					final double y = swy.prepareForInterpolationAndInterpolateI(tu, tv, false, false);

					// Stock bounding positions
					if (v_rect == 0 && u_rect == 0) {
						tileBoundingBox.x = (int) Math.floor(x);
						tileBoundingBox.y = (int) Math.floor(y);
						tileBoundingBox.width = (int) Math.ceil(x);
						tileBoundingBox.height = (int) Math.ceil(y);
					} else {
						tileBoundingBox.x = Math.min(tileBoundingBox.x, (int) Math.floor(x));
						tileBoundingBox.y = Math.min(tileBoundingBox.y, (int) Math.floor(y));
						tileBoundingBox.width = Math.max(tileBoundingBox.width, (int) Math.ceil(x));
						tileBoundingBox.height = Math.max(tileBoundingBox.height, (int) Math.ceil(y));
					}
				}
			}

			// Correct from positions to dimension
			tileBoundingBox.width = tileBoundingBox.width - tileBoundingBox.x;
			tileBoundingBox.height = tileBoundingBox.height - tileBoundingBox.y;
		}

	}

	public static Sequence applyTransformationToImage(String srcResultPath, String srcPath, String transformedSrcPath,
			String tgtPath, int intervals, double[][] cx, double[][] cy, Dimension registeredTgtDimension,
			DetailedProgressListener progressListener) {

		int srcChannels = getSequenceChannelCount(srcPath);
		DataType srcDataType = getSequenceDataType(srcPath);

		Dimension transformedSrcDimension = getSequenceSize(transformedSrcPath);
		int transformedSrcChannels = getSequenceChannelCount(transformedSrcPath);
		// System.out.println("src:" + srcPath);
		// System.out.println("channel in src:" + srcChannels);
		DataType transformedSrcDataType = getSequenceDataType(transformedSrcPath);
		Dimension tgtDimension = getSequenceSize(tgtPath);
		// int tgtChannels = getSequenceChannelCount(tgtPath);
		// DataType tgtDataType = getSequenceDataType(tgtPath);

		// Load B-Spline transformation model
		BSplineModel swx = new BSplineModel(cx, progressListener);
		BSplineModel swy = new BSplineModel(cy, progressListener);

		// Calculate tile's dimension and count
		Dimension tileDimension = new Dimension(500, 500);
		Dimension tgtTileCount = new Dimension((int) Math.ceil((double) tgtDimension.width / (double) tileDimension.width),
				(int) Math.ceil((double) tgtDimension.height / (double) tileDimension.height));
		int totalTgtTileCount = tgtTileCount.width * tgtTileCount.height;

		// Create full size result image
		Sequence resultSeq = new Sequence(new IcyBufferedImage(tgtTileCount.width * tileDimension.width,
				tgtTileCount.height * tileDimension.height, transformedSrcChannels, transformedSrcDataType));
		resultSeq.beginUpdate();
		double[][] resultData = Array2DUtil.arrayToDoubleArray(resultSeq.getDataXYC(0, 0), resultSeq.isSignedDataType());
		for (int c = 0; c < resultData.length; c++) {
			for (int y = 0; y < resultSeq.getHeight(); y++) {
				int yOff = y * resultSeq.getSizeX();
				for (int x = 0; x < resultSeq.getWidth(); x++) {
					resultData[c][x + yOff] = 100;
				}
			}
		}
		// Calculate Thread number
		int numProc = Runtime.getRuntime().availableProcessors();
		TileTransformProcessing[] threads = new TileTransformProcessing[numProc];

		// until all tiles have been treated
		int tileNo = 0;
		while (tileNo < totalTgtTileCount) {

			int thr;
			// treat as many tiles as the amount of available processors.
			for (thr = 0; thr < numProc && tileNo < totalTgtTileCount; thr++, tileNo++) {
				Rectangle tileRect = new Rectangle((tileNo % tgtTileCount.width) * tileDimension.width,
						(tileNo / tgtTileCount.width) * tileDimension.height, tileDimension.width, tileDimension.height);

				// // Get tile from source image
				// Rectangle srcRect = computeTransformSourceTileArea(intervals, swx,
				// swy, tileRect, srcDimension, tgtDimension,
				// registeredTgtDimension);
				// // System.out.println("tile no " + tileNo + "area: " + srcRect);
				// // System.out.println("tile no " + tileNo + "area: " + tileRect);
				// IcyBufferedImage srcTile = loadImageTile(srcPath, "img",
				// srcRect).getFirstImage();
				// // System.out.println("value "+ srcTile.getChannelBounds(0)[0] + " "
				// +
				// // srcTile.getChannelBounds(0)[1]);
				// BSplineModel[] srcModels = new BSplineModel[srcTile.getSizeC()];
				// for (int c = 0; c < srcModels.length; c++) {
				// srcModels[c] = new
				// BSplineModel(IcyBufferedImageUtil.extractChannel(srcTile, c), false,
				// 1);
				// srcModels[c].setPyramidDepth(0);
				// srcModels[c].startPyramids();
				// }
				// // Join threads
				// try {
				// for (int c = 0; c < srcModels.length; c++) {
				// srcModels[c].join();
				// }
				// } catch (InterruptedException e) {
				// System.out.println("Unexpected interruption exception " + e);
				// }
				// create thread to process tile
				threads[thr] = new TileTransformProcessing(intervals, swx, swy, tgtDimension, registeredTgtDimension, tileRect,
						tileNo, srcPath, srcChannels, srcDataType, transformedSrcPath, transformedSrcChannels,
						transformedSrcDataType, transformedSrcDimension, new Point(0, 0));
			}
			int usedThr = thr;

			for (thr = 0; thr < usedThr; thr++) {
				threads[thr].start();
			}

			// wait for tiles to be treated

			for (thr = 0; thr < usedThr; thr++) {
				TileTransformProcessing thread = threads[thr];
				try {
					thread.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
				// System.out.println("printing tile " + thread.tileRect);
				// write tiles to the full result sequence.
				double[][] resultTileData = Array2DUtil.arrayToDoubleArray(thread.transformedResultTile.getDataXYC(),
						thread.transformedResultTile.isSignedDataType());
				for (int y = 0, v = thread.tileRect.y; y < thread.tileRect.height; y++, v++) {
					int y_offset = y * thread.tileRect.width;
					int v_offset = v * resultSeq.getWidth();
					for (int x = 0, u = thread.tileRect.x; x < thread.tileRect.width; x++, u++) {
						for (int c = 0; c < resultData.length; c++) {
							resultData[c][u + v_offset] = resultTileData[c][x + y_offset];
						}
					}
				}
				threads[thr] = null;
			}

		}

		Array2DUtil.doubleArrayToSafeArray(resultData, resultSeq.getDataXYC(0, 0), resultSeq.isSignedDataType());
		resultSeq.dataChanged();
		resultSeq.endUpdate();
		return resultSeq;
	}

	private static Rectangle computeTransformSourceTileArea(int intervals, BSplineModel swx, BSplineModel swy,
			Rectangle rect, Dimension srcDimension, Dimension tgtDimension, Dimension registeredTgtDimension) {
		Rectangle area = new Rectangle();
		double minX = 0, maxX = 0, minY = 0, maxY = 0;
		boolean first = true;

		double factorX = (double) tgtDimension.width / (double) registeredTgtDimension.width;
		double factorY = (double) tgtDimension.height / (double) registeredTgtDimension.height;

		int tileEndHeight = rect.y + rect.height;
		int tileEndWidth = rect.x + rect.width;

		for (int v = rect.y - 2; v < tileEndHeight + 2; v++) {
			final double tv = (double) (v * intervals) / (double) (tgtDimension.height - 1) + 1.0F;

			for (int u = rect.x - 2; u < tileEndWidth + 2; u++) {
				final double tu = (double) (u * intervals) / (double) (tgtDimension.width - 1) + 1.0F;

				final double x = swx.prepareForInterpolationAndInterpolateI(tu, tv, false, false) * factorX;
				final double y = swy.prepareForInterpolationAndInterpolateI(tu, tv, false, false) * factorY;

				// Stock valid bounding positions
				if (x >= 0 && x < srcDimension.width && y >= 0 && y < srcDimension.height) {
					if (first) {
						first = false;
						minX = x;
						minY = y;
						maxX = x;
						maxY = y;
					} else {
						minX = Math.min(minX, x);
						minY = Math.min(minY, y);
						maxX = Math.max(maxX, x);
						maxY = Math.max(maxY, y);
					}
				}
			}
		}
		// Correct from positions to dimension
		area.x = first ? 0 : (int) Math.floor(minX);
		area.y = first ? 0 : (int) Math.floor(minY);
		area.width = first ? 0 : (int) Math.ceil(maxX - minX);
		area.height = first ? 0 : (int) Math.ceil(maxY - minY);
		return area;
	}

	private static class TileTransformProcessing extends Thread {
		/** B-spline deformation in x */
		private final BSplineModel swx;
		/** B-spline deformation in y */
		private final BSplineModel swy;
		/** number of intervals between B-spline coefficients */
		private final int intervals;
		/** dimension of the full target image */
		private final Dimension tgtFullDim;
		/** dimension of the registered target image */
		private final Dimension tgtRegisteredDim;
		/** target tile rectangle */
		private final Rectangle tileRect;
		/** result tile after processing */
		private IcyBufferedImage transformedResultTile;
		/** result source tile after processing */
		private IcyBufferedImage resultTile;
		/** result tile mask after processing */
		private IcyBufferedImage maskTile;
		/// ** number of the tile processed */
		private final int tileNo;
		// /** source image tile used for processing the target tile */
		// private final IcyBufferedImage srcTile;
		/** the source image path */
		private final String srcPath;
		/** the source tile channel count */
		private final int srcChannels;
		/** the source tile data type */
		private final DataType srcDataType;

		// /** models to interpolate colors */
		// private final BSplineModel[] srcModels;
		// /** source tile rectangle */
		// private final Rectangle srcTileRect;

		/** the transformed source path */
		private final String transformedSrcPath;
		/** the source tile channel count */
		private final int transformedSrcChannels;
		/** the source tile data type */
		private final DataType transformedSrcDataType;
		/** the source tile image dimension */
		private final Dimension transformedSrcDimension;

		/** the position in the full size image */
		private final Point fullSizePosition;

		// intervals, swx, swy, tgtDimension, registeredTgtDimension,tileRect/* ,
		// tileNo */, srcChannels, srcDataType, srcModels, srcRect
		public TileTransformProcessing(final int intervals, final BSplineModel swx, final BSplineModel swy,
				final Dimension tgtFullDim, final Dimension tgtRegisteredDim, final Rectangle tileRect, final int tileNo,
				String srcPath, final int srcChannels, final DataType srcDataType, String transformedSrcPath,
				final int transformedSrcChannels, final DataType transformedSrcDataType,
				final Dimension transformedSrcDimension, final Point fullSizePosition) {
			this.swx = swx;
			this.swy = swy;
			this.intervals = intervals;
			this.tgtFullDim = tgtFullDim;
			this.tgtRegisteredDim = tgtRegisteredDim;
			this.tileRect = tileRect;
			this.tileNo = tileNo;
			this.srcPath = srcPath;
			this.srcChannels = srcChannels;
			this.srcDataType = srcDataType;
			this.transformedSrcPath = transformedSrcPath;
			this.transformedSrcChannels = transformedSrcChannels;
			this.transformedSrcDataType = transformedSrcDataType;
			this.transformedSrcDimension = transformedSrcDimension;
			this.fullSizePosition = (fullSizePosition == null) ? new Point(0, 0) : fullSizePosition;
		}

		public int getTileNo() {
			return this.tileNo;
		}

		@Override
		public void run() {
			// Get tile from source image
			Rectangle srcRect = computeTransformSourceTileArea(intervals, swx, swy, tileRect, transformedSrcDimension,
					tgtFullDim, tgtRegisteredDim);
			// System.out.println("tile no " + tileNo + "area: " + srcRect);
			// System.out.println("tile no " + tileNo + "area: " + tileRect);
			IcyBufferedImage srcTile;
			IcyBufferedImage transformedSrcTile;
			srcRect.x += fullSizePosition.x;
			srcRect.y += fullSizePosition.y;
			if (srcRect.width * srcRect.height > 0) {
				srcTile = loadImageTile(srcPath, "img1", srcRect).getFirstImage();
				transformedSrcTile = loadImageTile(transformedSrcPath, "img", srcRect).getFirstImage();
			} else {
				transformedSrcTile = new IcyBufferedImage(2, 2, transformedSrcChannels, transformedSrcDataType);
				srcTile = new IcyBufferedImage(2, 2, srcChannels, transformedSrcDataType);
			}
			srcRect.x -= fullSizePosition.x;
			srcRect.y -= fullSizePosition.y;

			double[][] srcTileData = Array2DUtil.arrayToDoubleArray(srcTile.getDataXYC(), srcTile.isSignedDataType());
			double[][] transformedSrcTileData = Array2DUtil.arrayToDoubleArray(transformedSrcTile.getDataXYC(),
					transformedSrcTile.isSignedDataType());
			// compute color mode for interpolation in transformed source image
			/*
			 * BSplineModel[] transformedSrcModels = new
			 * BSplineModel[transformedSrcTile.getSizeC()]; for (int c = 0; c <
			 * transformedSrcModels.length; c++) { transformedSrcModels[c] = new
			 * BSplineModel(IcyBufferedImageUtil.extractChannel(transformedSrcTile,
			 * c), false, 1); transformedSrcModels[c].setPyramidDepth(0);
			 * transformedSrcModels[c].startPyramids(); } transformedSrcTile = null;
			 * // Join threads try { for (int c = 0; c < transformedSrcModels.length;
			 * c++) { transformedSrcModels[c].join(); } } catch (InterruptedException
			 * e) { System.out.println("Unexpected interruption exception " + e); } //
			 * compute color mode for interpolation in source image BSplineModel[]
			 * srcModels = new BSplineModel[srcTile.getSizeC()]; for (int c = 0; c <
			 * srcModels.length; c++) { srcModels[c] = new
			 * BSplineModel(IcyBufferedImageUtil.extractChannel(srcTile, c), false,
			 * 1); srcModels[c].setPyramidDepth(0); srcModels[c].startPyramids(); }
			 * srcTile = null; // Join threads try { for (int c = 0; c <
			 * srcModels.length; c++) { srcModels[c].join(); } } catch
			 * (InterruptedException e) {
			 * System.out.println("Unexpected interruption exception " + e); }
			 */
			double factorX = (double) tgtFullDim.width / (double) tgtRegisteredDim.width;
			double factorY = (double) tgtFullDim.height / (double) tgtRegisteredDim.height;

			// Compute warped coordinates
			int tileEndHeight = tileRect.y + tileRect.height;
			int tileEndWidth = tileRect.x + tileRect.width;

			transformedResultTile = new IcyBufferedImage(tileRect.width, tileRect.height, transformedSrcChannels,
					transformedSrcDataType);
			resultTile = new IcyBufferedImage(tileRect.width, tileRect.height, srcChannels, srcDataType);
			maskTile = new IcyBufferedImage(tileRect.width, tileRect.height, 1, DataType.UBYTE);

			transformedResultTile.beginUpdate();
			resultTile.beginUpdate();
			maskTile.beginUpdate();
			double[][] transformedResultData = Array2DUtil.arrayToDoubleArray(transformedResultTile.getDataXYC(),
					transformedResultTile.isSignedDataType());
			double[][] resultData = Array2DUtil.arrayToDoubleArray(resultTile.getDataXYC(), resultTile.isSignedDataType());
			byte[] maskData = maskTile.getDataXYAsByte(0);
			Rectangle srcLimits = new Rectangle();
			boolean firstLim = true;
			for (int v_rect = 0, v = tileRect.y; v < tileEndHeight; v++, v_rect++) {
				final int v_offset = v_rect * tileRect.width;
				final double tv = (double) (v * intervals) / (double) (tgtFullDim.height - 1) + 1.0F;

				for (int u_rect = 0, u = tileRect.x; u < tileEndWidth; u++, u_rect++) {
					final double tu = (double) (u * intervals) / (double) (tgtFullDim.width - 1) + 1.0F;

					final double x = swx.prepareForInterpolationAndInterpolateI(tu, tv, false, false) * factorX;
					final double y = swy.prepareForInterpolationAndInterpolateI(tu, tv, false, false) * factorY;

					double srcTileX = x - srcRect.x;
					double srcTileY = y - srcRect.y;

					srcTileX = Math.round(srcTileX);
					srcTileY = Math.round(srcTileY);

					if (firstLim) {
						firstLim = false;
						srcLimits.x = (int) srcTileX;
						srcLimits.y = (int) srcTileY;
						srcLimits.width = (int) srcTileX;
						srcLimits.height = (int) srcTileY;
					} else {
						srcLimits.x = Math.min(srcLimits.x, (int) srcTileX);
						srcLimits.y = Math.min(srcLimits.y, (int) srcTileY);
						srcLimits.width = Math.max(srcLimits.width, (int) srcTileX);
						srcLimits.height = Math.max(srcLimits.height, (int) srcTileY);
					}

					// Stock valid source tile positions
					if (srcTileX >= 0 && srcTileX < srcRect.width && srcTileY >= 0 && srcTileY < srcRect.height) {
						maskData[u_rect + v_offset] = (byte) DataType.UBYTE.getMaxValue();
						for (int c = 0; c < transformedResultTile.getSizeC(); c++) {
							transformedResultData[c][u_rect + v_offset] = transformedSrcTileData[c][(int) (srcTileX
									+ (srcTileY * transformedSrcTile.getSizeX()))];
							// transformedResultData[c][u_rect + v_offset] =
							// transformedSrcModels[c]
							// .prepareForInterpolationAndInterpolateI(srcTileX, srcTileY,
							// false, false);
						}
						for (int c = 0; c < resultTile.getSizeC(); c++) {
							resultData[c][u_rect + v_offset] = srcTileData[c][(int) (srcTileX + (srcTileY * srcTile.getSizeX()))];
							// resultData[c][u_rect + v_offset] =
							// srcModels[c].prepareForInterpolationAndInterpolateI(srcTileX,
							// srcTileY,
							// false, false);
						}
					} else {
						maskData[u_rect + v_offset] = (byte) 0;
						for (int c = 0; c < transformedResultTile.getSizeC(); c++) {
							transformedResultData[c][u_rect + v_offset] = 0;
						}
						for (int c = 0; c < resultTile.getSizeC(); c++) {
							resultData[c][u_rect + v_offset] = 0;
						}
					}
				}
			}
			// System.out.println("srcLimits:" + srcLimits);
			maskTile.dataChanged();
			maskTile.endUpdate();
			Array2DUtil.doubleArrayToSafeArray(transformedResultData, transformedResultTile.getDataXYC(),
					transformedResultTile.isSignedDataType());
			transformedResultTile.dataChanged();
			transformedResultTile.endUpdate();
			Array2DUtil.doubleArrayToSafeArray(resultData, resultTile.getDataXYC(), resultTile.isSignedDataType());
			resultTile.dataChanged();
			resultTile.endUpdate();
		}
	}

	public static Sequence getInterpolatedImage(Sequence srcSeq, DetailedProgressListener progressListener) {
		IcyBufferedImage ibi = srcSeq.getFirstImage();
		BSplineModel models[] = new BSplineModel[srcSeq.getSizeC()];

		for (int i = 0; i < models.length; i++) {
			models[i] = new BSplineModel(IcyBufferedImageUtil.extractChannel(ibi, i), false, 1, progressListener);
			models[i].setPyramidDepth(0);
			models[i].startPyramids();
		}
		for (int i = 0; i < models.length; i++) {
			try {
				models[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		Sequence result = new Sequence(
				new IcyBufferedImage(ibi.getWidth(), ibi.getHeight(), ibi.getSizeC(), ibi.getDataType_()));
		result.beginUpdate();
		double[][] resultData = Array2DUtil.arrayToDoubleArray(result.getDataXYC(0, 0), result.isSignedDataType());
		for (int y = 0; y < ibi.getSizeY(); y++) {
			int yOff = y * ibi.getSizeX();
			for (int x = 0; x < ibi.getSizeX(); x++) {
				for (int c = 0; c < resultData.length; c++) {
					resultData[c][x + yOff] = models[c].prepareForInterpolationAndInterpolateI(x, y, false, false);
				}
			}
		}
		Array2DUtil.doubleArrayToSafeArray(resultData, result.getDataXYC(0, 0), result.isSignedDataType());
		result.dataChanged();
		result.endUpdate();
		return result;
	}

	public static void applyAndSaveTransformationToBigImage(String srcResultPath, String transformedSrcResultPath,
			String srcPath, String transformedSrcPath, String tgtPath, int intervals, double[][] cx, double[][] cy,
			Dimension registeredTgtDimension, Rectangle resultTile, DetailedProgressListener progressListener)
			throws ServiceException, IOException, FormatException, InterruptedException, ExecutionException {
		notifyProgress(progressListener, 0d, "Transforming original image...");
		Dimension transformedSrcDimension = getSequenceSize(transformedSrcPath);
		int transformedSrcChannels = getSequenceChannelCount(transformedSrcPath);
		// System.out.println("src:" + srcPath);
		// System.out.println("channel in src:" + srcChannels);
		DataType transformedSrcDataType = getSequenceDataType(transformedSrcPath);
		int srcChannels = getSequenceChannelCount(srcPath);
		DataType srcDataType = getSequenceDataType(srcPath);

		Dimension tgtDimension = getSequenceSize(tgtPath);

		if (resultTile == null) {
			resultTile = new Rectangle(new Point(0, 0), tgtDimension);
		}

		resultTile = resultTile.intersection(new Rectangle(new Point(0, 0), tgtDimension));

		if (resultTile.isEmpty())
			return;

		// int tgtChannels = getSequenceChannelCount(tgtPath);
		// DataType tgtDataType = getSequenceDataType(tgtPath);

		// Load B-Spline transformation model
		BSplineModel swx = new BSplineModel(cx, progressListener);
		BSplineModel swy = new BSplineModel(cy, progressListener);

		// Calculate tile's dimension and count
		System.gc();
		long ram = Runtime.getRuntime().freeMemory();
		int numProc = Math.max(Runtime.getRuntime().availableProcessors() - 1, 1);
		System.out.println("Available memory: " + ram + " bytes, Available processors: " + numProc);
		// Divide by channel count
		ram /= transformedSrcChannels;
		// Divide by processor count
		ram /= numProc;
		double szMax = Math.sqrt(ram);
		// Divide by amount of elements treated at the same time
		// szMax /= 4;

		int tileSize = 256;
		while (szMax > tileSize * 2) {
			tileSize *= 2;
		}
		// tileSize = 500;

		System.out.println("Image size: " + transformedSrcDimension.width + "px*" + transformedSrcDimension.height
				+ "px. Tile size: " + tileSize);

		Dimension tileDimension = new Dimension(tileSize, tileSize);
		while (tileDimension.width > resultTile.width) {
			tileDimension.width /= 2;
		}
		while (tileDimension.height > resultTile.height) {
			tileDimension.height /= 2;
		}
		Dimension tgtTileCount = new Dimension((resultTile.width + tileDimension.width - 1) / tileDimension.width,
				(resultTile.height + tileDimension.height - 1) / tileDimension.height);
		int totalTgtTileCount = tgtTileCount.width * tgtTileCount.height;

		// Create full size result image
		BigImageSaver resultSaver = new BigImageSaver(new File(srcResultPath), resultTile.getSize(), srcChannels,
				srcDataType, tileDimension);
		BigImageSaver transformedResultSaver = new BigImageSaver(new File(transformedSrcResultPath), resultTile.getSize(),
				transformedSrcChannels, transformedSrcDataType, tileDimension);
		String maskResultPath = FilenameUtils.getFullPath(transformedSrcResultPath);
		maskResultPath += FilenameUtils.getBaseName(transformedSrcResultPath);
		maskResultPath += "_mask.";
		maskResultPath += FilenameUtils.getExtension(transformedSrcResultPath);
		BigImageSaver maskSaver = new BigImageSaver(new File(maskResultPath), resultTile.getSize(), 1, DataType.UBYTE,
				tileDimension);

		// Calculate Thread number
		int numThreads = Math.max((numProc + 1) / 2, 1);
		ExecutorService transformExecutor = Executors.newFixedThreadPool(numThreads);
		ExecutorCompletionService<TileTransformProcessing> transformCs = new ExecutorCompletionService<>(transformExecutor);

		// until all tiles have been treated
		int tileNo = 0;
		try {
			while (tileNo < totalTgtTileCount && !Thread.currentThread().isInterrupted()) {

				int thr;
				// treat as many tiles as the amount of available processors.
				for (thr = 0; thr < numThreads && tileNo < totalTgtTileCount
						&& !Thread.currentThread().isInterrupted(); thr++, tileNo++) {
					notifyProgress(progressListener, tileNo / (double) totalTgtTileCount,
							String.format("Processing tile %d of %d", tileNo + 1, totalTgtTileCount));

					Rectangle tileRect = new Rectangle((tileNo % tgtTileCount.width) * tileDimension.width,
							(tileNo / tgtTileCount.width) * tileDimension.height, tileDimension.width, tileDimension.height);
					if (tileRect.x + tileRect.width > resultTile.width) {
						tileRect.width -= tileRect.x + tileRect.width - resultTile.width;
					}
					if (tileRect.y + tileRect.height > resultTile.height) {
						tileRect.height -= tileRect.y + tileRect.height - resultTile.height;
					}

					// create thread to process tile
					TileTransformProcessing task = new TileTransformProcessing(intervals, swx, swy, resultTile.getSize(),
							registeredTgtDimension, tileRect, tileNo, srcPath, srcChannels, srcDataType, transformedSrcPath,
							transformedSrcChannels, transformedSrcDataType, transformedSrcDimension, resultTile.getLocation());
					transformCs.submit(task, task);
				}
				int usedThr = thr;

				// wait for tiles to be treated
				List<TileTransformProcessing> ls = new ArrayList<>();

				for (thr = 0; thr < usedThr; thr++) {
					Future<TileTransformProcessing> transFuture = transformCs.take();
					ls.add(transFuture.get());
				}
				ls.stream().sorted(Comparator.comparing(TileTransformProcessing::getTileNo)).forEach(task -> {
					TileSaveProcessing saveTask = new TileSaveProcessing(task.resultTile, task.transformedResultTile,
							task.maskTile, task.tileRect.getLocation(), resultSaver, transformedResultSaver, maskSaver);
					saveTask.run();
				});
			}
			notifyProgress(progressListener, 1d, "Processes all tiles");
			System.out.println("finished pool");
		} catch (ExecutionException e) {
			throw e;
		} finally {
			transformExecutor.shutdownNow();
			transformExecutor.awaitTermination(0, TimeUnit.SECONDS);
			resultSaver.closeWriter();
			transformedResultSaver.closeWriter();
			maskSaver.closeWriter();
		}

	}

	private static void notifyProgress(DetailedProgressListener progressListener, double progress, String message) {
		progressListener.notifyProgress(progress, message, null);
	}

	private static class TileSaveProcessing implements Runnable {

		/** result tile after processing */
		private IcyBufferedImage resultTile;
		/** transformed result tile after processing */
		private IcyBufferedImage transformedResultTile;
		/** result tile mask after processing */
		private IcyBufferedImage maskTile;
		/** position */
		private Point position;
		/** saver */
		BigImageSaver saver;
		/** saver */
		BigImageSaver transformedSaver;
		/** mask saver */
		BigImageSaver maskSaver;

		/**
		 * @param transformedResultTile
		 * @param position
		 * @param transformedSaver
		 */
		public TileSaveProcessing(IcyBufferedImage resultTile, IcyBufferedImage transformedResultTile,
				IcyBufferedImage maskTile, Point position, BigImageSaver saver, BigImageSaver transformedSaver,
				BigImageSaver maskSaver) {
			this.resultTile = resultTile;
			this.transformedResultTile = transformedResultTile;
			this.maskTile = maskTile;
			this.position = position;
			this.saver = saver;
			this.transformedSaver = transformedSaver;
			this.maskSaver = maskSaver;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Thread#run()
		 */
		@Override
		public void run() {
			try {
				saver.saveTile(resultTile, position);
				transformedSaver.saveTile(transformedResultTile, position);
				maskSaver.saveTile(maskTile, position);
			} catch (ServiceException | IOException | FormatException e) {
				e.printStackTrace();
			}
		}

	}

}
