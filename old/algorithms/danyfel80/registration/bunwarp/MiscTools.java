package algorithms.danyfel80.registration.bunwarp;

import java.awt.Point;
import java.awt.Rectangle;

import icy.common.listener.DetailedProgressListener;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array2DUtil;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class MiscTools {

	/**
	 * Put the image from an IcyBufferedImage into a double array.
	 *
	 * @param ibi
	 *          input, origin of the image
	 * @param image
	 *          output, the image in a double array
	 */
	public static void extractImage(IcyBufferedImage ibi, double[] outData) {
		int width = ibi.getWidth();
		int height = ibi.getHeight();
		int channels = ibi.getSizeC();
		IcyBufferedImage ibiD = IcyBufferedImageUtil.convertToType(ibi, DataType.DOUBLE, false);
		double[][] ibiDData = ibiD.getDataXYCAsDouble();
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				double val = 0;
				for (int c = 0; c < channels; c++) {
					val += ibiDData[c][x + y * width];
				}

				outData[x + y * width] = val / (double) channels;
			}
		}
	}

	public static IcyBufferedImage scale(IcyBufferedImage ibi, float scaleFactor) {
		return IcyBufferedImageUtil.scale(ibi, (int) Math.round(scaleFactor * (double) ibi.getWidth()),
				(int) Math.round(scaleFactor * (double) ibi.getHeight()));
	}

	/**
	 * Draw a line between two points. Bresenham's algorithm.
	 *
	 * @param canvas
	 *          canvas where we are painting
	 * @param x1
	 *          x- coordinate for first point
	 * @param y1
	 *          y- coordinate for first point
	 * @param x2
	 *          x- coordinate for second point
	 * @param y2
	 *          y- coordinate for second point
	 * @param color
	 *          line color
	 */
	static public void drawLine(double[][] canvas, int x1, int y1, int x2, int y2, double color) {
		int temp;
		int dy_neg = 1;
		int dx_neg = 1;
		int switch_x_y = 0;
		int neg_slope = 0;
		int tempx, tempy;
		int dx = x2 - x1;
		if (dx == 0)
			if (y1 > y2) {
				for (int n = y2; n <= y1; n++)
					drawPoint(canvas, n, x1, color);
				return;
			} else {
				for (int n = y1; n <= y2; n++)
					drawPoint(canvas, n, x1, color);
				return;
			}

		int dy = y2 - y1;
		if (dy == 0)
			if (x1 > x2) {
				for (int n = x2; n <= x1; n++)
					drawPoint(canvas, y1, n, color);
				return;
			} else {
				for (int n = x1; n <= x2; n++)
					drawPoint(canvas, y1, n, color);
				return;
			}

		float m = (float) dy / dx;

		if (m > 1 || m < -1) {
			temp = x1;
			x1 = y1;
			y1 = temp;
			temp = x2;
			x2 = y2;
			y2 = temp;
			dx = x2 - x1;
			dy = y2 - y1;
			m = (float) dy / dx;
			switch_x_y = 1;
		}

		if (x1 > x2) {
			temp = x1;
			x1 = x2;
			x2 = temp;
			temp = y1;
			y1 = y2;
			y2 = temp;
			dx = x2 - x1;
			dy = y2 - y1;
			m = (float) dy / dx;
		}

		if (m < 0) {
			if (dy < 0) {
				dy_neg = -1;
				dx_neg = 1;
			} else {
				dy_neg = 1;
				dx_neg = -1;
			}
			neg_slope = 1;
		}

		int d = 2 * (dy * dy_neg) - (dx * dx_neg);
		int incrH = 2 * dy * dy_neg;
		int incrHV = 2 * ((dy * dy_neg) - (dx * dx_neg));
		int x = x1;
		int y = y1;
		tempx = x;
		tempy = y;

		if (switch_x_y == 1) {
			temp = x;
			x = y;
			y = temp;
		}
		drawPoint(canvas, y, x, color);
		x = tempx;
		y = tempy;

		while (x < x2) {
			if (d <= 0) {
				x++;
				d += incrH;
			} else {
				d += incrHV;
				x++;
				if (neg_slope == 0)
					y++;
				else
					y--;
			}
			tempx = x;
			tempy = y;

			if (switch_x_y == 1) {
				temp = x;
				x = y;
				y = temp;
			}
			drawPoint(canvas, y, x, color);
			x = tempx;
			y = tempy;
		}
	}

	/**
	 * Plot a point in a canvas.
	 *
	 * @param canvas
	 *          canvas where we are painting
	 * @param x
	 *          x- coordinate for the point
	 * @param y
	 *          y- coordinate for the point
	 * @param color
	 *          point color
	 */
	static public void drawPoint(double[][] canvas, int y, int x, double color) {
		if (y < 0 || y >= canvas.length)
			return;
		if (x < 0 || x >= canvas[0].length)
			return;
		canvas[y][x] = color;
	}

	/**
	 * Draw an arrow between two points. The arrow head is in (x2,y2)
	 *
	 * @param canvas
	 *          canvas where we are painting
	 * @param x1
	 *          x- coordinate for the arrow origin
	 * @param y1
	 *          y- coordinate for the arrow origin
	 * @param x2
	 *          x- coordinate for the arrow head
	 * @param y2
	 *          y- coordinate for the arrow head
	 * @param color
	 *          arrow color
	 * @param arrow_size
	 *          arrow size
	 */
	public static void drawArrow(double[][] canvas, int x1, int y1, int x2, int y2, double color, int arrowSize) {
		drawLine(canvas, x1, y1, x2, y2, color);
		int arrow_size2 = 2 * arrowSize;

		// Do not draw the arrow_head if the arrow is very small
		if ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) < arrowSize * arrowSize)
			return;

		// Vertical arrow
		if (x2 == x1) {
			if (y2 > y1) {
				drawLine(canvas, x2, y2, x2 - arrowSize, y2 - arrow_size2, color);
				drawLine(canvas, x2, y2, x2 + arrowSize, y2 - arrow_size2, color);
			} else {
				drawLine(canvas, x2, y2, x2 - arrowSize, y2 + arrow_size2, color);
				drawLine(canvas, x2, y2, x2 + arrowSize, y2 + arrow_size2, color);
			}
		}

		// Horizontal arrow
		else if (y2 == y1) {
			if (x2 > x1) {
				drawLine(canvas, x2, y2, x2 - arrow_size2, y2 - arrowSize, color);
				drawLine(canvas, x2, y2, x2 - arrow_size2, y2 + arrowSize, color);
			} else {
				drawLine(canvas, x2, y2, x2 + arrow_size2, y2 - arrowSize, color);
				drawLine(canvas, x2, y2, x2 + arrow_size2, y2 + arrowSize, color);
			}
		}

		// Now we need to rotate the arrow head about the origin
		else {
			// Calculate the angle of rotation and adjust for the quadrant
			double t1 = Math.abs(new Integer(y2 - y1).doubleValue());
			double t2 = Math.abs(new Integer(x2 - x1).doubleValue());
			double theta = Math.atan(t1 / t2);
			if (x2 < x1) {
				if (y2 < y1)
					theta = Math.PI + theta;
				else
					theta = -(Math.PI + theta);
			} else if (x2 > x1 && y2 < y1)
				theta = 2 * Math.PI - theta;
			double cosTheta = Math.cos(theta);
			double sinTheta = Math.sin(theta);

			// Create the other points and translate the arrow to the origin
			Point p2 = new Point(-arrow_size2, -arrowSize);
			Point p3 = new Point(-arrow_size2, +arrowSize);

			// Rotate the points (without using matrices!)
			int x = new Long(Math.round((cosTheta * p2.x) - (sinTheta * p2.y))).intValue();
			p2.y = new Long(Math.round((sinTheta * p2.x) + (cosTheta * p2.y))).intValue();
			p2.x = x;
			x = new Long(Math.round((cosTheta * p3.x) - (sinTheta * p3.y))).intValue();
			p3.y = new Long(Math.round((sinTheta * p3.x) + (cosTheta * p3.y))).intValue();
			p3.x = x;

			// Translate back to desired location and add to polygon
			p2.translate(x2, y2);
			p3.translate(x2, y2);
			drawLine(canvas, x2, y2, p2.x, p2.y, color);
			drawLine(canvas, x2, y2, p3.x, p3.y, color);
		}
	}

	public static void applyTransformationToSourceMT(Sequence source, Sequence target, int intervals, double[][] cx,
			double[][] cy, DetailedProgressListener progressListener) {

		IcyBufferedImage result_imp = applyTransformationMT(source, target, intervals, cx, cy, progressListener);

		source.beginUpdate();
		source.setImage(0, 0, result_imp);
		source.dataChanged();
		source.endUpdate();
	}

	/**
	 * Apply a given B-spline transformation to the source (gray-scale) image. The
	 * result image is return. The target image is used to know the output size
	 * (Multi-thread version).
	 *
	 * @param sourceSeq
	 *          source image representation
	 * @param targetSeq
	 *          target image representation
	 * @param intervals
	 *          intervals in the deformation
	 * @param cx
	 *          x- B-spline coefficients
	 * @param cy
	 *          y- B-spline coefficients
	 * 
	 * @return result transformed image
	 */
	public static IcyBufferedImage applyTransformationMT(Sequence sourceSeq, Sequence targetSeq, int intervals,
			double[][] cx, double[][] cy, DetailedProgressListener progressListener) {
		final int targetHeight = targetSeq.getHeight();
		final int targetWidth = targetSeq.getWidth();

		// Compute the deformation
		// Set these coefficients to an interpolator
		BSplineModel swx = new BSplineModel(cx, progressListener);
		BSplineModel swy = new BSplineModel(cy, progressListener);

		// Compute the warped image

		BSplineModel[] sourceModels = new BSplineModel[sourceSeq.getSizeC()];
		for (int c = 0; c < sourceSeq.getSizeC(); c++) {
			sourceModels[c] = new BSplineModel(IcyBufferedImageUtil.extractChannel(sourceSeq.getFirstImage(), c), false, 1,
					progressListener);
			sourceModels[c].setPyramidDepth(0);
			sourceModels[c].startPyramids();
		}

		// Join threads
		try {
			for (int c = 0; c < sourceSeq.getSizeC(); c++) {
				sourceModels[c].join();
			}
		} catch (InterruptedException e) {
			System.out.println("Unexpected interruption exception " + e);
		}

		// Calculate warped RGB image
		IcyBufferedImage cp = new IcyBufferedImage(targetWidth, targetHeight, sourceSeq.getSizeC(),
				sourceSeq.getDataType_());

		// Check the number of processors in the computer
		int nproc = Runtime.getRuntime().availableProcessors();

		// We will use threads to display parts of the output image
		int blockHeight = targetHeight / nproc;
		if (targetHeight % 2 != 0)
			blockHeight++;

		int nThreads = nproc;

		Thread[] threads = new Thread[nThreads];
		Rectangle[] rects = new Rectangle[nThreads];
		IcyBufferedImage[] fpTiles = new IcyBufferedImage[nThreads];

		for (int i = 0; i < nThreads; i++) {
			// last block size is the rest of the window
			int y_start = i * blockHeight;

			if (nThreads - 1 == i)
				blockHeight = targetHeight - i * blockHeight;

			rects[i] = new Rectangle(0, y_start, targetWidth, blockHeight);

			// IJ.log("block = 0 " + (i*block_height) + " " + targetWidth + " " +
			// block_height );

			fpTiles[i] = new IcyBufferedImage(rects[i].width, rects[i].height, sourceSeq.getSizeC(),
					sourceSeq.getDataType_());

			threads[i] = new Thread(new ColorApplyTransformTile(swx, swy, sourceModels, targetWidth, targetHeight, intervals,
					rects[i], fpTiles[i]));
			threads[i].start();
		}

		for (int i = 0; i < nThreads; i++) {
			try {
				threads[i].join();
				threads[i] = null;
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		cp.beginUpdate();
		for (int i = 0; i < nThreads; i++) {
			for (int c = 0; c < cp.getSizeC(); c++) {
				cp.copyData(fpTiles[i], null, new Point(rects[i].x, rects[i].y), c, c);
			}
			fpTiles[i] = null;
			rects[i] = null;
		}

		cp.dataChanged();
		cp.endUpdate();
		return cp;

	}

	/**
	 * Class to apply transformation to color images in a concurrent way
	 * 
	 */
	private static class ColorApplyTransformTile implements Runnable {
		/** B-spline deformation in x */
		final BSplineModel swx;
		/** B-spline deformation in y */
		final BSplineModel swy;
		/** source image model */
		final BSplineModel[] sourceModels;
		/** target current width */
		final int targetCurrentWidth;
		/** target current height */
		final int targetCurrentHeight;
		/** number of intervals between B-spline coefficients */
		final int intervals;
		/** area of the image to be transformed */
		final Rectangle rect;
		/** resulting float processor for the red channel */
		final private IcyBufferedImage ibi;

		/**
		 * Constructor for color image transform
		 * 
		 * @param swx
		 *          B-spline deformation in x
		 * @param swy
		 *          B-spline deformation in y
		 * @param sourceR
		 *          red source image
		 * @param targetCurrentWidth
		 *          target current width
		 * @param targetCurrentHeight
		 *          target current height
		 * @param intervals
		 *          number of intervals between B-spline coefficients
		 * @param rect
		 *          area of the image to be transformed
		 * @param fpR
		 *          red channel processor to be updated
		 * @param fpG
		 *          green channel processor to be updated
		 * @param fpB
		 *          blue channel processor to be updated
		 */
		ColorApplyTransformTile(BSplineModel swx, BSplineModel swy, BSplineModel[] sourceModels, int targetCurrentWidth,
				int targetCurrentHeight, int intervals, Rectangle rect, IcyBufferedImage ibi) {
			this.swx = swx;
			this.swy = swy;
			this.sourceModels = sourceModels;
			this.targetCurrentWidth = targetCurrentWidth;
			this.targetCurrentHeight = targetCurrentHeight;
			this.intervals = intervals;
			this.rect = rect;
			this.ibi = ibi;
		}

		// ------------------------------------------------------------------
		/**
		 * Run method to update the intermediate window. Only the part defined by
		 * the rectangle will be updated (in this thread).
		 */
		public void run() {
			// Compute the warped image
			int auxTargetHeight = rect.y + rect.height;
			int auxTargetWidth = rect.x + rect.width;

			double[][] ibiArray = Array2DUtil.arrayToDoubleArray(ibi.getDataXYC(), ibi.isSignedDataType());

			final int sourceWidth = sourceModels[0].getWidth();
			final int sourceHeight = sourceModels[0].getHeight();
			Rectangle srcLimits = new Rectangle();
			boolean firstLim = true;
			for (int v_rect = 0, v = rect.y; v < auxTargetHeight; v++, v_rect++) {
				final int v_offset = v_rect * rect.width;
				final double tv = (double) (v * intervals) / (double) (targetCurrentHeight - 1) + 1.0F;

				for (int u_rect = 0, u = rect.x; u < auxTargetWidth; u++, u_rect++) {
					final double tu = (double) (u * intervals) / (double) (targetCurrentWidth - 1) + 1.0F;

					final double x = swx.prepareForInterpolationAndInterpolateI(tu, tv, false, false);
					final double y = swy.prepareForInterpolationAndInterpolateI(tu, tv, false, false);

					if (firstLim) {
						firstLim = false;
						srcLimits.x = (int) x;
						srcLimits.y = (int) y;
						srcLimits.width = (int) x;
						srcLimits.height = (int) y;
					} else {
						srcLimits.x = Math.min(srcLimits.x, (int) x);
						srcLimits.y = Math.min(srcLimits.y, (int) y);
						srcLimits.width = Math.max(srcLimits.width, (int) x);
						srcLimits.height = Math.max(srcLimits.height, (int) y);
					}

					if (x >= 0 && x < sourceWidth && y >= 0 && y < sourceHeight) {
						for (int c = 0; c < ibi.getSizeC(); c++) {
							ibiArray[c][u_rect + v_offset] = sourceModels[c].prepareForInterpolationAndInterpolateI(x, y, false,
									false);
						}
					} else {
						for (int c = 0; c < ibi.getSizeC(); c++) {
							ibiArray[c][u_rect + v_offset] = 0;
						}
					}

				}
			}
			System.out.println("limits: " + srcLimits);
			Array2DUtil.doubleArrayToSafeArray(ibiArray, ibi.getDataXYC(), ibi.isSignedDataType());
			ibi.dataChanged();
		} // end run method

	}

}
