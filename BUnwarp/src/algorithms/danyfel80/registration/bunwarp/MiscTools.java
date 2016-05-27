package algorithms.danyfel80.registration.bunwarp;

import java.awt.Point;

import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.type.DataType;

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
			p2.x = new Long(Math.round((sinTheta * p2.x) + (cosTheta * p2.y))).intValue();
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

}
