package danyfel80.registration.bspline.classic;

import java.awt.Point;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;

import algorithms.danyfel80.io.sequence.cursor.IcyBufferedImageCursor;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.image.IcyBufferedImageUtil.FilterType;
import icy.main.Icy;
import icy.sequence.Sequence;
import icy.type.DataType;

/**
 * Different tools for the bUnwarpJ interface.
 */
public class MiscTools {
	//------------------------------------------------------------------
	/**
	 * Scale an image with good quality in both up and down direction
	 */
	final static public IcyBufferedImage scale(final IcyBufferedImage source, final double scale) {
		if (scale == 1.0f)
			return IcyBufferedImageUtil.getCopy(source);
		else if (scale < 1.0f)
			return IcyBufferedImageUtil.scale(source, (int) (source.getWidth() * scale), (int) (source.getHeight() * scale),
					FilterType.BILINEAR);
		else {
			return IcyBufferedImageUtil.scale(source, (int) Math.round(source.getWidth() * scale),
					(int) Math.round(source.getHeight() * scale), FilterType.BILINEAR);
		}
	}
	//------------------------------------------------------------------

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
	static public void drawArrow(double[][] canvas, int x1, int y1, int x2, int y2, double color, int arrow_size) {
		drawLine(canvas, x1, y1, x2, y2, color);
		int arrow_size2 = 2 * arrow_size;

		// Do not draw the arrow_head if the arrow is very small
		if ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) < arrow_size * arrow_size)
			return;

		// Vertical arrow
		if (x2 == x1) {
			if (y2 > y1) {
				drawLine(canvas, x2, y2, x2 - arrow_size, y2 - arrow_size2, color);
				drawLine(canvas, x2, y2, x2 + arrow_size, y2 - arrow_size2, color);
			} else {
				drawLine(canvas, x2, y2, x2 - arrow_size, y2 + arrow_size2, color);
				drawLine(canvas, x2, y2, x2 + arrow_size, y2 + arrow_size2, color);
			}
		}

		// Horizontal arrow
		else if (y2 == y1) {
			if (x2 > x1) {
				drawLine(canvas, x2, y2, x2 - arrow_size2, y2 - arrow_size, color);
				drawLine(canvas, x2, y2, x2 - arrow_size2, y2 + arrow_size, color);
			} else {
				drawLine(canvas, x2, y2, x2 + arrow_size2, y2 - arrow_size, color);
				drawLine(canvas, x2, y2, x2 + arrow_size2, y2 + arrow_size, color);
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
			Point p2 = new Point(-arrow_size2, -arrow_size);
			Point p3 = new Point(-arrow_size2, +arrow_size);

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

	//------------------------------------------------------------------
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
					Point(canvas, n, x1, color);
				return;
			} else {
				for (int n = y1; n <= y2; n++)
					Point(canvas, n, x1, color);
				return;
			}

		int dy = y2 - y1;
		if (dy == 0)
			if (x1 > x2) {
				for (int n = x2; n <= x1; n++)
					Point(canvas, y1, n, color);
				return;
			} else {
				for (int n = x1; n <= x2; n++)
					Point(canvas, y1, n, color);
				return;
			}

		double m = (double) dy / dx;

		if (m > 1 || m < -1) {
			temp = x1;
			x1 = y1;
			y1 = temp;
			temp = x2;
			x2 = y2;
			y2 = temp;
			dx = x2 - x1;
			dy = y2 - y1;
			m = (double) dy / dx;
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
			m = (double) dy / dx;
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
		Point(canvas, y, x, color);
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
			Point(canvas, y, x, color);
			x = tempx;
			y = tempy;
		}
	}

	//------------------------------------------------------------------
	/**
	 * Put the image from an ImageProcessor into a double array.
	 *
	 * @param ip
	 *          input, origin of the image
	 * @param image
	 *          output, the image in a double array
	 */
	static public void extractImage(final IcyBufferedImage ip, double image[]) {
		if (ip.getSizeC() == 1) {
			IcyBufferedImage doubleImage = IcyBufferedImageUtil.convertToType(ip, DataType.DOUBLE, false);
			System.arraycopy(doubleImage.getDataXYAsDouble(0), 0, image, 0, image.length);
		} else {
			IcyBufferedImageCursor cursor = new IcyBufferedImageCursor(ip);
			for (int y = 0; y < ip.getHeight(); y++) {
				int yAccum = y * ip.getWidth();
				for (int x = 0; x < ip.getWidth(); x++) {
					double averageIntensity = 0d;
					for (int c = 0; c < ip.getSizeC(); c++) {
						averageIntensity += cursor.get(x, y, c);
					}
					averageIntensity /= ip.getSizeC();
					image[yAccum + x] = averageIntensity;
				}
			}

		}
	}

	//------------------------------------------------------------------
	/**
	 * Load a transformation from a file.
	 *
	 * @param filename
	 *          transformation file name
	 * @param cx
	 *          x- B-spline coefficients
	 * @param cy
	 *          y- B-spline coefficients
	 */
	public static void loadTransformation(String filename, final double[][] cx, final double[][] cy) {
		try {
			final FileReader fr = new FileReader(filename);
			final BufferedReader br = new BufferedReader(fr);
			String line;

			// Read number of intervals
			line = br.readLine();
			int lineN = 1;
			StringTokenizer st = new StringTokenizer(line, "=");
			if (st.countTokens() != 2) {
				br.close();
				fr.close();
				System.out.println("Line " + lineN + "+: Cannot read number of intervals");
				return;
			}
			st.nextToken();
			int intervals = Integer.valueOf(st.nextToken()).intValue();

			// Skip next 2 lines
			line = br.readLine();
			line = br.readLine();
			lineN += 2;

			// Read the cx coefficients
			for (int i = 0; i < intervals + 3; i++) {
				line = br.readLine();
				lineN++;
				st = new StringTokenizer(line);
				if (st.countTokens() != intervals + 3) {
					br.close();
					fr.close();
					System.out.println("Line " + lineN + ": Cannot read enough coefficients");
					return;
				}
				for (int j = 0; j < intervals + 3; j++)
					cx[i][j] = Double.valueOf(st.nextToken()).doubleValue();
			}

			// Skip next 2 lines
			line = br.readLine();
			line = br.readLine();
			lineN += 2;

			// Read the cy coefficients
			for (int i = 0; i < intervals + 3; i++) {
				line = br.readLine();
				lineN++;
				st = new StringTokenizer(line);
				if (st.countTokens() != intervals + 3) {
					br.close();
					fr.close();
					System.out.println("Line " + lineN + ": Cannot read enough coefficients");
					return;
				}
				for (int j = 0; j < intervals + 3; j++)
					cy[i][j] = Double.valueOf(st.nextToken()).doubleValue();
			}
			fr.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.err.println("File not found exception" + e);
			return;
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("IOException exception" + e);
			return;
		} catch (NumberFormatException e) {
			e.printStackTrace();
			System.err.println("Number format exception" + e);
			return;
		}
	}

	//------------------------------------------------------------------
	/**
	 * Save the elastic transformation.
	 *
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param cx
	 *          x- deformation coefficients
	 * @param cy
	 *          y- deformation coefficients
	 * @param filename
	 *          transformation file name
	 */
	public static void saveElasticTransformation(int intervals, double[][] cx, double[][] cy, String filename) {

		// Save the file
		try {
			final FileWriter fw = new FileWriter(filename);
			String aux;
			fw.write("Intervals=" + intervals + "\n\n");
			fw.write("X Coeffs -----------------------------------\n");
			for (int i = 0; i < intervals + 3; i++) {
				for (int j = 0; j < intervals + 3; j++) {
					aux = "" + cx[i][j];
					while (aux.length() < 21)
						aux = " " + aux;
					fw.write(aux + " ");
				}
				fw.write("\n");
			}
			fw.write("\n");
			fw.write("Y Coeffs -----------------------------------\n");
			for (int i = 0; i < intervals + 3; i++) {
				for (int j = 0; j < intervals + 3; j++) {
					aux = "" + cy[i][j];
					while (aux.length() < 21)
						aux = " " + aux;
					fw.write(aux + " ");
				}
				fw.write("\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("IOException exception" + e);
		} catch (SecurityException e) {
			e.printStackTrace();
			System.err.println("Security exception" + e);
		}
	}

	//------------------------------------------------------------------
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
	static public void Point(double[][] canvas, int y, int x, double color) {
		if (y < 0 || y >= canvas.length)
			return;
		if (x < 0 || x >= canvas[0].length)
			return;
		canvas[y][x] = color;
	}

	//------------------------------------------------------------------
	/**
	 * Show an image in a new bUnwarpJ window.
	 *
	 * @param title
	 *          image title
	 * @param array
	 *          image in a double array
	 */
	public static void showImage(final String title, final double[][] array) {
		int Ydim = array.length;
		int Xdim = array[0].length;

		final IcyBufferedImage fp = new IcyBufferedImage(Xdim, Ydim, 1, DataType.DOUBLE);
		IcyBufferedImageCursor fpCursor = new IcyBufferedImageCursor(fp);
		for (int i = 0; i < Ydim; i++)
			for (int j = 0; j < Xdim; j++)
				fpCursor.setSafe(j, i, 0, array[i][j]);
		fpCursor.commitChanges();
		final Sequence ip = new Sequence(title, fp);
		Icy.getMainInterface().addSequence(ip);
	} // end showImage
}
