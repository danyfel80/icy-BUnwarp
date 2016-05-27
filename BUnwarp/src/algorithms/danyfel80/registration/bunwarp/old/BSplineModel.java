package algorithms.danyfel80.registration.bunwarp.old;

import java.util.Stack;

import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.collection.array.Array2DUtil;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class BSplineModel extends Thread {

	private static final int minImageSize = 4;
	public static final int MAX_OUTPUT_SIZE = 1024;
	/**
	 * Image (after scaling)
	 */
	private Sequence seq;
	/**
	 * Original image, full-size without scaling
	 */
	private double[] originalImage;
	/**
	 * Working image at maximum resolution (after scaling)
	 */
	private double[] image;
	/**
	 * Subsampled output image
	 */
	private double[] subImage;
	/**
	 * Image spline coefficients
	 */
	private double[] coefficient;
	/**
	 * Subsampled output image B-spline coefficients
	 */
	private double[] subCoeffs;
	/**
	 * Sub-sampling factor at maximum resolution image (power of 2)
	 */
	private int maxSubsamplingFactor;
	/**
	 * Original image width (at full-resolution, without scaling)
	 */
	private int originalWidth;
	/**
	 * Original image height (at full-resolution, without scaling)
	 */
	private int originalHeight;
	/**
	 * Image width
	 */
	private int width;
	/**
	 * Image height
	 */
	private int height;
	/**
	 * Width of the subsampled output image
	 */
	private int subWidth;
	/**
	 * Height of the subsampled output image
	 */
	private int subHeight;
	/**
	 * Resolution pyramid depth
	 */
	private int pyramidDepth;
	/**
	 * Current pyramid depth
	 */
	private int currentDepth;
	/**
	 * Smallest image width
	 */
	private int smallestWidth;
	/**
	 * Smallest image height
	 */
	private int smallestHeight;

	// Flags
	/**
	 * Flag to check target image
	 */
	private boolean isTarget;
	/**
	 * Flag to check if coefficients are mirrored
	 */
	private boolean coefficientsAreMirrored;
	/**
	 * Flag to check if the output image is subsampled.
	 */
	private boolean isOutputSubsampled;

	// Indexes
	/**
	 * X index
	 */
	int[] xIndex;
	/**
	 * Y index
	 */
	int[] yIndex;

	// Weights of the spline
	/**
	 * X component of the spline weight
	 */
	double[] xWeight;
	/**
	 * Y component of the spline weight
	 */
	double[] yWeight;

	// Weights of the spline derivatives
	/**
	 * X component of the spline derivative weight
	 */
	double[] dxWeight;
	/**
	 * Y component of the spline derivative weight
	 */
	double[] dyWeight;

	// Weights of the spline second derivatives
	/**
	 * X component of the spline second derivative weight
	 */
	double[] d2xWeight;
	/**
	 * Y component of the spline second derivative weight
	 */
	double[] d2yWeight;

	/** Interpolation source (current or original) */
	private boolean fromCurrent;

	// Size of the image used for the interpolation
	/** width of the image used for the interpolation */
	private int widthToUse;
	/** height of the image used for the interpolation */
	private int heightToUse;

	// Stacks for the pyramid of images/coefficients
	/**
	 * stack of coefficients pyramid
	 */
	private final Stack<Object> cPyramid = new Stack<Object>();
	/**
	 * stack of image pyramid
	 */
	private final Stack<Object> imgPyramid = new Stack<Object>();

	// Current image (the size might be different from the original)
	/**
	 * current image (at the current resolution level)
	 */
	private double[] currentImage;
	/**
	 * current image spline coefficients
	 */
	private double[] currentCoefficient;
	/**
	 * current image width
	 */
	private int currentWidth;
	/**
	 * current image height
	 */
	private int currentHeight;

	// Some variables to speedup interpolation (precomputed)
	// All these information is set through prepareForInterpolation()
	// Indexes related
	/** precomputed x index */
	public int precXIndex[][];
	/** precomputed y index */
	public int precYIndex[][];
	// Weights of the splines related
	/** precomputed x component of the weight of the spline */
	private double precXWeight[][];
	/** precomputed y component of the weight of the spline */
	private double precYWeight[][];
	// Weights of the derivatives splines related
	/** precomputed x component of the weight of derivative spline */
	private double precDxWeight[][];
	/** precomputed y component of the weight of derivative spline */
	private double precDyWeight[][];
	// Weights of the second derivatives splines related
	/** precomputed x component of the weight of second derivative spline */
	private double precD2xWeight[][];
	/** precomputed y component of the weight of second derivative spline */
	private double precD2yWeight[][];

	/**
	 * Create image model: image and coefficient pyramid. When calling this
	 * constructor, the thread is not started, to do so, startPyramids must be
	 * called.
	 * 
	 * @param seq
	 * @param isTarget
	 * @param maxSubsamplingFactor
	 */
	public BSplineModel(Sequence seq, boolean isTarget, int maxSubsamplingFactor) {
		// Image
		this.seq = seq;

		// Image info
		this.isTarget = isTarget;
		this.maxSubsamplingFactor = maxSubsamplingFactor;
		width = seq.getWidth();
		height = seq.getHeight();
		coefficientsAreMirrored = true;

		// Speedup arrays
		xIndex = new int[4];
		yIndex = new int[4];
		xWeight = new double[4];
		yWeight = new double[4];
		dxWeight = new double[4];
		dyWeight = new double[4];
		d2xWeight = new double[4];
		d2yWeight = new double[4];
	}

	/**
	 * Initialize the B-spline model from a set of coefficients.
	 *
	 * @param c
	 *          Set of B-spline coefficients
	 */
	public BSplineModel(double[][] c) {
		// Get the size of the input array
		this.currentHeight = height = c.length;
		this.currentWidth = width = c[0].length;
		this.coefficientsAreMirrored = false;

		// Copy the array of coefficients
		coefficient = new double[height * width];
		int k = 0;
		for (int y = 0; y < height; y++, k += width)
			System.arraycopy(c[y], 0, coefficient, k, width);

		// Resize the speedup arrays
		xIndex = new int[4];
		yIndex = new int[4];
		xWeight = new double[4];
		yWeight = new double[4];
		dxWeight = new double[4];
		dyWeight = new double[4];
		d2xWeight = new double[4];
		d2yWeight = new double[4];
	}

	/**
	 * Initialize the B-spline model from a set of coefficients. The same as the
	 * previous function but now the coefficients are in a single row.
	 *
	 * @param c
	 *          Set of B-spline coefficients
	 * @param Ydim
	 *          Y-dimension of the set of coefficients
	 * @param Xdim
	 *          X-dimension of the set of coefficients
	 * @param offset
	 *          Offset of the beginning of the array with respect to the origin of
	 *          c
	 */
	public BSplineModel(final double[] c, final int Ydim, final int Xdim, final int offset) {
		// Get the size of the input array
		currentHeight = height = Ydim;
		currentWidth = width = Xdim;
		coefficientsAreMirrored = false;

		// Copy the array of coefficients
		coefficient = new double[height * width];
		System.arraycopy(c, offset, coefficient, 0, height * width);

		// Resize the speedup arrays
		xIndex = new int[4];
		yIndex = new int[4];
		xWeight = new double[4];
		yWeight = new double[4];
		dxWeight = new double[4];
		dyWeight = new double[4];
		d2xWeight = new double[4];
		d2yWeight = new double[4];
	}
	
	/**
	 * Start coefficient and image pyramids
	 */
	public void startPyramids() {
		subWidth = width;
		subHeight = height;

		// Output window must have a maximum size
		if (this.width > BSplineModel.MAX_OUTPUT_SIZE || this.height > BSplineModel.MAX_OUTPUT_SIZE) {
			isOutputSubsampled = true;
			// Calculate subsampled dimensions
			do {
				subWidth /= 2;
				subHeight /= 2;
			} while (subWidth > BSplineModel.MAX_OUTPUT_SIZE || subHeight > BSplineModel.MAX_OUTPUT_SIZE);
		} else {
			isOutputSubsampled = false;
		}
		// Initialize thread
		super.setDaemon(true);
		// And start it (call run)
		super.start();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Runnable#run()
	 */
	@Override
	public void run() {
		super.run();
		if (image == null && seq != null) {
			// Original image
			originalImage = new double[seq.getWidth() * seq.getHeight()];
			double[][] seqData = Array2DUtil.arrayToDoubleArray(seq.getDataXYC(0, 0), seq.isSignedDataType());
			for (int xy = 0; xy < originalImage.length; xy++) {
				double val = 0;
				for (int c = 0; c < seqData.length; c++) {
					val += seqData[c][xy];
				}
				originalImage[xy] = val / 3d;
			}
			originalHeight = height;
			originalWidth = width;

			// Copy the pixel array and scale if necessary
			if (maxSubsamplingFactor != 0) {
				final float scaleFactor = (float) (1.0f / maxSubsamplingFactor);
				seq = SequenceUtil.scale(seq, Math.round(seq.getWidth() * scaleFactor),
				    Math.round(seq.getHeight() * scaleFactor));
				width = seq.getWidth();
				height = seq.getHeight();
			}
			image = new double[seq.getWidth() * seq.getHeight()];
			seqData = Array2DUtil.arrayToDoubleArray(seq.getDataXYC(0, 0), seq.isSignedDataType());
			for (int xy = 0; xy < originalImage.length; xy++) {
				double val = 0;
				for (int c = 0; c < seqData.length; c++) {
					val += seqData[c][xy];
				}
				image[xy] = val / 3d;
			}

			// update sub-sampled output version information if necessary
			if (width <= subWidth) {
				subWidth = width;
				subHeight = height;
				subImage = image;
			}
		}
		coefficient = getBasicFromCardinal2D();

		if (coefficient != null) {
			buildCoefficientPyramid();
		} else {
			buildEmptyCoefficientPyramid();
		}

		if (isTarget || isOutputSubsampled) {
			buildImagePyramid();
		}

	}

	/**
	 * Get basic from cardinal: convert the 2D image from regular samples to
	 * standard B-spline coefficients.
	 * 
	 * @return array of standard B-spline coefficients
	 */
	private double[] getBasicFromCardinal2D() {
		if (image == null)
			return null;

		final double[] basic = new double[width * height];
		final double[] hLine = new double[width];
		final double[] vLine = new double[height];
		for (int y = 0; (y < height); y++) {
			extractRow(image, y, hLine);
			samplesToInterpolationCoefficient1D(hLine, 3, 0.0);
			putRow(basic, y, hLine);
		}
		for (int x = 0; (x < width); x++) {
			extractColumn(basic, width, x, vLine);
			samplesToInterpolationCoefficient1D(vLine, 3, 0.0);
			putColumn(basic, width, x, vLine);
		}
		return (basic);
	}

	/**
	 * Extract a row from the array
	 * 
	 * @param array
	 * @param y
	 *          Row position in the array
	 * @param row
	 *          Extracted row
	 */
	private void extractRow(double[] array, int y, double[] row) {
		y *= row.length;
		for (int x = 0; x < row.length; x++) {
			row[x] = array[x + y];
		}
	}

	/**
	 * Extract a column from the array.
	 * 
	 * @param array
	 * @param width
	 *          Width of the position of the column in the array
	 * @param x
	 *          Column position in the array
	 * @param col
	 *          Extracted column
	 */
	private void extractColumn(double[] array, int width, int x, double[] col) {
		for (int i = 0; (i < col.length); i++, x += width)
			col[i] = array[x];
	}

	private void samplesToInterpolationCoefficient1D(double[] c, int degree, double tolerance) {
		double[] z = new double[0];
		double lambda = 1.0;
		switch (degree) {
		case 3:
			z = new double[1];
			z[0] = Math.sqrt(3.0) - 2.0;
			break;
		case 7:
			z = new double[3];
			z[0] = -0.5352804307964381655424037816816460718339231523426924148812;
			z[1] = -0.122554615192326690515272264359357343605486549427295558490763;
			z[2] = -0.0091486948096082769285930216516478534156925639545994482648003;
			break;
		default:
		}
		// special case required by mirror boundaries
		if (c.length == 1) {
			return;
		}
		// compute the overall gain
		for (int k = 0; (k < z.length); k++) {
			lambda *= (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
		}
		// apply the gain
		for (int n = 0; (n < c.length); n++) {
			c[n] *= lambda;
		}
		// loop over all poles
		for (int k = 0; (k < z.length); k++) {
			// causal initialization
			c[0] = getInitialCausalCoefficientMirrorOffBounds(c, z[k], tolerance);
			// causal recursion
			for (int n = 1; (n < c.length); n++) {
				c[n] += z[k] * c[n - 1];
			}
			// anticausal initialization
			c[c.length - 1] = getInitialAntiCausalCoefficientMirrorOffBounds(c, z[k], tolerance);
			// anticausal recursion
			for (int n = c.length - 2; (0 <= n); n--) {
				c[n] = z[k] * (c[n + 1] - c[n]);
			}
		}
	}

	/**
	 * Get initial causal coefficients mirror off-bounds.
	 * 
	 * @param c
	 *          Coefficients
	 * @param z
	 * @param tolerance
	 * @return
	 */
	private double getInitialCausalCoefficientMirrorOffBounds(double[] c, double z, double tolerance) {
		double z1 = z;
		double zn = Math.pow(z, c.length);
		double sum = (1.0 + z) * (c[0] + zn * c[c.length - 1]);
		int horizon = c.length;
		if (0.0 < tolerance) {
			horizon = 2 + (int) (Math.log(tolerance) / Math.log(Math.abs(z)));
			horizon = (horizon < c.length) ? (horizon) : (c.length);
		}
		zn = zn * zn;
		for (int n = 1; (n < (horizon - 1)); n++) {
			z1 = z1 * z;
			zn = zn / z;
			sum = sum + (z1 + zn) * c[n];
		}
		return (sum / (1.0 - Math.pow(z, 2 * c.length)));
	}

	/**
	 * Get initial anti-causal coefficients mirror off-bounds.
	 * 
	 * @param c
	 *          Coefficients
	 * @param z
	 * @param tolerance
	 * @return bounds
	 */
	private double getInitialAntiCausalCoefficientMirrorOffBounds(double[] c, double z, double tolerance) {
		return (z * c[c.length - 1] / (z - 1.0));
	}

	/**
	 * Put a row in the array.
	 * 
	 * @param array
	 * @param y
	 *          Row position in the array
	 * @param row
	 *          Row to be put
	 */
	private void putRow(double[] array, int y, double[] row) {
		y *= row.length;
		for (int i = 0; (i < row.length); i++)
			array[y++] = (double) row[i];
	}

	/**
	 * Put a column in the array.
	 * 
	 * @param array
	 * @param width
	 *          Width of the position of the column in the array
	 * @param x
	 *          Column position in the array
	 * @param col
	 *          Column to be put
	 */
	private void putColumn(double[] array, int width, int x, double[] col) {
		for (int i = 0; (i < col.length); i++, x += width)
			array[x] = col[i];
	}

	/**
	 * Build the coefficients pyramid.
	 */
	private void buildCoefficientPyramid() {
		int fullWidth;
		int fullHeight;
		double[] fullDual = new double[width * height];
		int halfWidth = width;
		int halfHeight = height;
		basicToCardinal2D(coefficient, fullDual, width, height, 7);

		// We compute the coefficients pyramid
		for (int depth = 1; ((depth <= pyramidDepth) && (!super.isInterrupted())); depth++) {
			System.out.println("Building coefficients pyramid...");
			System.out.println((double) depth / pyramidDepth);
			fullWidth = halfWidth;
			fullHeight = halfHeight;
			halfWidth /= 2;
			halfHeight /= 2;

			// If the image is too small, we push the previous version of the
			// coefficients
			if (fullWidth <= BSplineModel.minImageSize || fullHeight <= BSplineModel.minImageSize) {
				if (this.isOutputSubsampled) {
					System.out.println("Coefficients pyramid " + fullWidth + "x" + fullHeight);
				}
				cPyramid.push(fullDual);
				cPyramid.push(new Integer(fullHeight));
				cPyramid.push(new Integer(fullWidth));
				halfWidth *= 2;
				halfHeight *= 2;
				continue;
			}

			// Otherwise, we reduce the coefficients by 2
			final double[] halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);
			final double[] halfCoefficient = getBasicFromCardinal2D(halfDual, halfWidth, halfHeight, 7);

			if (this.isOutputSubsampled) {
				System.out.println("Coefficients pyramid " + halfWidth + "x" + halfHeight);
			}
			cPyramid.push(halfCoefficient);
			cPyramid.push(new Integer(halfHeight));
			cPyramid.push(new Integer(halfWidth));

			fullDual = halfDual;

			// We store the coefficients of the corresponding subsampled
			// output if it exists.
			if (this.isOutputSubsampled && halfWidth == this.subWidth) {
				this.subCoeffs = halfCoefficient;
			}
		}
		smallestWidth = halfWidth;
		smallestHeight = halfHeight;
		currentDepth = pyramidDepth + 1;
	}

	/**
	 * Get half dual (2D).
	 * 
	 * @param fullDual
	 *          Full coefficients
	 * @param fullWidth
	 *          Full coefficients width
	 * @param fullHeight
	 *          Full coefficients height
	 * @return Half coefficients
	 */
	private double[] getHalfDual2D(double[] fullDual, int fullWidth, int fullHeight) {
		final int halfWidth = fullWidth / 2;
		final int halfHeight = fullHeight / 2;
		final double[] hLine = new double[fullWidth];
		final double[] hData = new double[halfWidth];
		final double[] vLine = new double[fullHeight];
		final double[] vData = new double[halfHeight];
		final double[] demiDual = new double[halfWidth * fullHeight];
		final double[] halfDual = new double[halfWidth * halfHeight];
		for (int y = 0; ((y < fullHeight) && (!super.isInterrupted())); y++) {
			extractRow(fullDual, y, hLine);
			reduceDual1D(hLine, hData);
			putRow(demiDual, y, hData);
		}
		for (int x = 0; ((x < halfWidth) && (!super.isInterrupted())); x++) {
			extractColumn(demiDual, halfWidth, x, vLine);
			reduceDual1D(vLine, vData);
			putColumn(halfDual, halfWidth, x, vData);
		}
		return (halfDual);
	}

	/**
	 * Reduces dual (1D).
	 * 
	 * @param c
	 * @param s
	 */
	private void reduceDual1D(double[] c, double[] s) {
		final double h[] = { 6.0 / 16.0, 4.0 / 16.0, 1.0 / 16.0 };
		if (2 <= s.length) {
			s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]);
			for (int i = 2, j = 1; (j < (s.length - 1)); i += 2, j++) {
				s[j] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1]) + h[2] * (c[i - 2] + c[i + 2]);
			}
			if (c.length == (2 * s.length)) {
				s[s.length - 1] = h[0] * c[c.length - 2] + h[1] * (c[c.length - 3] + c[c.length - 1])
				    + h[2] * (c[c.length - 4] + c[c.length - 1]);
			} else {
				s[s.length - 1] = h[0] * c[c.length - 3] + h[1] * (c[c.length - 4] + c[c.length - 2])
				    + h[2] * (c[c.length - 5] + c[c.length - 1]);
			}
		} else {
			switch (c.length) {
			case 3:
				s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]);
				break;
			case 2:
				s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + 2.0 * h[2] * c[1];
				break;
			default:
			}
		}
	}

	/**
	 * Get basic from cardinal (2D): convert a 2D signal from regular samples to
	 * standard B-spline coefficients.
	 * 
	 * @param cardinal
	 *          Sampled 2D signal
	 * @param width
	 *          Signal width
	 * @param height
	 *          Signal height
	 * @param degree
	 *          B-spline degree
	 * @return Array of standard B-spline coefficients
	 */
	private double[] getBasicFromCardinal2D(double[] cardinal, int width, int height, int degree) {
		final double[] basic = new double[width * height];
		final double[] hLine = new double[width];
		final double[] vLine = new double[height];
		for (int y = 0; ((y < height) && (!super.isInterrupted())); y++) {
			extractRow(cardinal, y, hLine);
			samplesToInterpolationCoefficient1D(hLine, degree, 0.0);
			putRow(basic, y, hLine);
		}
		for (int x = 0; ((x < width) && (!super.isInterrupted())); x++) {
			extractColumn(basic, width, x, vLine);
			samplesToInterpolationCoefficient1D(vLine, degree, 0.0);
			putColumn(basic, width, x, vLine);
		}
		return (basic);
	}

	/**
	 * Pass from basic to cardinal.
	 * 
	 * @param basic
	 *          Basic (standard B-splines) 2D array
	 * @param cardinal
	 *          Cardinal (sampled signal) 2D array
	 * @param width
	 * @param height
	 * @param degree
	 */
	private void basicToCardinal2D(double[] basic, double[] cardinal, int width, int height, int degree) {
		final double[] hLine = new double[width];
		final double[] vLine = new double[height];
		final double[] hData = new double[width];
		final double[] vData = new double[height];
		double[] h = null;
		switch (degree) {
		case 3:
			h = new double[2];
			h[0] = 2.0 / 3.0;
			h[1] = 1.0 / 6.0;
			break;
		case 7:
			h = new double[4];
			h[0] = 151.0 / 315.0;
			h[1] = 397.0 / 1680.0;
			h[2] = 1.0 / 42.0;
			h[3] = 1.0 / 5040.0;
			break;
		default:
			h = new double[1];
			h[0] = 1.0;
		}
		for (int y = 0; ((y < height) && (!super.isInterrupted())); y++) {
			extractRow(basic, y, hLine);
			symmetricFirMirrorOffBounds1D(h, hLine, hData);
			putRow(cardinal, y, hData);
		}
		for (int x = 0; ((x < width) && (!super.isInterrupted())); x++) {
			extractColumn(cardinal, width, x, vLine);
			symmetricFirMirrorOffBounds1D(h, vLine, vData);
			putColumn(cardinal, width, x, vData);
		}
	}

	/**
	 * Symmetric FIR filter with mirror off bounds (1D) conditions.
	 * 
	 * @param h
	 * @param c
	 * @param s
	 */
	private void symmetricFirMirrorOffBounds1D(double[] h, double[] c, double[] s) {
		switch (h.length) {
		case 2:
			if (2 <= c.length) {
				s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]);
				for (int i = 1; (i < (s.length - 1)); i++) {
					s[i] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1]);
				}
				s[s.length - 1] = h[0] * c[c.length - 1] + h[1] * (c[c.length - 2] + c[c.length - 1]);
			} else {
				s[0] = (h[0] + 2.0 * h[1]) * c[0];
			}
			break;
		case 4:
			if (6 <= c.length) {
				s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]) + h[3] * (c[2] + c[3]);
				s[1] = h[0] * c[1] + h[1] * (c[0] + c[2]) + h[2] * (c[0] + c[3]) + h[3] * (c[1] + c[4]);
				s[2] = h[0] * c[2] + h[1] * (c[1] + c[3]) + h[2] * (c[0] + c[4]) + h[3] * (c[0] + c[5]);
				for (int i = 3; (i < (s.length - 3)); i++) {
					s[i] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1]) + h[2] * (c[i - 2] + c[i + 2])
					    + h[3] * (c[i - 3] + c[i + 3]);
				}
				s[s.length - 3] = h[0] * c[c.length - 3] + h[1] * (c[c.length - 4] + c[c.length - 2])
				    + h[2] * (c[c.length - 5] + c[c.length - 1]) + h[3] * (c[c.length - 6] + c[c.length - 1]);
				s[s.length - 2] = h[0] * c[c.length - 2] + h[1] * (c[c.length - 3] + c[c.length - 1])
				    + h[2] * (c[c.length - 4] + c[c.length - 1]) + h[3] * (c[c.length - 5] + c[c.length - 2]);
				s[s.length - 1] = h[0] * c[c.length - 1] + h[1] * (c[c.length - 2] + c[c.length - 1])
				    + h[2] * (c[c.length - 3] + c[c.length - 2]) + h[3] * (c[c.length - 4] + c[c.length - 3]);
			} else {
				switch (c.length) {
				case 5:
					s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]) + h[3] * (c[2] + c[3]);
					s[1] = h[0] * c[1] + h[1] * (c[0] + c[2]) + h[2] * (c[0] + c[3]) + h[3] * (c[1] + c[4]);
					s[2] = h[0] * c[2] + h[1] * (c[1] + c[3]) + (h[2] + h[3]) * (c[0] + c[4]);
					s[3] = h[0] * c[3] + h[1] * (c[2] + c[4]) + h[2] * (c[1] + c[4]) + h[3] * (c[0] + c[3]);
					s[4] = h[0] * c[4] + h[1] * (c[3] + c[4]) + h[2] * (c[2] + c[3]) + h[3] * (c[1] + c[2]);
					break;
				case 4:
					s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]) + h[3] * (c[2] + c[3]);
					s[1] = h[0] * c[1] + h[1] * (c[0] + c[2]) + h[2] * (c[0] + c[3]) + h[3] * (c[1] + c[3]);
					s[2] = h[0] * c[2] + h[1] * (c[1] + c[3]) + h[2] * (c[0] + c[3]) + h[3] * (c[0] + c[2]);
					s[3] = h[0] * c[3] + h[1] * (c[2] + c[3]) + h[2] * (c[1] + c[2]) + h[3] * (c[0] + c[1]);
					break;
				case 3:
					s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]) + 2.0 * h[3] * c[2];
					s[1] = h[0] * c[1] + (h[1] + h[2]) * (c[0] + c[2]) + 2.0 * h[3] * c[1];
					s[2] = h[0] * c[2] + h[1] * (c[1] + c[2]) + h[2] * (c[0] + c[1]) + 2.0 * h[3] * c[0];
					break;
				case 2:
					s[0] = (h[0] + h[1] + h[3]) * c[0] + (h[1] + 2.0 * h[2] + h[3]) * c[1];
					s[1] = (h[0] + h[1] + h[3]) * c[1] + (h[1] + 2.0 * h[2] + h[3]) * c[0];
					break;
				case 1:
					s[0] = (h[0] + 2.0 * (h[1] + h[2] + h[3])) * c[0];
					break;
				default:
				}
			}
			break;
		default:
		}
	}

	/**
	 * Build an empty coefficient pyramid (for only-landmark registration).
	 */
	private void buildEmptyCoefficientPyramid() {
		int fullWidth;
		int fullHeight;

		int halfWidth = width;
		int halfHeight = height;
		final double[] fullDual = new double[] {};
		final double[] halfCoefficient = new double[] {};

		// We compute the coefficients pyramid
		for (int depth = 1; ((depth <= pyramidDepth) && (!super.isInterrupted())); depth++) {
			System.out.println("Building coefficients pyramid...");
			// IJ.showProgress((double) depth / pyramidDepth );
			fullWidth = halfWidth;
			fullHeight = halfHeight;
			halfWidth /= 2;
			halfHeight /= 2;

			// If the image is too small, we push the previous version of the
			// coefficients
			if (fullWidth <= BSplineModel.minImageSize || fullHeight <= BSplineModel.minImageSize) {
				if (this.isOutputSubsampled) {
					System.out.println("Coefficients pyramid " + fullWidth + "x" + fullHeight);
				}
				cPyramid.push(fullDual);
				cPyramid.push(new Integer(fullHeight));
				cPyramid.push(new Integer(fullWidth));
				halfWidth *= 2;
				halfHeight *= 2;
				continue;
			}

			// Otherwise, we reduce the coefficients by 2
			if (this.isOutputSubsampled) {
				System.out.println("Coefficients pyramid " + halfWidth + "x" + halfHeight);
			}
			cPyramid.push(halfCoefficient);
			cPyramid.push(new Integer(halfHeight));
			cPyramid.push(new Integer(halfWidth));

			// We store the coefficients of the corresponding subsampled
			// output if it exists.
			if (this.isOutputSubsampled && halfWidth == this.subWidth) {
				this.subCoeffs = halfCoefficient;
			}
		}
		smallestWidth = halfWidth;
		smallestHeight = halfHeight;
		currentDepth = pyramidDepth + 1;
	}

	/**
	 * Build the image pyramid
	 */
	private void buildImagePyramid() {
		int fullWidth;
		int fullHeight;
		double[] fullDual = new double[width * height];
		int halfWidth = width;
		int halfHeight = height;
		cardinalToDual2D(image, fullDual, width, height, 3);

		for (int depth = 1; depth <= pyramidDepth && !super.isInterrupted(); depth++) {
			System.out.println("Building image pyramid...");
			// IJ.showProgress((double) depth / pyramidDepth);

			fullWidth = halfWidth;
			fullHeight = halfHeight;

			halfWidth /= 2;
			halfHeight /= 2;

			if (fullWidth <= BSplineModel.minImageSize || fullHeight <= BSplineModel.minImageSize) {
				if (this.isOutputSubsampled) {
					System.out.println(" Image pyramid " + fullWidth + "x" + fullHeight);
				}
				imgPyramid.push(fullDual);
				imgPyramid.push(new Integer(fullHeight));
				imgPyramid.push(new Integer(fullWidth));
				halfWidth *= 2;
				halfHeight *= 2;
				continue;
			}

			final double[] halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);
			final double[] halfImage = new double[halfWidth * halfHeight];
			dualToCardinal2D(halfDual, halfImage, halfWidth, halfHeight, 3);

			if (this.isOutputSubsampled) {
				System.out.println(" Image pyramid " + halfWidth + "x" + halfHeight);
			}
			imgPyramid.push(halfImage);
			imgPyramid.push(new Integer(halfHeight));
			imgPyramid.push(new Integer(halfWidth));

			fullDual = halfDual;

			if (this.isOutputSubsampled && halfWidth == this.subWidth) {
				this.subImage = halfDual;
				// IJ.log("sub image set");
			}
		}

		// If the output sub-image has not been set yet, we keep reducing the image
		while (halfWidth > this.subWidth) {
			fullWidth = halfWidth;
			fullHeight = halfHeight;

			halfWidth /= 2;
			halfHeight /= 2;

			final double[] halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);
			final double[] halfImage = new double[halfWidth * halfHeight];
			dualToCardinal2D(halfDual, halfImage, halfWidth, halfHeight, 3);

			fullDual = halfDual;
			// We store the image version that matches the sub-sampled output (if
			// necessary)
			if (this.isOutputSubsampled && halfWidth == this.subWidth) {
				this.subImage = halfDual;
				// IJ.log("sub image set");
			}
		}
	}

	/**
	 * Passes from cardinal to dual (2D).
	 * 
	 * @param image2
	 * @param fullDual
	 * @param width
	 * @param height
	 * @param degree
	 */
	private void cardinalToDual2D(double[] cardinal, double[] dual, int width, int height, int degree) {
		basicToCardinal2D(getBasicFromCardinal2D(cardinal, width, height, degree), dual, width, height, 2 * degree + 1);
	}

	/**
	 * Pass from dual to cardinal (2D).
	 * 
	 * @param dual
	 * @param cardinal
	 * @param width
	 * @param height
	 * @param degree
	 */
	private void dualToCardinal2D(double[] dual, double[] cardinal, int width, int height, int degree) {
		basicToCardinal2D(getBasicFromCardinal2D(dual, width, height, 2 * degree + 1), cardinal, width, height, degree);
	}

	/**
	 * Sets the depth up to which the pyramids should be computed.
	 * 
	 * @param depth
	 */
	public void setPyramidDepth(int depth) {
		int proposedPyramidDepth = depth;

		int currentWidth = width;
		int currentHeight = height;
		int scale = 0;
		while (currentWidth >= minImageSize && currentHeight >= minImageSize) {
			currentWidth /= 2;
			currentHeight /= 2;
			scale++;
		}
		scale--;

		if (proposedPyramidDepth > scale) {
			proposedPyramidDepth = scale;
		}

		pyramidDepth = proposedPyramidDepth;
	}

	/**
	 * Pop one element from the coefficients and image pyramids.
	 */
	public void popFromPyramid() {
		// Pop coefficients
		if (cPyramid.isEmpty()) {
			currentWidth = width;
			currentHeight = height;
			currentCoefficient = coefficient;
		} else {
			currentWidth = ((Integer) cPyramid.pop()).intValue();
			currentHeight = ((Integer) cPyramid.pop()).intValue();
			currentCoefficient = (double[]) cPyramid.pop();
		}

		if (currentDepth > 0) {
			currentDepth--;
		}

		// Pop image
		if (isTarget && !imgPyramid.isEmpty()) {
			if (currentWidth != ((Integer) imgPyramid.pop()).intValue()) {
				System.out.println("I cannot understand");
			}
			if (currentHeight != ((Integer) imgPyramid.pop()).intValue()) {
				System.out.println("I cannot understand");
			}
			currentImage = (double[]) imgPyramid.pop();
		} else {
			currentImage = image;
		}
	}

	/**
	 * Prepare precomputations for a given image size. It calls
	 * prepareForInterpolation with ORIGINAL flag.
	 *
	 * @param Ydim
	 *          y- image dimension
	 * @param Xdim
	 *          x- image dimension
	 * @param intervals
	 *          intervals in the deformation
	 */
	public void precomputedPrepareForInterpolation(int Ydim, int Xdim, int intervals) {
		// Ask for memory
		precXIndex = new int[Xdim][4];
		precYIndex = new int[Ydim][4];
		precXWeight = new double[Xdim][4];
		precYWeight = new double[Ydim][4];
		precDxWeight = new double[Xdim][4];
		precDyWeight = new double[Ydim][4];
		precD2xWeight = new double[Xdim][4];
		precD2yWeight = new double[Ydim][4];

		boolean ORIGINAL = false;
		// Fill the precomputed weights and indexes for the Y axis
		for (int v = 0; v < Ydim; v++) {
			// Express the current point in Spline units
			final double tv = (double) (v * intervals) / (double) (Ydim - 1) + 1.0F;
			final double tu = 1.0F;

			// Compute all weights and indexes
			prepareForInterpolation(tu, tv, ORIGINAL);

			// Copy all values
			for (int k = 0; k < 4; k++) {
				precYIndex[v][k] = yIndex[k];
				precYWeight[v][k] = yWeight[k];
				precDyWeight[v][k] = dyWeight[k];
				precD2yWeight[v][k] = d2yWeight[k];
			}
		}

		// Fill the precomputed weights and indexes for the X axis
		for (int u = 0; u < Xdim; u++) {
			// Express the current point in Spline units
			final double tv = 1.0F;
			final double tu = (double) (u * intervals) / (double) (Xdim - 1) + 1.0F;

			// Compute all weights and indexes
			prepareForInterpolation(tu, tv, ORIGINAL);

			// Copy all values
			for (int k = 0; k < 4; k++) {
				precXIndex[u][k] = xIndex[k];
				precXWeight[u][k] = xWeight[k];
				precDxWeight[u][k] = dxWeight[k];
				precD2xWeight[u][k] = d2xWeight[k];
			}
		}
	}

	/**
	 * fromCurrent=true --> The interpolation is prepared to be done from the
	 * current image in the pyramid. fromCurrent=false --> The interpolation is
	 * prepared to be done from the original image.
	 *
	 * @param x
	 *          x- point coordinate
	 * @param y
	 *          y- point coordinate
	 * @param fromCurrent
	 *          flag to determine the image to do the interpolation from
	 */
	public void prepareForInterpolation(double x, double y, boolean fromCurrent) {
		// Remind this point for interpolation
		// this.x = x;
		// this.y = y;
		this.fromCurrent = fromCurrent;

		if (fromCurrent) {
			widthToUse = currentWidth;
			heightToUse = currentHeight;
		} else {
			widthToUse = width;
			heightToUse = height;
		}

		int ix = (int) x;
		int iy = (int) y;

		int twiceWidthToUse = 2 * widthToUse;
		int twiceHeightToUse = 2 * heightToUse;

		// Set X indexes
		// p is the index of the rightmost influencing spline
		int p = (0.0 <= x) ? (ix + 2) : (ix + 1);
		for (int k = 0; k < 4; p--, k++) {
			if (coefficientsAreMirrored) {
				int q = (p < 0) ? (-1 - p) : (p);
				if (twiceWidthToUse <= q)
					q -= twiceWidthToUse * (q / twiceWidthToUse);
				xIndex[k] = (widthToUse <= q) ? (twiceWidthToUse - 1 - q) : (q);
			} else
				xIndex[k] = (p < 0 || p >= widthToUse) ? (-1) : (p);
		}

		// Set Y indexes
		p = (0.0 <= y) ? (iy + 2) : (iy + 1);
		for (int k = 0; k < 4; p--, k++) {
			if (coefficientsAreMirrored) {
				int q = (p < 0) ? (-1 - p) : (p);
				if (twiceHeightToUse <= q)
					q -= twiceHeightToUse * (q / twiceHeightToUse);
				yIndex[k] = (heightToUse <= q) ? (twiceHeightToUse - 1 - q) : (q);
			} else
				yIndex[k] = (p < 0 || p >= heightToUse) ? (-1) : (p);
		}

		// Compute how much the sample depart from an integer position
		double ex = x - ((0.0 <= x) ? (ix) : (ix - 1));
		double ey = y - ((0.0 <= y) ? (iy) : (iy - 1));

		// Set X weights for the image and derivative interpolation
		double s = 1.0F - ex;
		dxWeight[0] = 0.5F * ex * ex;
		xWeight[0] = ex * dxWeight[0] / 3.0F; // Bspline03(x-ix-2)
		dxWeight[3] = -0.5F * s * s;
		xWeight[3] = s * dxWeight[3] / -3.0F; // Bspline03(x-ix+1)
		dxWeight[1] = 1.0F - 2.0F * dxWeight[0] + dxWeight[3];
		// xWeight[1] = 2.0F / 3.0F + (1.0F + ex) * dxWeight[3]; //
		// Bspline03(x-ix-1);
		xWeight[1] = MathTools.Bspline03(x - ix - 1);
		dxWeight[2] = 1.5F * ex * (ex - 4.0F / 3.0F);
		xWeight[2] = 2.0F / 3.0F - (2.0F - ex) * dxWeight[0]; // Bspline03(x-ix)

		d2xWeight[0] = ex;
		d2xWeight[1] = s - 2 * ex;
		d2xWeight[2] = ex - 2 * s;
		d2xWeight[3] = s;

		// Set Y weights for the image and derivative interpolation
		double t = 1.0F - ey;
		dyWeight[0] = 0.5F * ey * ey;
		yWeight[0] = ey * dyWeight[0] / 3.0F;
		dyWeight[3] = -0.5F * t * t;
		yWeight[3] = t * dyWeight[3] / -3.0F;
		dyWeight[1] = 1.0F - 2.0F * dyWeight[0] + dyWeight[3];
		yWeight[1] = 2.0F / 3.0F + (1.0F + ey) * dyWeight[3];
		dyWeight[2] = 1.5F * ey * (ey - 4.0F / 3.0F);
		yWeight[2] = 2.0F / 3.0F - (2.0F - ey) * dyWeight[0];

		d2yWeight[0] = ey;
		d2yWeight[1] = t - 2 * ey;
		d2yWeight[2] = ey - 2 * t;
		d2yWeight[3] = t;
	}

	/**
	 * Interpolate the image at a given point.
	 *
	 * @return image interpolation
	 */
	public double interpolateI() {
		// Only SplineDegree=3 is implemented
		double ival = 0.0F;
		for (int j = 0; j < 4; j++) {
			double s = 0.0F;
			int iy = yIndex[j];
			if (iy != -1) {
				int p = iy * widthToUse;
				for (int i = 0; i < 4; i++) {
					int ix = xIndex[i];
					if (ix != -1)
						if (fromCurrent)
							s += xWeight[i] * currentCoefficient[p + ix];
						else
							s += xWeight[i] * coefficient[p + ix];
				}
				ival += yWeight[j] * s;
			}
		}
		return ival;
	}

	/**
	 * Prepare for interpolation and interpolate
	 * 
	 * fromSub = true --> The interpolation is done from the subsampled version of
	 * the image else:
	 * 
	 * fromCurrent=true --> The interpolation is done from the current image in
	 * the pyramid. fromCurrent=false --> The interpolation is done from the
	 * original image.
	 *
	 * @param x
	 *          x- point coordinate
	 * @param y
	 *          y- point coordinate
	 * @param fromSub
	 *          flat to determine to do the interpolation from the subsampled
	 *          version of the image
	 * @param fromCurrent
	 *          flag to determine the image to do the interpolation from
	 *          interpolated value
	 */
	public double prepareForInterpolationAndInterpolateI(double x, double y, boolean fromSub, boolean fromCurrent) {
		int widthToUse = 0;
		int heightToUse = 0;
		final int[] xIndex = new int[4];
		final int[] yIndex = new int[4];
		final double[] xWeight = new double[4];
		final double[] dxWeight = new double[4];
		// double[] d2xWeight = new double[4];
		final double[] yWeight = new double[4];
		final double[] dyWeight = new double[4];
		// double[] d2yWeight = new double[4];

		if (fromSub && this.subCoeffs != null) {
			widthToUse = this.subWidth;
			heightToUse = this.subHeight;
		} else if (fromCurrent) {
			widthToUse = currentWidth;
			heightToUse = currentHeight;
		} else {
			widthToUse = width;
			heightToUse = height;
		}

		int ix = (int) x;
		int iy = (int) y;

		int twiceWidthToUse = 2 * widthToUse;
		int twiceHeightToUse = 2 * heightToUse;

		// Set X indexes
		// p is the index of the rightmost influencing spline
		int p = (0.0 <= x) ? (ix + 2) : (ix + 1);
		for (int k = 0; k < 4; p--, k++) {
			if (coefficientsAreMirrored) {
				int q = (p < 0) ? (-1 - p) : (p);
				if (twiceWidthToUse <= q)
					q -= twiceWidthToUse * (q / twiceWidthToUse);
				xIndex[k] = (widthToUse <= q) ? (twiceWidthToUse - 1 - q) : (q);
			} else
				xIndex[k] = (p < 0 || p >= widthToUse) ? (-1) : (p);
		}

		// Set Y indexes
		p = (0.0 <= y) ? (iy + 2) : (iy + 1);
		for (int k = 0; k < 4; p--, k++) {
			if (coefficientsAreMirrored) {
				int q = (p < 0) ? (-1 - p) : (p);
				if (twiceHeightToUse <= q)
					q -= twiceHeightToUse * (q / twiceHeightToUse);
				yIndex[k] = (heightToUse <= q) ? (twiceHeightToUse - 1 - q) : (q);
			} else
				yIndex[k] = (p < 0 || p >= heightToUse) ? (-1) : (p);
		}

		// Compute how much the sample depart from an integer position
		double ex = x - ((0.0 <= x) ? (ix) : (ix - 1));
		double ey = y - ((0.0 <= y) ? (iy) : (iy - 1));

		// Set X weights for the image and derivative interpolation
		double s = 1.0F - ex;
		dxWeight[0] = 0.5F * ex * ex;
		xWeight[0] = ex * dxWeight[0] / 3.0F; // Bspline03(x-ix-2)
		dxWeight[3] = -0.5F * s * s;
		xWeight[3] = s * dxWeight[3] / -3.0F; // Bspline03(x-ix+1)
		dxWeight[1] = 1.0F - 2.0F * dxWeight[0] + dxWeight[3];
		// xWeight[1] = 2.0F / 3.0F + (1.0F + ex) * dxWeight[3]; //
		// Bspline03(x-ix-1);
		xWeight[1] = MathTools.Bspline03(x - ix - 1);
		dxWeight[2] = 1.5F * ex * (ex - 4.0F / 3.0F);
		xWeight[2] = 2.0F / 3.0F - (2.0F - ex) * dxWeight[0]; // Bspline03(x-ix)

		// d2xWeight[0] = ex;
		// d2xWeight[1] = s-2*ex;
		// d2xWeight[2] = ex-2*s;
		// d2xWeight[3] = s;

		// Set Y weights for the image and derivative interpolation
		double t = 1.0F - ey;
		dyWeight[0] = 0.5F * ey * ey;
		yWeight[0] = ey * dyWeight[0] / 3.0F;
		dyWeight[3] = -0.5F * t * t;
		yWeight[3] = t * dyWeight[3] / -3.0F;
		dyWeight[1] = 1.0F - 2.0F * dyWeight[0] + dyWeight[3];
		yWeight[1] = 2.0F / 3.0F + (1.0F + ey) * dyWeight[3];
		dyWeight[2] = 1.5F * ey * (ey - 4.0F / 3.0F);
		yWeight[2] = 2.0F / 3.0F - (2.0F - ey) * dyWeight[0];

		// d2yWeight[0] = ey;
		// d2yWeight[1] = t-2*ey;
		// d2yWeight[2] = ey-2*t;
		// d2yWeight[3] = t;

		// Only SplineDegree=3 is implemented
		double ival = 0.0F;
		for (int j = 0; j < 4; j++) {
			s = 0.0F;
			iy = yIndex[j];
			if (iy != -1) {
				p = iy * widthToUse;
				for (int i = 0; i < 4; i++) {
					ix = xIndex[i];
					if (ix != -1) {
						if (fromSub && this.subCoeffs != null)
							s += xWeight[i] * this.subCoeffs[p + ix];
						else if (fromCurrent)
							s += xWeight[i] * currentCoefficient[p + ix];
						else
							s += xWeight[i] * coefficient[p + ix];
					}
				}
				ival += yWeight[j] * s;
			}
		}
		return ival;
	}

	/**
	 * Interpolate the X and Y derivatives of the image at a given point.
	 *
	 * @param D
	 *          output, X and Y derivatives of the image
	 * @param u
	 *          x-point coordinate
	 * @param v
	 *          y-point coordinate
	 */
	public void precomputedInterpolateD(double[] D, int u, int v) {
		// Only SplineDegree=3 is implemented
		D[0] = D[1] = 0.0F;
		for (int j = 0; j < 4; j++) {
			double sx = 0.0F, sy = 0.0F;
			int iy = precYIndex[v][j];
			if (iy != -1) {
				int p = iy * widthToUse;
				for (int i = 0; i < 4; i++) {
					int ix = precXIndex[u][i];
					if (ix != -1) {
						double c;
						if (fromCurrent)
							c = currentCoefficient[p + ix];
						else
							c = coefficient[p + ix];
						sx += precDxWeight[u][i] * c;
						sy += precXWeight[u][i] * c;
					}
				}
				D[0] += precYWeight[v][j] * sx;
				D[1] += precDyWeight[v][j] * sy;
			}
		}
	}

	/**
	 * Interpolate the XY, XX and YY derivatives of the image at a given point.
	 *
	 * @param D2
	 *          output, XY, XX and YY derivatives of the image
	 * @param u
	 *          x-point coordinate
	 * @param v
	 *          y-point coordinate
	 */
	public void precomputedInterpolateD2(double[] D2, int u, int v) {
		// Only SplineDegree=3 is implemented
		D2[0] = D2[1] = D2[2] = 0.0F;
		for (int j = 0; j < 4; j++) {
			double sxy = 0.0F, sxx = 0.0F, syy = 0.0F;
			int iy = precYIndex[v][j];
			if (iy != -1) {
				int p = iy * widthToUse;
				for (int i = 0; i < 4; i++) {
					int ix = precXIndex[u][i];
					if (ix != -1) {
						double c;
						if (fromCurrent)
							c = currentCoefficient[p + ix];
						else
							c = coefficient[p + ix];
						sxy += precDxWeight[u][i] * c;
						sxx += precD2xWeight[u][i] * c;
						syy += precXWeight[u][i] * c;
					}
				}
				D2[0] += precDyWeight[v][j] * sxy;
				D2[1] += precYWeight[v][j] * sxx;
				D2[2] += precD2yWeight[v][j] * syy;
			}
		}
	}

	/**
	 * Interpolate the image (or deformation) at a given point using the
	 * precomputed weights.
	 *
	 * @param u
	 *          x- point coordinate
	 * @param v
	 *          y- point coordinate
	 */
	public double precomputedInterpolateI(int u, int v) {
		// Only SplineDegree=3 is implemented
		double ival = 0.0F;
		for (int j = 0; j < 4; j++) {
			double s = 0.0F;
			int iy = precYIndex[v][j];
			if (iy != -1) {
				int p = iy * widthToUse;
				for (int i = 0; i < 4; i++) {
					int ix = precXIndex[u][i];
					if (ix != -1)
						if (fromCurrent)
							s += precXWeight[u][i] * currentCoefficient[p + ix];
						else
							s += precXWeight[u][i] * coefficient[p + ix];
				}
				ival += precYWeight[v][j] * s;
			}
		}
		return ival;
	}

	// ------------------------------------------------------------------
	/**
	 * There are two types of interpolation routines. Those that use precomputed
	 * weights and those that don't. An example of use of the ones without
	 * precomputation is the following:
	 * 
	 * <pre>
	 * 
	 * // Set of B-spline coefficients
	 * 
	 * double [][]c;
	 *
	 * // Set these coefficients to an interpolator
	 * 
	 * BSplineModel sw = new BSplineModel(c);
	 *
	 * // Compute the transformation mapping
	 * 
	 * for (int v = 0; v < ImageHeight; v++) { final double tv = (double)(v *
	 * intervals) / (double)(ImageHeight - 1) + 1.0F; for (int u = 0;
	 * u<ImageeWidth; u++) { final double tu = (double)(u * intervals) /
	 * (double)(ImageWidth - 1) + 1.0F; sw.prepareForInterpolation(tu, tv,
	 * ORIGINAL); interpolated_val[v][u] = sw.interpolateI(); } }
	 * 
	 * <pre>
	 */
	// ------------------------------------------------------------------

	/**
	 * Interpolate the X and Y derivatives of the image at a given point.
	 *
	 * @param D
	 *          output, interpolation the X and Y derivatives of the image
	 */
	public void interpolateD(double[] D) {
		// Only SplineDegree=3 is implemented
		D[0] = D[1] = 0.0F;
		for (int j = 0; j < 4; j++) {
			double sx = 0.0F, sy = 0.0F;
			int iy = yIndex[j];
			if (iy != -1) {
				int p = iy * widthToUse;
				for (int i = 0; i < 4; i++) {
					int ix = xIndex[i];
					if (ix != -1) {
						double c;
						if (fromCurrent)
							c = currentCoefficient[p + ix];
						else
							c = coefficient[p + ix];
						sx += dxWeight[i] * c;
						sy += xWeight[i] * c;
					}
				}
				D[0] += yWeight[j] * sx;
				D[1] += dyWeight[j] * sy;
			}
		}
	}

	/**
	 * @return The full-size image width.
	 */
	public int getWidth() {
		return width;
	}

	/**
	 * @return The full-size image height.
	 */
	public int getHeight() {
		return height;
	}

	/**
	 * @return The thread associated to this class.
	 */
	public Thread getThread() {
		return this;
	}

	/**
	 * Get subsampled output flag
	 * 
	 * @return true if the output needs to be subsampled
	 */
	public boolean isSubOutput() {
		return isOutputSubsampled;
	}

	/**
	 * @return Subsampled output image.
	 */
	public double[] getSubImage() {
		return subImage;
	}

	/**
	 * Get image (at the maximum resolution size determined by the scaling).
	 * 
	 * @return The less scaled image.
	 */
	public double[] getImage() {
		return image;
	}

	/**
	 * Get subsampled output width
	 * 
	 * @return subsampled output width
	 */
	public int getSubWidth() {
		return subWidth;
	}

	/**
	 * Get subsampled output height
	 * 
	 * @return subsampled output height
	 */
	public int getSubHeight() {
		return subHeight;
	}

	/**
	 * Get current image.
	 *
	 * @return the current image of the image/coefficients
	 */
	public double[] getCurrentImage() {
		return currentImage;
	}

	/**
	 * Get current height.
	 * 
	 * @return the current height of the image/coefficients
	 */
	public int getCurrentHeight() {
		return currentHeight;
	}

	/**
	 * Get current width.
	 * 
	 * @return the current width of the image/coefficients
	 */
	public int getCurrentWidth() {
		return currentWidth;
	}

	/**
	 * Get factor height.
	 * 
	 * @return the relationship between the current size of the image and the
	 *         original size
	 */
	public double getFactorHeight() {
		return (double) currentHeight / height;
	}

	/**
	 * Get factor width.
	 *
	 * @return the relationship between the current size of the image and the
	 *         original size.
	 */
	public double getFactorWidth() {
		return (double) currentWidth / width;
	}

	/**
	 * Get current depth.
	 * 
	 * @return the current depth of the image/coefficients
	 */
	public int getCurrentDepth() {
		return currentDepth;
	}

	/**
	 * Get original image.
	 * @return The original full-size image.
	 */
	public double[] getOriginalImage() {
		return originalImage;
	}
	
	/**
	 * Get original image width.
	 *
	 * @return the original full-size image width.
	 */
	public int getOriginalImageWidth() {
		return originalWidth;
	}

	/**
	 * Get original image height.
	 *
	 * @return the original full-size image height.
	 */
	public int getOriginalImageHeight() {
		return originalHeight;
	}

	/**
	 * Get sub-sampling factor at highest resolution
	 * 
	 * @return maxImageSubsamplingFactor
	 */
	public double getSubsamplingFactor() {
		return this.maxSubsamplingFactor;
	}

	/**
	 * Get width of precomputed vectors.
	 *
	 * @return the width of the precomputed vectors
	 */
	public int precomputedGetWidth() {
		return precYWeight.length;
	}

	/**
	 * Get precomputed weight of coefficient l,m
	 * 
	 * @param l
	 * @param m
	 * @param u
	 * @param v
	 * @return the weight of the coefficient l,m (prec_yWeight, prec_xWeight) in
	 *         the image interpolation
	 */
	public double precomputedGetWeightI(int l, int m, int u, int v) {
		return precYWeight[v][l] * precXWeight[u][m];
	}

	/**
	 * Get precomputed weight dx.
	 *
	 * @param l
	 * @param m
	 * @param u
	 * @param v
	 * @return the weight of the coefficient l,m (yIndex, xIndex) in the image
	 *         interpolation
	 */
	public double precomputedGetWeightDx(int l, int m, int u, int v) {
		return precYWeight[v][l] * precDxWeight[u][m];
	}

	/**
	 * Get precomputed weight dy
	 * 
	 * @param l
	 * @param m
	 * @param u
	 * @param v
	 * @return the weight of the coefficient l,m (prec_dyWeight, prec_xWeight) in
	 *         the image interpolation
	 */
	public double precomputedGetWeightDy(int l, int m, int u, int v) {
		return precDyWeight[v][l] * precXWeight[u][m];
	}

	/**
	 * Set spline coefficients. Copy coefficients to the model array.
	 *
	 * @param c
	 *          Set of B-spline coefficients
	 * @param Ydim
	 *          Y-dimension of the set of coefficients
	 * @param Xdim
	 *          X-dimension of the set of coefficients
	 * @param offset
	 *          Offset of the beginning of the array with respect to the origin of
	 *          c
	 */
	public void setCoefficients(final double[] c, final int Ydim, final int Xdim, final int offset) {
		System.arraycopy(c, offset, coefficient, 0, Ydim * Xdim);
	}

	/**
	 * Get image coefficient weight.
	 *
	 * @return the weight of the coefficient l,m (yWeight, xWeight) in the image
	 *         interpolation
	 */
	public double getWeightI(int l, int m) {
		return yWeight[l] * xWeight[m];
	}

	/**
	 * Prepare for interpolation and interpolate the image value and its
	 * derivatives
	 * 
	 * fromSub = true --> The interpolation is done from the subsampled
	 *                    version of the image
	 * else:
	 * 
	 * fromCurrent=true  --> The interpolation is done
	 *                       from the current image in the pyramid.
	 * fromCurrent=false --> The interpolation is done
	 *                       from the original image.
	 *
	 * @param x x- point coordinate
	 * @param y y- point coordinate
	 * @param D output, interpolation the X and Y derivatives of the image
	 * @param fromSub flat to determine to do the interpolation from the subsampled version of the image 
	 * @param fromCurrent flag to determine the image to do the interpolation from
	 * 		   interpolated value
	 */
	public double prepareForInterpolationAndInterpolateIAndD(double x,
			double y,
			double D[],
			boolean fromSub,
			boolean fromCurrent) {
		int widthToUse = 0;
		int heightToUse = 0;
		final int[] xIndex = new int[4];
		final int[] yIndex = new int[4];
		final double[] xWeight = new double[4];
		final double[] dxWeight = new double[4];
		//double[] d2xWeight = new double[4];
		final double[] yWeight = new double[4];
		final double[] dyWeight = new double[4];
		//double[] d2yWeight = new double[4];
		
		if (fromSub && this.subCoeffs != null)
		{
			widthToUse = this.subWidth;
			heightToUse = this.subHeight;
		}
		else if (fromCurrent)
		{
			widthToUse = currentWidth;
			heightToUse = currentHeight;
		}
		else 
		{
			widthToUse = width;
			heightToUse = height;
		}

		// integer x and y
		int ix = (int)x;
		int iy = (int)y;

		int twiceWidthToUse  = 2 * widthToUse;
		int twiceHeightToUse = 2 * heightToUse;

		// Set X indexes
		// p is the index of the rightmost influencing spline
		int p = (0.0 <= x) ? (ix + 2) : (ix + 1);
		for (int k = 0; k<4; p--, k++) {
			if (coefficientsAreMirrored) {
				int q = (p < 0) ? (-1 - p) : (p);
				if (twiceWidthToUse <= q) q -= twiceWidthToUse * (q / twiceWidthToUse);
				xIndex[k] = (widthToUse <= q) ? (twiceWidthToUse - 1 - q) : (q);
			} else
				xIndex[k] = (p<0 || p>=widthToUse) ? (-1):(p);
		}

		// Set Y indexes
		p = (0.0 <= y) ? (iy + 2) : (iy + 1);
		for (int k = 0; k<4; p--, k++) {
			if (coefficientsAreMirrored) {
				int q = (p < 0) ? (-1 - p) : (p);
				if (twiceHeightToUse <= q) q -= twiceHeightToUse * (q / twiceHeightToUse);
				yIndex[k] = (heightToUse <= q) ? (twiceHeightToUse - 1 - q) : (q);
			} else
				yIndex[k] = (p<0 || p>=heightToUse) ? (-1):(p);
		}

		// Compute how much the sample depart from an integer position
		double ex = x - ((0.0 <= x) ? (ix) : (ix - 1));
		double ey = y - ((0.0 <= y) ? (iy) : (iy - 1));

		// Set X weights for the image and derivative interpolation
		double s = 1.0F - ex;
		dxWeight[0] = 0.5F * ex * ex;
		xWeight[0]  = ex * dxWeight[0] / 3.0F; // Bspline03(x-ix-2)
		dxWeight[3] = -0.5F * s * s;
		xWeight[3]  = s * dxWeight[3] / -3.0F; // Bspline03(x-ix+1)
		dxWeight[1] = 1.0F - 2.0F * dxWeight[0] + dxWeight[3];
		//xWeight[1]  = 2.0F / 3.0F + (1.0F + ex) * dxWeight[3]; // Bspline03(x-ix-1);
		xWeight[1]  = MathTools.Bspline03(x-ix-1);
		dxWeight[2] = 1.5F * ex * (ex - 4.0F/ 3.0F);
		xWeight[2]  = 2.0F / 3.0F - (2.0F - ex) * dxWeight[0]; // Bspline03(x-ix)

		//d2xWeight[0] = ex;
		//d2xWeight[1] = s-2*ex;
		//d2xWeight[2] = ex-2*s;
		//d2xWeight[3] = s;

		// Set Y weights for the image and derivative interpolation
		double t = 1.0F - ey;
		dyWeight[0] = 0.5F * ey * ey;
		yWeight[0]  = ey * dyWeight[0] / 3.0F;
		dyWeight[3] = -0.5F * t * t;
		yWeight[3]  = t * dyWeight[3] / -3.0F;
		dyWeight[1] = 1.0F - 2.0F * dyWeight[0] + dyWeight[3];
		yWeight[1]  = 2.0F / 3.0F + (1.0F + ey) * dyWeight[3];
		dyWeight[2] = 1.5F * ey * (ey - 4.0F/ 3.0F);
		yWeight[2]  = 2.0F / 3.0F - (2.0F - ey) * dyWeight[0];

		//d2yWeight[0] = ey;
		//d2yWeight[1] = t-2*ey;
		//d2yWeight[2] = ey-2*t;
		//d2yWeight[3] = t;
		
		// Image value: Only SplineDegree=3 is implemented
		double ival=0.0F;
		for (int j = 0; j<4; j++) 
		{
			s = 0.0F;
			iy = yIndex[j];
			if (iy!=-1) 
			{
				p = iy*widthToUse;
				for (int i=0; i<4; i++) 
				{
					ix = xIndex[i];
					if (ix!=-1)
					{
						if (fromSub && this.subCoeffs != null)  
							s += xWeight[i] * this.subCoeffs[p + ix];
						else if (fromCurrent)
							s += xWeight[i] * currentCoefficient[p + ix];							
						else             
							s += xWeight[i] * coefficient[p + ix];
					}
				}
				ival+=yWeight[j] * s;
			}
		}
		
		// Derivatives: Only SplineDegree=3 is implemented
		D[0] = D[1] = 0.0F;
		for (int j = 0; j<4; j++) 
		{
			double sx = 0.0F, sy = 0.0F;
			iy = yIndex[j];
			if (iy!=-1) 
			{
				p = iy * widthToUse;
				for (int i=0; i<4; i++) 
				{
					ix = xIndex[i];
					if (ix!=-1) 
					{
						double c;
						if (fromSub && this.subCoeffs != null)  
							c = this.subCoeffs[p + ix];
						else if (fromCurrent) 
							c = currentCoefficient[p + ix];
						else             
							c = coefficient[p + ix];
						sx += dxWeight[i]*c;
						sy +=  xWeight[i]*c;
					}
				}
				D[0]+= yWeight[j] * sx;
				D[1]+=dyWeight[j] * sy;
			}
		}
		
		return ival;
	}


}
