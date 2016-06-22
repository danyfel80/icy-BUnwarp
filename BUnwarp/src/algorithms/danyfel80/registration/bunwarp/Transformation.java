package algorithms.danyfel80.registration.bunwarp;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.roi.ROI2D;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array2DUtil;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import plugins.danyfel80.registration.bunwarp.BUnwarp;
import plugins.kernel.roi.roi2d.ROI2DPoint;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class Transformation {

	/** float epsilon */
	private static final double FLT_EPSILON = (double) Float.intBitsToFloat((int) 0x33FFFFFF);
	/**
	 * pyramid flag to indicate the image information is taken from the pyramid
	 */
	private static final boolean PYRAMID = true;
	/**
	 * original flag to indicate the image information is taken from the original
	 * image
	 */
	private static final boolean ORIGINAL = false;

	/** degree of the B-splines involved in the transformation */
	private final int transformationSplineDegree = 3;

	// Some useful references
	/** reference to the first output image */
	private Sequence outputSeq1;
	/** reference to the second output image */
	private Sequence outputSeq2;
	/** pointer to the BUnwarp EzPlug */
	private BUnwarp plugin;

	// Images
	/** pointer to the source image representation */
	private Sequence sourceSeq;
	/** pointer to the target image representation */
	private Sequence targetSeq;
	/**
	 * pointer to the source image model (source image represented by B-splines)
	 */
	private BSplineModel sourceModel;
	/**
	 * pointer to the target image model (target image represented by B-splines)
	 */
	private BSplineModel targetModel;

	// Original buffered images
	/** initial source buffered image */
	private IcyBufferedImage originalSourceIBI;
	/** initial target buffered image */
	private IcyBufferedImage originalTargetIBI;

	// Landmarks
	/** pointer to the source landmarks */
	private List<ROI2DPoint> sourceLandmarks;
	/** pointer to the target landmarks */
	private List<ROI2DPoint> targetLandmarks;

	// Masks for the images
	/** pointer to the source mask */
	private ROI2D sourceMask;
	/** pointer to the target mask */
	private ROI2D targetMask;

	// Initial affine matrix pre-process
	/** percentage of shear correction in initial matrix */
	private double tweakShear = 0.0;
	/** percentage of scale correction in initial matrix */
	private double tweakScale = 0.0;
	/** percentage of anisotropy correction in initial matrix */
	private double tweakIso = 0.0;

	// Image size
	/** source image height */
	private int sourceHeight;
	/** source image width */
	private int sourceWidth;
	/** target image height */
	private int targetHeight;
	/** source image width */
	private int targetWidth;
	/** target image current height */
	private int targetCurrentHeight;
	/** target image current width */
	private int targetCurrentWidth;
	/** source image current height */
	private int sourceCurrentHeight;
	/** source image current width */
	private int sourceCurrentWidth;
	/** height factor in the target image */
	private double targetFactorHeight;
	/** width factor in the target image */
	private double targetFactorWidth;
	/** height factor in the source image */
	private double sourceFactorHeight;
	/** width factor in the source image */
	private double sourceFactorWidth;

	// Display variables
	/** direct similarity error for the current iteration */
	private double partialDirectSimilarityError;
	/** direct regularization error for the current iteration */
	private double partialDirectRegularizationError;
	/** direct landmarks error for the current iteration */
	private double partialDirectLandmarkError;
	/** direct consistency error for the current iteration */
	private double partialDirectConsitencyError;

	/** direct similarity error at the end of the registration */
	private double finalDirectSimilarityError;
	/** direct regularization error at the end of the registration */
	private double finalDirectRegularizationError;
	/** direct landmarks error at the end of the registration */
	private double finalDirectLandmarkError;
	/** direct consistency error at the end of the registration */
	private double finalDirectConsistencyError;

	/** inverse similarity error for the current iteration */
	private double partialInverseSimilarityError;
	/** inverse regularization error for the current iteration */
	private double partialInverseRegularizationError;
	/** inverse landmarks error for the current iteration */
	private double partialInverseLandmarkError;
	/** inverse consistency error for the current iteration */
	private double partialInverseConsitencyError;

	/** inverse similarity error at the end of the registration */
	private double finalInverseSimilarityError;
	/** inverse regularization error at the end of the registration */
	private double finalInverseRegularizationError;
	/** inverse landmarks error at the end of the registration */
	private double finalInverseLandmarkError;
	/** inverse consistency error at the end of the registration */
	private double finalInverseConsistencyError;

	// Transformation parameters
	/** minimum scale deformation */
	private int minScaleDeformation;
	/** maximum scale deformation */
	private int maxScaleDeformation;
	/** minimum scale image */
	private int minScaleImage;
	/** flag to specify the level of resolution in the output */
	private int outputLevel;
	/** flag to show the optimizer */
	private boolean showMarquardtOptim;
	/** divergence weight */
	private double divWeight;
	/** curl weight */
	private double curlWeight;
	/** landmark weight */
	private double landmarkWeight;
	/** weight for image similarity */
	private double imageWeight;
	/** weight for the deformations consistency */
	private double consistencyWeight;
	/** stopping threshold */
	private double stopThreshold;
	/** level of accuracy */
	private int accurateMode;

	// Transformation estimate
	/** number of intervals to place B-spline coefficients */
	private int intervals;
	/**
	 * x- B-spline coefficients keeping the transformation from source to target
	 */
	private double[][] cxSourceToTarget = null;
	/**
	 * y- B-spline coefficients keeping the transformation from source to target
	 */
	private double[][] cySourceToTarget = null;
	/**
	 * x- B-spline coefficients keeping the transformation from target to source
	 */
	private double[][] cxTargetToSource = null;
	/**
	 * y- B-spline coefficients keeping the transformation from target to source
	 */
	private double[][] cyTargetToSource = null;

	/** image model to interpolate the cxSourceToTarget coefficients */
	private BSplineModel swxSourceToTarget = null;
	/** image model to interpolate the cySourceToTarget coefficients */
	private BSplineModel swySourceToTarget = null;
	/** image model to interpolate the cxTargetToSource coefficients */
	private BSplineModel swxTargetToSource = null;
	/** image model to interpolate the cyTargetToSource coefficients */
	private BSplineModel swyTargetToSource = null;

	// Regularization temporary variables
	/** regularization P11 (source to target) matrix */
	private double[][] P11_SourceToTarget;
	/** regularization P22 (source to target) matrix */
	private double[][] P22_SourceToTarget;
	/** regularization P12 (source to target) matrix */
	private double[][] P12_SourceToTarget;

	/** regularization P11 (target to source) matrix */
	private double[][] P11_TargetToSource;
	/** regularization P22 (target to source) matrix */
	private double[][] P22_TargetToSource;
	/** regularization P12 (target to source) matrix */
	private double[][] P12_TargetToSource;

	/**
	 * Create an instance of Transformation.
	 * 
	 * @param sourceSeq
	 *          Image representation for the source
	 * @param targetSeq
	 *          Image representation for the target
	 * @param sourceModel
	 *          Source image model
	 * @param targetModel
	 *          Target image model
	 * @param sourceLandmarks
	 *          Landmarks in the source image
	 * @param targetLandmarks
	 *          Landmarks in the target image
	 * @param sourceMask
	 *          Source image mask
	 * @param targetMask
	 *          Target image mask
	 * @param minScaleDeformation
	 *          Minimum scale deformation
	 * @param maxScaleDeformation
	 *          Maximum scale deformation
	 * @param minScaleImage
	 *          Minimum image scale
	 * @param divWeight
	 *          Divergence weight
	 * @param curlWeight
	 *          Curl weight
	 * @param landmarkWeight
	 *          Landmark weight
	 * @param imageWeight
	 *          Weight for image similarity
	 * @param consistencyWeight
	 *          Weight for the deformations consistency
	 * @param stopThreshold
	 *          Stopping threshold
	 * @param outputLevel
	 *          Flag to specify the level of resolution in the output
	 * @param showMarquardtOptim
	 *          Flag to show the optimizer
	 * @param accurateMode
	 *          Level of accuracy
	 * @param outputSequence1
	 *          Pointer to the first output image
	 * @param outputSequence2
	 *          Pointer to the second output image
	 * @param plugin
	 *          Pointer to the BUnwarp EzPlug
	 */
	public Transformation(Sequence sourceSeq, Sequence targetSeq, BSplineModel sourceModel, BSplineModel targetModel,
	    List<ROI2DPoint> sourceLandmarks, List<ROI2DPoint> targetLandmarks, ROI2D sourceMask, ROI2D targetMask,
	    int minScaleDeformation, int maxScaleDeformation, int minScaleImage, double divWeight, double curlWeight,
	    double landmarkWeight, double imageWeight, double consistencyWeight, double stopThreshold, int outputLevel,
	    boolean showMarquardtOptim, int accurateMode, Sequence outputSequence1, Sequence outputSequence2,
	    BUnwarp plugin) {
		this.sourceSeq = sourceSeq;
		this.targetSeq = targetSeq;
		this.sourceModel = sourceModel;
		this.targetModel = targetModel;
		this.sourceLandmarks = sourceLandmarks;
		this.targetLandmarks = targetLandmarks;
		this.sourceMask = sourceMask;
		this.targetMask = targetMask;
		this.minScaleDeformation = minScaleDeformation;
		this.maxScaleDeformation = maxScaleDeformation;
		this.minScaleImage = minScaleDeformation;
		this.divWeight = divWeight;
		this.curlWeight = curlWeight;
		this.landmarkWeight = landmarkWeight;
		this.imageWeight = imageWeight;
		this.consistencyWeight = consistencyWeight;
		this.stopThreshold = stopThreshold;
		this.outputLevel = outputLevel;
		this.showMarquardtOptim = showMarquardtOptim;
		this.accurateMode = accurateMode;
		this.outputSeq1 = outputSequence1;
		this.outputSeq2 = outputSequence2;
		this.plugin = plugin;

		this.originalSourceIBI = sourceSeq.getFirstImage();
		this.originalTargetIBI = targetSeq.getFirstImage();

		this.sourceWidth = sourceModel.getWidth();
		this.sourceHeight = sourceModel.getHeight();
		this.targetWidth = targetModel.getWidth();
		this.targetHeight = targetModel.getHeight();
	}

	/**
	 * Unidirectional registration method. It applies unidirectional elastic
	 * registration to the selected source and target images.
	 */
	public void doUnidirectionalRegistration() {
		// This function can only be applied with splines of an odd order

		// Bring into consideration the image/coefficients at the smallest scale
		sourceModel.popFromPyramid();
		targetModel.popFromPyramid();

		targetCurrentHeight = targetModel.getCurrentHeight();
		targetCurrentWidth = targetModel.getCurrentWidth();

		targetFactorHeight = targetModel.getFactorHeight();
		targetFactorWidth = targetModel.getFactorWidth();

		sourceCurrentHeight = sourceModel.getCurrentHeight();
		sourceCurrentWidth = sourceModel.getCurrentWidth();

		sourceFactorHeight = sourceModel.getFactorHeight();
		sourceFactorWidth = sourceModel.getFactorWidth();

		// size correction factor
		int sizeCorrectionFactor = 0; // this.targetHeight / (1024 * (int)
		                              // Math.pow(2,
		                              // this.maxImageSubsamplingFactor));
		// System.out.println("Size correction factor = " + sizeCorrectionFactor);

		// Ask memory for the transformation coefficients
		intervals = (int) Math.pow(2, minScaleDeformation + sizeCorrectionFactor);

		cxTargetToSource = new double[intervals + 3][intervals + 3];
		cyTargetToSource = new double[intervals + 3][intervals + 3];

		// Build matrices for computing the regularization
		buildRegularizationTemporary(intervals, false);

		// Ask for memory for the residues
		final int K;
		if (targetLandmarks != null)
			K = targetLandmarks.size();
		else
			K = 0;
		double[] dxTargetToSource = new double[K];
		double[] dyTargetToSource = new double[K];
		computeInitialResidues(dxTargetToSource, dyTargetToSource, false);

		// Compute the affine transformation FROM THE TARGET TO THE SOURCE
		// coordinates
		// Notice that this matrix is independent of the scale (unless it was loaded
		// from
		// file), but the residues are not
		double[][] affineMatrix = null;
		affineMatrix = computeAffineMatrix(false);

		// MiscTools.printMatrix("source affine matrix", affineMatrix);

		// Incorporate the affine transformation into the spline coefficient
		for (int i = 0; i < intervals + 3; i++) {
			final double v = (double) ((i - 1) * (targetCurrentHeight - 1)) / (double) intervals;
			final double xv = affineMatrix[0][2] + affineMatrix[0][1] * v;
			final double yv = affineMatrix[1][2] + affineMatrix[1][1] * v;
			for (int j = 0; j < intervals + 3; j++) {
				final double u = (double) ((j - 1) * (targetCurrentWidth - 1)) / (double) intervals;
				cxTargetToSource[i][j] = xv + affineMatrix[0][0] * u;
				cyTargetToSource[i][j] = yv + affineMatrix[1][0] * u;
			}
		}

		// Now refine with the different scales
		int state; // state=-1 --> Finish
		// state= 0 --> Increase deformation detail
		// state= 1 --> Increase image detail
		// state= 2 --> Do nothing until the finest image scale
		if (minScaleDeformation == maxScaleDeformation)
			state = 1;
		else
			state = 0;
		int s = minScaleDeformation;
		int step = 0;
		computeTotalWorkload();

		while (state != -1) {
			int currentDepth = targetModel.getCurrentDepth();

			// Update the deformation coefficients only in states 0 and 1
			if (state == 0 || state == 1) {
				// Update the deformation coefficients with the error of the landmarks
				// The following conditional is now useless but it is there to allow
				// easy changes like applying the landmarks only in the coarsest
				// deformation
				if (s >= minScaleDeformation) {
					// Number of intervals at this scale and ask for memory
					intervals = (int) Math.pow(2, s + sizeCorrectionFactor);
					final double[][] newcxTargetToSource = new double[intervals + 3][intervals + 3];
					final double[][] newcyTargetToSource = new double[intervals + 3][intervals + 3];

					// Compute the coefficients at this scale
					boolean underconstrained = true;
					// FROM TARGET TO SOURCE.
					if (divWeight == 0 && curlWeight == 0)
						underconstrained = computeCoefficientsScale(intervals, dxTargetToSource, dyTargetToSource,
						    newcxTargetToSource, newcyTargetToSource, false);
					else
						underconstrained = computeCoefficientsScaleWithRegularization(intervals, dxTargetToSource, dyTargetToSource,
						    newcxTargetToSource, newcyTargetToSource, false);

					// Incorporate information from the previous scale
					if (!underconstrained || (step == 0 && landmarkWeight != 0)) {
						for (int i = 0; i < intervals + 3; i++)
							for (int j = 0; j < intervals + 3; j++) {
								cxTargetToSource[i][j] += newcxTargetToSource[i][j];
								cyTargetToSource[i][j] += newcyTargetToSource[i][j];
							}
					}

				}

				// Optimize deformation coefficients
				optimizeCoeffs(intervals, stopThreshold, cxTargetToSource, cyTargetToSource);

			}

			// Prepare for next iteration
			step++;
			switch (state) {
			case 0:
				// Finer details in the deformation
				if (s < maxScaleDeformation) {
					cxTargetToSource = propagateCoeffsToNextLevel(intervals, cxTargetToSource, 1);
					cyTargetToSource = propagateCoeffsToNextLevel(intervals, cyTargetToSource, 1);
					s++;
					intervals *= 2;

					// Prepare matrices for the regularization term
					buildRegularizationTemporary(intervals, false);

					if (currentDepth > minScaleImage)
						state = 1;
					else
						state = 0;
				} else if (currentDepth > minScaleImage)
					state = 1;
				else
					state = 2;
				break;
			case 1: // Finer details in the image, go on optimizing
			case 2: // Finer details in the image, do not optimize
				// Compute next state
				if (state == 1) {
					if (s == maxScaleDeformation && currentDepth == minScaleImage)
						state = 2;
					else if (s == maxScaleDeformation)
						state = 1;
					else
						state = 0;
				} else if (state == 2) {
					if (currentDepth == 0)
						state = -1;
					else
						state = 2;
				}

				// Pop another image and prepare the deformation
				if (currentDepth != 0) {
					double oldTargetCurrentHeight = targetCurrentHeight;
					double oldTargetCurrentWidth = targetCurrentWidth;

					sourceModel.popFromPyramid();
					targetModel.popFromPyramid();

					targetCurrentHeight = targetModel.getCurrentHeight();
					targetCurrentWidth = targetModel.getCurrentWidth();
					targetFactorHeight = targetModel.getFactorHeight();
					targetFactorWidth = targetModel.getFactorWidth();

					sourceCurrentHeight = sourceModel.getCurrentHeight();
					sourceCurrentWidth = sourceModel.getCurrentWidth();
					sourceFactorHeight = sourceModel.getFactorHeight();
					sourceFactorWidth = sourceModel.getFactorWidth();

					// Adapt the transformation to the new image size
					double targetFactorY = (targetCurrentHeight - 1) / (oldTargetCurrentHeight - 1);
					double targetFactorX = (targetCurrentWidth - 1) / (oldTargetCurrentWidth - 1);

					for (int i = 0; i < intervals + 3; i++)
						for (int j = 0; j < intervals + 3; j++) {
							cxTargetToSource[i][j] *= targetFactorX;
							cyTargetToSource[i][j] *= targetFactorY;
						}

					// Prepare matrices for the regularization term
					buildRegularizationTemporary(intervals, false);
				}
				break;
			}

			// In accurate_mode reduce the stopping threshold for the last iteration
			if ((state == 0 || state == 1) && s == maxScaleDeformation && currentDepth == minScaleImage + 1
			    && accurateMode == 1)
				stopThreshold /= 10;

		} // end while (state != -1).

		// Adapt coefficients if necessary
		if (sourceModel.getOriginalImageWidth() > this.targetCurrentWidth) {
			if (sourceModel.isSubOutput() || targetModel.isSubOutput())
				System.out.println("Adapting coefficients from " + this.sourceCurrentWidth + " to "
				    + this.originalSourceIBI.getWidth() + "...");
			// Adapt the transformation to the new image size
			double targetFactorY = (targetModel.getOriginalImageHeight() - 1) / (targetCurrentHeight - 1);
			double targetFactorX = (targetModel.getOriginalImageWidth() - 1) / (targetCurrentWidth - 1);

			for (int i = 0; i < intervals + 3; i++)
				for (int j = 0; j < intervals + 3; j++) {
					cxTargetToSource[i][j] *= targetFactorX;
					cyTargetToSource[i][j] *= targetFactorY;
				}
			this.targetCurrentHeight = targetModel.getOriginalImageHeight();
			this.targetCurrentWidth = targetModel.getOriginalImageWidth();
			this.sourceCurrentHeight = sourceModel.getOriginalImageHeight();
			this.sourceCurrentWidth = sourceModel.getOriginalImageWidth();
		}

		// Display final errors.
		if (this.outputLevel == 2) {
			if (this.imageWeight != 0) {
				System.out.println(" Optimal direct similarity error = " + this.finalDirectSimilarityError);
			}
			if (this.curlWeight != 0 || this.divWeight != 0) {
				System.out.println(" Optimal direct regularization error = " + this.finalDirectRegularizationError);
			}
			if (this.landmarkWeight != 0) {
				System.out.println(" Optimal direct landmark error = " + this.finalDirectLandmarkError);
			}
			if (this.consistencyWeight != 0) {
				System.out.println(" Optimal direct consistency error = " + this.finalDirectConsistencyError);
			}
		}

	}

	/**
	 * Build regularization temporary.
	 *
	 * @param intervals
	 *          intervals in the deformation
	 * @param bIsReverse
	 *          determines the transformation direction (source-target=TRUE or
	 *          target-source=FALSE)
	 */
	private void buildRegularizationTemporary(int intervals, boolean bIsReverse) {
		// M is the number of spline coefficients per row
		int M = intervals + 3;
		int M2 = M * M;

		// P11
		double[][] P11 = new double[M2][M2];
		if (bIsReverse)
			P11_SourceToTarget = P11;
		else
			P11_TargetToSource = P11;

		for (int i = 0; i < M2; i++)
			for (int j = 0; j < M2; j++)
				P11[i][j] = 0.0;
		build_Matrix_Rq1q2(intervals, divWeight, 2, 0, P11, bIsReverse);
		build_Matrix_Rq1q2(intervals, divWeight + curlWeight, 1, 1, P11, bIsReverse);
		build_Matrix_Rq1q2(intervals, curlWeight, 0, 2, P11, bIsReverse);

		// P22
		double[][] P22 = new double[M2][M2];
		if (bIsReverse)
			P22_SourceToTarget = P22;
		else
			P22_TargetToSource = P22;

		for (int i = 0; i < M2; i++)
			for (int j = 0; j < M2; j++)
				P22[i][j] = 0.0;
		build_Matrix_Rq1q2(intervals, divWeight, 0, 2, P22, bIsReverse);
		build_Matrix_Rq1q2(intervals, divWeight + curlWeight, 1, 1, P22, bIsReverse);
		build_Matrix_Rq1q2(intervals, curlWeight, 2, 0, P22, bIsReverse);

		// P12
		double[][] P12 = new double[M2][M2];
		if (bIsReverse)
			P12_SourceToTarget = P12;
		else
			P12_TargetToSource = P12;

		for (int i = 0; i < M2; i++)
			for (int j = 0; j < M2; j++)
				P12[i][j] = 0.0;
		build_Matrix_Rq1q2q3q4(intervals, 2 * divWeight, 2, 0, 1, 1, P12, bIsReverse);
		build_Matrix_Rq1q2q3q4(intervals, 2 * divWeight, 1, 1, 0, 2, P12, bIsReverse);
		build_Matrix_Rq1q2q3q4(intervals, -2 * curlWeight, 0, 2, 1, 1, P12, bIsReverse);
		build_Matrix_Rq1q2q3q4(intervals, -2 * curlWeight, 1, 1, 2, 0, P12, bIsReverse);
	}

	/**
	 * Build matrix Rq1q2.
	 */
	private void build_Matrix_Rq1q2(int intervals, double weight, int q1, int q2, double[][] R, boolean bIsReverse) {
		build_Matrix_Rq1q2q3q4(intervals, weight, q1, q2, q1, q2, R, bIsReverse);
	}

	/**
	 * Build matrix Rq1q2q3q4.
	 */
	private void build_Matrix_Rq1q2q3q4(int intervals, double weight, int q1, int q2, int q3, int q4, double[][] R,
	    boolean bIsReverse) {
		/*
		 * Let's define alpha_q as the q-th derivative of a B-spline
		 * 
		 * q n d B (x) alpha_q(x)= -------------- q dx
		 * 
		 * eta_q1q2(x,s1,s2)=integral_0^Xdim alpha_q1(x/h-s1) alpha_q2(x/h-s2)
		 * 
		 */
		double[][] etaq1q3 = new double[16][16];
		int Ydim = targetModel.getCurrentHeight();
		int Xdim = targetModel.getCurrentWidth();

		if (bIsReverse) {
			Ydim = sourceModel.getCurrentHeight();
			Xdim = sourceModel.getCurrentWidth();
		}

		build_Matrix_R_geteta(etaq1q3, q1, q3, Xdim, intervals);

		double[][] etaq2q4 = null;
		if (q2 != q1 || q4 != q3 || Ydim != Xdim) {
			etaq2q4 = new double[16][16];
			build_Matrix_R_geteta(etaq2q4, q2, q4, Ydim, intervals);
		} else
			etaq2q4 = etaq1q3;

		int M = intervals + 1;
		int Mp = intervals + 3;
		for (int l = -1; l <= M; l++)
			for (int k = -1; k <= M; k++)
				for (int n = -1; n <= M; n++)
					for (int m = -1; m <= M; m++) {
						int[] ip = new int[2];
						int[] jp = new int[2];
						boolean valid_i = build_Matrix_R_getetaIndex(l, n, intervals, ip);
						boolean valid_j = build_Matrix_R_getetaIndex(k, m, intervals, jp);
						if (valid_i && valid_j) {
							int mn = (n + 1) * Mp + (m + 1);
							int kl = (l + 1) * Mp + (k + 1);
							R[kl][mn] += weight * etaq1q3[jp[0]][jp[1]] * etaq2q4[ip[0]][ip[1]];
						}
					}
	}

	/**
	 * Build matrix R, get eta.
	 */
	private void build_Matrix_R_geteta(double[][] etaq1q2, int q1, int q2, int dim, int intervals) {
		boolean[][] done = new boolean[16][16];
		// Clear
		for (int i = 0; i < 16; i++)
			for (int j = 0; j < 16; j++) {
				etaq1q2[i][j] = 0;
				done[i][j] = false;
			}

		// Compute each integral we need
		int M = intervals + 1;
		double h = (double) dim / intervals;
		for (int ki1 = -1; ki1 <= M; ki1++)
			for (int ki2 = -1; ki2 <= M; ki2++) {
				int[] ip = new int[2];
				boolean valid_i = build_Matrix_R_getetaIndex(ki1, ki2, intervals, ip);
				if (valid_i && !done[ip[0]][ip[1]]) {
					etaq1q2[ip[0]][ip[1]] = build_Matrix_R_computeIntegral_aa(0, dim, ki1, ki2, h, q1, q2);
					done[ip[0]][ip[1]] = true;
				}
			}
	}

	/**
	 * Build matrix R, get eta index.
	 */
	private boolean build_Matrix_R_getetaIndex(int ki1, int ki2, int intervals, int[] ip) {
		ip[0] = 0;
		ip[1] = 0;

		// Determine the clipped inner limits of the intersection
		int kir = Math.min(intervals, Math.min(ki1, ki2) + 2);
		int kil = Math.max(0, Math.max(ki1, ki2) - 2);

		if (kil >= kir)
			return false;

		// Determine which are the pieces of the
		// function that lie in the intersection
		int two_i = 1;
		double ki;
		for (int i = 0; i <= 3; i++, two_i *= 2) {
			// First function
			ki = ki1 + i - 1.5; // Middle sample of the piece i
			if (kil <= ki && ki <= kir)
				ip[0] += two_i;

			// Second function
			ki = ki2 + i - 1.5; // Middle sample of the piece i
			if (kil <= ki && ki <= kir)
				ip[1] += two_i;
		}

		ip[0]--;
		ip[1]--;
		return true;
	}

	/**
	 * Compute the following integral
	 * <P>
	 * 
	 * <PRE>
	 * xF d^q1 3 x d^q2 3 x integral ----- B (--- - s1) ----- B (--- - s2) dx x0
	 * dx^q1 h dx^q2 h
	 * 
	 * <PRE>
	 */
	private double build_Matrix_R_computeIntegral_aa(double x0, double xF, double s1, double s2, double h, int q1,
	    int q2) {
		// Computes the following integral
		//
		// xF d^q1 3 x d^q2 3 x
		// integral ----- B (--- - s1) ----- B (--- - s2) dx
		// x0 dx^q1 h dx^q2 h

		// Form the spline coefficients
		double[][] C = new double[3][3];
		int[][] d = new int[3][3];
		double[][] s = new double[3][3];
		C[0][0] = 1;
		C[0][1] = 0;
		C[0][2] = 0;
		C[1][0] = 1;
		C[1][1] = -1;
		C[1][2] = 0;
		C[2][0] = 1;
		C[2][1] = -2;
		C[2][2] = 1;
		d[0][0] = 3;
		d[0][1] = 0;
		d[0][2] = 0;
		d[1][0] = 2;
		d[1][1] = 2;
		d[1][2] = 0;
		d[2][0] = 1;
		d[2][1] = 1;
		d[2][2] = 1;
		s[0][0] = 0;
		s[0][1] = 0;
		s[0][2] = 0;
		s[1][0] = -0.5;
		s[1][1] = 0.5;
		s[1][2] = 0;
		s[2][0] = 1;
		s[2][1] = 0;
		s[2][2] = -1;

		// Compute the integral
		double integral = 0;
		for (int k = 0; k < 3; k++) {
			double ck = C[q1][k];
			if (ck == 0)
				continue;
			for (int l = 0; l < 3; l++) {
				double cl = C[q2][l];
				if (cl == 0)
					continue;
				integral += ck * cl
				    * build_matrix_R_computeIntegral_BB(x0, xF, s1 + s[q1][k], s2 + s[q2][l], h, d[q1][k], d[q2][l]);
			}
		}
		return integral;
	}

	/**
	 * Compute the following integral
	 * 
	 * <PRE>
	 *           xF   n1  x          n2  x
	 *  integral     B  (--- - s1)  B  (--- - s2) dx
	 *           x0       h              h
	 * </PRE>
	 */
	private double build_matrix_R_computeIntegral_BB(double x0, double xF, double s1, double s2, double h, int n1,
	    int n2) {
		// Computes the following integral
		//
		// xF n1 x n2 x
		// integral B (--- - s1) B (--- - s2) dx
		// x0 h h

		// Change the variable so that the h disappears
		// X=x/h
		double xFp = xF / h;
		double x0p = x0 / h;

		// Form the spline coefficients
		double[] c1 = new double[n1 + 2];
		double fact_n1 = 1;
		for (int k = 2; k <= n1; k++)
			fact_n1 *= k;
		double sign = 1;
		for (int k = 0; k <= n1 + 1; k++, sign *= -1)
			c1[k] = sign * MathTools.nChooseK(n1 + 1, k) / fact_n1;

		double[] c2 = new double[n2 + 2];
		double fact_n2 = 1;
		for (int k = 2; k <= n2; k++)
			fact_n2 *= k;
		sign = 1;
		for (int k = 0; k <= n2 + 1; k++, sign *= -1)
			c2[k] = sign * MathTools.nChooseK(n2 + 1, k) / fact_n2;

		// Compute the integral
		double n1_2 = (double) ((n1 + 1)) / 2.0;
		double n2_2 = (double) ((n2 + 1)) / 2.0;
		double integral = 0;
		for (int k = 0; k <= n1 + 1; k++)
			for (int l = 0; l <= n2 + 1; l++) {
				integral += c1[k] * c2[l] * build_matrix_R_computeIntegral_xx(x0p, xFp, s1 + k - n1_2, s2 + l - n2_2, n1, n2);
			}
		return integral * h;
	}

	/**
	 * <P>
	 * 
	 * <PRE>
	 * Computation of the integral:
	 *             xF          q1       q2
	 *    integral       (x-s1)   (x-s2)     dx
	 *             x0          +        +
	 * </PRE>
	 */
	private double build_matrix_R_computeIntegral_xx(double x0, double xF, double s1, double s2, int q1, int q2) {
		// Computation of the integral
		// xF q1 q2
		// integral (x-s1) (x-s2) dx
		// x0 + +

		// Change of variable so that s1 is 0
		// X=x-s1 => x-s2=X-(s2-s1)
		double s2p = s2 - s1;
		double xFp = xF - s1;
		double x0p = x0 - s1;

		// Now integrate
		if (xFp < 0)
			return 0;

		// Move x0 to the first point where both integrals
		// are distinct from 0
		x0p = Math.max(x0p, Math.max(s2p, 0));
		if (x0p > xFp)
			return 0;

		// There is something to integrate
		// Evaluate the primitive at xF and x0
		double IxFp = 0;
		double Ix0p = 0;
		for (int k = 0; k <= q2; k++) {
			double aux = MathTools.nChooseK(q2, k) / (q1 + k + 1) * Math.pow(-s2p, q2 - k);
			IxFp += Math.pow(xFp, q1 + k + 1) * aux;
			Ix0p += Math.pow(x0p, q1 + k + 1) * aux;
		}

		return IxFp - Ix0p;
	}

	/**
	 * Compute the initial residues for the landmarks.
	 * <p>
	 * NOTE: The output vectors should be already resized
	 *
	 * @param dx
	 *          output, difference in x for each landmark
	 * @param dy
	 *          output, difference in y for each landmark
	 * @param bIsReverse
	 *          determines the transformation direction (source-target=TRUE or
	 *          target-source=FALSE). The output vectors should be already resized
	 */
	private void computeInitialResidues(final double[] dx, final double[] dy,

	    boolean bIsReverse) {

		// Auxiliary variables for registering in both directions.
		double auxFactorWidth = targetModel.getFactorWidth();
		double auxFactorHeight = targetModel.getFactorHeight();
		List<ROI2DPoint> auxSourcePh = this.sourceLandmarks;
		List<ROI2DPoint> auxTargetPh = this.targetLandmarks;

		if (bIsReverse) {
			auxFactorWidth = sourceModel.getFactorWidth();
			auxFactorHeight = sourceModel.getFactorHeight();
			auxSourcePh = this.targetLandmarks;
			auxTargetPh = this.sourceLandmarks;
		}

		List<Point2D> sourceVector = new ArrayList<>();
		if (auxSourcePh != null) {
			for (ROI2DPoint roi : auxSourcePh) {
				sourceVector.add(roi.getPoint());
			}
		}

		List<Point2D> targetVector = new ArrayList<>();
		if (auxTargetPh != null) {
			for (ROI2DPoint roi : auxTargetPh) {
				targetVector.add(roi.getPoint());
			}
		}
		// TODO check if K is correct
		int K = 0;

		if (auxTargetPh != null)
			K = auxTargetPh.size();

		for (int k = 0; k < K; k++) {
			final Point2D sourcePoint = sourceVector.get(k);
			final Point2D targetPoint = targetVector.get(k);
			dx[k] = auxFactorWidth * (sourcePoint.getX() - targetPoint.getX());
			dy[k] = auxFactorHeight * (sourcePoint.getY() - targetPoint.getY());
		}
	}

	/**
	 * Compute the affine matrix.
	 *
	 * @param bIsReverse
	 *          determines the transformation direction (source-target=TRUE or
	 *          target-source=FALSE)
	 */
	private double[][] computeAffineMatrix(boolean bIsReverse) {
		boolean adjust_size = false;

		final double[][] D = new double[3][3];
		final double[][] H = new double[3][3];
		final double[][] U = new double[3][3];
		final double[][] V = new double[3][3];
		final double[][] X = new double[2][3];
		final double[] W = new double[3];

		// Auxiliary variables to calculate inverse transformation
		List<ROI2DPoint> auxSourcePh = sourceLandmarks;
		List<ROI2DPoint> auxTargetPh = targetLandmarks;
		BSplineModel auxSource = sourceModel;
		BSplineModel auxTarget = targetModel;
		double auxFactorWidth = this.targetFactorWidth;
		double auxFactorHeight = this.targetFactorHeight;

		if (bIsReverse) {
			auxSourcePh = targetLandmarks;
			auxTargetPh = sourceLandmarks;
			auxSource = targetModel;
			auxTarget = sourceModel;
			auxFactorWidth = this.sourceFactorWidth;
			auxFactorHeight = this.sourceFactorHeight;
		}

		List<Point2D> sourceVector = new ArrayList<>();
		if (auxSourcePh != null) {
			for (ROI2DPoint roi : auxSourcePh) {
				sourceVector.add(roi.getPoint());
			}
		}

		List<Point2D> targetVector = new ArrayList<>();
		if (auxTargetPh != null) {
			for (ROI2DPoint roi : auxTargetPh) {
				targetVector.add(roi.getPoint());
			}
		}

		int removeLastPoint = 0;

		int n = targetVector.size();
		switch (n) {
		case 0:
			for (int i = 0; (i < 2); i++)
				for (int j = 0; (j < 3); j++)
					X[i][j] = 0.0;
			if (adjust_size) {
				// Make both images of the same size
				X[0][0] = (double) auxSource.getCurrentWidth() / auxTarget.getCurrentWidth();
				X[1][1] = (double) auxSource.getCurrentHeight() / auxTarget.getCurrentHeight();
			} else {
				// Make both images to be centered
				X[0][0] = X[1][1] = 1;
				X[0][2] = ((double) auxSource.getCurrentWidth() - auxTarget.getCurrentWidth()) / 2;
				X[1][2] = ((double) auxSource.getCurrentHeight() - auxTarget.getCurrentHeight()) / 2;
			}
			break;
		case 1:
			for (int i = 0; (i < 2); i++) {
				for (int j = 0; (j < 2); j++) {
					X[i][j] = (i == j) ? (1.0F) : (0.0F);
				}
			}
			X[0][2] = auxFactorWidth * (sourceVector.get(0).getX() - targetVector.get(0).getX());
			X[1][2] = auxFactorHeight * (sourceVector.get(0).getY() - targetVector.get(0).getY());
			break;
		case 2:
			final double x0 = auxFactorWidth * sourceVector.get(0).getX();
			final double y0 = auxFactorHeight * sourceVector.get(0).getY();
			final double x1 = auxFactorWidth * sourceVector.get(1).getX();
			final double y1 = auxFactorHeight * sourceVector.get(1).getY();
			final double u0 = auxFactorWidth * targetVector.get(0).getX();
			final double v0 = auxFactorHeight * targetVector.get(0).getY();
			final double u1 = auxFactorWidth * targetVector.get(1).getX();
			final double v1 = auxFactorHeight * targetVector.get(1).getY();
			sourceVector.add(new Point2D.Double((int) (x1 + y0 - y1), (int) (x1 + y1 - x0)));
			targetVector.add(new Point2D.Double((int) (u1 + v0 - v1), (int) (u1 + v1 - u0)));
			removeLastPoint = 1;
			n = 3;
		default:
			for (int i = 0; (i < 3); i++) {
				for (int j = 0; (j < 3); j++) {
					H[i][j] = 0.0F;
				}
			}
			for (int k = 0; (k < n); k++) {
				final Point2D sourcePoint = sourceVector.get(k);
				final Point2D targetPoint = targetVector.get(k);
				final double sx = auxFactorWidth * sourcePoint.getX();
				final double sy = auxFactorHeight * sourcePoint.getY();
				final double tx = auxFactorWidth * targetPoint.getX();
				final double ty = auxFactorHeight * targetPoint.getY();
				H[0][0] += tx * sx;
				H[0][1] += tx * sy;
				H[0][2] += tx;
				H[1][0] += ty * sx;
				H[1][1] += ty * sy;
				H[1][2] += ty;
				H[2][0] += sx;
				H[2][1] += sy;
				H[2][2] += 1.0F;
				D[0][0] += sx * sx;
				D[0][1] += sx * sy;
				D[0][2] += sx;
				D[1][0] += sy * sx;
				D[1][1] += sy * sy;
				D[1][2] += sy;
				D[2][0] += sx;
				D[2][1] += sy;
				D[2][2] += 1.0F;
			}
			MathTools.singularValueDecomposition(H, W, V);
			if ((Math.abs(W[0]) < FLT_EPSILON) || (Math.abs(W[1]) < FLT_EPSILON) || (Math.abs(W[2]) < FLT_EPSILON)) {
				return (computeRotationMatrix(bIsReverse));
			}
			for (int i = 0; (i < 3); i++) {
				for (int j = 0; (j < 3); j++) {
					V[i][j] /= W[j];
				}
			}
			for (int i = 0; (i < 3); i++) {
				for (int j = 0; (j < 3); j++) {
					U[i][j] = 0.0F;
					for (int k = 0; (k < 3); k++) {
						U[i][j] += D[i][k] * V[k][j];
					}
				}
			}
			for (int i = 0; (i < 2); i++) {
				for (int j = 0; (j < 3); j++) {
					X[i][j] = 0.0F;
					for (int k = 0; (k < 3); k++) {
						X[i][j] += U[i][k] * H[j][k];
					}
				}
			}
			break;
		}
		if (removeLastPoint > 0) {
			for (int i = 1; i <= removeLastPoint; i++) {
				sourceVector.remove(n - i);
				targetVector.remove(n - i);
			}
		}

		double centerX = (double) auxSource.getCurrentWidth();
		double centerY = (double) auxSource.getCurrentHeight();

		final AffineTransform matrix = new AffineTransform(X[0][0], X[1][0], X[0][1], X[1][1], X[0][2], X[1][2]);

		// MiscTools.printMatrix("source affine matrix before reg", X);

		// "Rigidize" matrix based on the user preferences
		regularizeMatrix(matrix, centerX, centerY);

		X[0][0] = matrix.getScaleX();
		X[0][1] = matrix.getShearX();

		X[1][0] = matrix.getShearY();
		X[1][1] = matrix.getScaleY();

		X[0][2] = matrix.getTranslateX();
		X[1][2] = matrix.getTranslateY();

		// MiscTools.printMatrix("source affine matrix after reg", X);

		return (X);
	}

	/**
	 * Compute the rotation matrix.
	 *
	 * @param bIsReverse
	 *          determines the transformation direction (source-target=TRUE or
	 *          target-source=FALSE)
	 * @return rotation matrix
	 */
	private double[][] computeRotationMatrix(boolean bIsReverse) {
		final double[][] X = new double[2][3];
		final double[][] H = new double[2][2];
		final double[][] V = new double[2][2];
		final double[] W = new double[2];

		double auxFactorWidth = this.targetFactorWidth;
		double auxFactorHeight = this.targetFactorHeight;
		List<ROI2DPoint> auxSourcePh = this.sourceLandmarks;
		List<ROI2DPoint> auxTargetPh = this.targetLandmarks;
		if (bIsReverse) {
			auxFactorWidth = this.sourceFactorWidth;
			auxFactorHeight = this.sourceFactorHeight;
			auxSourcePh = this.targetLandmarks;
			auxTargetPh = this.sourceLandmarks;
		}

		List<Point2D> sourceVector = new ArrayList<>();
		if (auxSourcePh != null) {
			for (ROI2DPoint roi : auxSourcePh) {
				sourceVector.add(roi.getPoint());
			}
		}

		List<Point2D> targetVector = new ArrayList<>();
		if (auxTargetPh != null) {
			for (ROI2DPoint roi : auxTargetPh) {
				targetVector.add(roi.getPoint());
			}
		}

		final int n = targetVector.size();
		switch (n) {
		case 0:
			for (int i = 0; (i < 2); i++) {
				for (int j = 0; (j < 3); j++) {
					X[i][j] = (i == j) ? (1.0F) : (0.0F);
				}
			}
			break;
		case 1:
			for (int i = 0; (i < 2); i++) {
				for (int j = 0; (j < 2); j++) {
					X[i][j] = (i == j) ? (1.0F) : (0.0F);
				}
			}
			X[0][2] = auxFactorWidth * (sourceVector.get(0).getX() - targetVector.get(0).getX());
			X[1][2] = auxFactorHeight * (sourceVector.get(0).getY() - targetVector.get(0).getY());
			break;
		default:
			double xTargetAverage = 0.0F;
			double yTargetAverage = 0.0F;

			for (int i = 0; (i < n); i++) {
				final Point2D p = targetVector.get(i);
				xTargetAverage += auxFactorWidth * p.getX();
				yTargetAverage += auxFactorHeight * p.getY();
			}

			xTargetAverage /= (double) n;
			yTargetAverage /= (double) n;

			final double[] xCenteredTarget = new double[n];
			final double[] yCenteredTarget = new double[n];

			for (int i = 0; (i < n); i++) {
				final Point2D p = targetVector.get(i);
				xCenteredTarget[i] = auxFactorWidth * p.getX() - xTargetAverage;
				yCenteredTarget[i] = auxFactorHeight * p.getY() - yTargetAverage;
			}

			double xSourceAverage = 0.0F;
			double ySourceAverage = 0.0F;

			for (int i = 0; (i < n); i++) {
				final Point2D p = sourceVector.get(i);
				xSourceAverage += auxFactorWidth * p.getX();
				ySourceAverage += auxFactorHeight * p.getY();
			}

			xSourceAverage /= (double) n;
			ySourceAverage /= (double) n;

			final double[] xCenteredSource = new double[n];
			final double[] yCenteredSource = new double[n];

			for (int i = 0; (i < n); i++) {
				final Point2D p = sourceVector.get(i);
				xCenteredSource[i] = auxFactorWidth * (double) p.getX() - xSourceAverage;
				yCenteredSource[i] = auxFactorHeight * (double) p.getY() - ySourceAverage;
			}

			for (int i = 0; (i < 2); i++) {
				for (int j = 0; (j < 2); j++) {
					H[i][j] = 0.0F;
				}
			}

			for (int k = 0; (k < n); k++) {
				H[0][0] += xCenteredTarget[k] * xCenteredSource[k];
				H[0][1] += xCenteredTarget[k] * yCenteredSource[k];
				H[1][0] += yCenteredTarget[k] * xCenteredSource[k];
				H[1][1] += yCenteredTarget[k] * yCenteredSource[k];
			}
			// COSS: Watch out that this H is the transpose of the one
			// defined in the text. That is why X=V*U^t is the inverse
			// of the rotation matrix.
			MathTools.singularValueDecomposition(H, W, V);
			if (((H[0][0] * H[1][1] - H[0][1] * H[1][0]) * (V[0][0] * V[1][1] - V[0][1] * V[1][0])) < 0.0F) {
				if (W[0] < W[1]) {
					V[0][0] *= -1.0F;
					V[1][0] *= -1.0F;
				} else {
					V[0][1] *= -1.0F;
					V[1][1] *= -1.0F;
				}
			}
			for (int i = 0; (i < 2); i++) {
				for (int j = 0; (j < 2); j++) {
					X[i][j] = 0.0F;
					for (int k = 0; (k < 2); k++) {
						X[i][j] += V[i][k] * H[j][k];
					}
				}
			}
			X[0][2] = xSourceAverage - X[0][0] * xTargetAverage - X[0][1] * yTargetAverage;
			X[1][2] = ySourceAverage - X[1][0] * xTargetAverage - X[1][1] * yTargetAverage;
			break;
		}
		return (X);
	}

	/**
	 * Regularize matrix to remove sharing, scaling, etc.
	 * 
	 * @param a
	 *          affine transform
	 * @param centerX
	 *          image center x- coordinate
	 * @param centerX
	 *          image center y- coordinate
	 */
	public void regularizeMatrix(AffineTransform a, double centerX, double centerY) {

		// Move to the center of the image
		a.translate(centerX, centerY);

		/*
		 * IJ.log(" A: " + a.getScaleX() + " " + a.getShearY() + " " + a.getShearX()
		 * + " " + a.getScaleY() + " " + a.getTranslateX() + " " + +
		 * a.getTranslateY() );
		 */

		// retrieves scaling, shearing, rotation and translation from an affine
		// transformation matrix A (which has translation values in the right
		// column)
		// by Daniel Berger for MIT-BCS Seung, April 19 2009

		// We assume that sheary=0
		// scalex=sqrt(A(1,1)*A(1,1)+A(2,1)*A(2,1));
		final double a11 = a.getScaleX();
		final double a21 = a.getShearY();
		final double scaleX = Math.sqrt(a11 * a11 + a21 * a21);
		// rotang=atan2(A(2,1)/scalex,A(1,1)/scalex);
		final double rotang = Math.atan2(a21 / scaleX, a11 / scaleX);

		// R=[[cos(-rotang) -sin(-rotang)];[sin(-rotang) cos(-rotang)]];

		// rotate back shearx and scaley
		// v=R*[A(1,2) A(2,2)]';
		final double a12 = a.getShearX();
		final double a22 = a.getScaleY();
		final double shearX = Math.cos(-rotang) * a12 - Math.sin(-rotang) * a22;
		final double scaleY = Math.sin(-rotang) * a12 + Math.cos(-rotang) * a22;

		// rotate back translation
		// v=R*[A(1,3) A(2,3)]';
		final double transX = Math.cos(-rotang) * a.getTranslateX() - Math.sin(-rotang) * a.getTranslateY();
		final double transY = Math.sin(-rotang) * a.getTranslateX() + Math.cos(-rotang) * a.getTranslateY();

		// TWEAK

		final double new_shearX = shearX * (1.0 - tweakShear);
		// final double new_shearY = 0; // shearY * (1.0 - tweakShear);

		final double avgScale = (scaleX + scaleY) / 2;
		final double aspectRatio = scaleX / scaleY;
		final double regAvgScale = avgScale * (1.0 - tweakScale) + 1.0 * tweakScale;
		final double regAspectRatio = aspectRatio * (1.0 - tweakIso) + 1.0 * tweakIso;

		// IJ.log("avgScale = " + avgScale + " aspectRatio = " + aspectRatio + "
		// regAvgScale = " + regAvgScale + " regAspectRatio = " + regAspectRatio);

		final double new_scaleY = (2.0 * regAvgScale) / (regAspectRatio + 1.0);
		final double new_scaleX = regAspectRatio * new_scaleY;

		final AffineTransform b = makeAffineMatrix(new_scaleX, new_scaleY, new_shearX, 0, rotang, transX, transY);

		// IJ.log("new_scaleX = " + new_scaleX + " new_scaleY = " + new_scaleY + "
		// new_shearX = " + new_shearX + " new_shearY = " + new_shearY);

		// Move back the center
		b.translate(-centerX, -centerY);

		a.setTransform(b);

	}

	/**
	 * Makes an affine transformation matrix from the given scale, shear, rotation
	 * and translation values if you want a uniquely retrievable matrix, give
	 * sheary=0
	 * 
	 * @param scalex
	 *          scaling in x
	 * @param scaley
	 *          scaling in y
	 * @param shearx
	 *          shearing in x
	 * @param sheary
	 *          shearing in y
	 * @param rotang
	 *          angle of rotation (in radians)
	 * @param transx
	 *          translation in x
	 * @param transy
	 *          translation in y
	 * @return affine transformation matrix
	 */
	public static AffineTransform makeAffineMatrix(final double scalex, final double scaley, final double shearx,
	    final double sheary, final double rotang, final double transx, final double transy) {
		/*
		 * %makes an affine transformation matrix from the given scale, shear,
		 * %rotation and translation values %if you want a uniquely retrievable
		 * matrix, give sheary=0 %by Daniel Berger for MIT-BCS Seung, April 19 2009
		 * 
		 * A=[[scalex shearx transx];[sheary scaley transy];[0 0 1]];
		 * A=[[cos(rotang) -sin(rotang) 0];[sin(rotang) cos(rotang) 0];[0 0 1]] * A;
		 */

		final double m00 = Math.cos(rotang) * scalex - Math.sin(rotang) * sheary;
		final double m01 = Math.cos(rotang) * shearx - Math.sin(rotang) * scaley;
		final double m02 = Math.cos(rotang) * transx - Math.sin(rotang) * transy;

		final double m10 = Math.sin(rotang) * scalex + Math.cos(rotang) * sheary;
		final double m11 = Math.sin(rotang) * shearx + Math.cos(rotang) * scaley;
		final double m12 = Math.sin(rotang) * transx + Math.cos(rotang) * transy;

		return new AffineTransform(m00, m10, m01, m11, m02, m12);
	}

	/**
	 * This code is an excerpt from doBidirectionalRegistration() to compute the
	 * exact number of steps.
	 */
	private void computeTotalWorkload() {
		// This code is an excerpt from doBidirectionalRegistration() to compute the
		// exact
		// number of steps

		// Now refine with the different scales
		int state; // state=-1 --> Finish
		// state= 0 --> Increase deformation detail
		// state= 1 --> Increase image detail
		// state= 2 --> Do nothing until the finest image scale
		if (minScaleDeformation == maxScaleDeformation)
			state = 1;
		else
			state = 0;
		int s = minScaleDeformation;
		int currentDepth = targetModel.getCurrentDepth();
		int workload = 0;
		while (state != -1) {
			// Update the deformation coefficients only in states 0 and 1
			if (state == 0 || state == 1) {
				// Optimize deformation coefficients
				if (imageWeight != 0)
					workload += 300 * (currentDepth + 1);
			}

			// Prepare for next iteration
			switch (state) {
			case 0:
				// Finer details in the deformation
				if (s < maxScaleDeformation) {
					s++;
					if (currentDepth > minScaleImage)
						state = 1;
					else
						state = 0;
				} else if (currentDepth > minScaleImage)
					state = 1;
				else
					state = 2;
				break;
			case 1: // Finer details in the image, go on optimizing
			case 2: // Finer details in the image, do not optimize
				// Compute next state
				if (state == 1) {
					if (s == maxScaleDeformation && currentDepth == minScaleImage)
						state = 2;
					else if (s == maxScaleDeformation)
						state = 1;
					else
						state = 0;
				} else if (state == 2) {
					if (currentDepth == 0)
						state = -1;
					else
						state = 2;
				}

				// Pop another image and prepare the deformation
				if (currentDepth != 0)
					currentDepth--;
				break;
			}
		}

		ProgressBar.resetProgressBar();
		ProgressBar.addWorkload(workload);
	}

	/**
	 * Compute the coefficients at this scale.
	 *
	 * @param intervals
	 *          input, number of intervals at this scale
	 * @param dx
	 *          input, x residue so far
	 * @param dy
	 *          input, y residue so far
	 * @param cx
	 *          output, x coefficients for splines
	 * @param cy
	 *          output, y coefficients for splines
	 * @param bIsReverse
	 *          determines the transformation direction (source-target=TRUE or
	 *          target-source=FALSE)
	 * @return under-constrained flag
	 */
	private boolean computeCoefficientsScale(final int intervals, // input, number
	                                                              // of intervals
	                                                              // at this scale
	    final double[] dx, // input, x residue so far
	    final double[] dy, // input, y residue so far
	    final double[][] cx, // output, x coefficients for splines
	    final double[][] cy, // output, y coefficients for splines
	    boolean bIsReverse) {

		List<ROI2DPoint> auxTargetPh = (bIsReverse) ? this.sourceLandmarks : this.targetLandmarks;

		int K = 0;
		if (auxTargetPh != null)
			K = auxTargetPh.size();
		boolean underconstrained = false;

		if (0 < K) {
			// Form the equation system Bc=d
			final double[][] B = new double[K][(intervals + 3) * (intervals + 3)];
			buildMatrixB(intervals, K, B, bIsReverse);

			// "Invert" the matrix B
			int Nunk = (intervals + 3) * (intervals + 3);
			double[][] iB = new double[Nunk][K];
			underconstrained = MathTools.invertMatrixSVD(K, Nunk, B, iB);

			// Now multiply iB times dx and dy respectively
			int ij = 0;
			for (int i = 0; i < intervals + 3; i++)
				for (int j = 0; j < intervals + 3; j++) {
					cx[i][j] = cy[i][j] = 0.0F;
					for (int k = 0; k < K; k++) {
						cx[i][j] += iB[ij][k] * dx[k];
						cy[i][j] += iB[ij][k] * dy[k];
					}
					ij++;
				}
		}
		return underconstrained;
	}

	/**
	 * Build the matrix for the landmark interpolation.
	 *
	 * @param intervals
	 *          Intervals in the deformation
	 * @param K
	 *          Number of landmarks
	 * @param B
	 *          System matrix of the landmark interpolation
	 * @param bIsReverse
	 *          determines the transformation direction (source-target=TRUE or
	 *          target-source=FALSE)
	 */
	private void buildMatrixB(int intervals, // Intervals in the deformation
	    int K, // Number of landmarks
	    double[][] B, // System matrix of the landmark interpolation
	    boolean bIsReverse) {

		// Auxiliary variables to calculate inverse transformation
		List<ROI2DPoint> auxTargetPh = this.targetLandmarks;
		double auxFactorWidth = this.targetFactorWidth;
		double auxFactorHeight = this.targetFactorHeight;

		if (bIsReverse) {
			auxTargetPh = this.sourceLandmarks;
			auxFactorWidth = this.sourceFactorWidth;
			auxFactorHeight = this.sourceFactorHeight;
		}

		List<Point2D> targetVector = new ArrayList<>();
		if (auxTargetPh != null) {
			for (ROI2DPoint roi : auxTargetPh) {
				targetVector.add(roi.getPoint());
			}
		}

		for (int k = 0; k < K; k++) {
			final Point2D targetPoint = targetVector.get(k);
			double x = auxFactorWidth * targetPoint.getX();
			double y = auxFactorHeight * targetPoint.getY();
			final double[] bx = xWeight(x, intervals, true, bIsReverse);
			final double[] by = yWeight(y, intervals, true, bIsReverse);
			for (int i = 0; i < intervals + 3; i++)
				for (int j = 0; j < intervals + 3; j++)
					B[k][(intervals + 3) * i + j] = by[i] * bx[j];
		}
	}

	/**
	 * Calculate the cubic B-spline x weight.
	 *
	 * @param x
	 *          x- value
	 * @param xIntervals
	 *          x- number of intervals
	 * @param extended
	 *          extended flat
	 * @param bIsReverse
	 *          flag to determine the transformation direction
	 *          (target-source=FALSE or source-target=TRUE)
	 * @return weights
	 */
	private double[] xWeight(final double x, final int xIntervals, final boolean extended, boolean bIsReverse) {
		int auxTargetCurrentWidth = (bIsReverse) ? this.sourceCurrentWidth : this.targetCurrentWidth;

		int length = xIntervals + 1;
		int j0 = 0, jF = xIntervals;
		if (extended) {
			length += 2;
			j0--;
			jF++;
		}
		final double[] b = new double[length];
		final double interX = (double) xIntervals / (double) (auxTargetCurrentWidth - 1);
		for (int j = j0; j <= jF; j++) {
			b[j - j0] = MathTools.Bspline03(x * interX - (double) j);
		}
		return (b);
	}

	/**
	 * Calculate the cubic B-spline y weight.
	 *
	 * @param y
	 *          y- value
	 * @param yIntervals
	 *          y- number of intervals
	 * @param extended
	 *          extended flat
	 * @param bIsReverse
	 *          flag to determine the transformation direction
	 *          (target-source=FALSE or source-target=TRUE)
	 * @return weights
	 */
	private double[] yWeight(final double y, final int yIntervals, final boolean extended, boolean bIsReverse) {

		int auxTargetCurrentHeight = (bIsReverse) ? this.sourceCurrentHeight : this.targetCurrentHeight;

		int length = yIntervals + 1;
		int i0 = 0, iF = yIntervals;
		if (extended) {
			length += 2;
			i0--;
			iF++;
		}
		final double[] b = new double[length];
		final double interY = (double) yIntervals / (double) (auxTargetCurrentHeight - 1);
		for (int i = i0; i <= iF; i++) {
			b[i - i0] = MathTools.Bspline03(y * interY - (double) i);
		}
		return (b);
	}

	/**
	 * Compute the coefficients scale with regularization.
	 *
	 * @param intervals
	 *          input, number of intervals at this scale
	 * @param dx
	 *          input, x residue so far
	 * @param dy
	 *          input, y residue so far
	 * @param cx
	 *          output, x coefficients for splines
	 * @param cy
	 *          output, y coefficients for splines
	 * @param bIsReverse
	 *          determines the transformation direction (source-target=TRUE or
	 *          target-source=FALSE)
	 * @return under-constrained flag
	 */
	private boolean computeCoefficientsScaleWithRegularization(final int intervals, final double[] dx, final double[] dy,
	    final double[][] cx, final double[][] cy, boolean bIsReverse) {

		double P11[][] = this.P11_TargetToSource;
		double P12[][] = this.P12_TargetToSource;
		double P22[][] = this.P22_TargetToSource;

		List<ROI2DPoint> auxTargetPh = this.targetLandmarks;
		if (bIsReverse) {
			auxTargetPh = this.sourceLandmarks;
			P11 = this.P11_SourceToTarget;
			P12 = this.P12_SourceToTarget;
			P22 = this.P22_SourceToTarget;
		}

		boolean underconstrained = true;
		int K = 0;
		if (auxTargetPh != null)
			K = auxTargetPh.size();

		if (0 < K) {
			// M is the number of spline coefficients per row
			int M = intervals + 3;
			int M2 = M * M;

			// Create A and b for the system Ac=b
			final double[][] A = new double[2 * M2][2 * M2];
			final double[] b = new double[2 * M2];
			for (int i = 0; i < 2 * M2; i++) {
				b[i] = 0.0;
				for (int j = 0; j < 2 * M2; j++)
					A[i][j] = 0.0;
			}

			// Get the matrix related to the landmarks
			final double[][] B = new double[K][M2];
			buildMatrixB(intervals, K, B, bIsReverse);

			// Fill the part of the equation system related to the landmarks
			// Compute 2 * B^t * B
			for (int i = 0; i < M2; i++) {
				for (int j = i; j < M2; j++) {
					double bitbj = 0; // bi^t * bj, i.e., column i x column j
					for (int l = 0; l < K; l++)
						bitbj += B[l][i] * B[l][j];
					bitbj *= 2;
					// int ij=i*M2+j;
					A[M2 + i][M2 + j] = A[M2 + j][M2 + i] = A[i][j] = A[j][i] = bitbj;
				}
			}

			// Compute 2 * B^t * [dx dy]
			for (int i = 0; i < M2; i++) {
				double bitdx = 0;
				double bitdy = 0;
				for (int l = 0; l < K; l++) {
					bitdx += B[l][i] * dx[l];
					bitdy += B[l][i] * dy[l];
				}
				bitdx *= 2;
				bitdy *= 2;
				b[i] = bitdx;
				b[M2 + i] = bitdy;
			}

			// Get the matrices associated to the regularization
			// Copy P11 symmetrized to the equation system
			for (int i = 0; i < M2; i++)
				for (int j = 0; j < M2; j++) {
					double aux = P11[i][j];
					A[i][j] += aux;
					A[j][i] += aux;
				}

			// Copy P22 symmetrized to the equation system
			for (int i = 0; i < M2; i++)
				for (int j = 0; j < M2; j++) {
					double aux = P22[i][j];
					A[M2 + i][M2 + j] += aux;
					A[M2 + j][M2 + i] += aux;
				}

			// Copy P12 and P12^t to their respective places
			for (int i = 0; i < M2; i++)
				for (int j = 0; j < M2; j++) {
					A[i][M2 + j] = P12[i][j]; // P12
					A[M2 + i][j] = P12[j][i]; // P12^t
				}

			// Now solve the system
			// Invert the matrix A
			double[][] iA = new double[2 * M2][2 * M2];
			underconstrained = MathTools.invertMatrixSVD(2 * M2, 2 * M2, A, iA);

			// Now multiply iB times b and distribute in cx and cy
			int ij = 0;
			for (int i = 0; i < intervals + 3; i++)
				for (int j = 0; j < intervals + 3; j++) {
					cx[i][j] = cy[i][j] = 0.0F;
					for (int l = 0; l < 2 * M2; l++) {
						cx[i][j] += iA[ij][l] * b[l];
						cy[i][j] += iA[M2 + ij][l] * b[l];
					}
					ij++;
				}
		}
		return underconstrained;
	}

	/**
	 * Optimize the B-spline coefficients (bidirectional method).
	 *
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param thChangef
	 * @param cxTargetToSource
	 *          x- B-spline coefficients storing the target to source deformation
	 * @param cyTargetToSource
	 *          y- B-spline coefficients storing the target to source deformation
	 * @param cxSourceToTarget
	 *          x- B-spline coefficients storing the source to target deformation
	 * @param cySourceToTarget
	 *          y- B-spline coefficients storing the source to target deformation
	 * @return energy function value
	 */
	private double optimizeCoeffs(int intervals, double thChangef, double[][] cxTargetToSource,
	    double[][] cyTargetToSource, double[][] cxSourceToTarget, double[][] cySourceToTarget) {
		if (plugin != null && plugin.isPluginInterrumped())
			return 0.0;

		if (sourceModel.isSubOutput()) {
			System.out.println(" -----\n Intervals = " + intervals + "x" + intervals);
			System.out.println(" Source Image Size = " + this.sourceCurrentWidth + "x" + this.sourceCurrentHeight);
		}

		final double TINY = FLT_EPSILON;
		final double EPS = 3.0e-8F;
		final double FIRSTLAMBDA = 1;
		final int MAXITER_OPTIMCOEFF = 300;
		final int CUMULATIVE_SIZE = 5;

		int int3 = intervals + 3;
		int halfM = 2 * int3 * int3;
		int quarterM = halfM / 2;
		int threeQuarterM = quarterM * 3;
		int M = halfM * 2;

		double rescuedf, f;
		double[] x = new double[M];
		double[] rescuedx = new double[M];
		double[] diffx = new double[M];
		double[] rescuedgrad = new double[M];
		double[] grad = new double[M];
		double[] diffgrad = new double[M];
		double[] Hdx = new double[M];
		double[] rescuedhess = new double[M * M];
		double[] hess = new double[M * M];
		// double []safehess = new double [M*M];
		double[] proposedHess = new double[M * M];
		boolean[] optimize = new boolean[M];
		int i, j, p, iter = 1;
		boolean skip_update, ill_hessian;
		double improvementx = (double) Math.sqrt(TINY), lambda = FIRSTLAMBDA, max_normx, distx, aux, gmax;
		double fac, fae, dgdx, dxHdx, sumdiffg, sumdiffx;

		CumulativeQueue lastBest = new CumulativeQueue(CUMULATIVE_SIZE);

		for (i = 0; i < M; i++)
			optimize[i] = true;

		/* Form the vector with the current guess for the optimization */
		for (i = 0, p = 0; i < intervals + 3; i++)
			for (j = 0; j < intervals + 3; j++, p++) {
				x[p] = cxTargetToSource[i][j];
				x[quarterM + p] = cxSourceToTarget[i][j];

				x[halfM + p] = cyTargetToSource[i][j];
				x[threeQuarterM + p] = cySourceToTarget[i][j];
			}

		/* Prepare the precomputed weights for interpolation */
		this.swxTargetToSource = new BSplineModel(x, intervals + 3, intervals + 3, 0);
		this.swyTargetToSource = new BSplineModel(x, intervals + 3, intervals + 3, halfM);
		this.swxTargetToSource.precomputedPrepareForInterpolation(targetModel.getCurrentHeight(),
		    targetModel.getCurrentWidth(), intervals);
		this.swyTargetToSource.precomputedPrepareForInterpolation(targetModel.getCurrentHeight(),
		    targetModel.getCurrentWidth(), intervals);

		this.swxSourceToTarget = new BSplineModel(x, intervals + 3, intervals + 3, quarterM);
		this.swySourceToTarget = new BSplineModel(x, intervals + 3, intervals + 3, threeQuarterM);
		this.swxSourceToTarget.precomputedPrepareForInterpolation(sourceModel.getCurrentHeight(),
		    sourceModel.getCurrentWidth(), intervals);
		this.swySourceToTarget.precomputedPrepareForInterpolation(sourceModel.getCurrentHeight(),
		    sourceModel.getCurrentWidth(), intervals);

		/* First computation of the energy */
		f = energyFunction(x, intervals, grad, false, false);

		if (showMarquardtOptim)
			System.out.println("f(1)=" + f);

		/*
		 * Initially the hessian is the identity matrix multiplied by the first
		 * function value
		 */
		for (i = 0, p = 0; i < M; i++)
			for (j = 0; j < M; j++, p++)
				if (i == j)
					hess[p] = 1.0F;
				else
					hess[p] = 0.0F;

		rescuedf = f;
		for (i = 0, p = 0; i < M; i++) {
			rescuedx[i] = x[i];
			rescuedgrad[i] = grad[i];
			for (j = 0; j < M; j++, p++)
				rescuedhess[p] = hess[p];
		}

		// Maximum iteration number
		int maxiter = MAXITER_OPTIMCOEFF * (sourceModel.getCurrentDepth() + 1);

		ProgressBar.stepProgressBar();

		int last_successful_iter = 0;

		boolean stop = plugin != null && plugin.isPluginInterrumped();

		while (iter < maxiter && !stop) {
			/* Compute new x ------------------------------------------------- */
			Marquardt_it(x, optimize, grad, hess, lambda);

			/* Stopping criteria --------------------------------------------- */
			/* Compute difference with the previous iteration */
			max_normx = improvementx = 0;
			for (i = 0; i < M; i++) {
				diffx[i] = x[i] - rescuedx[i];
				distx = Math.abs(diffx[i]);
				improvementx += distx * distx;
				aux = Math.abs(rescuedx[i]) < Math.abs(x[i]) ? x[i] : rescuedx[i];
				max_normx += aux * aux;
			}

			if (TINY < max_normx)
				improvementx = improvementx / max_normx;

			improvementx = (double) Math.sqrt(Math.sqrt(improvementx));

			/*
			 * If there is no change with respect to the old geometry then finish the
			 * iterations
			 */
			if (improvementx < Math.sqrt(TINY))
				break;

			/* Estimate the new function value -------------------------------- */
			f = energyFunction(x, intervals, grad, false, false);
			iter++;
			if (showMarquardtOptim)
				System.out.println("f(" + iter + ")=" + f + " lambda=" + lambda);
			ProgressBar.stepProgressBar();

			/* Update lambda -------------------------------------------------- */
			if (rescuedf > f) {
				// We save the last energy terms values in order to be displayed.
				this.finalDirectConsistencyError = this.partialDirectConsitencyError;
				this.finalDirectSimilarityError = this.partialDirectSimilarityError;
				this.finalDirectRegularizationError = this.partialDirectRegularizationError;
				this.finalDirectLandmarkError = this.partialDirectLandmarkError;

				this.finalInverseConsistencyError = this.partialInverseConsitencyError;
				this.finalInverseSimilarityError = this.partialInverseSimilarityError;
				this.finalInverseRegularizationError = this.partialInverseRegularizationError;
				this.finalInverseLandmarkError = this.partialInverseLandmarkError;

				/* Check if the improvement is only residual */
				lastBest.push_back(rescuedf - f);
				if (lastBest.currentSize() == CUMULATIVE_SIZE && lastBest.getSum() / f < thChangef)
					break;

				/*
				 * If we have improved then estimate the hessian, update the geometry,
				 * and decrease the lambda
				 */
				/* Estimate the hessian ....................................... */
				if (showMarquardtOptim)
					System.out.println("  Accepted");
				if ((last_successful_iter++ % 10) == 0 && outputLevel > -1)
					updateOutputs(x, intervals);

				/* Estimate the difference between gradients */
				for (i = 0; i < M; i++)
					diffgrad[i] = grad[i] - rescuedgrad[i];

				/* Multiply this difference by the current inverse of the hessian */
				for (i = 0, p = 0; i < M; i++) {
					Hdx[i] = 0.0F;
					for (j = 0; j < M; j++, p++)
						Hdx[i] += hess[p] * diffx[j];
				}

				/* Calculate dot products for the denominators ................ */
				dgdx = dxHdx = sumdiffg = sumdiffx = 0.0F;
				skip_update = true;
				for (i = 0; i < M; i++) {
					dgdx += diffgrad[i] * diffx[i];
					dxHdx += diffx[i] * Hdx[i];
					sumdiffg += diffgrad[i] * diffgrad[i];
					sumdiffx += diffx[i] * diffx[i];
					if (Math.abs(grad[i]) >= Math.abs(rescuedgrad[i]))
						gmax = Math.abs(grad[i]);
					else
						gmax = Math.abs(rescuedgrad[i]);
					if (gmax != 0 && Math.abs(diffgrad[i] - Hdx[i]) > Math.sqrt(EPS) * gmax)
						skip_update = false;
				}

				/* Update hessian ............................................. */
				/* Skip if fac not sufficiently positive */
				if (dgdx > Math.sqrt(EPS * sumdiffg * sumdiffx) && !skip_update) {
					fae = 1.0F / dxHdx;
					fac = 1.0F / dgdx;

					/* Update the hessian after BFGS formula */
					for (i = 0, p = 0; i < M; i++)
						for (j = 0; j < M; j++, p++) {
							if (i <= j)
								proposedHess[p] = hess[p] + fac * diffgrad[i] * diffgrad[j] - fae * (Hdx[i] * Hdx[j]);
							else
								proposedHess[p] = proposedHess[j * M + i];
						}

					ill_hessian = false;
					if (!ill_hessian) {
						for (i = 0, p = 0; i < M; i++)
							for (j = 0; j < M; j++, p++)
								hess[p] = proposedHess[p];
					} else if (showMarquardtOptim)
						System.out.println("Hessian cannot be safely updated, ill-conditioned");

				} else if (showMarquardtOptim)
					System.out.println("Hessian cannot be safely updated");

				/* Update geometry and lambda ................................. */
				rescuedf = f;
				for (i = 0, p = 0; i < M; i++) {
					rescuedx[i] = x[i];
					rescuedgrad[i] = grad[i];
					for (j = 0; j < M; j++, p++)
						rescuedhess[p] = hess[p];
				}
				if (1e-4 < lambda)
					lambda = lambda / 10;
			} else {
				/*
				 * else, if it is worse, then recover the last geometry and increase
				 * lambda, saturate lambda with FIRSTLAMBDA
				 */
				for (i = 0, p = 0; i < M; i++) {
					x[i] = rescuedx[i];
					grad[i] = rescuedgrad[i];
					for (j = 0; j < M; j++, p++)
						hess[p] = rescuedhess[p];
				}
				if (lambda < 1.0 / TINY)
					lambda *= 10;
				else
					break;
				if (lambda < FIRSTLAMBDA)
					lambda = FIRSTLAMBDA;
			}

			stop = plugin != null && plugin.isPluginInterrumped();
		}

		// Copy the values back to the input arrays
		for (i = 0, p = 0; i < intervals + 3; i++)
			for (j = 0; j < intervals + 3; j++, p++) {
				cxTargetToSource[i][j] = x[p];
				cxSourceToTarget[i][j] = x[quarterM + p];

				cyTargetToSource[i][j] = x[halfM + p];
				cySourceToTarget[i][j] = x[threeQuarterM + p];
			}

		ProgressBar.skipProgressBar(maxiter - iter);
		return f;
	}

	/**
	 * Energy function to be minimized by the optimizer in the bidirectional case.
	 *
	 * @param c
	 *          Input: Deformation coefficients
	 * @param intervals
	 *          Input: Number of intervals for the deformation
	 * @param grad
	 *          Output: Gradient of the function
	 * @param only_image
	 *          Input: if true, only the image term is considered and not the
	 *          regularization
	 * @param show_error
	 *          Input: if true, an image is shown with the error
	 * @return value of the energy function for these deformation coefficients
	 */
	private double energyFunction(final double[] c, final int intervals, double[] grad, final boolean only_image,
	    final boolean show_error) {
		final int M = c.length / 2;
		final int halfM = M / 2;
		double[] x1 = new double[M];
		double[] auxGrad1 = new double[M];
		double[] auxGrad2 = new double[M];

		for (int i = 0, p = 0; i < halfM; i++, p++) {
			x1[p] = c[i];
			x1[p + halfM] = c[i + M];
		}

		// Source to Target evaluation (Similarity + Landmarks + Regularization)
		// double f = evaluateSimilarity(x1, intervals, auxGrad1, only_image,
		// show_error, false);
		double f = evaluateSimilarityMultiThread(x1, intervals, auxGrad1, only_image, false);

		double[] x2 = new double[M];
		for (int i = halfM, p = 0; i < M; i++, p++) {
			x2[p] = c[i];
			x2[p + halfM] = c[i + M];
		}

		// Target to Source evaluation (Similarity + Landmarks + Regularization)
		// f += evaluateSimilarity(x2, intervals, auxGrad2, only_image, show_error,
		// true);
		f += evaluateSimilarityMultiThread(x2, intervals, auxGrad2, only_image, true);

		// Gradient composition.
		for (int i = 0, p = 0; i < halfM; i++, p++) {
			grad[p] = auxGrad1[i];
			grad[p + halfM] = auxGrad2[i];
			grad[p + M] = auxGrad1[i + halfM];
			grad[p + M + halfM] = auxGrad2[i + halfM];
		}

		double f_consistency = 0;

		// Consistency term
		if (this.consistencyWeight != 0) {
			// Consistency gradient.
			double[] vgradcons = new double[grad.length];

			// f_consistency = evaluateConsistency(intervals, vgradcons);
			f_consistency = evaluateConsistencyMultiThread(intervals, vgradcons);

			// Update gradient.
			for (int i = 0; i < grad.length; i++)
				grad[i] += vgradcons[i];
		}

		return f + f_consistency;

	}

	/**
	 * Evaluate the similarity between the source and the target images but also
	 * the transformation regularization and and landmarks energy term if
	 * necessary. Multi-threading version.
	 *
	 * @param c
	 *          Input: Deformation coefficients
	 * @param intervals
	 *          Input: Number of intervals for the deformation
	 * @param grad
	 *          Output: Gradient of the similarity
	 * @param only_image
	 *          Input: if true, only the image term is considered and not the
	 *          regularization
	 * @param bIsReverse
	 *          Input: flag to determine the transformation direction
	 *          (target-source=FALSE or source-target=TRUE)
	 * @return images similarity value
	 */

	private double evaluateSimilarityMultiThread(final double[] c, final int intervals, double[] grad,
	    final boolean only_image, boolean bIsReverse) {

		// Auxiliary variables for changing from source to target and inversely
		final BSplineModel auxTarget = (!bIsReverse) ? targetModel : sourceModel;
		final BSplineModel auxSource = (!bIsReverse) ? sourceModel : targetModel;

		final ROI2D auxTargetMsk = (!bIsReverse) ? targetMask : sourceMask;
		final ROI2D auxSourceMsk = (!bIsReverse) ? sourceMask : targetMask;

		final List<ROI2DPoint> auxTargetPh = (!bIsReverse) ? targetLandmarks : sourceLandmarks;
		final List<ROI2DPoint> auxSourcePh = (!bIsReverse) ? sourceLandmarks : targetLandmarks;

		final BSplineModel swx = (!bIsReverse) ? swxTargetToSource : swxSourceToTarget;
		final BSplineModel swy = (!bIsReverse) ? swyTargetToSource : swySourceToTarget;

		final double auxFactorWidth = (!bIsReverse) ? this.targetModel.getFactorWidth() : this.sourceFactorWidth;
		final double auxFactorHeight = (!bIsReverse) ? this.targetModel.getFactorHeight() : this.sourceFactorHeight;

		final double P11[][] = (!bIsReverse) ? this.P11_TargetToSource : this.P11_SourceToTarget;
		final double P12[][] = (!bIsReverse) ? this.P12_TargetToSource : this.P12_SourceToTarget;
		final double P22[][] = (!bIsReverse) ? this.P22_TargetToSource : this.P12_SourceToTarget;

		final int auxTargetCurrentWidth = (!bIsReverse) ? this.targetCurrentWidth : this.sourceCurrentWidth;
		final int auxTargetCurrentHeight = (!bIsReverse) ? this.targetCurrentHeight : this.sourceCurrentHeight;

		final int cYdim = intervals + 3;
		final int cXdim = cYdim;
		final int Nk = cYdim * cXdim;
		final int twiceNk = 2 * Nk;

		final double[] vgradreg = new double[grad.length];
		final double[] vgradland = new double[grad.length];

		// Set the transformation coefficients to the interpolator
		swx.setCoefficients(c, cYdim, cXdim, 0);
		swy.setCoefficients(c, cYdim, cXdim, Nk);

		// Initialize gradient
		for (int k = 0; k < twiceNk; k++)
			vgradreg[k] = vgradland[k] = grad[k] = 0.0F;

		// Estimate the similarity and gradient between both images
		double imageSimilarity = 0.0;

		// final int Ydim = auxTarget.getCurrentHeight();
		// final int Xdim = auxTarget.getCurrentWidth();

		// Image similarity calculated in a concurrent way
		if (imageWeight != 0) {
			// Check the number of processors in the computer
			final int nproc = Runtime.getRuntime().availableProcessors();

			// We will use threads to calculate the similarity of the different
			// parts of the target and source image
			int block_height = auxTargetCurrentHeight / nproc;
			if (auxTargetCurrentHeight % 2 != 0)
				block_height++;

			// We use as many threads as processors
			final int nThreads = nproc;

			Thread[] threads = new Thread[nThreads];
			Rectangle[] rects = new Rectangle[nThreads];

			// Every thread will provide the corresponding similarity value and
			// gradient
			final double[][] grad_thread = new double[nThreads][grad.length];
			// Result array:
			// First result is the partial image similarity and second the number of
			// pixels
			final double[][] result = new double[nThreads][2];
			// Number of processed pixels (taking into account the masks)
			int n = 0;

			for (int i = 0; i < nThreads; i++) {
				// Last block goes to the end of the window
				int y_start = i * block_height;
				if (nThreads - 1 == i)
					block_height = auxTargetCurrentHeight - i * block_height;

				// Corresponding rectangle
				rects[i] = new Rectangle(0, y_start, auxTargetCurrentWidth, block_height);

				// Create threads and start them.
				threads[i] = new Thread(new EvaluateSimilarityTile(auxTarget, auxSource, auxTargetMsk, auxSourceMsk, swx, swy,
				    auxFactorWidth, auxFactorHeight, intervals, grad_thread[i], result[i], rects[i]));
				threads[i].start();
			}

			// Wait for the threads to finish
			for (int i = 0; i < nThreads; i++) {
				try {
					threads[i].join();
					threads[i] = null;
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}

			// Accumulate results
			for (int i = 0; i < nThreads; i++) {
				imageSimilarity += result[i][0];
				n += result[i][1];
			}

			// Average image similarity
			imageSimilarity /= n;
			// Average gradients
			for (int i = 0; i < nThreads; i++) {
				for (int j = 0; j < grad.length; j++)
					grad[j] += (grad_thread[i][j] / n);
			}

		}

		// Compute regularization term
		// ..............................................
		double regularization = 0.0;
		if (!only_image) {
			for (int i = 0; i < Nk; i++)
				for (int j = 0; j < Nk; j++) {
					regularization += c[i] * P11[i][j] * c[j] + // c1^t P11 c1
					    c[Nk + i] * P22[i][j] * c[Nk + j] + // c2^t P22 c2
					    c[i] * P12[i][j] * c[Nk + j];// c1^t P12 c2
					vgradreg[i] += 2 * P11[i][j] * c[j]; // 2 P11 c1
					vgradreg[Nk + i] += 2 * P22[i][j] * c[Nk + j]; // 2 P22 c2
					vgradreg[i] += P12[i][j] * c[Nk + j]; // P12 c2
					vgradreg[Nk + i] += P12[j][i] * c[j]; // P12^t c1
				}
			regularization *= 1.0 / (auxTargetCurrentHeight * auxTargetCurrentWidth);
			for (int k = 0; k < twiceNk; k++)
				vgradreg[k] *= 1.0 / (auxTargetCurrentHeight * auxTargetCurrentWidth);
		}

		// Compute landmark error and derivative ...............................
		// Get the list of landmarks
		double landmarkError = 0.0;
		int K = 0;
		if (auxTargetPh != null)
			K = auxTargetPh.size();

		if (landmarkWeight != 0) {
			List<Point2D> sourceVector = new ArrayList<>();
			if (auxSourcePh != null) {
				for (ROI2DPoint roi : auxSourcePh) {
					sourceVector.add(roi.getPoint());
				}
			}

			List<Point2D> targetVector = new ArrayList<>();
			if (auxTargetPh != null) {
				for (ROI2DPoint roi : auxTargetPh) {
					targetVector.add(roi.getPoint());
				}
			}

			for (int kp = 0; kp < K; kp++) {
				// Get the landmark coordinate in the target image
				final Point2D sourcePoint = sourceVector.get(kp);
				final Point2D targetPoint = targetVector.get(kp);
				double u = auxFactorWidth * targetPoint.getX();
				double v = auxFactorHeight * targetPoint.getY();

				// Express it in "spline" units
				double tu = (double) (u * intervals) / (double) (auxTargetCurrentWidth - 1) + 1.0F;
				double tv = (double) (v * intervals) / (double) (auxTargetCurrentHeight - 1) + 1.0F;

				// Transform this coordinate to the source image
				swx.prepareForInterpolation(tu, tv, false);
				double x = swx.interpolateI();
				swy.prepareForInterpolation(tu, tv, false);
				double y = swy.interpolateI();

				// Substract the result from the residual
				double dx = auxFactorWidth * sourcePoint.getX() - x;
				double dy = auxFactorHeight * sourcePoint.getY() - y;

				// Add to landmark error
				landmarkError += dx * dx + dy * dy;

				// Compute the derivative with respect to all the c coefficients
				for (int l = 0; l < 4; l++)
					for (int m = 0; m < 4; m++) {
						if (swx.yIndex[l] == -1 || swx.xIndex[m] == -1)
							continue;
						int k = swx.yIndex[l] * cYdim + swx.xIndex[m];

						// There's also a multiplication by 2 that I will do later
						// Derivative related to X deformation
						vgradland[k] -= dx * swx.getWeightI(l, m);

						// Derivative related to Y deformation
						vgradland[k + Nk] -= dy * swy.getWeightI(l, m);
					}
			}
		}

		if (K != 0) {
			landmarkError *= landmarkWeight / K;
			double aux = 2.0 * landmarkWeight / K;
			// This is the 2 coming from the derivative
			// computation that I would do at the end
			for (int k = 0; k < twiceNk; k++)
				vgradland[k] *= aux;
		}
		if (only_image)
			landmarkError = 0;

		// Finish computations
		// .............................................................
		// Add all gradient terms (similarity + regularization + landmarks)
		for (int k = 0; k < twiceNk; k++)
			grad[k] += vgradreg[k] + vgradland[k];

		if (showMarquardtOptim) {
			String s = bIsReverse ? new String("(t-s)") : new String("(s-t)");
			if (imageWeight != 0) {
				System.out.println("    Image          error " + s + ": " + imageSimilarity);
				if (bIsReverse)
					this.partialInverseSimilarityError = imageSimilarity;
				else
					this.partialDirectSimilarityError = imageSimilarity;

			}
			if (landmarkWeight != 0) {
				System.out.println("    Landmark       error " + s + ": " + landmarkError);
				if (bIsReverse)
					this.partialInverseLandmarkError = landmarkError;
				else
					this.partialDirectLandmarkError = landmarkError;
			}
			if (divWeight != 0 || curlWeight != 0) {
				System.out.println("    Regularization error " + s + ": " + regularization);
				if (bIsReverse)
					this.partialInverseRegularizationError = regularization;
				else
					this.partialDirectRegularizationError = regularization;

			}
		}
		return imageSimilarity + landmarkError + regularization;
	}

	/**
	 * Class to run concurrent similarity evaluation
	 * 
	 */
	private class EvaluateSimilarityTile implements Runnable {
		// Fields
		/** current target image */
		final BSplineModel auxTarget;
		/** current source image */
		final BSplineModel auxSource;
		/** target mask */
		final ROI2D auxTargetMsk;
		/** source mask */
		final ROI2D auxSourceMsk;
		/** B-spline deformation in x */
		final BSplineModel swx;
		/** B-spline deformation in y */
		final BSplineModel swy;
		/** factor width */
		final double auxFactorWidth;
		/** factor height */
		final double auxFactorHeight;
		/** number of intervals between B-spline coefficients */
		final int intervals;
		/** similarity gradient */
		final double[] grad;
		/**
		 * evaluation results: image similarity value for the current rectangle and
		 * number of pixels that have been evaluated
		 */
		final double[] result;
		/** rectangle containing the area of the image to be evaluated */
		final Rectangle rect;

		/**
		 * Evaluate similarity tile constructor
		 * 
		 * @param auxTarget
		 *          current target image
		 * @param auxSource
		 *          current source image
		 * @param auxTargetMsk
		 *          target mask
		 * @param auxSourceMsk
		 *          source mask
		 * @param swx
		 *          B-spline deformation in x
		 * @param swy
		 *          B-spline deformation in y
		 * @param auxFactorWidth
		 *          factor width
		 * @param auxFactorHeight
		 *          factor height
		 * @param intervals
		 *          number of intervals between B-spline coefficients
		 * @param grad
		 *          similarity gradient (output)
		 * @param result
		 *          output results: image similarity value for the current rectangle
		 *          and number of pixels that have been evaluated
		 * @param rect
		 *          rectangle containing the area of the image to be evaluated
		 */
		EvaluateSimilarityTile(final BSplineModel auxTarget, final BSplineModel auxSource, final ROI2D auxTargetMsk,
		    final ROI2D auxSourceMsk, final BSplineModel swx, final BSplineModel swy, final double auxFactorWidth,
		    final double auxFactorHeight, final int intervals, final double[] grad, final double[] result,
		    final Rectangle rect) {
			this.auxTarget = auxTarget;
			this.auxSource = auxSource;

			this.auxTargetMsk = auxTargetMsk;
			this.auxSourceMsk = auxSourceMsk;

			this.swx = swx;
			this.swy = swy;

			this.auxFactorWidth = auxFactorWidth;
			this.auxFactorHeight = auxFactorHeight;

			this.intervals = intervals;

			this.grad = grad;

			this.result = result;

			this.rect = rect;
		}

		// ------------------------------------------------------------------
		/**
		 * Run method to evaluate the similarity of source and target images. Only
		 * the part defined by the rectangle will be evaluated.
		 */

		public void run() {
			final int cYdim = intervals + 3;
			final int cXdim = cYdim;
			final int Nk = cYdim * cXdim;
			final int twiceNk = 2 * Nk;

			double imageSimilarity = 0.0;

			// The rectangle marks the area of the image to be treated.
			int uv = rect.y * rect.width + rect.x;
			final int Ydim = rect.y + rect.height;
			final int Xdim = rect.x + rect.width;

			// Loop over all points in the source image (rectangle)
			int n = 0;

			final double[] I1D = new double[2]; // Space for the first derivatives of
			                                    // I1

			final double[] targetCurrentImage = auxTarget.getCurrentImage();

			for (int v = rect.y; v < Ydim; v++) {
				for (int u = rect.x; u < Xdim; u++, uv++) {
					// Compute image term
					// .....................................................

					// Check if this point is in the target mask
					if (auxTargetMsk == null || auxTargetMsk.contains(u / auxFactorWidth, v / auxFactorHeight)) {
						// Compute value in the source image
						final double I2 = targetCurrentImage[uv];

						// Compute the position of this point in the target
						double x = swx.precomputedInterpolateI(u, v);
						double y = swy.precomputedInterpolateI(u, v);

						// Check if this point is in the source mask
						if (auxSourceMsk == null || auxSourceMsk.contains(x / auxFactorWidth, y / auxFactorHeight)) {
							// Compute the value of the target at that point
							final double I1 = auxSource.prepareForInterpolationAndInterpolateIAndD(x, y, I1D, false, PYRAMID);

							final double I1dx = I1D[0], I1dy = I1D[1];

							final double error = I2 - I1;
							final double error2 = error * error;
							imageSimilarity += error2;

							// Compute the derivative with respect to all the c coefficients
							// Cost of the derivatives = 16*(3 mults + 2 sums)
							// Current cost = 359 mults + 346 sums
							for (int l = 0; l < 4; l++)
								for (int m = 0; m < 4; m++) {
									if (swx.prec_yIndex[v][l] == -1 || swx.prec_xIndex[u][m] == -1)
										continue;

									// Note: It's the same to take the indexes and weightI from
									// swx than from swy
									double weightI = swx.precomputedGetWeightI(l, m, u, v);

									int k = swx.prec_yIndex[v][l] * cYdim + swx.prec_xIndex[u][m];

									// Compute partial result
									// There's also a multiplication by 2 that I will
									// do later
									double aux = -error * weightI;

									// Derivative related to X deformation
									grad[k] += aux * I1dx;

									// Derivative related to Y deformation
									grad[k + Nk] += aux * I1dy;
								}
							n++; // Another point has been successfully evaluated
						}
					}
				}
			}

			// Average the image related terms (now i do the 1/n outside)
			if (n != 0) {
				imageSimilarity *= imageWeight;
				double aux = imageWeight * 2.0; // This is the 2 coming from the
				// derivative that I would do later
				for (int k = 0; k < twiceNk; k++)
					grad[k] *= aux;
			} else
				imageSimilarity = 1 / FLT_EPSILON;

			// Set result (image similarity value for the current rectangle
			// and number of pixels that have been evaluated)
			this.result[0] = imageSimilarity;
			this.result[1] = n;

		} // end run method

	}

	/**
	 * In this function the system (H+lambda*Diag(H))*update=gradient is solved
	 * for update. H is the hessian of the function f, gradient is the gradient of
	 * the function f, Diag(H) is a matrix with the diagonal of H.
	 */
	private void Marquardt_it(double[] x, boolean[] optimize, double[] gradient, double[] Hessian, double lambda) {
		// final double TINY = FLT_EPSILON;
		final int M = x.length;

		// Find the threshold for the most important components
		double[] sortedgradient = new double[M];
		for (int i = 0; i < M; i++)
			sortedgradient[i] = Math.abs(gradient[i]);
		Arrays.sort(sortedgradient);

		double largestGradient = sortedgradient[M - 1];

		// We set the threshold gradient at 9% of the largest value.
		double gradient_th = 0.09 * largestGradient;

		// We count the number of values over the threshold.
		int Mused = 0;
		for (int i = 0; i < M; i++)
			if (sortedgradient[i] >= gradient_th)
				Mused++;

		double[][] u = new double[Mused][Mused];
		// double [][] v = null; //new double [Mused][Mused];
		// double [] w = null; //new double [Mused];
		double[] g = new double[Mused];
		double[] update = new double[Mused];
		boolean[] optimizep = new boolean[M];

		System.arraycopy(optimize, 0, optimizep, 0, M);

		lambda += 1.0F;

		int m = 0, i;

		// Take the Mused components with big gradients
		for (i = 0; i < M; i++)
			if (optimizep[i] && Math.abs(gradient[i]) >= gradient_th) {
				m++;
				if (m == Mused)
					break;
			} else
				optimizep[i] = false;
		// Set the rest to 0
		for (i = i + 1; i < M; i++)
			optimizep[i] = false;

		// Gradient descent
		// for (int i=0; i<M; i++) if (optimizep[i]) x[i]-=0.01*gradient[i];
		// if (true) return;

		/*
		 * u will be a copy of the Hessian where we take only those components
		 * corresponding to variables being optimized
		 */
		int kr = 0, iw = 0;
		for (int ir = 0; ir < M; kr = kr + M, ir++) {
			if (optimizep[ir]) {
				int jw = 0;
				for (int jr = 0; jr < M; jr++)
					if (optimizep[jr])
						u[iw][jw++] = Hessian[kr + jr];
				g[iw] = gradient[ir];
				u[iw][iw] *= lambda;
				iw++;
			}
		}

		// Solve he equation system
		/* SVD u=u*w*v^t */
		update = MathTools.linearLeastSquares(u, g);
		if (update == null) {
			System.out.println("Error when calculating linear least square solution...");
			return;
		}

		/* x = x - update */
		kr = 0;
		for (int kw = 0; kw < M; kw++)
			if (optimizep[kw])
				x[kw] -= update[kr++];

	}

	/**
	 * Method to update both current outputs (source-target and target-source).
	 *
	 * @param c
	 *          B-spline coefficients
	 * @param intervals
	 *          number of intervals in the deformation
	 */
	private void updateOutputs(final double[] c, int intervals) {
		final int M = c.length / 2;
		final int halfM = M / 2;
		double[] x1 = new double[M];

		for (int i = 0, p = 0; i < halfM; i++, p++) {
			x1[p] = c[i];
			x1[p + halfM] = c[i + M];
		}

		double[] x2 = new double[M];
		for (int i = halfM, p = 0; i < M; i++, p++) {
			x2[p] = c[i];
			x2[p + halfM] = c[i + M];
		}
		// Updates.
		updateCurrentOutput(x1, intervals, false);
		updateCurrentOutput(x2, intervals, true);
	}

	/**
	 * Method to update a current output (multi-thread).
	 *
	 * @param c
	 *          B-spline coefficients
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param bIsReverse
	 *          flag to decide the deformation direction (source-target,
	 *          target-source)
	 */
	private void updateCurrentOutput(final double[] c, int intervals, boolean bIsReverse) {
		// Set the coefficients to an interpolator
		int cYdim = intervals + 3;
		int cXdim = cYdim;
		int Nk = cYdim * cXdim;

		BSplineModel auxTarget = targetModel;
		BSplineModel auxSource = sourceModel;
		ROI2D auxTargetMsk = targetMask;
		ROI2D auxSourceMsk = sourceMask;
		BSplineModel swx = swxTargetToSource;
		BSplineModel swy = swyTargetToSource;
		int auxTargetWidth = this.targetWidth;
		int auxTargetHeight = this.targetHeight;
		int auxTargetCurrentWidth = this.targetCurrentWidth;
		int auxTargetCurrentHeight = this.targetCurrentHeight;
		int auxSourceWidth = this.sourceWidth;
		int auxSourceHeight = this.sourceHeight;
		Sequence auxSourceSeq = this.sourceSeq;
		Sequence outputSeq = this.outputSeq1;
		double auxFactorWidth = this.targetFactorWidth;
		double auxFactorHeight = this.targetFactorHeight;
		double subFactorWidth = targetModel.isSubOutput() ? (targetModel.getWidth() / targetModel.getSubWidth()) : 1;
		double subFactorHeight = targetModel.isSubOutput() ? (targetModel.getHeight() / targetModel.getSubHeight()) : 1;

		// Change if necessary
		if (bIsReverse) {
			auxTarget = sourceModel;
			auxSource = targetModel;
			auxTargetMsk = sourceMask;
			auxSourceMsk = targetMask;
			swx = swxSourceToTarget;
			swy = swySourceToTarget;
			auxTargetWidth = this.sourceWidth;
			auxTargetHeight = this.sourceHeight;
			auxTargetCurrentWidth = sourceCurrentWidth;
			auxTargetCurrentHeight = sourceCurrentHeight;
			auxSourceWidth = this.targetWidth;
			auxSourceHeight = this.targetHeight;
			auxSourceSeq = this.targetSeq;
			outputSeq = this.outputSeq2;
			auxFactorWidth = this.sourceFactorWidth;
			auxFactorHeight = this.sourceFactorHeight;
			subFactorWidth = sourceModel.isSubOutput() ? (sourceModel.getWidth() / sourceModel.getSubWidth()) : 1;
			subFactorHeight = sourceModel.isSubOutput() ? (sourceModel.getHeight() / sourceModel.getSubHeight()) : 1;
		}

		swx.setCoefficients(c, cYdim, cXdim, 0);
		swy.setCoefficients(c, cYdim, cXdim, Nk);

		// Compute the deformed image
		outputSeq.beginUpdate();
		IcyBufferedImage ibi = outputSeq.getImage(0, 0);

		int uv = 0;

		// Check the number of processors in the computer
		int nproc = Runtime.getRuntime().availableProcessors();

		// We will use threads to display parts of the output image
		int block_height = auxTargetHeight / ((int) subFactorHeight * nproc);
		if (auxTargetHeight % 2 != 0)
			block_height++;

		int nThreads = nproc; /*
		                       * (nproc > 1) ? (nproc / 2) : 1; if
		                       * (this.accurate_mode == MainDialog.MONO_MODE)
		                       * nThreads *= 2;
		                       */

		Thread[] threads = new Thread[nThreads];
		Rectangle[] rects = new Rectangle[nThreads];
		IcyBufferedImage[] ibi_tile = new IcyBufferedImage[nThreads];
		for (int i = 0; i < nThreads; i++) {
			// Last block goes to the end of the window
			int x_start = i * block_height;
			if (nThreads - 1 == i)
				block_height = auxTargetHeight / (int) subFactorHeight - i * block_height;
			/*
			 * IJ.log("block height " + block_height); IJ.log("Update : 0 " + " "+
			 * x_start +" " + (auxTargetWidth / (int)subFactorWidth) + " " +
			 * block_height); IJ.log("auxFactorWidth = " + auxFactorWidth +
			 * " auxFactorHeight = " + auxFactorHeight);
			 */
			rects[i] = new Rectangle(0, x_start, auxTargetWidth / (int) subFactorWidth, block_height);

			ibi_tile[i] = new IcyBufferedImage(rects[i].width, rects[i].height, 1, DataType.FLOAT);

			threads[i] = new Thread(new OutputTileMaker(swx, swy, auxSource, auxTarget, auxSourceMsk, auxTargetMsk,
			    auxFactorWidth * subFactorWidth, auxFactorHeight * subFactorHeight, auxTargetCurrentHeight,
			    auxTargetCurrentWidth, rects[i], ibi_tile[i]));
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

		for (int i = 0; i < nThreads; i++) {
			ibi.copyData(ibi_tile[i], null, new Point(rects[i].x, rects[i].y));
			ibi_tile[i] = null;
			rects[i] = null;
		}

		outputSeq.setImage(0, 0, ibi);
		outputSeq.endUpdate();

		// Draw the grid on the target image ...............................

		// We take the values from the original image
		auxTargetHeight = bIsReverse ? this.originalSourceIBI.getHeight() : this.originalTargetIBI.getHeight();
		auxTargetWidth = bIsReverse ? this.originalSourceIBI.getWidth() : this.originalTargetIBI.getWidth();
		auxSourceHeight = bIsReverse ? this.originalTargetIBI.getHeight() : this.originalSourceIBI.getHeight();
		auxSourceWidth = bIsReverse ? this.originalTargetIBI.getWidth() : this.originalSourceIBI.getWidth();

		auxFactorWidth = (double) auxTargetCurrentWidth / auxTargetWidth;
		auxFactorHeight = (double) auxTargetCurrentHeight / auxTargetHeight;

		// Some initialization
		int stepv = Math.min(Math.max(10, auxTargetHeight / 15), 60);
		int stepu = Math.min(Math.max(10, auxTargetWidth / 15), 60);
		final double transformedImage[][] = new double[auxSourceHeight][auxSourceWidth];
		double grid_colour = -1e-10;
		uv = 0;
		for (int v = 0; v < auxSourceHeight; v++)
			for (int u = 0; u < auxSourceWidth; u++, uv++) {
				transformedImage[v][u] = auxSource.getOriginalImage()[uv];
				if (transformedImage[v][u] > grid_colour)
					grid_colour = transformedImage[v][u];
			}

		// Draw grid
		for (int v = 0; v < auxTargetHeight + stepv; v += stepv)
			for (int u = 0; u < auxTargetWidth + stepu; u += stepu) {
				// down_u/v are the coordinates in the current image
				double down_u = u * auxFactorWidth;
				double down_v = v * auxFactorHeight;

				// tv,tu are the corresponding coordinates in the interpolator
				final double tv = (double) (down_v * intervals) / (double) (auxTargetCurrentHeight - 1) + 1.0F;
				final double tu = (double) (down_u * intervals) / (double) (auxTargetCurrentWidth - 1) + 1.0F;

				// x,y are the coordinates after the transformation
				swx.prepareForInterpolation(tu, tv, ORIGINAL);
				double x = swx.interpolateI();
				swy.prepareForInterpolation(tu, tv, ORIGINAL);
				double y = swy.interpolateI();

				// up_x, up_y are the transformed coordinates in the original image
				double up_x = x / auxFactorWidth;
				double up_y = y / auxFactorHeight;

				// Draw horizontal line
				int uh = u + stepu;
				if (uh < auxTargetWidth + stepu) {
					final double down_uh = uh * auxFactorWidth;
					final double tuh = (double) (down_uh * intervals) / (double) (auxTargetCurrentWidth - 1) + 1.0F;
					swx.prepareForInterpolation(tuh, tv, ORIGINAL);
					final double xh = swx.interpolateI();
					swy.prepareForInterpolation(tuh, tv, ORIGINAL);
					final double yh = swy.interpolateI();
					final double up_xh = xh / auxFactorWidth;
					final double up_yh = yh / auxFactorHeight;
					MiscTools.drawLine(transformedImage, (int) Math.round(up_x), (int) Math.round(up_y), (int) Math.round(up_xh),
					    (int) Math.round(up_yh), grid_colour);
				}

				// Draw vertical line
				int vv = v + stepv;
				if (vv < auxTargetHeight + stepv) {
					double down_vv = vv * auxFactorHeight;
					final double tvv = (double) (down_vv * intervals) / (double) (auxTargetCurrentHeight - 1) + 1.0F;
					swx.prepareForInterpolation(tu, tvv, ORIGINAL);
					double xv = swx.interpolateI();
					swy.prepareForInterpolation(tu, tvv, ORIGINAL);
					double yv = swy.interpolateI();
					double up_xv = xv / auxFactorWidth;
					double up_yv = yv / auxFactorHeight;
					MiscTools.drawLine(transformedImage, (int) Math.round(up_x), (int) Math.round(up_y), (int) Math.round(up_xv),
					    (int) Math.round(up_yv), grid_colour);
				}
			}

		// Update the target image plus
		IcyBufferedImage ibig = new IcyBufferedImage(auxSourceWidth, auxSourceHeight, 1, DataType.DOUBLE);
		double[] ibigData = ibig.getDataXYAsDouble(0);
		for (int v = 0; v < auxSourceHeight; v++)
			for (int u = 0; u < auxSourceWidth; u++)
				ibigData[u + v * auxSourceWidth] = transformedImage[v][u];
		// ibig.putPixelValue(u, v, transformedImage[v][u]);

		ibig.dataChanged();

		auxSourceSeq.beginUpdate();
		auxSourceSeq.setImage(0, 0, ibig);
		auxSourceSeq.endUpdate();
	}

	/* ------------------------------------------------------------------------ */
	/**
	 * Class to run concurrent tile windows updaters (for intermediate results)
	 * 
	 */
	private class OutputTileMaker implements Runnable {
		final BSplineModel swx;
		final BSplineModel swy;
		final BSplineModel auxSource;
		final BSplineModel auxTarget;
		final ROI2D auxTargetMsk;
		final ROI2D auxSourceMsk;
		final double auxFactorWidth;
		final double auxFactorHeight;
		final int auxTargetCurrentHeight;
		final int auxTargetCurrentWidth;
		final Rectangle rect;
		final private IcyBufferedImage fp;

		// ------------------------------------------------------------------
		/**
		 * Output tile maker constructor
		 * 
		 * @param swx
		 *          B-spline interpolator for transformation in x-
		 * @param swy
		 *          B-spline interpolator for transformation in y-
		 * @param auxSource
		 *          source image
		 * @param auxTarget
		 *          target image
		 * @param auxSourceMsk
		 *          source mask
		 * @param auxTargetMsk
		 *          target mask
		 * @param auxFactorWidth
		 *          width factor
		 * @param auxFactorHeight
		 *          height factor
		 * @param auxTargetCurrentHeight
		 *          current target height
		 * @param auxTargetCurrentWidth
		 *          current target width
		 * @param rect
		 *          retangle with the coordinates of the output image to be updated
		 * @param fp
		 *          processor to be updated
		 */
		OutputTileMaker(final BSplineModel swx, final BSplineModel swy, final BSplineModel auxSource,
		    final BSplineModel auxTarget, final ROI2D auxSourceMsk, final ROI2D auxTargetMsk, final double auxFactorWidth,
		    final double auxFactorHeight, final int auxTargetCurrentHeight, final int auxTargetCurrentWidth,
		    final Rectangle rect, final IcyBufferedImage fp) {
			this.swx = swx;
			this.swy = swy;
			this.auxSource = auxSource;
			this.auxTarget = auxTarget;
			this.auxTargetMsk = auxTargetMsk;
			this.auxSourceMsk = auxSourceMsk;
			this.auxFactorWidth = auxFactorWidth;
			this.auxFactorHeight = auxFactorHeight;
			this.auxTargetCurrentWidth = auxTargetCurrentWidth;
			this.auxTargetCurrentHeight = auxTargetCurrentHeight;
			this.rect = rect;
			this.fp = fp;
		}

		// ------------------------------------------------------------------
		/**
		 * Run method to update the intermediate window. Only the part defined by
		 * the rectangle will be updated (in this thread).
		 */
		public void run() {
			int uv = rect.y * rect.width + rect.x;
			int auxTargetHeight = rect.y + rect.height;
			int auxTargetWidth = rect.x + rect.width;

			// Subsampling (output) factors
			// Note: we get them from original image, since they will be used in the
			// masks and the
			// masks store the information relative to the original sizes (without
			// scaling).
			final int subFactorT = this.auxTarget.getOriginalImageWidth() / this.auxTarget.getSubWidth();
			final int subFactorS = this.auxSource.getOriginalImageWidth() / this.auxSource.getSubWidth();

			final boolean fromSubT = (auxTarget.isSubOutput());
			final boolean fromSubS = (auxSource.isSubOutput());

			final double[] tImage = fromSubT ? auxTarget.getSubImage() : auxTarget.getImage();
			final float[] f_array = fp.getDataXYAsFloat(0);

			for (int v_rect = 0, v = rect.y; v < auxTargetHeight; v++, v_rect++) {
				final int v_offset = v_rect * rect.width;

				for (int u_rect = 0, u = rect.x; u < auxTargetWidth; u++, uv++, u_rect++) {
					if (auxTargetMsk == null || auxTargetMsk.contains(u * subFactorT, v * subFactorT)) {
						double down_u = u * auxFactorWidth;
						double down_v = v * auxFactorHeight;
						final double tv = (double) (down_v * intervals) / (double) (auxTargetCurrentHeight - 1) + 1.0F;
						final double tu = (double) (down_u * intervals) / (double) (auxTargetCurrentWidth - 1) + 1.0F;
						double x = swx.prepareForInterpolationAndInterpolateI(tu, tv, fromSubT, ORIGINAL);
						double y = swy.prepareForInterpolationAndInterpolateI(tu, tv, fromSubT, ORIGINAL);
						double up_x = x / auxFactorWidth;
						double up_y = y / auxFactorHeight;
						if (auxSourceMsk == null || auxSourceMsk.contains(up_x * subFactorS, up_y * subFactorS)) {
							double sourceValue = auxSource.prepareForInterpolationAndInterpolateI(up_x, up_y, fromSubS, ORIGINAL);
							// fp.putPixelValue(u_rect, v_rect, tImage[uv] - sourceValue);
							f_array[u_rect + v_offset] = (float) (tImage[uv] - sourceValue);
						} else
						  // fp.putPixelValue(u_rect, v_rect, 0);
						  f_array[u_rect + v_offset] = 0;
					} else
					  // fp.putPixelValue(u_rect, v_rect, 0);
					  f_array[u_rect + v_offset] = 0;
				}
			}

		} // end run
	}

	/**
	 * Calculate the geometric error between the source-target and target-source
	 * deformations. The corresponding coefficients are assumed to be at
	 * swxTargetToSource, swyTargetToSource, swxSourceToTarget and
	 * swySourceToTarget. Multi-thread version.
	 *
	 * @param intervals
	 *          Input: Number of intervals for the deformation
	 * @param grad
	 *          Output: Gradient of the function
	 * 
	 * @return geometric error between the source-target and target-source
	 *         deformations.
	 */

	private double evaluateConsistencyMultiThread(final int intervals, double[] grad) {
		// Consistency values
		double f_direct = 0.0;
		double f_inverse = 0.0;

		// Check the number of processors in the computer
		final int nproc = Runtime.getRuntime().availableProcessors();

		// We will use threads to calculate the similarity of the different
		// parts of the target and source image
		int block_height_target = this.targetCurrentHeight / nproc;
		if (this.targetCurrentHeight % 2 != 0)
			block_height_target++;

		int block_height_source = this.sourceCurrentHeight / nproc;
		if (this.sourceCurrentHeight % 2 != 0)
			block_height_source++;

		// We use as many threads as processors
		final int nThreads = nproc;

		Thread[] threads = new Thread[nThreads];
		Rectangle[] rect_target = new Rectangle[nThreads];
		Rectangle[] rect_source = new Rectangle[nThreads];

		// Every thread will provide the corresponding consistency value and
		// gradient
		final double[][] grad_direct = new double[nThreads][grad.length];
		final double[][] grad_inverse = new double[nThreads][grad.length];
		// Result array:
		// First result is the direct partial consistency, second the number of
		// pixels (direct),
		// third the inverse partical consistency and fourth the number of pixels
		// (inverse)
		final double[][] result = new double[nThreads][4];
		// Number of processed pixels (taking into account the masks)
		int n_direct = 0;
		int n_inverse = 0;

		for (int i = 0; i < nThreads; i++) {
			// Last block goes to the end of the window
			int y_start_target = i * block_height_target;
			int y_start_source = i * block_height_source;
			if (nThreads - 1 == i) {
				block_height_target = this.targetCurrentHeight - i * block_height_target;
				block_height_source = this.sourceCurrentHeight - i * block_height_source;
			}

			// Corresponding rectangles
			rect_target[i] = new Rectangle(0, y_start_target, this.targetCurrentHeight, block_height_target);
			rect_source[i] = new Rectangle(0, y_start_source, this.sourceCurrentHeight, block_height_source);

			// Create threads and start them.
			threads[i] = new Thread(new EvaluateConsistencyTile(this, grad_direct[i], grad_inverse[i], result[i],
			    rect_target[i], rect_source[i]));
			threads[i].start();
		}

		// Wait for the threads to finish
		for (int i = 0; i < nThreads; i++) {
			try {
				threads[i].join();
				threads[i] = null;
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		// Accumulate results
		for (int i = 0; i < nThreads; i++) {
			f_direct += result[i][0];
			n_direct += result[i][1];
			f_inverse += result[i][2];
			n_inverse += result[i][3];
		}

		// Average consistency
		f_direct /= n_direct;
		f_inverse /= n_inverse;

		// Average and combine gradients
		for (int i = 0; i < nThreads; i++) {
			for (int j = 0; j < grad.length; j++)
				grad[j] += (grad_direct[i][j] / n_direct) + (grad_inverse[i][j] / n_inverse);
		}

		this.partialDirectConsitencyError = this.consistencyWeight * f_direct;
		this.partialInverseConsitencyError = this.consistencyWeight * f_inverse;

		double consistencyDirectError = (n_direct == 0) ? 1.0 / FLT_EPSILON : (this.consistencyWeight * f_direct);
		double consistencyInverseError = (n_inverse == 0) ? 1.0 / FLT_EPSILON : (this.consistencyWeight * f_inverse);

		if (showMarquardtOptim) {
			System.out.println("    Consistency Error (s-t): " + consistencyDirectError);
			System.out.println("    Consistency Error (t-s): " + consistencyInverseError);
		}

		if (n_direct == 0 || n_inverse == 0)
			return 1 / FLT_EPSILON;
		return (this.consistencyWeight * (f_direct + f_inverse));
	}

	/**
	 * Class to run concurrent consistency evaluation
	 * 
	 */
	private class EvaluateConsistencyTile implements Runnable {

		/** transformation object, it contains all the registration information */
		final Transformation transf;
		/** output direct gradient array */
		final double[] grad_direct;
		/** output inverse gradient array */
		final double[] grad_inverse;
		/**
		 * direct and inverse consistency error values and number of pixels (f_dir,
		 * n_dir, f_inv, n_inv)
		 */
		final double[] result;
		/** rectangle marking the target area to be evaluated */
		final Rectangle rect_target;
		/** rectangle marking the source area to be evaluated */
		final Rectangle rect_source;

		/**
		 * Evaluate consistency tile constructor
		 * 
		 * @param transf
		 *          transformation object, it contains all the registration
		 *          information
		 * @param grad_direct
		 *          output direct gradient array
		 * @param grad_inverse
		 *          output inverse gradient array
		 * @param result
		 *          direct and inverse consistency error values and number of pixels
		 *          (f_dir, n_dir, f_inv, n_inv)
		 * @param rect_target
		 *          rectangle marking the target area to be evaluated
		 * @param rect_source
		 *          rectangle marking the source area to be evaluated
		 */
		EvaluateConsistencyTile(Transformation transf, double[] grad_direct, double[] grad_inverse, double[] result,
		    Rectangle rect_target, Rectangle rect_source) {
			this.transf = transf;

			this.grad_direct = grad_direct;
			this.grad_inverse = grad_inverse;

			this.result = result;

			this.rect_target = rect_target;
			this.rect_source = rect_source;
		}

		/**
		 * Run method to evaluate the transformation consistency between source and
		 * target images. Only the part defined by the rectangle will be evaluated.
		 */
		public void run() {
			final int cYdim = this.transf.intervals + 3;
			final int cXdim = cYdim;
			final int Nk = cYdim * cXdim;
			final int twiceNk = 2 * Nk;

			// Initialize gradient
			for (int k = 0; k < this.grad_direct.length; k++)
				this.grad_direct[k] = 0.0F;

			// The target rectangle marks the area of the target image to be treated.
			final int YdimT = rect_target.y + rect_target.height;
			final int XdimT = rect_target.x + rect_target.width;

			// Compute the deformation
			// Set these coefficients to an interpolator
			final BSplineModel swx_direct = this.transf.swxTargetToSource;
			final BSplineModel swy_direct = this.transf.swyTargetToSource;

			final BSplineModel swx_inverse = this.transf.swxSourceToTarget;
			final BSplineModel swy_inverse = this.transf.swySourceToTarget;

			// *********** Compute the geometric error and gradient (DIRECT)
			// ***********
			double f_direct = 0;
			int n_direct = 0;

			for (int v = rect_target.y; v < YdimT; v++)
				for (int u = rect_target.x; u < XdimT; u++) {
					// Check if this point is in the target mask
					if (this.transf.targetMask == null || this.transf.targetMask.contains(u / this.transf.targetFactorWidth,
					    v / this.transf.targetFactorHeight)) {

						final int x = (int) Math.round(swx_direct.precomputedInterpolateI(u, v));
						final int y = (int) Math.round(swy_direct.precomputedInterpolateI(u, v));

						if (x >= 0 && x < this.transf.sourceCurrentWidth && y >= 0 && y < this.transf.sourceCurrentHeight) {
							final double x2 = swx_inverse.precomputedInterpolateI(x, y);
							final double y2 = swy_inverse.precomputedInterpolateI(x, y);
							double aux1 = u - x2;
							double aux2 = v - y2;

							f_direct += aux1 * aux1 + aux2 * aux2;

							// Compute the derivative with respect to all the c coefficients
							// Derivatives from direct coefficients.
							for (int l = 0; l < 4; l++)
								for (int m = 0; m < 4; m++) {
									if (swx_direct.prec_yIndex[v][l] == -1 || swx_direct.prec_xIndex[u][m] == -1)
										continue;

									double dddx = swx_direct.precomputedGetWeightI(l, m, u, v);
									double dixx = swx_inverse.precomputedGetWeightDx(l, m, x, y);
									double diyy = swy_inverse.precomputedGetWeightDy(l, m, x, y);

									double weightIx = (dixx + diyy) * dddx;

									double dddy = swy_direct.precomputedGetWeightI(l, m, u, v);
									double dixy = swx_inverse.precomputedGetWeightDy(l, m, x, y);
									double diyx = swy_inverse.precomputedGetWeightDx(l, m, x, y);

									double weightIy = (diyx + dixy) * dddy;

									int k = swx_direct.prec_yIndex[v][l] * cYdim + swx_direct.prec_xIndex[u][m];

									// Derivative related to X deformation
									this.grad_direct[k] += -aux1 * weightIx;

									// Derivative related to Y deformation
									this.grad_direct[k + twiceNk] += -aux2 * weightIy;
								}

							// Derivatives from inverse coefficients.
							for (int l = 0; l < 4; l++)
								for (int m = 0; m < 4; m++) {
									// d inverse(direct(x)) / d c_inverse
									if (swx_inverse.prec_yIndex[y][l] == -1 || swx_inverse.prec_xIndex[x][m] == -1)
										continue;

									double weightI = swx_inverse.precomputedGetWeightI(l, m, x, y);

									int k = swx_inverse.prec_yIndex[y][l] * cYdim + swx_inverse.prec_xIndex[x][m];

									// Derivative related to X deformation
									this.grad_direct[k + Nk] += -aux1 * weightI;

									// Derivative related to Y deformation
									this.grad_direct[k + Nk + twiceNk] += -aux2 * weightI;
								}

							n_direct++; // Another point has been successfully evaluated
						}
					} // end if mask.
				}

			if (n_direct != 0) {
				// Average the image related terms
				double aux = consistencyWeight * 2.0; // This is the 2 coming from the
				// derivative that I would do later
				for (int k = 0; k < grad_direct.length; k++)
					grad_direct[k] *= aux;
			}

			// Inverse gradient
			// Initialize
			for (int k = 0; k < this.grad_inverse.length; k++)
				this.grad_inverse[k] = 0.0F;

			// *********** Compute the geometric error and gradient (INVERSE)
			// ***********
			// The source rectangle marks the area of the source image to be treated.
			final int YdimS = rect_source.y + rect_source.height;
			final int XdimS = rect_source.x + rect_source.width;

			double f_inverse = 0;
			int n_inverse = 0;
			for (int v = rect_source.y; v < YdimS; v++)
				for (int u = rect_source.x; u < XdimS; u++) {
					// Check if this point is in the target mask
					if (this.transf.sourceMask == null || this.transf.sourceMask.contains(u / this.transf.sourceFactorWidth,
					    v / this.transf.sourceFactorHeight)) {
						final int x = (int) Math.round(swx_inverse.precomputedInterpolateI(u, v));
						final int y = (int) Math.round(swy_inverse.precomputedInterpolateI(u, v));

						if (x >= 0 && x < this.transf.targetCurrentWidth && y >= 0 && y < this.transf.targetCurrentHeight) {
							final double x2 = swx_direct.precomputedInterpolateI(x, y);
							final double y2 = swy_direct.precomputedInterpolateI(x, y);
							double aux1 = u - x2;
							double aux2 = v - y2;

							f_inverse += aux1 * aux1 + aux2 * aux2;

							// Compute the derivative with respect to all the c coefficients
							// Derivatives from direct coefficients.
							for (int l = 0; l < 4; l++)
								for (int m = 0; m < 4; m++) {
									// d direct(inverse(x)) / d c_direct
									if (swx_direct.prec_yIndex[y][l] == -1 || swx_direct.prec_xIndex[x][m] == -1)
										continue;

									double weightI = swx_direct.precomputedGetWeightI(l, m, x, y);

									int k = swx_direct.prec_yIndex[y][l] * cYdim + swx_direct.prec_xIndex[x][m];

									// Derivative related to X deformation
									this.grad_inverse[k] += -aux1 * weightI;

									// Derivative related to Y deformation
									this.grad_inverse[k + twiceNk] += -aux2 * weightI;
								}
							// Derivatives from inverse coefficients.
							for (int l = 0; l < 4; l++)
								for (int m = 0; m < 4; m++) {
									if (swx_inverse.prec_yIndex[v][l] == -1 || swx_inverse.prec_xIndex[u][m] == -1)
										continue;

									double diix = swx_inverse.precomputedGetWeightI(l, m, u, v);
									double ddxx = swx_direct.precomputedGetWeightDx(l, m, x, y);
									double ddyy = swy_direct.precomputedGetWeightDy(l, m, x, y);

									double weightIx = (ddxx + ddyy) * diix;

									double diiy = swy_inverse.precomputedGetWeightI(l, m, u, v);
									double ddxy = swx_direct.precomputedGetWeightDy(l, m, x, y);
									double ddyx = swy_direct.precomputedGetWeightDx(l, m, x, y);

									double weightIy = (ddyx + ddxy) * diiy;

									int k = swx_inverse.prec_yIndex[v][l] * cYdim + swx_inverse.prec_xIndex[u][m];

									// Derivative related to X deformation
									this.grad_inverse[k + Nk] += -aux1 * weightIx;

									// Derivative related to Y deformation
									this.grad_inverse[k + Nk + twiceNk] += -aux2 * weightIy;
								}

							n_inverse++; // Another point has been successfully evaluated
						}
					} // end if mask
				} // end inverse geometric error calculation

			if (n_inverse != 0) {
				// Average the image related terms
				double aux = consistencyWeight * 2.0; // This is the 2 coming from the
				// derivative that I would do later
				for (int k = 0; k < this.grad_inverse.length; k++)
					this.grad_inverse[k] *= aux;
			}

			// Save results
			this.result[0] = f_direct;
			this.result[1] = n_direct;
			this.result[2] = f_inverse;
			this.result[3] = n_inverse;

		} // end run method

	}

	/**
	 * Optimize the B-spline coefficients (unidirectional method).
	 *
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param thChangef
	 * @param cxTargetToSource
	 *          x- B-spline coefficients storing the target to source deformation
	 * @param cyTargetToSource
	 *          y- B-spline coefficients storing the target to source deformation
	 * @return energy function value
	 */
	private double optimizeCoeffs(int intervals, double thChangef, double[][] cxTargetToSource,
	    double[][] cyTargetToSource) {
		if (sourceModel.isSubOutput()) {
			System.out.println(" -----\n Intervals = " + intervals + "x" + intervals);
			System.out.println(" Source Image Size = " + this.sourceCurrentWidth + "x" + this.sourceCurrentHeight);
		}

		if (plugin != null && plugin.isPluginInterrumped())
			return 0.0;

		final double TINY = FLT_EPSILON;
		final double EPS = 3.0e-8F;
		final double FIRSTLAMBDA = 1;
		final int MAXITER_OPTIMCOEFF = 300;
		final int CUMULATIVE_SIZE = 5;

		int int3 = intervals + 3;
		int halfM = int3 * int3;
		int M = halfM * 2;

		double rescuedf, f;
		double[] x = new double[M];
		double[] rescuedx = new double[M];
		double[] diffx = new double[M];
		double[] rescuedgrad = new double[M];
		double[] grad = new double[M];
		double[] diffgrad = new double[M];
		double[] Hdx = new double[M];
		double[] rescuedhess = new double[M * M];
		double[] hess = new double[M * M];
		double[] proposedHess = new double[M * M];
		boolean[] optimize = new boolean[M];
		int i, j, p, iter = 1;
		boolean skip_update, ill_hessian;
		double improvementx = (double) Math.sqrt(TINY), lambda = FIRSTLAMBDA, max_normx, distx, aux, gmax;
		double fac, fae, dgdx, dxHdx, sumdiffg, sumdiffx;

		CumulativeQueue lastBest = new CumulativeQueue(CUMULATIVE_SIZE);

		for (i = 0; i < M; i++)
			optimize[i] = true;

		/* Form the vector with the current guess for the optimization */
		for (i = 0, p = 0; i < intervals + 3; i++)
			for (j = 0; j < intervals + 3; j++, p++) {
				x[p] = cxTargetToSource[i][j];
				x[halfM + p] = cyTargetToSource[i][j];
			}

		/* Prepare the precomputed weights for interpolation */
		this.swxTargetToSource = new BSplineModel(x, intervals + 3, intervals + 3, 0);
		this.swyTargetToSource = new BSplineModel(x, intervals + 3, intervals + 3, halfM);
		this.swxTargetToSource.precomputedPrepareForInterpolation(targetModel.getCurrentHeight(),
		    targetModel.getCurrentWidth(), intervals);
		this.swyTargetToSource.precomputedPrepareForInterpolation(targetModel.getCurrentHeight(),
		    targetModel.getCurrentWidth(), intervals);

		// First computation of the energy (similarity + landmarks + regularization)
		// f = evaluateSimilarity(x, intervals, grad, false, false, false);
		f = evaluateSimilarityMultiThread(x, intervals, grad, false, false);

		if (showMarquardtOptim)
			System.out.println("f(1)=" + f);

		/*
		 * Initially the hessian is the identity matrix multiplied by the first
		 * function value
		 */
		for (i = 0, p = 0; i < M; i++)
			for (j = 0; j < M; j++, p++)
				if (i == j)
					hess[p] = 1.0F;
				else
					hess[p] = 0.0F;

		rescuedf = f;
		for (i = 0, p = 0; i < M; i++) {
			rescuedx[i] = x[i];
			rescuedgrad[i] = grad[i];
			for (j = 0; j < M; j++, p++)
				rescuedhess[p] = hess[p];
		}

		// Maximum iteration number
		int maxiter = MAXITER_OPTIMCOEFF * (sourceModel.getCurrentDepth() + 1);

		ProgressBar.stepProgressBar();

		int last_successful_iter = 0;

		boolean stop = plugin != null && plugin.isPluginInterrumped();

		while (iter < maxiter && !stop) {
			/* Compute new x ------------------------------------------------- */
			Marquardt_it(x, optimize, grad, hess, lambda);

			/* Stopping criteria --------------------------------------------- */
			/* Compute difference with the previous iteration */
			max_normx = improvementx = 0;
			for (i = 0; i < M; i++) {
				diffx[i] = x[i] - rescuedx[i];
				distx = Math.abs(diffx[i]);
				improvementx += distx * distx;
				aux = Math.abs(rescuedx[i]) < Math.abs(x[i]) ? x[i] : rescuedx[i];
				max_normx += aux * aux;
			}

			if (TINY < max_normx)
				improvementx = improvementx / max_normx;

			improvementx = (double) Math.sqrt(Math.sqrt(improvementx));

			/*
			 * If there is no change with respect to the old geometry then finish the
			 * iterations
			 */
			if (improvementx < Math.sqrt(TINY))
				break;

			/* Estimate the new function value -------------------------------- */
			// f = evaluateSimilarity(x, intervals, grad, false, false, false);
			f = evaluateSimilarityMultiThread(x, intervals, grad, false, false);

			iter++;
			if (showMarquardtOptim)
				System.out.println("f(" + iter + ")=" + f + " lambda=" + lambda);
			ProgressBar.stepProgressBar();

			/* Update lambda -------------------------------------------------- */
			if (rescuedf > f) {
				// We save the last energy terms values in order to be displayed.
				this.finalDirectConsistencyError = this.partialDirectConsitencyError;
				this.finalDirectSimilarityError = this.partialDirectSimilarityError;
				this.finalDirectRegularizationError = this.partialDirectRegularizationError;
				this.finalDirectLandmarkError = this.partialDirectLandmarkError;

				/* Check if the improvement is only residual */
				lastBest.push_back(rescuedf - f);
				if (lastBest.currentSize() == CUMULATIVE_SIZE && lastBest.getSum() / f < thChangef)
					break;

				/*
				 * If we have improved then estimate the hessian, update the geometry,
				 * and decrease the lambda
				 */
				/* Estimate the hessian ....................................... */
				if (showMarquardtOptim)
					System.out.println("  Accepted");
				if ((last_successful_iter++ % 10) == 0 && outputLevel > -1)
					updateCurrentOutput(x, intervals, false);

				/* Estimate the difference between gradients */
				for (i = 0; i < M; i++)
					diffgrad[i] = grad[i] - rescuedgrad[i];

				/* Multiply this difference by the current inverse of the hessian */
				for (i = 0, p = 0; i < M; i++) {
					Hdx[i] = 0.0F;
					for (j = 0; j < M; j++, p++)
						Hdx[i] += hess[p] * diffx[j];
				}

				/* Calculate dot products for the denominators ................ */
				dgdx = dxHdx = sumdiffg = sumdiffx = 0.0F;
				skip_update = true;
				for (i = 0; i < M; i++) {
					dgdx += diffgrad[i] * diffx[i];
					dxHdx += diffx[i] * Hdx[i];
					sumdiffg += diffgrad[i] * diffgrad[i];
					sumdiffx += diffx[i] * diffx[i];
					if (Math.abs(grad[i]) >= Math.abs(rescuedgrad[i]))
						gmax = Math.abs(grad[i]);
					else
						gmax = Math.abs(rescuedgrad[i]);
					if (gmax != 0 && Math.abs(diffgrad[i] - Hdx[i]) > Math.sqrt(EPS) * gmax)
						skip_update = false;
				}

				/* Update hessian ............................................. */
				/* Skip if fac not sufficiently positive */
				if (dgdx > Math.sqrt(EPS * sumdiffg * sumdiffx) && !skip_update) {
					fae = 1.0F / dxHdx;
					fac = 1.0F / dgdx;

					/* Update the hessian after BFGS formula */
					for (i = 0, p = 0; i < M; i++)
						for (j = 0; j < M; j++, p++) {
							if (i <= j)
								proposedHess[p] = hess[p] + fac * diffgrad[i] * diffgrad[j] - fae * (Hdx[i] * Hdx[j]);
							else
								proposedHess[p] = proposedHess[j * M + i];
						}

					ill_hessian = false;
					if (!ill_hessian) {
						for (i = 0, p = 0; i < M; i++)
							for (j = 0; j < M; j++, p++)
								hess[p] = proposedHess[p];
					} else if (showMarquardtOptim)
						System.out.println("Hessian cannot be safely updated, ill-conditioned");

				} else if (showMarquardtOptim)
					System.out.println("Hessian cannot be safely updated");

				/* Update geometry and lambda ................................. */
				rescuedf = f;
				for (i = 0, p = 0; i < M; i++) {
					rescuedx[i] = x[i];
					rescuedgrad[i] = grad[i];
					for (j = 0; j < M; j++, p++)
						rescuedhess[p] = hess[p];
				}
				if (1e-4 < lambda)
					lambda = lambda / 10;
			} else {
				/*
				 * else, if it is worse, then recover the last geometry and increase
				 * lambda, saturate lambda with FIRSTLAMBDA
				 */
				for (i = 0, p = 0; i < M; i++) {
					x[i] = rescuedx[i];
					grad[i] = rescuedgrad[i];
					for (j = 0; j < M; j++, p++)
						hess[p] = rescuedhess[p];
				}
				if (lambda < 1.0 / TINY)
					lambda *= 10;
				else
					break;
				if (lambda < FIRSTLAMBDA)
					lambda = FIRSTLAMBDA;
			}

			stop = plugin != null && plugin.isPluginInterrumped();
		}

		// Copy the values back to the input arrays
		for (i = 0, p = 0; i < intervals + 3; i++)
			for (j = 0; j < intervals + 3; j++, p++) {
				cxTargetToSource[i][j] = x[p];

				cyTargetToSource[i][j] = x[halfM + p];
			}

		ProgressBar.skipProgressBar(maxiter - iter);
		return f;
	}

	/**
	 * Propagate deformation coefficients to the next level.
	 *
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param c
	 *          B-spline coefficients
	 * @param expansionFactor
	 *          due to the change of size in the represented image
	 * @return propagated coefficients
	 */
	private double[][] propagateCoeffsToNextLevel(int intervals, final double[][] c, double expansionFactor) {
		// Expand the coefficients for the next scale
		intervals *= 2;
		double[][] cs_expand = new double[intervals + 7][intervals + 7];

		// Upsample
		for (int i = 0; i < intervals + 7; i++)
			for (int j = 0; j < intervals + 7; j++) {
				// If it is not in an even sample then set it to 0
				if (i % 2 == 0 || j % 2 == 0)
					cs_expand[i][j] = 0.0F;
				else {
					// Now look for this sample in the coarser level
					int ipc = (i - 1) / 2;
					int jpc = (j - 1) / 2;
					cs_expand[i][j] = c[ipc][jpc];
				}
			}

		// Define the FIR filter
		double[][] u2n = new double[4][];
		u2n[0] = null;
		u2n[1] = new double[3];
		u2n[1][0] = 0.5F;
		u2n[1][1] = 1.0F;
		u2n[1][2] = 0.5F;
		u2n[2] = null;
		u2n[3] = new double[5];
		u2n[3][0] = 0.125F;
		u2n[3][1] = 0.5F;
		u2n[3][2] = 0.75F;
		u2n[3][3] = 0.5F;
		u2n[3][4] = 0.125F;
		int[] half_length_u2n = { 0, 1, 0, 2 };
		int kh = half_length_u2n[transformationSplineDegree];

		// Apply the u2n filter to rows
		double[][] cs_expand_aux = new double[intervals + 7][intervals + 7];

		for (int i = 1; i < intervals + 7; i += 2)
			for (int j = 0; j < intervals + 7; j++) {
				cs_expand_aux[i][j] = 0.0F;
				for (int k = -kh; k <= kh; k++)
					if (j + k >= 0 && j + k <= intervals + 6)
						cs_expand_aux[i][j] += u2n[transformationSplineDegree][k + kh] * cs_expand[i][j + k];
			}

		// Apply the u2n filter to columns
		for (int i = 0; i < intervals + 7; i++)
			for (int j = 0; j < intervals + 7; j++) {
				cs_expand[i][j] = 0.0F;
				for (int k = -kh; k <= kh; k++)
					if (i + k >= 0 && i + k <= intervals + 6)
						cs_expand[i][j] += u2n[transformationSplineDegree][k + kh] * cs_expand_aux[i + k][j];
			}

		// Copy the central coefficients to c
		double[][] newc = new double[intervals + 3][intervals + 3];
		for (int i = 0; i < intervals + 3; i++)
			for (int j = 0; j < intervals + 3; j++)
				newc[i][j] = cs_expand[i + 2][j + 2] * expansionFactor;

		// Return the new set of coefficients
		return newc;
	}

	/**
	 * Show the direct transformation results(multi-thread version).
	 */
	public void showDirectResults() {
		showTransformationMultiThread(intervals, cxTargetToSource, cyTargetToSource, false);
	}

	/**
	 * Show the transformation (multi-thread version).
	 *
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param cx
	 *          x- deformation coefficients
	 * @param cy
	 *          y- deformation coefficients
	 * @param bIsReverse
	 *          flag to determine the transformation direction
	 *          (target-source=FALSE or source-target=TRUE)
	 */
	private void showTransformationMultiThread(final int intervals, final double[][] cx, // Input,
	                                                                                     // spline
	                                                                                     // coefficients
	    final double[][] cy, boolean bIsReverse) {

		Sequence outputSeq = (!bIsReverse) ? this.outputSeq1 : this.outputSeq2;

		// Calculate tranformation results
		ProgressBar.setProgressBarMessage("Calculating result window...");
		Sequence result_imp = applyTransformationMultiThread(intervals, cx, cy, bIsReverse);

		plugin.removeSequence(outputSeq);

		outputSeq = result_imp;

		plugin.addSequence(outputSeq);

	}

	/**
	 * Apply the final transformation (multi-thread version).
	 *
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param cx
	 *          x- deformation coefficients
	 * @param cy
	 *          y- deformation coefficients
	 * @param bIsReverse
	 *          flag to determine the transformation direction
	 *          (target-source=FALSE or source-target=TRUE)
	 * 
	 * @return output images (depending on the output level)
	 */
	private Sequence applyTransformationMultiThread(final int intervals, final double[][] cx, // Input,
	                                                                                          // spline
	                                                                                          // coefficients
	    final double[][] cy, boolean bIsReverse) {
		// BSplineModel auxTarget = targetModel;
		// BSplineModel auxSource = sourceModel;
		ROI2D auxTargetMsk = targetMask;
		ROI2D auxSourceMsk = sourceMask;
		int auxTargetWidth = this.originalTargetIBI.getWidth();
		int auxTargetHeight = this.originalTargetIBI.getHeight();
		IcyBufferedImage originalIBI = this.originalSourceIBI;

		// Change if necessary
		if (bIsReverse) {
			// auxTarget = sourceModel;
			// auxSource = targetModel;
			auxTargetMsk = sourceMask;
			auxSourceMsk = targetMask;
			auxTargetWidth = this.originalSourceIBI.getWidth();
			auxTargetHeight = this.originalSourceIBI.getHeight();
			originalIBI = this.originalTargetIBI;
		}

		final Sequence is = new Sequence();
		final String s = bIsReverse ? new String("Target") : new String("Source");

		// Create transformation B-spline models
		BSplineModel swx = new BSplineModel(cx);
		BSplineModel swy = new BSplineModel(cy);

		// We compute the deformation (transformation_x and transformation_y) on the
		// fly

		// /* GRAY SCALE IMAGES */
		// if(!(originalIP instanceof ColorProcessor))
		// {
		// final FloatProcessor fp = new FloatProcessor(auxTargetWidth,
		// auxTargetHeight);
		// final FloatProcessor fp_mask = new FloatProcessor(auxTargetWidth,
		// auxTargetHeight);
		// final FloatProcessor fp_target = new FloatProcessor(auxTargetWidth,
		// auxTargetHeight, auxTarget.getOriginalImage());
		//
		//
		// // take original processor if necessary
		// if(auxSource.getOriginalImageWidth() > auxSource.getWidth())
		// {
		// auxSource = new BSplineModel( originalIP, false, 1);
		// auxSource.setPyramidDepth(0);
		// auxSource.startPyramids();
		//
		// // Join thread
		// try {
		// auxSource.getThread().join();
		// } catch (InterruptedException e) {
		// IJ.error("Unexpected interruption exception " + e);
		// }
		// }
		//
		// // Check the number of processors in the computer
		// int nproc = Runtime.getRuntime().availableProcessors();
		//
		// // We will use threads to display parts of the output image
		// int block_height = auxTargetHeight / nproc;
		// if (auxTargetHeight % 2 != 0)
		// block_height++;
		//
		//
		// int nThreads = nproc; /*(nproc > 1) ? (nproc / 2) : 1;
		// if (this.accurate_mode == MainDialog.MONO_MODE)
		// nThreads *= 2;*/
		//
		//
		// Thread[] threads = new Thread[nThreads];
		// Rectangle[] rects = new Rectangle[nThreads];
		// FloatProcessor[] fp_tile = new FloatProcessor[nThreads];
		// FloatProcessor[] fp_mask_tile = new FloatProcessor[nThreads];
		//
		// for (int i=0; i<nThreads; i++)
		// {
		// // last block size is the rest of the window
		// int y_start = i*block_height;
		//
		// if (nThreads-1 == i)
		// block_height = auxTargetHeight - i*block_height;
		//
		// rects[i] = new Rectangle(0, y_start, auxTargetWidth, block_height);
		//
		// //IJ.log("block = 0 " + (i*block_height) + " " + auxTargetWidth + " " +
		// block_height );
		//
		// fp_tile[i] = new FloatProcessor(rects[i].width, rects[i].height);
		// fp_mask_tile[i] = new FloatProcessor(rects[i].width, rects[i].height);
		//
		// threads[i] = new Thread(new GrayscaleResultTileMaker(swx, swy, auxSource,
		// auxTargetWidth, auxTargetHeight,
		// auxTargetMsk, auxSourceMsk,
		// rects[i], fp_tile[i], fp_mask_tile[i]));
		// threads[i].start();
		// }
		//
		// for (int i=0; i<nThreads; i++)
		// {
		// try {
		// threads[i].join();
		// threads[i] = null;
		// } catch (InterruptedException e) {
		// e.printStackTrace();
		// }
		// }
		//
		// for (int i=0; i<nThreads; i++)
		// {
		// fp.insert(fp_tile[i], rects[i].x, rects[i].y);
		// fp_tile[i] = null;
		// fp_mask.insert(fp_mask_tile[i], rects[i].x, rects[i].y);
		// fp_mask_tile[i] = null;
		// rects[i] = null;
		// }
		//
		// fp.resetMinAndMax();
		//
		// // Add slices to result stack
		// is.addSlice("Registered " + s + " Image", fp);
		// //if (outputLevel > -1)
		// is.addSlice("Target Image", fp_target);
		// //if (outputLevel > -1)
		// is.addSlice("Warped Source Mask",fp_mask);
		// }
		// else /* COLOR IMAGES */
		{
			BSplineModel[] sourceModels = new BSplineModel[originalIBI.getSizeC()];
			for (int c = 0; c < originalIBI.getSizeC(); c++) {
				sourceModels[c] = new BSplineModel(IcyBufferedImageUtil.extractChannel(originalIBI, c), false, 1);
				sourceModels[c].setPyramidDepth(0);
				sourceModels[c].startPyramids();
			}

			// Join threads
			try {
				for (int c = 0; c < originalIBI.getSizeC(); c++) {
					sourceModels[c].join();
				}
			} catch (InterruptedException e) {
				System.err.println("Unexpected interruption exception " + e);
			}

			// Calculate warped RGB image
			IcyBufferedImage ibi = new IcyBufferedImage(auxTargetWidth, auxTargetHeight, originalIBI.getSizeC(),
			    originalIBI.getDataType_());
			IcyBufferedImage ibi_mask = new IcyBufferedImage(auxTargetWidth, auxTargetHeight, originalIBI.getSizeC(),
			    originalIBI.getDataType_());

			// Check the number of processors in the computer
			int nproc = Runtime.getRuntime().availableProcessors();

			// We will use threads to display parts of the output image
			int block_height = auxTargetHeight / nproc;
			if (auxTargetHeight % 2 != 0)
				block_height++;

			int nThreads = nproc; /*
			                       * (nproc > 1) ? (nproc / 2) : 1; if
			                       * (this.accurate_mode == MainDialog.MONO_MODE)
			                       * nThreads *= 2;
			                       */

			Thread[] threads = new Thread[nThreads];
			Rectangle[] rects = new Rectangle[nThreads];
			IcyBufferedImage[] ibi_tile = new IcyBufferedImage[nThreads];
			IcyBufferedImage[] ibi_mask_tile = new IcyBufferedImage[nThreads];

			for (int i = 0; i < nThreads; i++) {
				// last block size is the rest of the window
				int y_start = i * block_height;

				if (nThreads - 1 == i)
					block_height = auxTargetHeight - i * block_height;

				rects[i] = new Rectangle(0, y_start, auxTargetWidth, block_height);

				// IJ.log("block = 0 " + (i*block_height) + " " + auxTargetWidth + " " +
				// block_height );

				ibi_tile[i] = new IcyBufferedImage(rects[i].width, rects[i].height, originalIBI.getSizeC(),
				    originalIBI.getDataType_());
				ibi_mask_tile[i] = new IcyBufferedImage(rects[i].width, rects[i].height, originalIBI.getSizeC(),
				    originalIBI.getDataType_());

				threads[i] = new Thread(new ColorResultTileMaker(swx, swy, sourceModels, auxTargetWidth, auxTargetHeight,
				    auxTargetMsk, auxSourceMsk, rects[i], ibi_tile[i], ibi_mask_tile[i]));
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

			for (int i = 0; i < nThreads; i++) {
				ibi.copyData(ibi_tile[i], null, new Point(rects[i].x, rects[i].y));
				ibi_tile[i] = null;

				ibi_mask.copyData(ibi_mask_tile[i], null, new Point(rects[i].x, rects[i].y));
				ibi_mask_tile[i] = null;
				rects[i] = null;
			}

			ibi.dataChanged();

			// Add slices to result stack
			is.beginUpdate();
			is.setName("Registered " + s + " Image");
			is.addImage(ibi);
			// if (outputLevel > -1)
			is.addImage(bIsReverse ? this.originalSourceIBI : this.originalTargetIBI);
			// if (outputLevel > -1)
			is.addImage(ibi_mask);

		} // end caculate warped color image

		if (outputLevel == 2) {
			computeDeformationVectors(intervals, cx, cy, is, bIsReverse);
			computeDeformationGrid(intervals, cx, cy, is, bIsReverse);
		}

		is.endUpdate();
		return is;
	}

	/**
	 * Class to run concurrent tile windows for final results (color)
	 * 
	 */
	private class ColorResultTileMaker implements Runnable {
		final BSplineModel swx;
		final BSplineModel swy;
		final BSplineModel[] sourceModels;
		final int auxTargetCurrentWidth;
		final int auxTargetCurrentHeight;
		final ROI2D auxTargetMsk;
		final ROI2D auxSourceMsk;
		final Rectangle rect;
		final private IcyBufferedImage ibiTile;
		final private IcyBufferedImage ibiMaskTile;

		// ------------------------------------------------------------------
		/**
		 * Color result tile maker constructor
		 * 
		 * @param swx
		 *          B-spline interpolator for transformation in x-
		 * @param swy
		 *          B-spline interpolator for transformation in y-
		 * @param sourceR
		 *          red source image
		 * @param sourceG
		 *          green source image
		 * @param sourceB
		 *          blue source image
		 * @param auxTargetCurrentWidth
		 *          current target height
		 * @param auxTargetCurrentHeight
		 *          current target height
		 * @param auxTargetMsk
		 *          target mask
		 * @param auxSourceMsk
		 *          source mask
		 * @param rect
		 *          retangle with the coordinates of the output image to be updated
		 * @param fpR
		 *          red channel processor to be updated
		 * @param fpG
		 *          green channel processor to be updated
		 * @param fpB
		 *          blue channel processor to be updated
		 * @param cp_mask
		 *          mask color processor to be updated
		 */
		ColorResultTileMaker(BSplineModel swx, BSplineModel swy, BSplineModel[] sourceModels, int auxTargetCurrentWidth,
		    int auxTargetCurrentHeight, ROI2D auxTargetMsk, ROI2D auxSourceMsk, Rectangle rect, IcyBufferedImage fpB,
		    IcyBufferedImage cp_mask) {
			this.swx = swx;
			this.swy = swy;
			this.sourceModels = sourceModels;
			this.auxTargetCurrentWidth = auxTargetCurrentWidth;
			this.auxTargetCurrentHeight = auxTargetCurrentHeight;
			this.auxTargetMsk = auxTargetMsk;
			this.auxSourceMsk = auxSourceMsk;
			this.rect = rect;
			this.ibiTile = fpB;
			this.ibiMaskTile = cp_mask;
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

			double[][] ibiData = Array2DUtil.arrayToDoubleArray(ibiTile.getDataXYC(), ibiTile.isSignedDataType());
			double[][] ibiMaskData = Array2DUtil.arrayToDoubleArray(ibiMaskTile.getDataXYC(), ibiMaskTile.isSignedDataType());

			for (int v_rect = 0, v = rect.y; v < auxTargetHeight; v++, v_rect++) {
				final int v_offset = v_rect * rect.width;
				final double tv = (double) (v * intervals) / (double) (auxTargetCurrentHeight - 1) + 1.0F;

				for (int u_rect = 0, u = rect.x; u < auxTargetWidth; u++, u_rect++) {

					final double tu = (double) (u * intervals) / (double) (auxTargetCurrentWidth - 1) + 1.0F;
					final double transformation_x_v_u = swx.prepareForInterpolationAndInterpolateI(tu, tv, false, ORIGINAL);
					final double transformation_y_v_u = swy.prepareForInterpolationAndInterpolateI(tu, tv, false, ORIGINAL);

					if (auxTargetMsk != null && !auxTargetMsk.contains(u, v)) {
						for (int c = 0; c < ibiTile.getSizeC(); c++) {
							ibiData[c][u_rect + v_offset] = 0;
							ibiMaskData[c][u_rect + v_offset] = 0;
						}
					} else {

						final double x = transformation_x_v_u;
						final double y = transformation_y_v_u;
						if (auxSourceMsk == null || auxSourceMsk.contains(x, y)) {
							for (int c = 0; c < ibiTile.getSizeC(); c++) {
								ibiData[c][u_rect + v_offset] = sourceModels[c].prepareForInterpolationAndInterpolateI(x, y, false,
								    ORIGINAL);
								ibiMaskData[c][u_rect + v_offset] = ibiMaskTile.getDataTypeMax();
							}
						} else {
							for (int c = 0; c < ibiTile.getSizeC(); c++) {
								ibiData[c][u_rect + v_offset] = 0;
								ibiMaskData[c][u_rect + v_offset] = 0;
							}
						}
					}
				}
			}

			Array2DUtil.doubleArrayToSafeArray(ibiData, ibiTile.getDataXYC(), ibiTile.isSignedDataType());
			Array2DUtil.doubleArrayToSafeArray(ibiMaskData, ibiMaskTile.getDataXYC(), ibiMaskTile.isSignedDataType());
			/*
			 * (new ImagePlus("Red", fpR)).show(); (new ImagePlus("Green",
			 * fpG)).show(); (new ImagePlus("Blue", fpB)).show();
			 */

		} /* end run method */

	}

	/**
	 * Compute and draw the final deformation vectors.
	 *
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param cx
	 *          x- deformation coefficients
	 * @param cy
	 *          y- deformation coefficients
	 * @param is
	 *          image stack where we want to show the deformation vectors
	 * @param bIsReverse
	 *          flag to determine the transformation direction
	 *          (target-source=FALSE or source-target=TRUE)
	 */
	private void computeDeformationVectors(int intervals, double[][] cx, double[][] cy, Sequence is, boolean bIsReverse) {
		// Auxiliar variables for changing from source to target and inversely
		ROI2D auxTargetMsk = this.targetMask;
		ROI2D auxSourceMsk = this.sourceMask;
		int auxTargetCurrentHeight = this.targetCurrentHeight;
		int auxTargetCurrentWidth = this.targetCurrentWidth;

		// Change if necessary
		if (bIsReverse) {
			auxTargetMsk = this.sourceMask;
			auxSourceMsk = this.targetMask;
			auxTargetCurrentHeight = this.sourceCurrentHeight;
			auxTargetCurrentWidth = this.sourceCurrentWidth;
		}

		// Initialize output image
		int stepv = Math.min(Math.max(10, auxTargetCurrentHeight / 15), 30);
		int stepu = Math.min(Math.max(10, auxTargetCurrentWidth / 15), 30);

		final double transformedImage[][] = new double[auxTargetCurrentHeight][auxTargetCurrentWidth];

		for (int v = 0; v < auxTargetCurrentHeight; v++)
			for (int u = 0; u < auxTargetCurrentWidth; u++)
				transformedImage[v][u] = is.getDataTypeMax();

		// Ask for memory for the transformation
		double[][] transformation_x = new double[auxTargetCurrentHeight][auxTargetCurrentWidth];
		double[][] transformation_y = new double[auxTargetCurrentHeight][auxTargetCurrentWidth];

		// Compute the deformation
		computeDeformation(intervals, cx, cy, transformation_x, transformation_y, bIsReverse);

		// Show shift field ........................................
		// Show deformation vectors
		for (int v = 0; v < auxTargetCurrentHeight; v += stepv)
			for (int u = 0; u < auxTargetCurrentWidth; u += stepu)
				if (auxTargetMsk == null || auxTargetMsk.contains(u, v)) {
					final double x = transformation_x[v][u];
					final double y = transformation_y[v][u];
					if (auxSourceMsk == null || auxSourceMsk.contains(x, y))
						MiscTools.drawArrow(transformedImage, u, v, (int) Math.round(x), (int) Math.round(y), 0, 2);
				}

		// Set it to the image stack
		IcyBufferedImage fp = new IcyBufferedImage(auxTargetCurrentWidth, auxTargetCurrentHeight, is.getSizeC(),
		    is.getDataType_());
		double[][] fpData = Array2DUtil.arrayToDoubleArray(fp.getDataXYC(), fp.isSignedDataType());
		for (int v = 0; v < auxTargetCurrentHeight; v++)
			for (int u = 0; u < auxTargetCurrentWidth; u++)
				for (int c = 0; c < is.getSizeC(); c++)
					fpData[c][u + v * auxTargetCurrentWidth] = transformedImage[v][u];

		Array2DUtil.doubleArrayToSafeArray(fpData, fp.getDataXYC(), fp.isSignedDataType());
		fp.dataChanged();
		is.addImage(fp);
	}

	/**
	 * Compute the deformation.
	 *
	 * @param intervals
	 *          input, number of intervals
	 * @param cx
	 *          input, X B-spline coefficients
	 * @param cy
	 *          input, Y B-spline coefficients
	 * @param transformation_x
	 *          output, X transformation map
	 * @param transformation_y
	 *          output, Y transformation map
	 * @param bIsReverse
	 *          determines the transformation direction (source-target=TRUE or
	 *          target-source=FALSE)
	 */
	private void computeDeformation(final int intervals, final double[][] cx, final double[][] cy,
	    final double[][] transformation_x, final double[][] transformation_y, boolean bIsReverse) {

		int auxTargetCurrentHeight = this.targetCurrentHeight;
		int auxTargetCurrentWidth = this.targetCurrentWidth;

		if (bIsReverse) {
			auxTargetCurrentHeight = this.sourceCurrentHeight;
			auxTargetCurrentWidth = this.sourceCurrentWidth;
		}

		/*
		 * // Set these coefficients to an interpolator BSplineModel swx = new
		 * BSplineModel(cx); BSplineModel swy = new BSplineModel(cy);
		 * 
		 * 
		 * // Compute the transformation mapping for (int v=0;
		 * v<auxTargetCurrentHeight; v++) { final double tv = (double)(v *
		 * intervals) / (double)(auxTargetCurrentHeight - 1) + 1.0F; for (int u = 0;
		 * u<auxTargetCurrentWidth; u++) { final double tu = (double)(u * intervals)
		 * / (double)(auxTargetCurrentWidth - 1) + 1.0F;
		 * swx.prepareForInterpolation(tu, tv, ORIGINAL); transformation_x[v][u] =
		 * swx.interpolateI(); swy.prepareForInterpolation(tu, tv, ORIGINAL);
		 * transformation_y[v][u] = swy.interpolateI(); } }
		 * 
		 */

		Thread x_thread = new Thread(
		    new ConcurrentDeformation(cx, auxTargetCurrentHeight, auxTargetCurrentWidth, transformation_x, intervals));

		Thread y_thread = new Thread(
		    new ConcurrentDeformation(cy, auxTargetCurrentHeight, auxTargetCurrentWidth, transformation_y, intervals));

		x_thread.start();
		y_thread.start();

		try {
			x_thread.join();
			y_thread.join();
			x_thread = null;
			y_thread = null;
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Class to concurrently calculate the two deformation mapping tables
	 * 
	 */
	private class ConcurrentDeformation extends Thread {
		final double[][] c;
		final int auxTargetCurrentHeight;
		final int auxTargetCurrentWidth;
		final double[][] transformation;
		final int intervals;

		ConcurrentDeformation(double[][] c, int auxTargetCurrentHeight, int auxTargetCurrentWidth,
		    double[][] transformation, int intervals) {
			this.c = c;
			this.auxTargetCurrentWidth = auxTargetCurrentWidth;
			this.auxTargetCurrentHeight = auxTargetCurrentHeight;
			this.transformation = transformation;
			this.intervals = intervals;
		}

		// ------------------------------------------------------------------
		/**
		 * Run method to calculate the corresponding X or Y transformation table.
		 */
		public void run() {
			// Set these coefficients to an interpolator
			BSplineModel sw = new BSplineModel(c);

			// Compute the transformation mapping
			for (int v = 0; v < auxTargetCurrentHeight; v++) {
				final double tv = (double) (v * intervals) / (double) (auxTargetCurrentHeight - 1) + 1.0F;
				for (int u = 0; u < auxTargetCurrentWidth; u++) {
					final double tu = (double) (u * intervals) / (double) (auxTargetCurrentWidth - 1) + 1.0F;
					transformation[v][u] = sw.prepareForInterpolationAndInterpolateI(tu, tv, false, ORIGINAL);
				}
			}
		} /* end run */
	}

	/**
	 * Compute and draw the final deformation grid.
	 *
	 * @param intervals
	 *          number of intervals in the deformation
	 * @param cx
	 *          x- deformation coefficients
	 * @param cy
	 *          y- deformation coefficients
	 * @param is
	 *          image stack where we want to show the deformation grid
	 * @param bIsReverse
	 *          flag to determine the transformation direction
	 *          (target-source=FALSE or source-target=TRUE)
	 */
	private void computeDeformationGrid(int intervals, double[][] cx, double[][] cy, Sequence is, boolean bIsReverse) {

		int auxTargetCurrentHeight = this.targetCurrentHeight;
		int auxTargetCurrentWidth = this.targetCurrentWidth;

		if (bIsReverse) {
			auxTargetCurrentHeight = sourceCurrentHeight;
			auxTargetCurrentWidth = sourceCurrentWidth;
		}
		// Initialize output image
		int stepv = Math.min(Math.max(10, auxTargetCurrentHeight / 15), 30);
		int stepu = Math.min(Math.max(10, auxTargetCurrentWidth / 15), 30);
		final double transformedImage[][] = new double[auxTargetCurrentHeight][auxTargetCurrentWidth];
		for (int v = 0; v < auxTargetCurrentHeight; v++)
			for (int u = 0; u < auxTargetCurrentWidth; u++)
				transformedImage[v][u] = is.getDataTypeMax();

		// Ask for memory for the transformation
		double[][] transformation_x = new double[auxTargetCurrentHeight][auxTargetCurrentWidth];
		double[][] transformation_y = new double[auxTargetCurrentHeight][auxTargetCurrentWidth];

		// Compute the deformation
		computeDeformation(intervals, cx, cy, transformation_x, transformation_y, bIsReverse);

		// Show deformed grid ........................................
		// Show deformation vectors
		for (int v = 0; v < auxTargetCurrentHeight; v += stepv)
			for (int u = 0; u < auxTargetCurrentWidth; u += stepu) {
				final double x = transformation_x[v][u];
				final double y = transformation_y[v][u];
				// Draw horizontal line
				int uh = u + stepu;
				if (uh < auxTargetCurrentWidth) {
					final double xh = transformation_x[v][uh];
					final double yh = transformation_y[v][uh];
					MiscTools.drawLine(transformedImage, (int) Math.round(x), (int) Math.round(y), (int) Math.round(xh),
					    (int) Math.round(yh), 0);
				}

				// Draw vertical line
				int vv = v + stepv;
				if (vv < auxTargetCurrentHeight) {
					final double xv = transformation_x[vv][u];
					final double yv = transformation_y[vv][u];
					MiscTools.drawLine(transformedImage, (int) Math.round(x), (int) Math.round(y), (int) Math.round(xv),
					    (int) Math.round(yv), 0);
				}
			}

		// Set it to the image stack
		IcyBufferedImage fp = new IcyBufferedImage(auxTargetCurrentWidth, auxTargetCurrentHeight, is.getSizeC(),
		    is.getDataType_());
		double[][] fpData = Array2DUtil.arrayToDoubleArray(fp.getDataXYC(), fp.isSignedDataType());
		for (int v = 0; v < auxTargetCurrentHeight; v++)
			for (int u = 0; u < auxTargetCurrentWidth; u++)
				for (int c = 0; c < fp.getSizeC(); c++)
					fpData[c][u + v * auxTargetCurrentWidth] = transformedImage[v][u];

		Array2DUtil.doubleArrayToSafeArray(fpData, fp.getDataXYC(), fp.isSignedDataType());
		fp.dataChanged();
		is.addImage(fp);
	}

	/**
	 * Registration method. It applies the consistent and elastic registration
	 * algorithm to the selected source and target images.
	 */
	public void doBidirectionalRegistration() {
		// This function can only be applied with splines of an odd order

		// Bring into consideration the image/coefficients at the smallest scale
		sourceModel.popFromPyramid();
		targetModel.popFromPyramid();

		// size correction factor
		int sizeCorrectionFactor = 0; // this.targetHeight / (1024 * (int)
		                              // Math.pow(2,
		                              // this.maxImageSubsamplingFactor));
		// System.out.println("Size correction factor = " + sizeCorrectionFactor);

		targetCurrentHeight = targetModel.getCurrentHeight();
		targetCurrentWidth = targetModel.getCurrentWidth();

		targetFactorHeight = targetModel.getFactorHeight();
		targetFactorWidth = targetModel.getFactorWidth();

		sourceCurrentHeight = sourceModel.getCurrentHeight();
		sourceCurrentWidth = sourceModel.getCurrentWidth();

		sourceFactorHeight = sourceModel.getFactorHeight();
		sourceFactorWidth = sourceModel.getFactorWidth();

		// Ask memory for the transformation coefficients
		intervals = (int) Math.pow(2, minScaleDeformation + sizeCorrectionFactor);

		cxTargetToSource = new double[intervals + 3][intervals + 3];
		cyTargetToSource = new double[intervals + 3][intervals + 3];

		// Build matrices for computing the regularization
		buildRegularizationTemporary(intervals, false);
		buildRegularizationTemporary(intervals, true);

		// Ask for memory for the residues
		final int K;
		if (targetLandmarks != null)
			K = targetLandmarks.size();
		else
			K = 0;
		double[] dxTargetToSource = new double[K];
		double[] dyTargetToSource = new double[K];
		computeInitialResidues(dxTargetToSource, dyTargetToSource, false);
		computeInitialResidues(dxTargetToSource, dyTargetToSource, false);

		// Compute the affine transformation FROM THE TARGET TO THE SOURCE
		// coordinates
		// Notice that this matrix is independent of the scale (unless it was loaded
		// from
		// file), but the residues are not
		double[][] affineMatrix = null;
		// NOTE: after version 1.1 the landmarks are always used to calculate
		// an initial affine transformation (whether the landmarks weight is 0 or
		// not).
		affineMatrix = computeAffineMatrix(false);

		// Incorporate the affine transformation into the spline coefficient
		for (int i = 0; i < intervals + 3; i++) {
			final double v = (double) ((i - 1) * (targetCurrentHeight - 1)) / (double) intervals;
			final double xv = affineMatrix[0][2] + affineMatrix[0][1] * v;
			final double yv = affineMatrix[1][2] + affineMatrix[1][1] * v;
			for (int j = 0; j < intervals + 3; j++) {
				final double u = (double) ((j - 1) * (targetCurrentWidth - 1)) / (double) intervals;
				cxTargetToSource[i][j] = xv + affineMatrix[0][0] * u;
				cyTargetToSource[i][j] = yv + affineMatrix[1][0] * u;
			}
		}

		// Compute the affine transformation FROM THE SOURCE TO THE TARGET
		// coordinates
		// Notice again that this matrix is independent of the scale, but the
		// residues are not
		// Ask for memory for the residues
		final int K2;
		if (sourceLandmarks != null)
			K2 = sourceLandmarks.size();
		else
			K2 = 0;
		double[] dxSourceToTarget = new double[K2];
		double[] dySourceToTarget = new double[K2];
		computeInitialResidues(dxSourceToTarget, dySourceToTarget, true);
		computeInitialResidues(dxSourceToTarget, dySourceToTarget, true);

		cxSourceToTarget = new double[intervals + 3][intervals + 3];
		cySourceToTarget = new double[intervals + 3][intervals + 3];

		// NOTE: after version 1.1 the landmarks are always used to calculate
		// an initial affine transformation (whether the landmarks weight is 0 or
		// not).
		affineMatrix = computeAffineMatrix(true);

		// Incorporate the affine transformation into the spline coefficient
		for (int i = 0; i < intervals + 3; i++) {
			final double v = (double) ((i - 1) * (sourceCurrentHeight - 1)) / (double) intervals;
			final double xv = affineMatrix[0][2] + affineMatrix[0][1] * v;
			final double yv = affineMatrix[1][2] + affineMatrix[1][1] * v;
			for (int j = 0; j < intervals + 3; j++) {
				final double u = (double) ((j - 1) * (sourceCurrentWidth - 1)) / (double) intervals;
				cxSourceToTarget[i][j] = xv + affineMatrix[0][0] * u;
				cySourceToTarget[i][j] = yv + affineMatrix[1][0] * u;
			}
		}

		// Now refine with the different scales
		int state; // state=-1 --> Finish
		// state= 0 --> Increase deformation detail
		// state= 1 --> Increase image detail
		// state= 2 --> Do nothing until the finest image scale
		if (minScaleDeformation == maxScaleDeformation)
			state = 1;
		else
			state = 0;
		int s = minScaleDeformation;
		int step = 0;
		computeTotalWorkload();

		while (state != -1) {
			int currentDepth = targetModel.getCurrentDepth();

			// Update the deformation coefficients only in states 0 and 1
			if (state == 0 || state == 1) {
				// Update the deformation coefficients with the error of the landmarks
				// The following conditional is now useless but it is there to allow
				// easy changes like applying the landmarks only in the coarsest
				// deformation
				if (s >= minScaleDeformation) {
					// Number of intervals at this scale and ask for memory
					intervals = (int) Math.pow(2, s + sizeCorrectionFactor);
					final double[][] newcxTargetToSource = new double[intervals + 3][intervals + 3];
					final double[][] newcyTargetToSource = new double[intervals + 3][intervals + 3];

					final double[][] newcxSourceToTarget = new double[intervals + 3][intervals + 3];
					final double[][] newcySourceToTarget = new double[intervals + 3][intervals + 3];

					// Compute the coefficients at this scale
					boolean underconstrained = true;
					// FROM TARGET TO SOURCE.
					if (divWeight == 0 && curlWeight == 0)
						underconstrained = computeCoefficientsScale(intervals, dxTargetToSource, dyTargetToSource,
						    newcxTargetToSource, newcyTargetToSource, false);
					else
						underconstrained = computeCoefficientsScaleWithRegularization(intervals, dxTargetToSource, dyTargetToSource,
						    newcxTargetToSource, newcyTargetToSource, false);

					// Incorporate information from the previous scale
					if (!underconstrained || (step == 0 && landmarkWeight != 0)) {
						for (int i = 0; i < intervals + 3; i++)
							for (int j = 0; j < intervals + 3; j++) {
								cxTargetToSource[i][j] += newcxTargetToSource[i][j];
								cyTargetToSource[i][j] += newcyTargetToSource[i][j];
							}
					}

					// FROM SOURCE TO TARGET.
					underconstrained = true;
					if (divWeight == 0 && curlWeight == 0)
						underconstrained = computeCoefficientsScale(intervals, dxSourceToTarget, dySourceToTarget,
						    newcxSourceToTarget, newcySourceToTarget, true);
					else
						underconstrained = computeCoefficientsScaleWithRegularization(intervals, dxSourceToTarget, dySourceToTarget,
						    newcxSourceToTarget, newcySourceToTarget, true);

					// Incorporate information from the previous scale
					if (!underconstrained || (step == 0 && landmarkWeight != 0)) {
						for (int i = 0; i < intervals + 3; i++)
							for (int j = 0; j < intervals + 3; j++) {
								cxSourceToTarget[i][j] += newcxSourceToTarget[i][j];
								cySourceToTarget[i][j] += newcySourceToTarget[i][j];
							}
					}
				}

				// Optimize deformation coefficients
				// if (imageWeight!=0)
				optimizeCoeffs(intervals, stopThreshold, cxTargetToSource, cyTargetToSource, cxSourceToTarget,
				    cySourceToTarget);
			}

			// Prepare for next iteration
			step++;
			switch (state) {
			case 0:
				// Finer details in the deformation
				if (s < maxScaleDeformation) {
					cxTargetToSource = propagateCoeffsToNextLevel(intervals, cxTargetToSource, 1);
					cyTargetToSource = propagateCoeffsToNextLevel(intervals, cyTargetToSource, 1);
					cxSourceToTarget = propagateCoeffsToNextLevel(intervals, cxSourceToTarget, 1);
					cySourceToTarget = propagateCoeffsToNextLevel(intervals, cySourceToTarget, 1);
					s++;
					intervals *= 2;

					// Prepare matrices for the regularization term
					buildRegularizationTemporary(intervals, false);
					buildRegularizationTemporary(intervals, true);

					if (currentDepth > minScaleImage)
						state = 1;
					else
						state = 0;
				} else if (currentDepth > minScaleImage)
					state = 1;
				else
					state = 2;
				break;
			case 1: // Finer details in the image, go on optimizing
			case 2: // Finer details in the image, do not optimize
				// Compute next state
				if (state == 1) {
					if (s == maxScaleDeformation && currentDepth == minScaleImage)
						state = 2;
					else if (s == maxScaleDeformation)
						state = 1;
					else
						state = 0;
				} else if (state == 2) {
					if (currentDepth == 0)
						state = -1;
					else
						state = 2;
				}

				// Pop another image and prepare the deformation
				if (currentDepth != 0) {
					double oldTargetCurrentHeight = targetCurrentHeight;
					double oldTargetCurrentWidth = targetCurrentWidth;
					double oldSourceCurrentHeight = sourceCurrentHeight;
					double oldSourceCurrentWidth = sourceCurrentWidth;

					sourceModel.popFromPyramid();
					targetModel.popFromPyramid();

					targetCurrentHeight = targetModel.getCurrentHeight();
					targetCurrentWidth = targetModel.getCurrentWidth();
					targetFactorHeight = targetModel.getFactorHeight();
					targetFactorWidth = targetModel.getFactorWidth();

					sourceCurrentHeight = sourceModel.getCurrentHeight();
					sourceCurrentWidth = sourceModel.getCurrentWidth();
					sourceFactorHeight = sourceModel.getFactorHeight();
					sourceFactorWidth = sourceModel.getFactorWidth();

					// Adapt the transformation to the new image size
					double targetFactorY = (targetCurrentHeight - 1) / (oldTargetCurrentHeight - 1);
					double targetFactorX = (targetCurrentWidth - 1) / (oldTargetCurrentWidth - 1);
					double sourceFactorY = (sourceCurrentHeight - 1) / (oldSourceCurrentHeight - 1);
					double sourceFactorX = (sourceCurrentWidth - 1) / (oldSourceCurrentWidth - 1);

					for (int i = 0; i < intervals + 3; i++)
						for (int j = 0; j < intervals + 3; j++) {
							cxTargetToSource[i][j] *= targetFactorX;
							cyTargetToSource[i][j] *= targetFactorY;
							cxSourceToTarget[i][j] *= sourceFactorX;
							cySourceToTarget[i][j] *= sourceFactorY;
						}

					// Prepare matrices for the regularization term
					buildRegularizationTemporary(intervals, false);
					buildRegularizationTemporary(intervals, true);
				}
				break;
			}

			// In accurate_mode reduce the stopping threshold for the last iteration
			if ((state == 0 || state == 1) && s == maxScaleDeformation && currentDepth == minScaleImage + 1
			    && accurateMode == 1)
				stopThreshold /= 10;

		} // end while (state != -1).

		// Adapt coefficients if necessary
		if (sourceModel.getOriginalImageWidth() > this.sourceCurrentWidth) {
			if (sourceModel.isSubOutput() || targetModel.isSubOutput())
				System.out.println("Adapting coefficients from " + this.sourceCurrentWidth + " to "
				    + sourceModel.getOriginalImageWidth() + "...");
			// Adapt the transformation to the new image size
			double targetFactorY = (targetModel.getOriginalImageHeight() - 1) / (targetCurrentHeight - 1);
			double targetFactorX = (targetModel.getOriginalImageWidth() - 1) / (targetCurrentWidth - 1);
			double sourceFactorY = (sourceModel.getOriginalImageHeight() - 1) / (sourceCurrentHeight - 1);
			double sourceFactorX = (sourceModel.getOriginalImageWidth() - 1) / (sourceCurrentWidth - 1);

			for (int i = 0; i < intervals + 3; i++)
				for (int j = 0; j < intervals + 3; j++) {
					cxTargetToSource[i][j] *= targetFactorX;
					cyTargetToSource[i][j] *= targetFactorY;
					cxSourceToTarget[i][j] *= sourceFactorX;
					cySourceToTarget[i][j] *= sourceFactorY;
				}
			this.targetCurrentHeight = targetModel.getOriginalImageHeight();
			this.targetCurrentWidth = targetModel.getOriginalImageWidth();
			this.sourceCurrentHeight = sourceModel.getOriginalImageHeight();
			this.sourceCurrentWidth = sourceModel.getOriginalImageWidth();
		}

		// Display final errors.
		if (this.outputLevel == 2) {
			if (this.imageWeight != 0) {
				System.out.println(" Optimal direct similarity error = " + this.finalDirectSimilarityError);
				System.out.println(" Optimal inverse similarity error = " + this.finalInverseSimilarityError);
			}
			if (this.curlWeight != 0 || this.divWeight != 0) {
				System.out.println(" Optimal direct regularization error = " + this.finalDirectRegularizationError);
				System.out.println(" Optimal inverse regularization error = " + this.finalInverseRegularizationError);
			}
			if (this.landmarkWeight != 0) {
				System.out.println(" Optimal direct landmark error = " + this.finalDirectLandmarkError);
				System.out.println(" Optimal inverse landmark error = " + this.finalInverseLandmarkError);
			}
			if (this.consistencyWeight != 0) {
				System.out.println(" Optimal direct consistency error = " + this.finalDirectConsistencyError);
				System.out.println(" Optimal inverse consistency error = " + this.finalInverseConsistencyError);
			}
		}

	}

	/**
	 * Show the inverse transformation results (multi-thread version).
	 */
	public void showInverseResults() {
		showTransformationMultiThread(intervals, cxSourceToTarget, cySourceToTarget, true);
	}

	public void getRegisteredSource(Sequence srcTgtSeq) {
		applyTransformationMT(srcTgtSeq, targetSeq, intervals, cxTargetToSource, cyTargetToSource);
	}

	public void getRegisteredTarget(Sequence tgtTgtSeq) {
		applyTransformationMT(tgtTgtSeq, sourceSeq, intervals, cxSourceToTarget, cySourceToTarget);
	}

	private void applyTransformationMT(Sequence source, Sequence target, int intervals, double[][] cx, double[][] cy) {
		// Apply transformation
		MiscTools.applyTransformationToSourceMT(source, target, intervals, cx, cy);
	}

	public Sequence getRegisteredSource(String srcResultPath, String srcPath, String transformedSrcPath, String tgtPath) {
		return BigImageTools.applyTransformationToImage(srcResultPath, srcPath, transformedSrcPath, tgtPath, intervals,
		    cxTargetToSource, cyTargetToSource, new Dimension(targetWidth, targetHeight));
	}

	public Sequence getRegisteredTarget(String tgtResultPath, String srcPath, String tgtPath, String transformedTgtPath) {
		return BigImageTools.applyTransformationToImage(tgtResultPath, tgtPath, transformedTgtPath, srcPath, intervals,
		    cxSourceToTarget, cySourceToTarget, new Dimension(sourceWidth, sourceHeight));
	}

	public double[][] getCxSourceToTarget() {
		return cxSourceToTarget;
	}

	public double[][] getCySourceToTarget() {
		return cySourceToTarget;
	}

	public double[][] getCxTargetToSource() {
		return cxTargetToSource;
	}

	public double[][] getCyTargetToSource() {
		return cyTargetToSource;
	}

	public int getIntervals() {
		return intervals;
	}

	public void saveBigRegisteredSource(String srcResultPath, String transformedSrcResultPath, String srcPath,
	    String transformedSrcPath, String tgtPath, Rectangle tile)
	    throws ServiceException, IOException, FormatException, InterruptedException {
		BigImageTools.applyAndSaveTransformationToBigImage(srcResultPath, transformedSrcResultPath, srcPath,
		    transformedSrcPath, tgtPath, intervals, cxTargetToSource, cyTargetToSource,
		    new Dimension(targetWidth, targetHeight), plugin, tile);
	}

	public void saveBigRegisteredTarget(String tgtResultPath, String transformedTgtResultPath, String tgtPath,
	    String transformedTgtPath, String srcPath, Rectangle tile)
	    throws ServiceException, IOException, FormatException, InterruptedException {
		BigImageTools.applyAndSaveTransformationToBigImage(tgtResultPath, transformedTgtResultPath, tgtPath,
		    transformedTgtPath, srcPath, intervals, cxSourceToTarget, cySourceToTarget,
		    new Dimension(targetWidth, targetHeight), plugin, tile);
	}


}
