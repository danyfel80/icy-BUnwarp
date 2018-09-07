package algorithms.danyfel80.registration.bunwarp;

import java.awt.Rectangle;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;

import icy.common.listener.DetailedProgressListener;
import icy.image.IcyBufferedImage;
import icy.main.Icy;
import icy.roi.ROI2D;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
//import plugins.danyfel80.registration.bunwarp.BUnwarpSimple;
import plugins.kernel.roi.roi2d.ROI2DPoint;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class BUnwarpper implements Runnable {
	// Images
	/** image representation for the source */
	private Sequence sourceSeq;
	/** image representation for the target */
	private Sequence targetSeq;
	/** source image model */
	private BSplineModel sourceModel;
	/** target image model */
	private BSplineModel targetModel;

	// Landmarks
	/** point handler for the landmarks in the source image */
	private List<ROI2DPoint> sourceLandmarks;
	/** point handler for the landmarks in the target image */
	private List<ROI2DPoint> targetLandmarks;

	// Masks for the images
	/** source image mask */
	private ROI2D sourceMask;
	/** target image mask */
	private ROI2D targetMask;

	// Transformation parameters
	/** maximum image subsampling factor */
	private int maxImageSubsamplingFactor;
	/** minimum scale deformation */
	private int minScaleDeformation;
	/** maximum scale deformation */
	private int maxScaleDeformation;
	/** minimum image scale */
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
	/** maximum depth for the image pyramid */
	private int imagePyramidDepth;

	/** warp transformation */
	private Transformation warp;

	private DetailedProgressListener progressListener;

	/*
	 * List<ROI2DPoint> srcLandmarks, List<ROI2DPoint> tgtLandmarks, ROI2DPolygon
	 * srcMask, ROI2DPolygon tgtMask, RegistrationModeEnum mode, Integer
	 * maxSubsamplingFactor, MinimumScaleDeformationEnum minScaleDef,
	 * MaximumScaleDeformationEnum maxScaleDef, Double divWeight, Double
	 * curlWeight, Double landmarkWeight, Double consistencyWeight, Double
	 * imageWeight, Double stopThreshold, Boolean showProcess, EzPlug plugin
	 */
	public BUnwarpper(final Sequence sourceSequence, final Sequence targetSequence,
			final List<ROI2DPoint> sourceLandmarks, final List<ROI2DPoint> targetLandmarks, final ROI2D sourceMask,
			final ROI2D targetMask, final int maxImageSubsamplingFactor, final int minScaleDeformation,
			final int maxScaleDeformation, final int minScaleImage, final double divWeight, final double curlWeight,
			final double landmarkWeight, final double imageWeight, final double consistencyWeight, final double stopThreshold,
			final int outputLevel, final boolean showMarquardtOptim, final int accurateMode,
			DetailedProgressListener progressListener) {
		this.sourceSeq = sourceSequence;
		this.targetSeq = targetSequence;
		this.sourceLandmarks = sourceLandmarks;
		this.targetLandmarks = targetLandmarks;
		this.sourceMask = sourceMask;
		this.targetMask = targetMask;
		this.maxImageSubsamplingFactor = maxImageSubsamplingFactor;
		this.minScaleDeformation = minScaleDeformation;
		this.maxScaleDeformation = maxScaleDeformation;
		this.minScaleImage = minScaleImage;
		this.divWeight = divWeight;
		this.curlWeight = curlWeight;
		this.landmarkWeight = landmarkWeight;
		this.imageWeight = imageWeight;
		this.consistencyWeight = consistencyWeight;
		this.stopThreshold = stopThreshold;
		this.outputLevel = outputLevel;
		this.showMarquardtOptim = showMarquardtOptim;
		this.accurateMode = accurateMode;

		setProgressListener(progressListener);

		createSourceImage(this.accurateMode < RegistrationModeEnum.MONO.getNumber());
		createTargetImage();
	}

	public void setProgressListener(DetailedProgressListener listener) {
		this.progressListener = listener;
	}

	protected void notifyProgress(double progress, String message) {
		if (this.progressListener != null) {
			this.progressListener.notifyProgress(progress, message, null);
		}
	}

	private void createSourceImage(boolean isReverse) {
		sourceModel = new BSplineModel(sourceSeq.getFirstImage(), isReverse, (int) Math.pow(2, maxImageSubsamplingFactor),
				(progress, message, data) -> {
					notifyProgress(progress, "Creating source image model: ");
					return false;
				});
		computeImagePyramidDepth();
		sourceModel.setPyramidDepth(imagePyramidDepth + minScaleImage);
	}

	private void computeImagePyramidDepth() {
		imagePyramidDepth = maxScaleDeformation - minScaleDeformation + 1;
	}

	private void createTargetImage() {
		targetModel = new BSplineModel(targetSeq.getFirstImage(), true, (int) Math.pow(2, maxImageSubsamplingFactor),
				(progress, message, data) -> {
					notifyProgress(progress, "Creating target image model: ");
					return false;
				});
		computeImagePyramidDepth();
		targetModel.setPyramidDepth(imagePyramidDepth + minScaleImage);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Thread#run()
	 */
	@Override
	public void run() {

		// Start pyramids
		notifyProgress(0.001d, "Starting image pyramids...");

		if (targetModel.getWidth() > BSplineModel.MAX_OUTPUT_SIZE || targetModel.getHeight() > BSplineModel.MAX_OUTPUT_SIZE
				|| sourceModel.getWidth() > BSplineModel.MAX_OUTPUT_SIZE
				|| sourceModel.getHeight() > BSplineModel.MAX_OUTPUT_SIZE)
			System.out.println("Starting image pyramids...");

		sourceModel.startPyramids();
		targetModel.startPyramids();

		try {
			sourceModel.join();
			targetModel.join();
		} catch (InterruptedException e) {
			System.err.println("Unhandled interruption: " + e);
		}

		// Create output image (source-target)
		final Sequence[] outputSeqs = initializeOutputSeqs();

		// If mono mode, reset consistency weight
		if (this.accurateMode == RegistrationModeEnum.MONO.getNumber())
			this.consistencyWeight = 0.0;

		// IJ.log("FinalAction: maxImageSubsamplingFactor = " +
		// maxImageSubsamplingFactor);

		// Prepare registration parameters
		Sequence srcSeq = SequenceUtil.getCopy(sourceSeq);
		Sequence tgtSeq = SequenceUtil.getCopy(targetSeq);
		if (!Icy.getMainInterface().isHeadLess()) {
			Icy.getMainInterface().addSequence(srcSeq);
			Icy.getMainInterface().addSequence(tgtSeq);
		}
		warp = new Transformation(srcSeq, tgtSeq, sourceModel, targetModel, sourceLandmarks, targetLandmarks, sourceMask,
				targetMask, minScaleDeformation, maxScaleDeformation, minScaleImage, divWeight, curlWeight, landmarkWeight,
				imageWeight, consistencyWeight, stopThreshold, outputLevel, showMarquardtOptim, accurateMode, outputSeqs[0],
				outputSeqs[1], (progress, message, data) -> {
					notifyProgress(progress, "Transformation: " + message);
					return false;
				});

		// Perform the registration
		notifyProgress(0.01d, "Registering...");

		long start = System.currentTimeMillis(); // start timing

		if (accurateMode == RegistrationModeEnum.MONO.getNumber()) {
			// Do unidirectional registration
			warp.doUnidirectionalRegistration();
			// Show results.
			if (sourceModel.isSubOutput()) {
				System.out.println("Calculating final transformed source image");
			}
			warp.showDirectResults();

		} else {
			// Do bidirectional registration
			warp.doBidirectionalRegistration();
			if (sourceModel.isSubOutput()) {
				System.out.println("Calculating final transformed source image");
			}
			warp.showDirectResults();
			if (targetModel.isSubOutput()) {
				System.out.println("Calculating final transformed target image");
			}
			warp.showInverseResults();
		}

		System.out.println("Intervals: " + warp.getIntervals());

		long stop = System.currentTimeMillis(); // stop timing
		if (outputLevel == 2)
			System.out.println("\nRegistration time: " + (stop - start) + "ms"); // print

		if (!Icy.getMainInterface().isHeadLess()) {
			Icy.getMainInterface().closeSequence(srcSeq);
			Icy.getMainInterface().closeSequence(tgtSeq);
		}
	}

	private Sequence[] initializeOutputSeqs() {
		int Xdimt = targetModel.getWidth();
		int Ydimt = targetModel.getHeight();
		int Xdims = sourceModel.getWidth();
		int Ydims = sourceModel.getHeight();
		double[] tImage = targetModel.isSubOutput() ? targetModel.getSubImage() : targetModel.getImage();
		double[] sImage = sourceModel.isSubOutput() ? sourceModel.getSubImage() : sourceModel.getImage();
		int sSubFactorX = 1;
		int sSubFactorY = 1;
		int tSubFactorX = 1;
		int tSubFactorY = 1;
		Sequence[] outputSeqs = new Sequence[2];

		String extraTitleS = "";
		String extraTitleT = "";

		if (targetModel.isSubOutput() || sourceModel.isSubOutput())
			System.out.println("Initializing output windows...");

		// If the output (difference) images are subsampled (because they were
		// larger than the maximum size), update variables.
		if (targetModel.isSubOutput()) {
			tSubFactorX = Xdimt / targetModel.getSubWidth();
			tSubFactorY = Ydimt / targetModel.getSubHeight();
			extraTitleT = " (Subsampled)";
			Xdimt = targetModel.getSubWidth();
			Ydimt = targetModel.getSubHeight();
		}

		if (sourceModel.isSubOutput()) {
			sSubFactorX = Xdims / sourceModel.getSubWidth();
			sSubFactorY = Ydims / sourceModel.getSubHeight();
			extraTitleS = " (Subsampled)";
			Xdims = sourceModel.getSubWidth();
			Ydims = sourceModel.getSubHeight();
		}

		// Float processor for the output source-target image.
		final IcyBufferedImage ibi = new IcyBufferedImage(Xdimt, Ydimt, 1, DataType.FLOAT);
		float[] ibiData = ibi.getDataXYAsFloat(0);

		for (int i = 0; i < Ydimt; i++) {
			final int i_offset_t = i * Xdimt;
			final int i_offset_s = i * Xdims;
			final int i_s_sub = i * sSubFactorY;
			final int i_t_sub = i * tSubFactorY;

			for (int j = 0; j < Xdimt; j++) {

				if ((sourceMask == null || targetMask == null
						|| (sourceMask.contains(j * sSubFactorX, i_s_sub) && targetMask.contains(j * tSubFactorX, i_t_sub)))
						&& j < Xdims && i < Ydims)
					ibiData[j + i_offset_t] = (float) (tImage[i_offset_t + j] - sImage[i_offset_s + j]);
				else {
					ibiData[j + i_offset_t] = 0;
				}

			}
		}
		ibi.dataChanged();

		final Sequence seq1 = new Sequence("Output Source-Target" + extraTitleS, ibi);
		if (!Icy.getMainInterface().isHeadLess()) {
			Icy.getMainInterface().addSequence(seq1);
		}

		outputSeqs[0] = seq1;

		// Create output image (target-source) if necessary

		if (this.accurateMode != RegistrationModeEnum.MONO.getNumber()) {
			final IcyBufferedImage ibi2 = new IcyBufferedImage(Xdims, Ydims, 1, DataType.FLOAT);
			float[] ibi2Data = ibi2.getDataXYAsFloat(0);

			for (int i = 0; i < Ydims; i++) {
				int i_offset_t = i * Xdimt;
				int i_offset_s = i * Xdims;
				int i_s_sub = i * sSubFactorY;
				int i_t_sub = i * tSubFactorY;

				for (int j = 0; j < Xdims; j++)
					if ((targetMask == null || sourceMask != null
							|| (targetMask.contains(j * tSubFactorX, i_t_sub) && sourceMask.contains(j * sSubFactorX, i_s_sub)))
							&& i < Ydimt && j < Xdimt)
						ibi2Data[j + i_offset_s] = (float) (sImage[i_offset_s + j] - tImage[i_offset_t + j]);
					else
						ibi2Data[j + i_offset_s] = 0;
			}
			ibi2.dataChanged();

			final Sequence seq2 = new Sequence("Output Target-Source" + extraTitleT, ibi2);
			if (!Icy.getMainInterface().isHeadLess()) {
				Icy.getMainInterface().addSequence(seq2);
			}

			outputSeqs[1] = seq2;
		} else
			outputSeqs[1] = null;

		return outputSeqs;
	}

	public void getRegisteredSource(Sequence srcTgtSeq) {
		warp.getRegisteredSource(srcTgtSeq);
	}

	public void getRegisteredTarget(Sequence tgtTgtSeq) {
		warp.getRegisteredTarget(tgtTgtSeq);
	}

	/**
	 * Computes the registered image in the specified result file path using
	 * srcPath as the transformed image.
	 * 
	 * @param srcResultPath
	 *          Path of the result image
	 * @param srcPath
	 *          Path of the image to be transformed
	 * @param tgtPath
	 *          Path of the base image
	 * @return The resulting sequence
	 */
	public Sequence getRegisteredSource(String srcResultPath, String srcPath, String transformedSrcPath, String tgtPath) {
		return warp.getRegisteredSource(srcResultPath, srcPath, transformedSrcPath, tgtPath);
	}

	/**
	 * Computes the registered image in the specified result file path using
	 * tgtPath as the transformed image.
	 * 
	 * @param tgtResultPath
	 *          Path of the result image
	 * @param srcPath
	 *          Path of the image to be transformed
	 * @param tgtPath
	 *          Path of the base image
	 * @return The resulting sequence
	 */
	public Sequence getRegisteredTarget(String tgtResultPath, String srcPath, String tgtPath, String transformedTgtPath) {
		return warp.getRegisteredTarget(tgtResultPath, srcPath, tgtPath, transformedTgtPath);
	}

	public double[][] getCxSourceToTarget() {
		return warp.getCxSourceToTarget();
	}

	public double[][] getCySourceToTarget() {
		return warp.getCySourceToTarget();
	}

	public double[][] getCxTargetToSource() {
		return warp.getCxTargetToSource();
	}

	public double[][] getCyTargetToSource() {
		return warp.getCyTargetToSource();
	}

	public int getIntervals() {
		return warp.getIntervals();
	}

	public void saveBigRegisteredSource(String srcResultPath, String transformedSrcResultPath, String srcPath,
			String transformedSrcPath, String tgtPath, Rectangle tile)
			throws ServiceException, IOException, FormatException, InterruptedException, ExecutionException {
		warp.saveBigRegisteredSource(srcResultPath, transformedSrcResultPath, srcPath, transformedSrcPath, tgtPath, tile);
	}

	public void saveBigRegisteredTarget(String tgtResultPath, String transformedTgtResultPath, String tgtPath,
			String transformedTgtPath, String srcPath, Rectangle tile)
			throws ServiceException, IOException, FormatException, InterruptedException, ExecutionException {
		warp.saveBigRegisteredTarget(tgtResultPath, transformedTgtResultPath, tgtPath, transformedTgtPath, srcPath, tile);
	}

}
