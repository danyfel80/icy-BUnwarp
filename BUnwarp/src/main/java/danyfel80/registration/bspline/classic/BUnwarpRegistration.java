package danyfel80.registration.bspline.classic;

import java.awt.geom.Point2D;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import algorithms.danyfel80.io.sequence.cursor.IcyBufferedImageCursor;
import icy.common.listener.DetailedProgressListener;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import plugins.kernel.roi.roi2d.ROI2DArea;

public class BUnwarpRegistration {
	/* Internal Variables */
	private Sequence sourceSequence;
	private Sequence targetSequence;
	private ROI2DArea sourceMask;
	private ROI2DArea targetMask;
	private List<Point2D> sourceLandmarks;
	private List<Point2D> targetLandmarks;
	private RegistrationMode registrationMode;
	private int initialSubsampleFactor;
	private Sequence transformedSourceSequence;
	private Sequence transformedTargetSequence;
	private DeformationScale initialDeformationScale;
	private DeformationScale finalDeformationScale;
	private double divWeight;
	private double curlWeight;
	private double landmarkWeight;
	private double imageWeight;
	private double consistencyWeight;
	private double stopThreshold;
	private boolean showProcess;
	private Set<DetailedProgressListener> progressListeners;
	private Set<DetailedProgressListener> progressOutputListeners;

	private IcyBufferedImage originalSourceImage;
	private IcyBufferedImage originalTargetImage;
	private BSplineModel sourceBSplineModel;
	private BSplineModel targetBSplineModel;

	private Sequence sourceProgressOutput;
	private Sequence targetProgressOutput;
	private Transformation transformation;
	private Sequence directResult;
	private Sequence indirectResult;

	/**
	 * @return The sourceSequence.
	 */
	public Sequence getSourceSequence() {
		return sourceSequence;
	}

	/**
	 * @param sourceSequence
	 *          The sourceSequence to set.
	 */
	public void setSourceSequence(Sequence sourceSequence) {
		this.sourceSequence = sourceSequence;
	}

	/**
	 * @return The targetSequence.
	 */
	public Sequence getTargetSequence() {
		return targetSequence;
	}

	/**
	 * @param targetSequence
	 *          The targetSequence to set.
	 */
	public void setTargetSequence(Sequence targetSequence) {
		this.targetSequence = targetSequence;
	}

	/**
	 * @return The source mask.
	 */
	public ROI2DArea getSourceMask() {
		return sourceMask;
	}

	/**
	 * @param mask
	 *          The source mask.
	 */
	public void setSourceMask(ROI2DArea mask) {
		this.sourceMask = mask;
	}

	/**
	 * @return The target mask.
	 */
	public ROI2DArea getTargetMask() {
		return targetMask;
	}

	/**
	 * @param mask
	 *          The target mask
	 */
	public void setTargetMask(ROI2DArea mask) {
		this.targetMask = mask;
	}

	/**
	 * @return The sourceLandmarks.
	 */
	public List<Point2D> getSourceLandmarks() {
		return sourceLandmarks;
	}

	/**
	 * @param sourceLandmarks
	 *          The sourceLandmarks to set.
	 */
	public void setSourceLandmarks(List<Point2D> sourceLandmarks) {
		this.sourceLandmarks = sourceLandmarks;
	}

	/**
	 * @return The targetLandmarks.
	 */
	public List<Point2D> getTargetLandmarks() {
		return targetLandmarks;
	}

	/**
	 * @param targetLandmarks
	 *          The targetLandmarks to set.
	 */
	public void setTargetLandmarks(List<Point2D> targetLandmarks) {
		this.targetLandmarks = targetLandmarks;
	}

	/**
	 * @return The registrationMode.
	 */
	public RegistrationMode getRegistrationMode() {
		return registrationMode;
	}

	/**
	 * @param registrationMode
	 *          The registrationMode to set.
	 */
	public void setRegistrationMode(RegistrationMode registrationMode) {
		this.registrationMode = registrationMode;
	}

	/**
	 * @return The initialSubsampleFactor.
	 */
	public int getInitialSubsampleFactor() {
		return initialSubsampleFactor;
	}

	/**
	 * @param initialSubsampleFactor
	 *          The initialSubsampleFactor to set.
	 */
	public void setInitialSubsampleFactor(int initialSubsampleFactor) {
		this.initialSubsampleFactor = initialSubsampleFactor;
	}

	/**
	 * @return The transformedSourceSequence.
	 */
	public Sequence getTransformedSourceSequence() {
		return transformedSourceSequence;
	}

	/**
	 * @param transformedSourceSequence
	 *          The transformedSourceSequence to set.
	 */
	public void setTransformedSourceSequence(Sequence transformedSourceSequence) {
		this.transformedSourceSequence = transformedSourceSequence;
	}

	/**
	 * @return The transformedTargetSequence.
	 */
	public Sequence getTransformedTargetSequence() {
		return transformedTargetSequence;
	}

	/**
	 * @param transformedTargetSequence
	 *          The transformedTargetSequence to set.
	 */
	public void setTransformedTargetSequence(Sequence transformedTargetSequence) {
		this.transformedTargetSequence = transformedTargetSequence;
	}

	/**
	 * @return The initialDeformationScale.
	 */
	public DeformationScale getInitialDeformationScale() {
		return initialDeformationScale;
	}

	/**
	 * @param initialDeformationScale
	 *          The initialDeformationScale to set.
	 */
	public void setInitialDeformationScale(DeformationScale initialDeformationScale) {
		this.initialDeformationScale = initialDeformationScale;
	}

	/**
	 * @return The finalDeformationScale.
	 */
	public DeformationScale getFinalDeformationScale() {
		return finalDeformationScale;
	}

	/**
	 * @param finalDeformationScale
	 *          The finalDeformationScale to set.
	 */
	public void setFinalDeformationScale(DeformationScale finalDeformationScale) {
		this.finalDeformationScale = finalDeformationScale;
	}

	/**
	 * @return The divWeight.
	 */
	public double getDivWeight() {
		return divWeight;
	}

	/**
	 * @param divWeight
	 *          The divWeight to set.
	 */
	public void setDivWeight(double divWeight) {
		this.divWeight = divWeight;
	}

	/**
	 * @return The curlWeight.
	 */
	public double getCurlWeight() {
		return curlWeight;
	}

	/**
	 * @param curlWeight
	 *          The curlWeight to set.
	 */
	public void setCurlWeight(double curlWeight) {
		this.curlWeight = curlWeight;
	}

	/**
	 * @return The landmarkWeight.
	 */
	public double getLandmarkWeight() {
		return landmarkWeight;
	}

	/**
	 * @param landmarkWeight
	 *          The landmarkWeight to set.
	 */
	public void setLandmarkWeight(double landmarkWeight) {
		this.landmarkWeight = landmarkWeight;
	}

	/**
	 * @return The imageWeight.
	 */
	public double getImageWeight() {
		return imageWeight;
	}

	/**
	 * @param imageWeight
	 *          The imageWeight to set.
	 */
	public void setImageWeight(double imageWeight) {
		this.imageWeight = imageWeight;
	}

	/**
	 * @return The consistencyWeight.
	 */
	public double getConsistencyWeight() {
		return consistencyWeight;
	}

	/**
	 * @param consistencyWeight
	 *          The consistencyWeight to set.
	 */
	public void setConsistencyWeight(double consistencyWeight) {
		this.consistencyWeight = consistencyWeight;
	}

	/**
	 * @return The stopThreshold.
	 */
	public double getStopThreshold() {
		return stopThreshold;
	}

	/**
	 * @param stopThreshold
	 *          The stopThreshold to set.
	 */
	public void setStopThreshold(double stopThreshold) {
		this.stopThreshold = stopThreshold;
	}

	/**
	 * @return The showProcess.
	 */
	public boolean isShowProcess() {
		return showProcess;
	}

	/**
	 * @param showProcess
	 *          The showProcess to set.
	 */
	public void setShowProcess(boolean showProcess) {
		this.showProcess = showProcess;
	}

	public void addProgressListener(DetailedProgressListener listener) {
		this.progressListeners.add(listener);
	}

	public void removeProgressListener(DetailedProgressListener listener) {
		this.progressListeners.remove(listener);
	}

	private void notifyProgress(double progress, String message) {
		for (DetailedProgressListener listener: progressListeners) {
			listener.notifyProgress(progress, message, null);
		}
	}

	public void addProgressOutputListener(DetailedProgressListener listener) {
		this.progressOutputListeners.add(listener);
	}

	public void removeProgressOutputListener(DetailedProgressListener listener) {
		this.progressOutputListeners.remove(listener);
	}

	public Sequence getDirectResult() {
		return directResult;
	}

	public Sequence getIndirectResult() {
		return indirectResult;
	}

	public BUnwarpRegistration() {
		progressListeners = new HashSet<>();
		progressOutputListeners = new HashSet<>();
	}

	public void compute() throws InterruptedException {
		copyOriginalImages();
		createImageBSplineModels();
		notifyProgress(Double.NaN, "Starting image pyramids...");
		computePyramids();
		if (isShowProcess()) {
			initializeProgressOutput();
		}
		adjustConsistencyWeight();
		createTransformation();
		computeTransformation();
		restablishOriginalImages();
	}

	private void copyOriginalImages() {
		this.originalSourceImage = sourceSequence.getFirstImage();
		this.originalTargetImage = targetSequence.getFirstImage();
	}

	private void createImageBSplineModels() {
		this.sourceBSplineModel = new BSplineModel(sourceSequence.getFirstImage(),
				getRegistrationMode() != RegistrationMode.MONO, (int) Math.pow(2, getInitialSubsampleFactor()));
		sourceBSplineModel.setPyramidDepth(getPyramidDepth());

		this.targetBSplineModel = new BSplineModel(targetSequence.getFirstImage(), true,
				(int) Math.pow(2, getInitialSubsampleFactor()));
		targetBSplineModel.setPyramidDepth(getPyramidDepth());
	}

	private int getPyramidDepth() {
		return getFinalDeformationScale().getNumber() - getInitialDeformationScale().getNumber();
	}

	private void computePyramids() throws InterruptedException {
		sourceBSplineModel.startPyramids();
		targetBSplineModel.startPyramids();
		sourceBSplineModel.getThread().join();
		targetBSplineModel.getThread().join();
	}

	private void initializeProgressOutput() {
		initializeSourceProgressOutput();
		initializeTargetProgressOutput();

	}

	private void initializeSourceProgressOutput() {

		int Ydimt = targetBSplineModel.getHeight();
		int Xdimt = targetBSplineModel.getWidth();
		int Xdims = sourceBSplineModel.getWidth();
		int Ydims = sourceBSplineModel.getHeight();
		double[] tImage = targetBSplineModel.isSubOutput()? targetBSplineModel.getSubImage(): targetBSplineModel.getImage();
		double[] sImage = sourceBSplineModel.isSubOutput()? sourceBSplineModel.getSubImage(): sourceBSplineModel.getImage();
		int sSubFactorX = 1;
		int sSubFactorY = 1;
		int tSubFactorX = 1;
		int tSubFactorY = 1;

		String extraTitleS = "";

		if (targetBSplineModel.isSubOutput() || sourceBSplineModel.isSubOutput())
			notifyProgress(Double.NaN, "Initializing source -> target window...");

		// If the output (difference) images are subsampled (because they were
		// larger than the maximum size), update variables.
		if (targetBSplineModel.isSubOutput()) {
			tSubFactorX = Xdimt / targetBSplineModel.getSubWidth();
			tSubFactorY = Ydimt / targetBSplineModel.getSubHeight();
			Xdimt = targetBSplineModel.getSubWidth();
			Ydimt = targetBSplineModel.getSubHeight();
		}

		if (sourceBSplineModel.isSubOutput()) {
			sSubFactorX = Xdims / sourceBSplineModel.getSubWidth();
			sSubFactorY = Ydims / sourceBSplineModel.getSubHeight();
			extraTitleS = " (Subsampled)";
			Xdims = sourceBSplineModel.getSubWidth();
			Ydims = sourceBSplineModel.getSubHeight();
		}

		// Float processor for the output source-target image.
		final IcyBufferedImage fp = new IcyBufferedImage(Xdimt, Ydimt, 1, DataType.DOUBLE);
		IcyBufferedImageCursor fpCursor = new IcyBufferedImageCursor(fp);

		for (int i = 0; i < Ydimt; i++) {
			final int i_offset_t = i * Xdimt;
			final int i_offset_s = i * Xdims;
			final int i_s_sub = i * sSubFactorY;
			final int i_t_sub = i * tSubFactorY;

			for (int j = 0; j < Xdimt; j++) {

				if (sourceMask.contains(j * sSubFactorX, i_s_sub) && targetMask.contains(j * tSubFactorX, i_t_sub) && j < Xdims
						&& i < Ydims)
					fpCursor.setSafe(j, i, 0, tImage[i_offset_t + j] - sImage[i_offset_s + j]);
				else {
					fpCursor.setSafe(j, i, 0, 0d);
				}

			}
		}
		fp.updateChannelsBounds();

		sourceProgressOutput = new Sequence("Output Source-Target" + extraTitleS, fp);
		notifySourceProgressOutputAvailable(sourceProgressOutput);
	}

	private void notifySourceProgressOutputAvailable(Sequence sourceOutput) {
		progressOutputListeners.forEach(l -> l.notifyProgress(0, "", sourceOutput));
	}

	private void initializeTargetProgressOutput() {
		int Ydimt = targetBSplineModel.getHeight();
		int Xdimt = targetBSplineModel.getWidth();
		int Xdims = sourceBSplineModel.getWidth();
		int Ydims = sourceBSplineModel.getHeight();
		double[] tImage = targetBSplineModel.isSubOutput()? targetBSplineModel.getSubImage(): targetBSplineModel.getImage();
		double[] sImage = sourceBSplineModel.isSubOutput()? sourceBSplineModel.getSubImage(): sourceBSplineModel.getImage();
		int sSubFactorX = 1;
		int sSubFactorY = 1;
		int tSubFactorX = 1;
		int tSubFactorY = 1;

		String extraTitleT = "";

		if (targetBSplineModel.isSubOutput() || sourceBSplineModel.isSubOutput())
			notifyProgress(Double.NaN, "Initializing target -> source window...");

		// If the output (difference) images are subsampled (because they were
		// larger than the maximum size), update variables.
		if (targetBSplineModel.isSubOutput()) {
			tSubFactorX = Xdimt / targetBSplineModel.getSubWidth();
			tSubFactorY = Ydimt / targetBSplineModel.getSubHeight();
			extraTitleT = " (Subsampled)";
			Xdimt = targetBSplineModel.getSubWidth();
			Ydimt = targetBSplineModel.getSubHeight();
		}

		if (sourceBSplineModel.isSubOutput()) {
			sSubFactorX = Xdims / sourceBSplineModel.getSubWidth();
			sSubFactorY = Ydims / sourceBSplineModel.getSubHeight();
			Xdims = sourceBSplineModel.getSubWidth();
			Ydims = sourceBSplineModel.getSubHeight();
		}

		if (getRegistrationMode() != RegistrationMode.MONO) {
			final IcyBufferedImage fp2 = new IcyBufferedImage(Xdims, Ydims, 1, DataType.DOUBLE);
			IcyBufferedImageCursor fp2Cursor = new IcyBufferedImageCursor(fp2);

			for (int i = 0; i < Ydims; i++) {
				int i_offset_t = i * Xdimt;
				int i_offset_s = i * Xdims;
				int i_s_sub = i * sSubFactorY;
				int i_t_sub = i * tSubFactorY;

				for (int j = 0; j < Xdims; j++)
					if (targetMask.contains(j * tSubFactorX, i_t_sub) && sourceMask.contains(j * sSubFactorX, i_s_sub)
							&& i < Ydimt && j < Xdimt)
						fp2Cursor.setSafe(j, i, 0, sImage[i_offset_s + j] - tImage[i_offset_t + j]);
					else
						fp2Cursor.setSafe(j, i, 0, 0d);
			}
			fp2Cursor.commitChanges();

			targetProgressOutput = new Sequence("Output Target-Source" + extraTitleT, fp2);
			notifyTargetProgressOutputAvailable(targetProgressOutput);
		}
	}

	private void notifyTargetProgressOutputAvailable(Sequence targetOutput) {
		progressOutputListeners.forEach(l -> l.notifyProgress(0, "", targetOutput));
	}

	private void adjustConsistencyWeight() {
		if (registrationMode == RegistrationMode.MONO) {
			setConsistencyWeight(0d);
		}
	}

	private void createTransformation() {
		transformation = new Transformation(getSourceSequence(), getTargetSequence(), sourceBSplineModel,
				targetBSplineModel, getSourceLandmarks(), getTargetLandmarks(), getSourceMask(), getTargetMask(), null, null,
				getInitialDeformationScale().getNumber(), getFinalDeformationScale().getNumber(), 0, getDivWeight(),
				getCurlWeight(), getLandmarkWeight(), getImageWeight(), getConsistencyWeight(), getStopThreshold(),
				(isShowProcess()? 2: 0), isShowProcess(), getRegistrationMode().getNumber(), "", "", sourceProgressOutput,
				targetProgressOutput, originalSourceImage, originalTargetImage);
		transformation.addProgressListener((double progress, String message, Object data) -> {
			notifyProgress(progress, message);
			return true;
		});
	}

	public Transformation getTransformation() {
		return transformation;
	}

	private void computeTransformation() throws InterruptedException {
		notifyProgress(Double.NaN, "Registering...");
		long registrationStartTime = System.currentTimeMillis();
		if (getRegistrationMode() == RegistrationMode.MONO) {
			transformation.doUnidirectionalRegistration();
		} else {
			transformation.doBidirectionalRegistration();
		}

		if (sourceBSplineModel.isSubOutput()) {
			notifyProgress(Double.NaN, "Calculating final transformed source image");
		}
		directResult = transformation.getDirectResults();

		if (getRegistrationMode() != RegistrationMode.MONO) {
			if (targetBSplineModel.isSubOutput()) {
				notifyProgress(Double.NaN, "Calculating final transformed target image");
			}
			indirectResult = transformation.getInverseResults();
		}
		long registrationEndTime = System.currentTimeMillis();
		System.out.format("Registration time %d milliseconds.\n", registrationEndTime - registrationStartTime);

	}

	private void restablishOriginalImages() {
		getSourceSequence().setImage(0, 0, originalSourceImage);
		getTargetSequence().setImage(0, 0, originalTargetImage);
	}

}
