/**
 * 
 */
package algorithms.danyfel80.registration.bunwarp.old;

import java.awt.geom.Point2D;
import java.util.List;

import algorithms.danyfel80.registration.bunwarp.MaximumScaleDeformationEnum;
import algorithms.danyfel80.registration.bunwarp.MinimumScaleDeformationEnum;
import algorithms.danyfel80.registration.bunwarp.RegistrationModeEnum;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import ij.IJ;
import ij.ImagePlus;
import loci.plugins.in.MainDialog;
import plugins.adufour.ezplug.EzPlug;
import plugins.kernel.roi.roi2d.ROI2DPoint;
import plugins.kernel.roi.roi2d.ROI2DPolygon;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class BUnwarpper extends Thread {

	// Input
	private Sequence srcSeq;
	private Sequence tgtSeq;

	private RegistrationModeEnum mode;
	private int maxSubsamplingFactor;
	private MinimumScaleDeformationEnum minScaleDef;
	private MaximumScaleDeformationEnum maxScaleDef;
	private double divWeight;
	private double curlWeight;
	private double landmarkWeight;
	private double consistencyWeight;
	private double imageWeight;
	private double stopThreshold;
	private boolean showProcess;
	private EzPlug plugin;

	private List<ROI2DPoint> srcLandmarks;
	private List<ROI2DPoint> tgtLandmarks;

	private ROI2DPolygon srcMask;
	private ROI2DPolygon tgtMask;
	
	private double[][] srcAffineMatrix;
	private double[][] tgtAffineMatrix;

	// Output

	// Internal
	private BSplineModel srcModel = null;
	private BSplineModel tgtModel = null;
	private int imagePyramidDepth;
	private int minImgScale = 0;

	public BUnwarpper(Sequence srcSeq, Sequence tgtSeq, List<ROI2DPoint> srcLandmarks, List<ROI2DPoint> tgtLandmarks,
	    ROI2DPolygon srcMask, ROI2DPolygon tgtMask, RegistrationModeEnum mode, Integer maxSubsamplingFactor,
	    MinimumScaleDeformationEnum minScaleDef, MaximumScaleDeformationEnum maxScaleDef, Double divWeight,
	    Double curlWeight, Double landmarkWeight, Double consistencyWeight, Double imageWeight, Double stopThreshold,
	    Boolean showProcess, EzPlug plugin) {
		this.srcSeq = srcSeq;
		this.tgtSeq = tgtSeq;
		this.srcLandmarks = srcLandmarks;
		this.tgtLandmarks = tgtLandmarks;
		this.srcMask = srcMask;
		this.tgtMask = tgtMask;
		this.mode = mode;
		this.maxSubsamplingFactor = maxSubsamplingFactor;
		this.minScaleDef = minScaleDef;
		this.maxScaleDef = maxScaleDef;
		this.divWeight = divWeight;
		this.curlWeight = curlWeight;
		this.landmarkWeight = landmarkWeight;
		this.consistencyWeight = consistencyWeight;
		this.imageWeight = imageWeight;
		this.stopThreshold = stopThreshold;
		this.showProcess = showProcess;
		this.plugin = plugin;
		createSourceImage(mode);
		createTargetImage();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Thread#run()
	 */
	@Override
	public void run() {
		super.run();

		// Start pyramids
		System.out.println("Starting image pyramids...");
		if (tgtModel.getWidth() > BSplineModel.MAX_OUTPUT_SIZE || tgtModel.getHeight() > BSplineModel.MAX_OUTPUT_SIZE
		    || srcModel.getWidth() > BSplineModel.MAX_OUTPUT_SIZE || srcModel.getHeight() > BSplineModel.MAX_OUTPUT_SIZE) {
			System.out.println("Starting image pyramids...");
		}
		srcModel.startPyramids();
		tgtModel.startPyramids();

		// Wait for the pyramids to be done
		try {
			srcModel.getThread().join();
			tgtModel.getThread().join();
		} catch (InterruptedException e) {
			System.err.println("Unexpected interruption exception" + e);
		}

		// Create output image (source-target)
		final Sequence[] outputSeqs = initializeOutputSeqs();

		// If mono mode, reset consistency weight
		if (this.mode == RegistrationModeEnum.MONO)
			this.consistencyWeight = 0.0;

		// Prepare registration parameters
		final Transformation warp = new Transformation(srcSeq, tgtSeq, srcModel, tgtModel, srcLandmarks, tgtLandmarks,
		    srcMask, tgtMask, srcAffineMatrix, tgtAffineMatrix, minScaleDef.getNumber(), maxScaleDef.getNumber(), minImgScale, divWeight,
		    curlWeight, landmarkWeight, imageWeight, consistencyWeight, stopThreshold, ((showProcess)? 2: 1), showProcess, mode, "", "",
		    outputSeqs[0], outputSeqs[1], plugin);

		// Perform the registration
		System.out.println("Registering...");

		long start = System.currentTimeMillis(); // start timing

		if (this.mode == RegistrationModeEnum.MONO) {
			// Do unidirectional registration
			warp.doUnidirectionalRegistration();
			// Save transformation
			warp.saveDirectTransformation();
			// Show results.
			if (srcModel.isSubOutput()) {
				IJ.log("Calculating final transformed source image");
			}
			warp.showDirectResults();
		} else {
			// Do bidirectional registration
			warp.doBidirectionalRegistration();
			// Save transformations
			warp.saveDirectTransformation();
			warp.saveInverseTransformation();
			// Show results.
			if (srcModel.isSubOutput()) {
				IJ.log("Calculating final transformed source image");
			}
			warp.showDirectResults();
			if (tgtModel.isSubOutput()) {
				IJ.log("Calculating final transformed target image");
			}
			warp.showInverseResults();
		}

		long stop = System.currentTimeMillis(); // stop timing
		// Print execution time
		if (showProcess)
			System.out.println("\nRegistration time: " + (stop - start) + "ms");
	}

	/**
	 * @return
	 */
	private Sequence[] initializeOutputSeqs() {
		int Ydimt = tgtModel.getHeight();
		int Xdimt = tgtModel.getWidth();
		int Xdims = srcModel.getWidth();
		int Ydims = srcModel.getHeight();
		double[] tImage = tgtModel.isSubOutput() ? tgtModel.getSubImage() : tgtModel.getImage();
		double[] sImage = srcModel.isSubOutput() ? srcModel.getSubImage() : srcModel.getImage();
		int sSubFactorX = 1;
		int sSubFactorY = 1;
		int tSubFactorX = 1;
		int tSubFactorY = 1;
		Sequence[] outputSeq = new Sequence[2];

		String extraTitleS = "";
		String extraTitleT = "";

		if (tgtModel.isSubOutput() || srcModel.isSubOutput()) {
			System.out.println("Initializing output windows...");
		}

		// If the output (difference) images are subsampled (because they were
		// larger than the maximum size), update variables.
		if (tgtModel.isSubOutput()) {
			tSubFactorX = Xdimt / tgtModel.getSubWidth();
			tSubFactorY = Ydimt / tgtModel.getSubHeight();
			extraTitleT = " (Subsampled)";
			Xdimt = tgtModel.getSubWidth();
			Ydimt = tgtModel.getSubHeight();
		}

		if (srcModel.isSubOutput()) {
			sSubFactorX = Xdims / srcModel.getSubWidth();
			sSubFactorY = Ydims / srcModel.getSubHeight();
			extraTitleS = " (Subsampled)";
			Xdims = srcModel.getSubWidth();
			Ydims = srcModel.getSubHeight();
		}

		// Float processor for the output source-target image.
		final Sequence stSeq = new Sequence(new IcyBufferedImage(Xdimt, Ydimt, 1, DataType.FLOAT));
		stSeq.beginUpdate();
		float[] stData = stSeq.getDataXYAsFloat(0, 0, 0);

		for (int i = 0; i < Ydimt; i++) {
			final int i_offset_t = i * Xdimt;
			final int i_offset_s = i * Xdims;
			final int i_s_sub = i * sSubFactorY;
			final int i_t_sub = i * tSubFactorY;

			for (int j = 0; j < Xdimt; j++) {

				if ((srcMask == null || (srcMask.contains(new Point2D.Double(j * sSubFactorX, i_s_sub))
				    && tgtMask.contains(new Point2D.Double(j * tSubFactorX, i_t_sub)))) && j < Xdims && i < Ydims)
					stData[j + i_offset_t] = (float) (tImage[i_offset_t + j] - sImage[i_offset_s + j]);
				else {
					stData[j + i_offset_t] = 0;
				}

			}
		}
		stSeq.dataChanged();
		stSeq.endUpdate();

		stSeq.setName("Output Source-Target" + extraTitleS);
		plugin.addSequence(stSeq);

		outputSeq[0] = stSeq;

		// Create output image (target-source) if necessary

		if (this.mode != RegistrationModeEnum.MONO) {
			final Sequence tsSeq = new Sequence(new IcyBufferedImage(Xdims, Ydims, 1, DataType.FLOAT));
			tsSeq.beginUpdate();
			float[] tsData = tsSeq.getDataXYAsFloat(0, 0, 0);

			for (int i = 0; i < Ydims; i++) {
				int i_offset_t = i * Xdimt;
				int i_offset_s = i * Xdims;
				int i_s_sub = i * sSubFactorY;
				int i_t_sub = i * tSubFactorY;

				for (int j = 0; j < Xdims; j++)
					if ((srcMask == null || (tgtMask.contains(new Point2D.Double(j * tSubFactorX, i_t_sub))
					    && srcMask.contains(new Point2D.Double(j * sSubFactorX, i_s_sub)))) && i < Ydimt && j < Xdimt)
						tsData[j + i_offset_s] = (float) (sImage[i_offset_s + j] - tImage[i_offset_t + j]);
					else
						tsData[j + i_offset_s] = 0;
			}
			tsSeq.dataChanged();
			tsSeq.endUpdate();

			tsSeq.setName("Output Target-Source" + extraTitleT);
			plugin.addSequence(tsSeq);
			outputSeq[1] = tsSeq;
		} else
			outputSeq[1] = null;

		return outputSeq;
	}

	/**
	 * Creates the source image, i.e. initialize the B-spline model for the source
	 * image. The resolution pyramid is not started here, but in startPyramids.
	 * 
	 * @param mode
	 *          The registration mode. If mode == MONO, the registration is
	 *          performed only in one direction
	 */
	private void createSourceImage(RegistrationModeEnum mode) {
		srcModel = new BSplineModel(srcSeq, mode != RegistrationModeEnum.MONO, (int) Math.pow(2, maxSubsamplingFactor));
		computeImagePyramidDepth();
		srcModel.setPyramidDepth(imagePyramidDepth + minImgScale);
	}

	/**
	 * Creates the target image, i.e. initialize the B-spline model for the target
	 * image. The resolution pyramid is not started here, but in startPyramids.
	 */
	private void createTargetImage() {
		tgtModel = new BSplineModel(tgtSeq, true, (int) Math.pow(2, maxSubsamplingFactor));
		computeImagePyramidDepth();
		tgtModel.setPyramidDepth(imagePyramidDepth + minImgScale);
	}

	/**
	 * Compute the depth of the image resolution pyramid.
	 */
	private void computeImagePyramidDepth() {
		this.imagePyramidDepth = maxScaleDef.getNumber() - minScaleDef.getNumber() + 1;
	}

	public Sequence getDirectTransform() {
		// TODO Auto-generated method stub
		return null;
	}

	public Sequence getInverseTransform() {
		// TODO Auto-generated method stub
		return null;
	}

}
