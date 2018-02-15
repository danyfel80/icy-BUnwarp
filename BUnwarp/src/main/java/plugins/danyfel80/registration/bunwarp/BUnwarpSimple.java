package plugins.danyfel80.registration.bunwarp;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import algorithms.danyfel80.registration.bunwarp.BUnwarpper;
import algorithms.danyfel80.registration.bunwarp.MaximumScaleDeformationEnum;
import algorithms.danyfel80.registration.bunwarp.MinimumScaleDeformationEnum;
import algorithms.danyfel80.registration.bunwarp.RegistrationModeEnum;
import icy.gui.dialog.MessageDialog;
import icy.image.IcyBufferedImage;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarEnum;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.vars.lang.Var;
import plugins.adufour.vars.lang.VarInteger;
import plugins.adufour.vars.lang.VarSequence;
import plugins.kernel.roi.roi2d.ROI2DPoint;
import plugins.kernel.roi.roi2d.ROI2DPolygon;

/**
 * BUnwarp Plugin based on ImageJ plugin "bUnwarpJ" by Ignacio
 * Arganda-Carreras and Jan Kybic which extends previous UnwarpJ project by
 * Carlos Oscar Sanchez Sorzano.
 * Adapted by Daniel Gonzalez Obando
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class BUnwarpSimple extends BUnwarp {

	// Input variables
	// Images
	// - Source image
	EzVarSequence inSrcSeq = new EzVarSequence("Source");
	// - Target image
	EzVarSequence inTgtSeq = new EzVarSequence("Target");

	// Parameters
	// - Registration mode
	EzVarEnum<RegistrationModeEnum> inMode = new EzVarEnum<>("Mode", RegistrationModeEnum.values(),
			RegistrationModeEnum.ACCURATE);
	// - Subsampling factor
	EzVarInteger inSubsampleFactor = new EzVarInteger("Image Subsampling Factor", 0, 0, 7, 1);

	// - Advanced Parameters
	// - Source transformation target image
	EzVarSequence inSrcTgtSeq = new EzVarSequence("Transformation Source");
	// - Target transformation target image
	EzVarSequence inTgtTgtSeq = new EzVarSequence("Transformation Target");

	// - Initial deformation
	EzVarEnum<MinimumScaleDeformationEnum> inIniDef = new EzVarEnum<>("Initial deformation",
			MinimumScaleDeformationEnum.values(), MinimumScaleDeformationEnum.VERY_COARSE);
	// - Final deformation
	EzVarEnum<MaximumScaleDeformationEnum> inFnlDef = new EzVarEnum<>("Final Deformation",
			MaximumScaleDeformationEnum.values(), MaximumScaleDeformationEnum.FINE);

	// Weights
	// - Divergence Weight
	EzVarDouble inDivWeight = new EzVarDouble("Divergence Weight");
	// - Curl Weight
	EzVarDouble inCurlWeight = new EzVarDouble("Curl Weight");
	// - Landmark Weight
	EzVarDouble inLandmarkWeight = new EzVarDouble("Landmark Weight");
	// - Image Weight
	EzVarDouble inImageWeight = new EzVarDouble("Image Weight");
	// - Consistency Weight
	EzVarDouble inConsistencyWeight = new EzVarDouble("Consistency Weight");

	EzGroup weightsGroup = new EzGroup("Weights", inDivWeight, inCurlWeight, inLandmarkWeight, inImageWeight,
			inConsistencyWeight);

	// - Stop threshold
	EzVarDouble inStopThreshold = new EzVarDouble("Stop Threshold");
	// - Show process
	EzVarBoolean inShowProcess = new EzVarBoolean("Show Process", false);

	EzGroup	outputSequenceGroup	= new EzGroup("Transformed Output", inSrcTgtSeq, inTgtTgtSeq);
	EzGroup	advancedParamsGroup	= new EzGroup("Advanced Parameters", inIniDef, inFnlDef, outputSequenceGroup,
			weightsGroup, inStopThreshold, inShowProcess);

	// Output variables

	VarSequence	outSrcSeq	= new VarSequence("Source Registered", (Sequence) null);
	VarSequence	outTgtSeq	= new VarSequence("Target Registered", (Sequence) null);

	Var<double[][]>	outCxSourceToTarget	= new Var<double[][]>("Cx Source to Target", new double[1][1]);
	Var<double[][]>	outCySourceToTarget	= new Var<double[][]>("Cy Source to Target", new double[1][1]);
	Var<double[][]>	outCxTargetToSource	= new Var<double[][]>("Cx Target to Source", new double[1][1]);
	Var<double[][]>	outCyTargetToSource	= new Var<double[][]>("Cy Target to Source", new double[1][1]);
	VarInteger			outIntervals				= new VarInteger("Transform intervals", 0);

	// Internal variables
	Sequence	srcSeq;
	Sequence	tgtSeq;
	Sequence	srcTgtSeq;
	Sequence	tgtTgtSeq;

	IcyBufferedImage	originalSrcIBI;
	IcyBufferedImage	originalTgtIBI;

	BUnwarpper	bu;
	Thread			but;

	/*
	 * (non-Javadoc)
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.
	 * VarList)
	 */
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add(inSrcSeq.name, inSrcSeq.getVariable());
		inputMap.add(inTgtSeq.name, inTgtSeq.getVariable());
		inputMap.add(inSrcTgtSeq.name, inSrcTgtSeq.getVariable());
		inputMap.add(inTgtTgtSeq.name, inTgtTgtSeq.getVariable());
		inputMap.add(inMode.name, inMode.getVariable());
		inputMap.add(inSubsampleFactor.name, inSubsampleFactor.getVariable());
		inputMap.add(inIniDef.name, inIniDef.getVariable());
		inputMap.add(inFnlDef.name, inFnlDef.getVariable());
		inputMap.add(inDivWeight.name, inDivWeight.getVariable());
		inputMap.add(inCurlWeight.name, inCurlWeight.getVariable());
		inputMap.add(inLandmarkWeight.name, inLandmarkWeight.getVariable());
		inputMap.add(inImageWeight.name, inImageWeight.getVariable());
		inputMap.add(inConsistencyWeight.name, inConsistencyWeight.getVariable());
		inputMap.add(inStopThreshold.name, inStopThreshold.getVariable());
		inputMap.add(inShowProcess.name, inShowProcess.getVariable());

		inDivWeight.setValue(0d);
		inCurlWeight.setValue(0d);
		inLandmarkWeight.setValue(0d);
		inImageWeight.setValue(1d);
		inConsistencyWeight.setValue(10d);
		inStopThreshold.setValue(1e-2);
	}

	/*
	 * (non-Javadoc)
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util
	 * .VarList)
	 */
	@Override
	public void declareOutput(VarList outputMap) {
		outputMap.add(outSrcSeq.getName(), outSrcSeq);
		outputMap.add(outTgtSeq.getName(), outTgtSeq);

		outputMap.add(outCxSourceToTarget.getName(), outCxSourceToTarget);
		outputMap.add(outCySourceToTarget.getName(), outCySourceToTarget);
		outputMap.add(outCxTargetToSource.getName(), outCxTargetToSource);
		outputMap.add(outCyTargetToSource.getName(), outCyTargetToSource);
		outputMap.add(outIntervals.getName(), outIntervals);
	}

	/*
	 * (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#initialize()
	 */
	@Override
	protected void initialize() {
		inSrcSeq.setToolTipText("Source(floating) sequence used to perform the registration.");
		inTgtSeq.setToolTipText("Target(fixed) sequence used to perform the registration.");
		inMode.setToolTipText(
				"Mode of interpolation: Mono uses source -> target transformation. Fast or Accurate use source <-> target transformation.");
		inSubsampleFactor
				.setToolTipText("Level of subsampling of the source and target sequences to perform the registration.");

		inSrcTgtSeq.setToolTipText("Sequence used to apply source transformation.");
		inTgtTgtSeq.setToolTipText("Sequence used to apply target transformation.");

		inIniDef.setToolTipText("Sets the initial transformation detail.");
		inFnlDef.setToolTipText("Sets the final transformation detail.");

		inDivWeight.setToolTipText(
				"Weight related to the divergence of the tensors in the transformation. Higher value means result will have less divergence.");
		inCurlWeight.setToolTipText(
				"Weight related to the curl of the tensors in the transformation. Higher value means result will have less curl.");
		inLandmarkWeight.setToolTipText(
				"Weight related to landmarks present on the sequence. Higher value means landmarks have more impact on the result. Landmarks must be ROI2DPoints in the sequence.");
		inImageWeight.setToolTipText(
				"Weight related to image intensities. Higher value means image intensities will have more impact on the result.");
		inConsistencyWeight.setToolTipText(
				"When the mode is set to Fast or Accurate, this weight represents the similarity constraint on the s->t and t->s transformations. The higher the value, the more similar the transformations will be.");
		inStopThreshold.setToolTipText(
				"This is the optimization stop criteria. When the optimization changes the transformation less than the given value, the process ends and the result is shown.");

		inShowProcess
				.setToolTipText("If checked, more details of the transformation will be shown at the end of the procedure.");

		addEzComponent(inSrcSeq);
		addEzComponent(inTgtSeq);
		// addEzComponent(inSrcTgtSeq);
		// addEzComponent(inTgtTgtSeq);
		addEzComponent(inMode);
		addEzComponent(inSubsampleFactor);
		outputSequenceGroup.setFoldedState(true);
		weightsGroup.setFoldedState(true);
		advancedParamsGroup.setFoldedState(true);
		addEzComponent(advancedParamsGroup);

		inDivWeight.setValue(0d);
		inCurlWeight.setValue(0d);
		inLandmarkWeight.setValue(0d);
		inImageWeight.setValue(1d);
		inConsistencyWeight.setValue(10d);
		inStopThreshold.setValue(1e-2);

		inMode.addVarChangeListener(new EzVarListener<RegistrationModeEnum>() {

			@Override
			public void variableChanged(EzVar<RegistrationModeEnum> source, RegistrationModeEnum newValue) {
				inConsistencyWeight.setEnabled(newValue != RegistrationModeEnum.MONO);
			}
		});
		inMode.setValue(RegistrationModeEnum.MONO);
	}

	/*
	 * (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#execute()
	 */
	@Override
	protected void execute() {
		this.isPluginInterrupted = false;
		if (validateInput() != 0) {
			return;
		}
		srcSeq = inSrcSeq.getValue();
		tgtSeq = inTgtSeq.getValue();
		srcTgtSeq = inSrcTgtSeq.getValue();
		tgtTgtSeq = inTgtTgtSeq.getValue();
		if (srcTgtSeq == null) srcTgtSeq = srcSeq;
		if (tgtTgtSeq == null) tgtTgtSeq = tgtSeq;

		originalSrcIBI = srcSeq.getFirstImage();
		originalTgtIBI = tgtSeq.getFirstImage();

		List<? extends ROI> srcLandmarks = srcSeq.getROIs(ROI2DPoint.class);
		List<? extends ROI> tgtLandmarks = tgtSeq.getROIs(ROI2DPoint.class);

		Comparator<ROI> comp = new Comparator<ROI>() {
			@Override
			public int compare(ROI o1, ROI o2) {
				return o1.getName().compareTo(o2.getName());
			}
		};

		srcLandmarks.sort(comp);
		tgtLandmarks.sort(comp);

		ROI2DPolygon srcMask = null;
		ROI2DPolygon tgtMask = null;
		if (srcSeq.getROICount(ROI2DPolygon.class) > 0) {
			srcMask = (ROI2DPolygon) srcSeq.getROIs(ROI2DPolygon.class).get(0);
		}
		if (tgtSeq.getROICount(ROI2DPolygon.class) > 0) {
			tgtMask = (ROI2DPolygon) tgtSeq.getROIs(ROI2DPolygon.class).get(0);
		}
		if (srcMask == null) {
			List<Point2D> pts = new ArrayList<>();
			pts.add(new Point2D.Double(0, 0));
			pts.add(new Point2D.Double(0, srcSeq.getHeight()));
			pts.add(new Point2D.Double(srcSeq.getWidth(), srcSeq.getHeight()));
			pts.add(new Point2D.Double(srcSeq.getWidth(), 0));
			srcMask = new ROI2DPolygon(pts);
		}

		if (tgtMask == null) {
			List<Point2D> pts = new ArrayList<>();
			pts.add(new Point2D.Double(0, 0));
			pts.add(new Point2D.Double(0, tgtSeq.getHeight()));
			pts.add(new Point2D.Double(tgtSeq.getWidth(), tgtSeq.getHeight()));
			pts.add(new Point2D.Double(tgtSeq.getWidth(), 0));
			tgtMask = new ROI2DPolygon(pts);
		}

		@SuppressWarnings("unchecked")
		BUnwarpper buLocal = new BUnwarpper(srcSeq, tgtSeq, (List<ROI2DPoint>) srcLandmarks,
				(List<ROI2DPoint>) tgtLandmarks, srcMask, tgtMask, inSubsampleFactor.getValue(),
				inIniDef.getValue().getNumber(), inFnlDef.getValue().getNumber(), 0, inDivWeight.getValue(),
				inCurlWeight.getValue(), inLandmarkWeight.getValue(), inImageWeight.getValue(), inConsistencyWeight.getValue(),
				inStopThreshold.getValue(), inShowProcess.getValue() ? 2 : 1, inShowProcess.getValue(),
				inMode.getValue().getNumber(), this);
		bu = buLocal;
		but = new Thread(bu);
		but.start();
		try {
			but.join();
			but = null;
		} catch (InterruptedException e) {
			System.err.println("Thread interrupted: " + e.getMessage());
		}
		if (!this.isPluginInterrupted) {
			Sequence srcTgtCopySeq = SequenceUtil.getCopy(srcTgtSeq);
			bu.getRegisteredSource(srcTgtCopySeq);
			outSrcSeq.setValue(srcTgtSeq);
			addSequence(srcTgtCopySeq);

			outCxSourceToTarget.setValue(bu.getCxSourceToTarget());
			outCySourceToTarget.setValue(bu.getCxSourceToTarget());
			outIntervals.setValue(bu.getIntervals());

			if (inMode.getValue() != RegistrationModeEnum.MONO) {
				Sequence tgtTgtCopySeq = SequenceUtil.getCopy(tgtTgtSeq);
				bu.getRegisteredTarget(tgtTgtCopySeq);
				outTgtSeq.setValue(tgtTgtCopySeq);
				addSequence(tgtTgtCopySeq);

				outCxTargetToSource.setValue(bu.getCxSourceToTarget());
				outCyTargetToSource.setValue(bu.getCxSourceToTarget());
			}
		}

	}

	/**
	 * Validate plugin input variables
	 * 
	 * @return 0 if input is valid, else a positive number with the error code.
	 */
	private int validateInput() {
		if (inSrcSeq.getValue() == null || inTgtSeq.getValue() == null) {
			MessageDialog.showDialog("Error", "Please choose two valid images", MessageDialog.ERROR_MESSAGE);
			return 1;
		}
		if (inSrcSeq.getValue().getROICount(ROI2DPolygon.class) > 1
				|| inTgtSeq.getValue().getROICount(ROI2DPolygon.class) > 1) {
			MessageDialog.showDialog("Error", "Please define a single mask for each input sequence.",
					MessageDialog.ERROR_MESSAGE);
			return 2;
		}
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#stopExecution()
	 */
	@Override
	public void stopExecution() {
		isPluginInterrupted = true;
		if (but != null && but.isAlive()) {
			try {
				but.join();
				but = null;
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#clean()
	 */
	@Override
	public void clean() {}

	/*
	 * (non-Javadoc)
	 * @see plugins.danyfel80.registration.bunwarp.BUnwarp#restoreAll()
	 */
	@Override
	public void restoreAll() {
		ungrayInputImages();
		if (getUI() != null) {
			getUI().setProgressBarMessage("");
			getUI().setProgressBarValue(0);
		}
		Runtime.getRuntime().gc();
	}

	/**
	 * Restore original input sequences.
	 */
	private void ungrayInputImages() {
		srcSeq.setImage(0, 0, originalSrcIBI);
		tgtSeq.setImage(0, 0, originalTgtIBI);
	};

}
