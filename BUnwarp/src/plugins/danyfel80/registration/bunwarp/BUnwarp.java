package plugins.danyfel80.registration.bunwarp;

import algorithms.danyfel80.registration.bunwarp.BUnwarpper;
import algorithms.danyfel80.registration.bunwarp.MaximumScaleDeformationEnum;
import algorithms.danyfel80.registration.bunwarp.MinimumScaleDeformationEnum;
import algorithms.danyfel80.registration.bunwarp.RegistrationModeEnum;
import icy.gui.dialog.MessageDialog;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarEnum;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class BUnwarp extends EzPlug implements Block, EzStoppable {

	// Input variables
	// - Source image
	EzVarSequence inSrcSeq = new EzVarSequence("Source");
	// - Target image
	EzVarSequence inTgtSeq = new EzVarSequence("Target");
	// - Registration mode
	EzVarEnum<RegistrationModeEnum> inMode = new EzVarEnum<>("Mode", RegistrationModeEnum.values(),
	    RegistrationModeEnum.ACCURATE);
	// - Subsampling factor
	EzVarInteger inSubsampleFactor = new EzVarInteger("Image Subsampling Factor", 0, 0, 7, 1);
	// - Advanced Parameters
	// - Initial deformation
	EzVarEnum<MinimumScaleDeformationEnum> inIniDef = new EzVarEnum<>("Initial deformation",
	    MinimumScaleDeformationEnum.values(), MinimumScaleDeformationEnum.VERY_COARSE);
	// - Final deformation
	EzVarEnum<MaximumScaleDeformationEnum> inFnlDef = new EzVarEnum<>("Final Deformation",
	    MaximumScaleDeformationEnum.values(), MaximumScaleDeformationEnum.FINE);
	// - Weights
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

	EzGroup advancedParamsGroup = new EzGroup("Advanced Parameters", inIniDef, inFnlDef, weightsGroup, inStopThreshold,
	    inShowProcess);

	// Output variables

	// Internal variables

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.
	 * VarList)
	 */
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add(inSrcSeq.name, inSrcSeq.getVariable());
		inputMap.add(inTgtSeq.name, inTgtSeq.getVariable());
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
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util
	 * .VarList)
	 */
	@Override
	public void declareOutput(VarList outputMap) {
		// TODO Auto-generated method stub

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#initialize()
	 */
	@Override
	protected void initialize() {
		addEzComponent(inSrcSeq);
		addEzComponent(inTgtSeq);
		addEzComponent(inMode);
		addEzComponent(advancedParamsGroup);

		inDivWeight.setValue(0.1);
		inCurlWeight.setValue(0.1);
		inLandmarkWeight.setValue(0d);
		inImageWeight.setValue(1d);
		inConsistencyWeight.setValue(10d);

		inMode.addVarChangeListener(new EzVarListener<RegistrationModeEnum>() {

			@Override
			public void variableChanged(EzVar<RegistrationModeEnum> source, RegistrationModeEnum newValue) {
				inConsistencyWeight.setEnabled(newValue != RegistrationModeEnum.MONO);
			}
		});
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#execute()
	 */
	@Override
	protected void execute() {
		if (validateInput() != 0) {
			return;
		}
		BUnwarpper bu = new BUnwarpper(inSrcSeq.getValue(), inTgtSeq.getValue(), inMode.getValue(),
		    inSubsampleFactor.getValue(), inIniDef.getValue(), inFnlDef.getValue(), inDivWeight.getValue(),
		    inCurlWeight.getValue(), inLandmarkWeight.getValue(), inConsistencyWeight.getValue(), inImageWeight.getValue(),
		    inStopThreshold.getValue(), inShowProcess.getValue(), this);
		bu.start();
		try {
			bu.join();
			addSequence(bu.getDirectTransform());
			if (inMode.getValue() != RegistrationModeEnum.MONO) {
				addSequence(bu.getInverseTransform());
			}
		} catch (InterruptedException e) {
			System.err.println("Thread interrupted: " + e.getMessage());
		}
		
	}

	private int validateInput() {
		if (inSrcSeq.getValue() == null || inTgtSeq.getValue() == null) {
			MessageDialog.showDialog("Error", "Please choose two valid images", MessageDialog.ERROR_MESSAGE);
			return 1;
		}
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#stopExecution()
	 */
	@Override
	public void stopExecution() {
		// TODO Auto-generated method stub

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#clean()
	 */
	@Override
	public void clean() {
		// TODO Auto-generated method stub

	}

}
