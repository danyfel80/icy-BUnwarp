package plugins.danyfel80.registration.elastic;

import java.io.File;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import danyfel80.registration.bspline.big.BUnwarpRegistrationBig;
import danyfel80.registration.bspline.big.BUnwarpRegistrationBigException;
import danyfel80.registration.bspline.classic.DeformationScale;
import danyfel80.registration.bspline.classic.RegistrationMode;
import icy.gui.main.MainInterfaceBatch;
import icy.main.Icy;
import icy.sequence.Sequence;
import icy.system.IcyHandledException;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarEnum;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.vars.lang.VarBoolean;
import plugins.adufour.vars.lang.VarROIArray;
import plugins.kernel.roi.roi2d.ROI2DPoint;

public class BigBUnwarp extends EzPlug implements EzStoppable, Block {
	/* Input Variables */
	private EzVarFile varSourceFile;
	private EzVarFile varTargetFile;
	private VarROIArray varSourceLandmarks;
	private VarROIArray varTargetLandmarks;

	private EzVarEnum<RegistrationMode> varRegistrationMode;
	private EzVarInteger varInitialSubsampleFactor;

	private EzVarFile varTransformedSourceFile;
	private EzVarFile varTransformedTargetFile;

	private EzVarEnum<DeformationScale> varInitialDeformationScale;
	private EzVarEnum<DeformationScale> varFinalDeformationScale;

	private EzVarDouble varDivWeight;
	private EzVarDouble varCurlWeight;
	private EzVarDouble varLandmarkWeight;
	private EzVarDouble varImageWeight;
	private EzVarDouble varConsistencyWeight;

	private EzVarDouble varStopThreshold;
	private EzVarBoolean varShowProcess;

	private EzGroup advancedParamsGroup;

	private BUnwarpRegistrationBig registration;
	private List<Sequence> progressSequences;

	private EzVarFile varResultSourceFile;
	private EzVarFile varResultTargetFile;
	private VarBoolean varEndSignalVar;

	@Override
	protected void initialize() {
		initializeInputVariables();
		initializeInputListeners();
		initializeDefaultInputValues();
		initializeVariableGroups();
		addInputComponents();
	}

	private void initializeInputVariables() {
		varSourceFile = new EzVarFile("Source image file", null);
		varTargetFile = new EzVarFile("Target image file", null);
		varSourceLandmarks = new VarROIArray("Source image landmarks");
		varTargetLandmarks = new VarROIArray("Target image landmarks");
		varRegistrationMode = new EzVarEnum<>("Registration mode", RegistrationMode.values());
		varInitialSubsampleFactor = new EzVarInteger("Starting image subsampling factor", 0, 0, 7, 1);

		varTransformedSourceFile = new EzVarFile("Transformed source image file", null);
		varTransformedTargetFile = new EzVarFile("Transformed target image file", null);

		varInitialDeformationScale = new EzVarEnum<>("Initial deformation size", Arrays.stream(DeformationScale.values())
				.limit(DeformationScale.values().length - 1).toArray(DeformationScale[]::new));
		varFinalDeformationScale = new EzVarEnum<>("Final deformation size", DeformationScale.values());

		varDivWeight = new EzVarDouble("Divergence weight");
		varCurlWeight = new EzVarDouble("Curl weight");
		varLandmarkWeight = new EzVarDouble("Landmark weight");
		varImageWeight = new EzVarDouble("Image weight");
		varConsistencyWeight = new EzVarDouble("Consistency weight");

		varStopThreshold = new EzVarDouble("Stop threshold");
		varShowProcess = new EzVarBoolean("Show process", false);

		varResultSourceFile = new EzVarFile("Result source image file", null);
		varResultTargetFile = new EzVarFile("Result target image file", null);
	}

	private void initializeInputListeners() {
		varRegistrationMode.addVarChangeListener(
				(sourceVar, newMode) -> varConsistencyWeight.setEnabled(newMode != RegistrationMode.MONO));

		varSourceFile.addVarChangeListener((sourceVar, newFile) -> varTransformedSourceFile.setValue(newFile));
		varTargetFile.addVarChangeListener((sourceVar, newFile) -> varTransformedTargetFile.setValue(newFile));
	}

	private void initializeDefaultInputValues() {
		varRegistrationMode.setValue(RegistrationMode.MONO);
		varInitialSubsampleFactor.setValue(0);
		varInitialDeformationScale.setValue(DeformationScale.VERY_COARSE);
		varFinalDeformationScale.setValue(DeformationScale.FINE);
		varDivWeight.setValue(0d);
		varCurlWeight.setValue(0d);
		varLandmarkWeight.setValue(0d);
		varImageWeight.setValue(1d);
		varConsistencyWeight.setValue(10d);
		varStopThreshold.setValue(1e-2);
		varShowProcess.setValue(false);
	}

	private void initializeVariableGroups() {
		EzGroup transformedSequenceGroup = new EzGroup("Images to be transformed", varTransformedSourceFile,
				varTransformedTargetFile);
		EzGroup weightsGroup = new EzGroup("Weights", varDivWeight, varCurlWeight, varLandmarkWeight, varImageWeight,
				varConsistencyWeight);
		advancedParamsGroup = new EzGroup("Advanced Parameters", varInitialDeformationScale, varFinalDeformationScale,
				transformedSequenceGroup, weightsGroup, varStopThreshold, varShowProcess);
		transformedSequenceGroup.setFoldedState(true);
		weightsGroup.setFoldedState(true);
		advancedParamsGroup.setFoldedState(true);
	}

	private void addInputComponents() {
		addEzComponent(varSourceFile);
		addEzComponent(varTargetFile);
		addEzComponent(varRegistrationMode);
		addEzComponent(varInitialSubsampleFactor);
		addEzComponent(advancedParamsGroup);
	}

	@Override
	public void declareInput(VarList inputMap) {
		initializeInputVariables();
		initializeDefaultInputValues();
		addInputVariablesToMap(inputMap);
	}

	private void addInputVariablesToMap(VarList inputMap) {
		inputMap.add(varSourceFile.name, varSourceFile.getVariable());
		inputMap.add(varTargetFile.name, varTargetFile.getVariable());
		inputMap.add(varSourceLandmarks.getName(), varSourceLandmarks);
		inputMap.add(varTargetLandmarks.getName(), varTargetLandmarks);
		inputMap.add(varRegistrationMode.name, varRegistrationMode.getVariable());
		inputMap.add(varInitialSubsampleFactor.name, varInitialSubsampleFactor.getVariable());
		inputMap.add(varInitialDeformationScale.name, varInitialDeformationScale.getVariable());
		inputMap.add(varFinalDeformationScale.name, varFinalDeformationScale.getVariable());
		inputMap.add(varTransformedSourceFile.name, varTransformedSourceFile.getVariable());
		inputMap.add(varTransformedTargetFile.name, varTransformedTargetFile.getVariable());
		inputMap.add(varDivWeight.name, varDivWeight.getVariable());
		inputMap.add(varCurlWeight.name, varCurlWeight.getVariable());
		inputMap.add(varLandmarkWeight.name, varLandmarkWeight.getVariable());
		inputMap.add(varImageWeight.name, varImageWeight.getVariable());
		inputMap.add(varConsistencyWeight.name, varConsistencyWeight.getVariable());
		inputMap.add(varStopThreshold.name, varStopThreshold.getVariable());
		inputMap.add(varShowProcess.name, varShowProcess.getVariable());
		inputMap.add(varResultSourceFile.name, varResultSourceFile.getVariable());
		inputMap.add(varResultTargetFile.name, varResultTargetFile.getVariable());
	}

	@Override
	public void declareOutput(VarList outputMap) {
		varEndSignalVar = new VarBoolean("Signal", false);
	}

	@Override
	protected void execute() {
		notifyProgress(Double.NaN, "Starting process...");
		createRegistration();
		readParameters();
		try {
			computeRegistration();
		} catch (InterruptedException e) {
			throw new IcyHandledException("Registration interrupted", e);
		} catch (BUnwarpRegistrationBigException e) {
			e.printStackTrace();
			throw new IcyHandledException("Could not perform registration", e);
		}
	}

	private void createRegistration() {
		registration = new BUnwarpRegistrationBig();
	}

	private void notifyProgress(double progress, String message) {
		if (!isHeadLess()) {
			getUI().setProgressBarMessage(message);
			getUI().setProgressBarValue(progress);
		}
	}

	private void readParameters() {
		this.registration.setSourceFile(varSourceFile.getValue(true));
		this.registration.setTargetFile(varTargetFile.getValue(true));
		if (varSourceLandmarks.getValue() != null && varSourceLandmarks.getValue().length > 0
				&& varTargetLandmarks.getValue() != null && varTargetLandmarks.getValue().length > 0) {

			List<ROI2DPoint> sourceLandmarks = Arrays.stream(varSourceLandmarks.getValue())
					.filter(roi -> roi instanceof ROI2DPoint).map(roi -> (ROI2DPoint) roi).collect(Collectors.toList());
			List<ROI2DPoint> targetLandmarks = Arrays.stream(varTargetLandmarks.getValue())
					.filter(roi -> roi instanceof ROI2DPoint).map(roi -> (ROI2DPoint) roi).collect(Collectors.toList());

			if (sourceLandmarks.size() == targetLandmarks.size()) {
				this.registration.setSourceLandmarkROIs(sourceLandmarks);
				this.registration.setTargetLandmarkROIs(targetLandmarks);
			}

		}

		this.registration.setRegistrationMode(varRegistrationMode.getValue(true));
		this.registration.setInitialSubsampleFactor(varInitialSubsampleFactor.getValue(true));
		this.registration.setInitialDeformationScale(varInitialDeformationScale.getValue(true));
		this.registration.setFinalDeformationScale(varFinalDeformationScale.getValue(true));
		this.registration.setDivWeight(varDivWeight.getValue(true));
		this.registration.setCurlWeight(varCurlWeight.getValue(true));
		this.registration.setLandmarkWeight(varLandmarkWeight.getValue(true));
		this.registration.setImageWeight(varImageWeight.getValue(true));
		if (registration.getRegistrationMode() != RegistrationMode.MONO)
			this.registration.setConsistencyWeight(varConsistencyWeight.getValue(true));
		this.registration.setStopThreshold(varStopThreshold.getValue(true));
		this.registration.setShowProcess(varShowProcess.getValue(true));

		this.registration.setTransformedSourceFile(varTransformedSourceFile.getValue(true));
		this.registration.setTransformedTargetFile(varTransformedTargetFile.getValue(true));

		if (varResultSourceFile.getValue() == null) {
			this.registration.setTransformedSourceOutputFile(getOutputFileName(varTransformedSourceFile.getValue(true)));
		} else {
			this.registration.setTransformedSourceOutputFile(varResultSourceFile.getValue(true));
		}
		if (varResultTargetFile.getValue() == null) {
			this.registration.setTransformedTargetOutputFile(getOutputFileName(varTransformedTargetFile.getValue(true)));
		} else {
			this.registration.setTransformedTargetOutputFile(varResultTargetFile.getValue(true));
		}
	}

	private File getOutputFileName(File inputFile) {
		Path inputFilePath = inputFile.toPath();
		String inputFileName = inputFilePath.getFileName().toString();
		int indexOfExtension = inputFileName.lastIndexOf('.');
		if (indexOfExtension >= 1) {
			inputFileName = inputFileName.substring(0, indexOfExtension);
		}
		inputFileName += "_Warped.tif";
		return inputFilePath.resolveSibling(inputFileName).toFile();
	}

	private void computeRegistration() throws InterruptedException, BUnwarpRegistrationBigException {
		addProgressListeners();
		if (isHeadLess())
			varEndSignalVar.setValue(false);
		try {
			registration.compute();
		} finally {
			cleanProgress();
		}
		if (isHeadLess())
			varEndSignalVar.setValue(true);
	}

	private void addProgressListeners() {
		if (!isHeadLess()) {
			registration.addProgressListener((double progress, String message, Object data) -> {
				getUI().setProgressBarMessage(message);
				getUI().setProgressBarValue(progress);
				return true;
			});
		} else {
			final DecimalFormat formatter = new java.text.DecimalFormat("000.##");
			registration.addProgressListener((double progress, String message, Object data) -> {
				System.out.format("Registering %s%%: %s\n", formatter.format(progress * 100d), message);
				return true;
			});
		}
		progressSequences = new ArrayList<>();
		if (varShowProcess.getValue() && !(Icy.getMainInterface() instanceof MainInterfaceBatch)) {
			registration.addProgressOutputListener((double progress, String message, Object data) -> {
				progressSequences.add((Sequence) data);
				addSequence((Sequence) data);
				return true;
			});
		}
	}

	private void cleanProgress() {
		progressSequences.forEach(s -> {
			removeSequence(s);
		});
		if (!isHeadLess()) {
			getUI().setProgressBarValue(0);
		}
	}

	@Override
	public void clean() {}

}
