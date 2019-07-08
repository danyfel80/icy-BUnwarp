package plugins.danyfel80.registration.elastic;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import danyfel80.registration.bspline.classic.BUnwarpRegistration;
import danyfel80.registration.bspline.classic.DeformationScale;
import danyfel80.registration.bspline.classic.RegistrationMode;
import danyfel80.registration.bspline.classic.Transformation;
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
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.vars.lang.Var;
import plugins.adufour.vars.lang.VarSequence;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DPoint;

public class SimpleBUnwarp extends EzPlug implements EzStoppable, Block {
	/* Input Variables */
	private EzVarSequence varSourceSequence;
	private EzVarSequence varTargetSequence;
	private Var<List<ROI2DPoint>> varSourceLandmarks;
	private Var<List<ROI2DPoint>> varTargetLandmarks;
	private Var<ROI2DArea> varSourceMask;
	private Var<ROI2DArea> varTargetMask;
	private EzVarEnum<RegistrationMode> varRegistrationMode;
	private EzVarInteger varInitialSubsampleFactor;

	private EzVarSequence varTransformedSourceSequence;
	private EzVarSequence varTransformedTargetSequence;

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

	private BUnwarpRegistration registration;
	private List<Sequence> progressSequences;

	/* Output Variables */
	private VarSequence varResultSourceSequence;
	private VarSequence varResultTargetSequence;
	private Var<Transformation> varResultTransformation;

	@Override
	protected void initialize() {
		initializeInputVariables();
		initializeInputListeners();
		initializeDefaultInputValues();
		initializeVariableGroups();
		addInputComponents();
	}

	private void initializeInputVariables() {
		varSourceSequence = new EzVarSequence("Source image");
		varTargetSequence = new EzVarSequence("Target image");
		varSourceLandmarks = new Var<List<ROI2DPoint>>("Source landmarks", new ArrayList<>());
		varTargetLandmarks = new Var<List<ROI2DPoint>>("Target landmarks", new ArrayList<>());
		varSourceMask = new Var<ROI2DArea>("Source mask", new ROI2DArea());
		varTargetMask = new Var<ROI2DArea>("Target mask", new ROI2DArea());
		varRegistrationMode = new EzVarEnum<>("Registration mode", RegistrationMode.values());
		varInitialSubsampleFactor = new EzVarInteger("Starting image subsampling factor", 0, 0, 7, 1);

		varTransformedSourceSequence = new EzVarSequence("Transformed source image");
		varTransformedTargetSequence = new EzVarSequence("Transformed target image");

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
	}

	private void initializeInputListeners() {
		varSourceSequence
				.addVarChangeListener((sourceVar, newSequence) -> varTransformedSourceSequence.setValue(newSequence));
		varTargetSequence
				.addVarChangeListener((sourceVar, newSequence) -> varTransformedTargetSequence.setValue(newSequence));

		varRegistrationMode.addVarChangeListener(
				(sourceVar, newMode) -> varConsistencyWeight.setEnabled(newMode != RegistrationMode.MONO));
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
		EzGroup transformedSequenceGroup = new EzGroup("Images to be transformed", varTransformedSourceSequence,
				varTransformedTargetSequence);
		EzGroup weightsGroup = new EzGroup("Weights", varDivWeight, varCurlWeight, varLandmarkWeight, varImageWeight,
				varConsistencyWeight);
		advancedParamsGroup = new EzGroup("Advanced Parameters", varInitialDeformationScale, varFinalDeformationScale,
				transformedSequenceGroup, weightsGroup, varStopThreshold, varShowProcess);
		transformedSequenceGroup.setFoldedState(true);
		weightsGroup.setFoldedState(true);
		advancedParamsGroup.setFoldedState(true);
	}

	private void addInputComponents() {
		addEzComponent(varSourceSequence);
		addEzComponent(varTargetSequence);
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
		inputMap.add(varSourceSequence.name, varSourceSequence.getVariable());
		inputMap.add(varTargetSequence.name, varTargetSequence.getVariable());
		inputMap.add(varSourceLandmarks.getName(), varSourceLandmarks);
		inputMap.add(varTargetLandmarks.getName(), varTargetLandmarks);
		inputMap.add(varSourceMask.getName(), varSourceMask);
		inputMap.add(varTargetMask.getName(), varTargetMask);
		inputMap.add(varRegistrationMode.name, varRegistrationMode.getVariable());
		inputMap.add(varInitialSubsampleFactor.name, varInitialSubsampleFactor.getVariable());
		inputMap.add(varInitialDeformationScale.name, varInitialDeformationScale.getVariable());
		inputMap.add(varFinalDeformationScale.name, varFinalDeformationScale.getVariable());
		inputMap.add(varTransformedSourceSequence.name, varTransformedSourceSequence.getVariable());
		inputMap.add(varTransformedTargetSequence.name, varTransformedTargetSequence.getVariable());
		inputMap.add(varDivWeight.name, varDivWeight.getVariable());
		inputMap.add(varCurlWeight.name, varCurlWeight.getVariable());
		inputMap.add(varLandmarkWeight.name, varLandmarkWeight.getVariable());
		inputMap.add(varImageWeight.name, varImageWeight.getVariable());
		inputMap.add(varConsistencyWeight.name, varConsistencyWeight.getVariable());
		inputMap.add(varStopThreshold.name, varStopThreshold.getVariable());
		inputMap.add(varShowProcess.name, varShowProcess.getVariable());
	}

	@Override
	public void declareOutput(VarList outputMap) {
		initializeOutputVariables();
		addOutputVariablesToMap(outputMap);
	}

	private void initializeOutputVariables() {
		varResultSourceSequence = new VarSequence("Result Source-Target", null);
		varResultTargetSequence = new VarSequence("Result Target-Source", null);
		varResultTransformation = new Var<>("Transformation model", Transformation.class);
	}

	private void addOutputVariablesToMap(VarList outputMap) {
		outputMap.add(varResultSourceSequence.getName(), varResultSourceSequence);
		outputMap.add(varResultTargetSequence.getName(), varResultTargetSequence);
		outputMap.add(varResultTransformation.getName(), varResultTransformation);
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
		}
		showResults();
	}

	private void createRegistration() {
		registration = new BUnwarpRegistration();
	}

	private void notifyProgress(double progress, String message) {
		if (!isHeadLess()) {
			getUI().setProgressBarMessage(message);
			getUI().setProgressBarValue(progress);
		}
	}

	private void readParameters() {
		this.registration.setSourceSequence(varSourceSequence.getValue(true));
		this.registration.setTargetSequence(varTargetSequence.getValue(true));
		this.registration.setRegistrationMode(varRegistrationMode.getValue(true));
		this.registration.setInitialSubsampleFactor(varInitialSubsampleFactor.getValue(true));
		this.registration.setTransformedSourceSequence(varTransformedSourceSequence.getValue(true));
		this.registration.setTransformedTargetSequence(varTransformedTargetSequence.getValue(true));
		this.registration.setInitialDeformationScale(varInitialDeformationScale.getValue(true));
		this.registration.setFinalDeformationScale(varFinalDeformationScale.getValue(true));
		this.registration.setDivWeight(varDivWeight.getValue(true));
		this.registration.setCurlWeight(varCurlWeight.getValue(true));
		this.registration.setLandmarkWeight(varLandmarkWeight.getValue(true));
		this.registration.setImageWeight(varImageWeight.getValue(true));
		this.registration.setConsistencyWeight(varConsistencyWeight.getValue(true));
		this.registration.setStopThreshold(varStopThreshold.getValue(true));
		this.registration.setShowProcess(varShowProcess.getValue(true));

		if (varSourceLandmarks.getValue() == null || varSourceLandmarks.getValue().isEmpty()
				|| varTargetLandmarks.getValue() == null || varTargetLandmarks.getValue().isEmpty()) {
			extractLandmarks();
		}
		this.registration.setSourceLandmarks(
				varSourceLandmarks.getValue(true).stream().map(l -> l.getPoint()).collect(Collectors.toList()));
		this.registration.setTargetLandmarks(
				varTargetLandmarks.getValue(true).stream().map(l -> l.getPoint()).collect(Collectors.toList()));

		if (varSourceMask.getValue() == null || varSourceMask.getValue().isEmpty() || varTargetMask.getValue() == null
				|| varTargetMask.getValue().isEmpty()) {
			extractMasks();
		}

		this.registration.setSourceMask(varSourceMask.getValue(true));
		this.registration.setTargetMask(varTargetMask.getValue(true));
	}

	private void extractLandmarks() {
		varSourceLandmarks.setValue(this.registration.getSourceSequence().getROIs(ROI2DPoint.class, false));
		varTargetLandmarks.setValue(this.registration.getTargetSequence().getROIs(ROI2DPoint.class, false));
		organizeLandmarks();
	}

	private void organizeLandmarks() {
		if (varSourceLandmarks.getValue().size() != varTargetLandmarks.getValue().size())
			throw new RuntimeException(String.format("Source landmarks(%d) and target landmarks(%d) have different size.",
					varSourceLandmarks.getValue().size(), varTargetLandmarks.getValue().size()));
		varSourceLandmarks.setValue(varSourceLandmarks.getValue().stream().sorted(Comparator.comparing(ROI2DPoint::getName))
				.collect(Collectors.toList()));
		varTargetLandmarks.setValue(varTargetLandmarks.getValue().stream().sorted(Comparator.comparing(ROI2DPoint::getName))
				.collect(Collectors.toList()));
		for (Iterator<ROI2DPoint> itSource = varSourceLandmarks.getValue().iterator(), itTarget = varTargetLandmarks
				.getValue().iterator(); itSource.hasNext();) {
			ROI2DPoint sourceLandmark = itSource.next();
			ROI2DPoint targetLandmark = itTarget.next();
			if (!sourceLandmark.getName().equals(targetLandmark.getName()))
				throw new RuntimeException("No corresponding landmark for " + sourceLandmark.getName());
		}
	}

	private void extractMasks() {
		varSourceMask.setValue(new ROI2DArea());
		varSourceMask.getValue().addRect(0, 0, registration.getSourceSequence().getWidth(),
				registration.getSourceSequence().getHeight());
		varTargetMask.setValue(new ROI2DArea());
		varTargetMask.getValue().addRect(0, 0, registration.getTargetSequence().getWidth(),
				registration.getTargetSequence().getHeight());
	}

	private void computeRegistration() throws InterruptedException {
		addProgressListeners();
		try {
			registration.compute();
		} finally {
			cleanProgress();
		}
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

	private void showResults() {
		Sequence transformedInputSource = registration.getDirectResult();
		transformedInputSource.setName(registration.getSourceSequence().getName() + "_Warped");
		Sequence transformedSource = registration.getTransformedDirectResult();

		if (!isHeadLess()) {
			addSequence(transformedInputSource);
			addSequence(transformedSource);
		} else {
			varResultSourceSequence.setValue(transformedSource);
			varResultTransformation.setValue(registration.getTransformation());
		}
		if (registration.getRegistrationMode() != RegistrationMode.MONO) {
			Sequence transformedInputTarget = registration.getIndirectResult();
			transformedInputTarget.setName(registration.getTargetSequence().getName() + "_Warped");
			Sequence transformedTarget = registration.getTransformedIndirectResult();
			if (!isHeadLess()) {
				addSequence(transformedInputTarget);
				addSequence(transformedTarget);
			} else {
				varResultTargetSequence.setValue(transformedTarget);
			}
		}
	}

	@Override
	public void clean() {}

}
