package plugins.danyfel80.registration.bunwarp;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.io.FilenameUtils;

import algorithms.danyfel80.bigimage.BigImageLoader;
import algorithms.danyfel80.registration.bunwarp.BUnwarpper;
import algorithms.danyfel80.registration.bunwarp.MaximumScaleDeformationEnum;
import algorithms.danyfel80.registration.bunwarp.MinimumScaleDeformationEnum;
import algorithms.danyfel80.registration.bunwarp.ProgressBar;
import algorithms.danyfel80.registration.bunwarp.RegistrationModeEnum;
import icy.common.exception.UnsupportedFormatException;
import icy.gui.dialog.MessageDialog;
import icy.image.IcyBufferedImage;
import icy.roi.ROI;
import icy.sequence.Sequence;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarDoubleArrayNative;
import plugins.adufour.ezplug.EzVarEnum;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.kernel.importer.LociImporterPlugin;
import plugins.kernel.roi.roi2d.ROI2DPoint;
import plugins.kernel.roi.roi2d.ROI2DPolygon;

/**
 * BUnwarp plugin for big images. It uses image downsampling for registering and
 * tile processing for transform reconstruction.
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageBUnwarp extends BUnwarp {

	// - Source image file path
	EzVarFile inSrcFile = new EzVarFile("Source file", "");
	// - Target image file path
	EzVarFile inTgtFile = new EzVarFile("Target file", "");
	// - Source transformation image file path
	EzVarFile inSrcResultFile = new EzVarFile("Transformed source file", "");
	// - Target transformation image file path
	EzVarFile inTgtResultFile = new EzVarFile("Transformed target file", "");

	// Parameters
	// - Registration mode
	EzVarEnum<RegistrationModeEnum> inMode = new EzVarEnum<>("Mode", RegistrationModeEnum.values(),
	    RegistrationModeEnum.ACCURATE);
	double[][] scales = {{1.0, 0.16}};
	EzVarDoubleArrayNative inUsedScales = new EzVarDoubleArrayNative("Registration scales", scales, true);
	// - Subsampling factor
	EzVarInteger inSubsampleFactor = new EzVarInteger("Image Subsampling Factor", 0, 0, 7, 1);
	// - Advanced Parameters
	// - Initial deformation
	EzVarEnum<MinimumScaleDeformationEnum> inIniDef = new EzVarEnum<>("Initial deformation",
	    MinimumScaleDeformationEnum.values(), MinimumScaleDeformationEnum.VERY_COARSE);
	// - Final deformation
	EzVarEnum<MaximumScaleDeformationEnum> inFnlDef = new EzVarEnum<>("Final Deformation",
	    MaximumScaleDeformationEnum.values(), MaximumScaleDeformationEnum.VERY_FINE);

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

	EzGroup advancedParamsGroup = new EzGroup("Advanced Parameters", inIniDef, inFnlDef, weightsGroup, inStopThreshold,
	    inShowProcess);

	// Internal variables
	BUnwarpper bu;

	Sequence srcSeq;
	Sequence tgtSeq;
	Sequence srcDefSeq;
	Sequence tgtDefSeq;

	IcyBufferedImage originalSrcIBI;
	IcyBufferedImage originalTgtIBI;

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.
	 * VarList)
	 */
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add(inSrcFile.name, inSrcFile.getVariable());
		inputMap.add(inTgtFile.name, inTgtFile.getVariable());
		inputMap.add(inSrcResultFile.name, inSrcResultFile.getVariable());
		inputMap.add(inTgtResultFile.name, inTgtResultFile.getVariable());
		inputMap.add(inMode.name, inMode.getVariable());
		inputMap.add(inUsedScales.name, inUsedScales.getVariable());
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
	 * 
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util
	 * .VarList)
	 */
	@Override
	public void declareOutput(VarList outputMap) {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#initialize()
	 */
	@Override
	protected void initialize() {
		addEzComponent(inSrcFile);
		addEzComponent(inTgtFile);
		addEzComponent(inSrcResultFile);
		addEzComponent(inTgtResultFile);
		addEzComponent(inMode);
		addEzComponent(inUsedScales);
		addEzComponent(inSubsampleFactor);
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
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#execute()
	 */
	@Override
	protected void execute() {
		this.isPluginInterrupted = false;
		if (validateInput() != 0) {
			return;
		}

		String srcPath = inSrcFile.getValue().getPath();
		String tgtPath = inTgtFile.getValue().getPath();
		System.out.println(srcPath);

		String transformedSrcPath;
		String transformedTgtPath;
		if (inSrcResultFile.getValue() == null) {
			transformedSrcPath = inSrcFile.getValue().getPath();
		} else {
			transformedSrcPath = inSrcResultFile.getValue().getPath();
		}
		
		if (inTgtResultFile.getValue() == null) {
			transformedTgtPath = inTgtFile.getValue().getPath();
		} else {
			transformedTgtPath = inTgtResultFile.getValue().getPath();
		}
		
		String srcResultPath;
		String tgtResultPath;

		srcResultPath = FilenameUtils.getFullPath(srcPath);
		srcResultPath += FilenameUtils.getBaseName(srcPath);
		srcResultPath += "_BUnwarp.";
		srcResultPath += FilenameUtils.getExtension(srcPath);
		
		tgtResultPath = FilenameUtils.getFullPath(tgtPath);
		tgtResultPath += FilenameUtils.getBaseName(tgtPath);
		tgtResultPath += "_BUnwarp.";
		tgtResultPath += FilenameUtils.getExtension(tgtPath);
		
		long startTime = System.nanoTime();
		
		ProgressBar.setProgressBarMessage("Loading source image...");
		
//		srcSeq = BigImageTools.loadSubsampledSequence(inSrcFile.getValue().getPath(), inSrcFile.getValue().getName());
//		Runtime.getRuntime().gc();
//		tgtSeq = BigImageTools.loadSubsampledSequence(inTgtFile.getValue().getPath(), inTgtFile.getValue().getName());
//		Runtime.getRuntime().gc();
		
		try {
			srcSeq = BigImageLoader.loadDownsampledImage(inSrcFile.getValue().getPath(), null, 1000, 1000);
			tgtSeq = BigImageLoader.loadDownsampledImage(inTgtFile.getValue().getPath(), null, 1000, 1000);
		} catch (UnsupportedFormatException | IOException e1) {
			e1.printStackTrace();
			return;
		}
		
		addSequence(srcSeq);
		addSequence(tgtSeq);

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
		bu.start();
		try {
			bu.join();
		} catch (InterruptedException e) {
			System.err.println("Thread interrupted: " + e.getMessage());
		}
		if (!this.isPluginInterrupted) {
			Sequence result = bu.getRegisteredSource(srcResultPath, transformedSrcPath, tgtPath);
			addSequence(result);
			if (inMode.getValue() != RegistrationModeEnum.MONO) {
				Sequence result1 = bu.getRegisteredTarget(tgtResultPath, transformedTgtPath, srcPath);
				addSequence(result1);
			}
		}
		long endTime = System.nanoTime();
		long totalTime = (endTime - startTime);
		System.out.println(String.format("Done (%d millisecs)", totalTime / 1000000));
	}

	/**
	 * Validate plugin input variables
	 * 
	 * @return 0 if input is valid, else a positive number with the error code.
	 */
	private int validateInput() {
		@SuppressWarnings("resource")
		LociImporterPlugin p = new LociImporterPlugin();
		if (!p.acceptFile(inSrcFile.getValue().getPath())) {
			MessageDialog.showDialog("Error", "Invalid source file.", MessageDialog.ERROR_MESSAGE);
			return 1;
		}
		if (!p.acceptFile(inTgtFile.getValue().getPath())) {
			MessageDialog.showDialog("Error", "Invalid target file.", MessageDialog.ERROR_MESSAGE);
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
		isPluginInterrupted = true;
		if (bu != null && bu.isAlive()) {
			try {
				bu.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#clean()
	 */
	@Override
	public void clean() {
	}

	/*
	 * (non-Javadoc)
	 * 
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
