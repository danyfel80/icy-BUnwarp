package plugins.danyfel80.bigimage.thresholding;

import java.io.File;

import algorithms.danyfel80.bigimage.thresholding.BigImageThresholder;
import icy.gui.dialog.MessageDialog;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarFile;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageThresholding extends EzPlug implements EzStoppable {

	EzVarFile inInputFile = new EzVarFile("Input file", "");
	EzVarFile inOutputFile = new EzVarFile("Output file", "");
	EzVarDouble inMinValue = new EzVarDouble("Minimum value", 0, 0, 255, Double.MIN_VALUE);
	EzVarDouble inMaxValue = new EzVarDouble("Maximum value", 1, 0, 255, Double.MIN_VALUE);

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#initialize()
	 */
	@Override
	protected void initialize() {
		addEzComponent(inInputFile);
		addEzComponent(inOutputFile);
		addEzComponent(inMinValue);
		addEzComponent(inMaxValue);
	}

	BigImageThresholder thresholder;

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
		System.out.println("valid input");
		File inputFile = inInputFile.getValue();
		File outputFile = inOutputFile.getValue();
		double minValue = inMinValue.getValue();
		double maxValue = inMaxValue.getValue();

		
		thresholder = new BigImageThresholder(inputFile, outputFile, minValue, maxValue);
		System.out.println("Thresholder created");
		thresholder.execute();
		System.out.println("Thresholder executed");
	}

	private int validateInput() {
		if (inInputFile.getValue() == null || !inInputFile.getValue().exists()) {
			MessageDialog.showDialog("Error", "Select a valid input file. (Error: -1)", MessageDialog.ERROR_MESSAGE);
			return -1;
		}
		if (inOutputFile.getValue() == null || inOutputFile.getValue().getParentFile() == null
		    || !inOutputFile.getValue().getParentFile().exists()) {
			MessageDialog.showDialog("Error", "Select a valid output file. (Error: -2)", MessageDialog.ERROR_MESSAGE);
			return -2;
		}
		if (inMinValue.getValue() < 0 || inMinValue.getValue() > 255) {
			MessageDialog.showDialog("Error", "Select a valid min value between 0 and 255. (Error: -3)",
			    MessageDialog.ERROR_MESSAGE);
			return -3;
		}
		if (inMaxValue.getValue() < 0 || inMaxValue.getValue() > 255) {
			MessageDialog.showDialog("Error", "Select a valid max value between 0 and 255. (Error: -4)",
			    MessageDialog.ERROR_MESSAGE);
			return -4;
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
		if (thresholder != null) {
			thresholder.stopExecution();
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

}
