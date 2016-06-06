/**
 * 
 */
package plugins.danyfel80.bigimage;

import java.io.IOException;

import algorithms.danyfel80.bigimage.BigImageLoader;
import icy.common.exception.UnsupportedFormatException;
import icy.sequence.Sequence;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class LoadBigImage extends EzPlug implements Block {

	private EzVarFile inFile = new EzVarFile("Image path", "");
	private EzVarInteger inMaxWidth = new EzVarInteger("Max width");
	private EzVarInteger inMaxHeight = new EzVarInteger("Max height");

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.
	 * VarList)
	 */
	@Override
	public void declareInput(VarList inputMap) {
		// TODO Auto-generated method stub

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
	 * @see plugins.adufour.ezplug.EzPlug#clean()
	 */
	@Override
	public void clean() {
		// TODO Auto-generated method stub

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#execute()
	 */
	@Override
	protected void execute() {
		String path = inFile.getValue().getAbsolutePath();
		int maxWidth = inMaxWidth.getValue();
		int maxHeight = inMaxHeight.getValue();

		Sequence s;
		try {
			long startTime = System.nanoTime();
			s = BigImageLoader.loadDownsampledImage(path, maxWidth, maxHeight);
			long endTime = System.nanoTime();
			addSequence(s);
			System.out.println("Loaded in " + ((endTime - startTime) / 1000000) + "msecs.");

		} catch (UnsupportedFormatException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see plugins.adufour.ezplug.EzPlug#initialize()
	 */
	@Override
	protected void initialize() {
		addEzComponent(inFile);
		addEzComponent(inMaxWidth);
		addEzComponent(inMaxHeight);
	}

}
