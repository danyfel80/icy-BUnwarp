/**
 * 
 */
package plugins.danyfel80.bigimage;

import java.awt.Rectangle;
import java.io.IOException;

import algorithms.danyfel80.bigimage.BigImageLoader;
import icy.common.exception.UnsupportedFormatException;
import icy.sequence.Sequence;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.vars.lang.VarBoolean;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class LoadBigImage extends EzPlug implements Block {

	private EzVarFile inFile = new EzVarFile("Image path", "");
	// Downsampling
	private EzVarInteger inMaxWidth = new EzVarInteger("Max width");
	private EzVarInteger inMaxHeight = new EzVarInteger("Max height");

	// Tiling
	private EzVarBoolean inIsTiled = new EzVarBoolean("Load Tile", false);
	private EzVarInteger inTileX = new EzVarInteger("Tile x");
	private EzVarInteger inTileY = new EzVarInteger("Tile y");
	private EzVarInteger inTileW = new EzVarInteger("Tile width");
	private EzVarInteger inTileH = new EzVarInteger("Tile height");

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.
	 * VarList)
	 */
	@Override
	public void declareInput(VarList inputMap) {
		inputMap.add(inFile.name, inFile.getVariable());
		inputMap.add(inMaxWidth.name, inMaxWidth.getVariable());
		inputMap.add(inMaxHeight.name, inMaxHeight.getVariable());
		inputMap.add(inIsTiled.name, inIsTiled.getVariable());
		inputMap.add(inTileX.name, inTileX.getVariable());
		inputMap.add(inTileY.name, inTileY.getVariable());
		inputMap.add(inTileW.name, inTileW.getVariable());
		inputMap.add(inTileH.name, inTileH.getVariable());
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
		boolean isTiled = inIsTiled.getValue();
		int tileX = inTileX.getValue();
		int tileY = inTileY.getValue();
		int tileW = inTileW.getValue();
		int tileH = inTileH.getValue();

		Sequence s;
		try {
			long startTime = System.nanoTime();
			BigImageLoader.setPluginGUI(this.getUI());
			s = BigImageLoader.loadDownsampledImage(path, isTiled? new Rectangle(tileX, tileY, tileW, tileH): null, maxWidth, maxHeight);
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

		EzGroup downsamplingGroup = new EzGroup("Downsampling", inMaxWidth, inMaxHeight);
		addEzComponent(downsamplingGroup);

		inIsTiled.addVarChangeListener(new EzVarListener<Boolean>() {
			@Override
			public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
				inTileX.setVisible(newValue);
				inTileY.setVisible(newValue);
				inTileW.setVisible(newValue);
				inTileH.setVisible(newValue);
			}
		});
		EzGroup tileGroup = new EzGroup("Tiles", inIsTiled, inTileX, inTileY, inTileW, inTileH);
		addEzComponent(tileGroup);
		inIsTiled.setValue(false);
	}

}
