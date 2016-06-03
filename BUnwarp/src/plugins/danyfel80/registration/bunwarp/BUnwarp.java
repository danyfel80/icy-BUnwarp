/**
 * 
 */
package plugins.danyfel80.registration.bunwarp;

import plugins.adufour.blocks.lang.Block;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;

/**
 * Abstract Class for BUnwarp plugin
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public abstract class BUnwarp extends EzPlug implements Block, EzStoppable {

	protected boolean isPluginInterrupted;
	
	public boolean isPluginInterrumped() {
		return this.isPluginInterrupted;
	}

	public abstract void restoreAll();

}
