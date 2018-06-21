/**
 * 
 */
package plugins.danyfel80.registration.bunwarp;

import plugins.adufour.blocks.lang.Block;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;

/**
 * Abstract Class for BUnwarp plugin
 * 
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public abstract class BUnwarp extends EzPlug implements Block, EzStoppable {

	protected void notifyProgress(double progress, String message) {
		if (!isHeadLess()) {
			getUI().setProgressBarMessage(message);
			getUI().setProgressBarValue(progress);
		}
	};

	public abstract void restoreAll();

	@Override
	public void clean() {
		// Nothing to do here
	}

}
