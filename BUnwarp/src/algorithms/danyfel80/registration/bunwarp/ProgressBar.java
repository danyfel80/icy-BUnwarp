/**
 * 
 */
package algorithms.danyfel80.registration.bunwarp;

import plugins.adufour.ezplug.EzPlug;

/**
 * This class implements the interactions when dealing with the plugin progress
 * bar.
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class ProgressBar {
	
	private static final long TIME_QUANTUM = 50L;

  private static volatile long lastTime = System.currentTimeMillis();
  private static volatile int completed = 0;
  private static volatile int workload = 0;
  private static EzPlug plugin;
  
  public static synchronized void setPlugin(EzPlug plugin) {
  	ProgressBar.plugin = plugin;
  }
  
//  public static synchronized EzPlug getPlugin() {
//  	return ProgressBar.plugin;
//  }
  
  /**
   * Extend the amount of work to perform by <code>batch</code>.
   *
   * @param batch Additional amount of work that need be performed.
   */
  public static synchronized void addWorkload (final int batch)
  {
     workload += batch;
  } /* end addWorkload */

  /**
   * Erase the progress bar and cancel pending operations.
   */
  public static synchronized void resetProgressBar ()
  {
     final long timeStamp = System.currentTimeMillis();
     if ((timeStamp - lastTime) < TIME_QUANTUM) {
        try {
           Thread.sleep(TIME_QUANTUM - timeStamp + lastTime);
        } catch (InterruptedException e) {
           System.err.println("Unexpected interruption exception" + e);
        }
     }
     lastTime = timeStamp;
     completed = 0;
     workload = 0;
     if (plugin != null && plugin.getUI() != null) {
    	 plugin.getUI().setProgressBarValue(1.0);
     }
  }

  /**
   * Perform <code>stride</code> operations at once.
   *
   * @param stride Amount of work that is skipped.
   */
  public static synchronized void skipProgressBar (final int stride)
  {
     completed += stride - 1;
     stepProgressBar();
  }

  /**
   * Perform <code>1</code> operation unit.
   */
  public static synchronized void stepProgressBar ()
  {
     final long timeStamp = System.currentTimeMillis();
     completed = completed + 1;
     if ((TIME_QUANTUM <= (timeStamp - lastTime)) | (completed == workload)) {
        lastTime = timeStamp;
        if (plugin != null && plugin.getUI() != null) {
        	plugin.getUI().setProgressBarValue((double)completed / (double)workload);
        }
     }
  }

  /**
   * Acknowledge that <code>batch</code> work has been performed.
   *
   * @param batch Completed amount of work.
   */
  public static synchronized void workloadDone (final int batch)
  {
     workload -= batch;
     completed -= batch;
  }
  
  public static synchronized void setProgressBarMessage (final String message) {
  	if (plugin != null && plugin.getUI() != null)
  		plugin.getUI().setProgressBarMessage(message);
  }
  
  public static synchronized void setProgressBarValue (final double value) {
  	if (plugin != null && plugin.getUI() != null)
  		plugin.getUI().setProgressBarValue(value);
  }
  
}
