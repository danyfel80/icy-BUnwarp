/**
 * 
 */
package algorithms.danyfel80.registration.bunwarp;

import icy.sequence.Sequence;
import plugins.adufour.ezplug.EzPlug;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class BUnwarpper extends Thread {
	
	// Input
	private Sequence srcSeq;
	private Sequence tgtSeq;
	private RegistrationModeEnum mode;
	private int subSampleFactor;
	private MinimumScaleDeformationEnum minScaleDef;
	private MaximumScaleDeformationEnum maxScaleDef;
	private double divWeight;
	private double curlWeight;
	private double landmarkWeight;
	private double consistencyWeight;
	private double imageWeight;
	private double stopThreshold;
	private boolean showProcess;
	private EzPlug plugin;

	// Output
	
	// Internal
	
	public BUnwarpper(Sequence srcSeq, Sequence tgtSeq, RegistrationModeEnum mode, Integer subSampleFactor,
	    MinimumScaleDeformationEnum minScaleDef, MaximumScaleDeformationEnum maxScaleDef, Double divWeight, Double curlWeight,
	    Double landmarkWeight, Double consistencyWeight, Double imageWeight, Double stopThreshold, Boolean showProcess, EzPlug plugin) {
		this.srcSeq = srcSeq;
		this.tgtSeq = tgtSeq;
		this.mode = mode;
		this.subSampleFactor = subSampleFactor;
		this.minScaleDef = minScaleDef;
		this.maxScaleDef = maxScaleDef;
		this.divWeight = divWeight;
		this.curlWeight = curlWeight;
		this.landmarkWeight = landmarkWeight;
		this.consistencyWeight = consistencyWeight;
		this.imageWeight = imageWeight;
		this.stopThreshold = stopThreshold;
		this.showProcess = showProcess;
		this.plugin = plugin;
	}

	/* (non-Javadoc)
	 * @see java.lang.Thread#run()
	 */
	@Override
	public void run() {
		super.run();
	}

	public Sequence getDirectTransform() {
		// TODO Auto-generated method stub
		return null;
	}

	public Sequence getInverseTransform() {
		// TODO Auto-generated method stub
		return null;
	}

	
	
}
