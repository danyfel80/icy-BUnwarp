package algorithms.danyfel80.registration.bunwarp;

import java.util.List;

import icy.roi.ROI2D;
import plugins.danyfel80.registration.bunwarp.BUnwarp;
import plugins.kernel.roi.roi2d.ROI2DPoint;

/**
 * Class to perform BUnwarp on big images
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigBUnwarpper extends Thread{
	private String srcPath;
	private String tgtPath;
	private String transformedSrcPath;
	private String transformedTgtPath;
	private String srcResultPath;
	private String tgtResultPath;
	
	private List<ROI2DPoint> srcLandmarks;
	private List<ROI2DPoint> tgtLandmarks;
	
	private ROI2D srcMask;
	private ROI2D tgtMask;
	
	private int subsampleFactor;
	private double[] usedScales;
	private int initialDeformation;
	private int finalDeformation;
  
	private double divWeight;
	private double curlWeight;
	private double landmarkWeight;
	private double imageWeight;
	private double consistencyWeight;
  
	private double stopThreshold;
  
	private boolean showProcess;
	private int mode;
	private BUnwarp plugin;
	
	/**
	 * Constructor
	 * @param srcPath Source image path. Used for registering.
	 * @param tgtPath Target image path. Used for registering.
	 * @param transformedSrcPath Transformed source image path. Used to get final results.
	 * @param transformedTgtPath Transformed target image path. Used to get final results.
	 * @param srcResultPath Final result source image path. 
	 * @param tgtResultPath Final result target image path.
	 * @param srcLandmarks Source image landmarks.
	 * @param tgtLandmarks Source image landmarks.
	 * @param srcMask Source mask. Indicates the area to register.
	 * @param tgtMask Target mask. Indicates the area to register.
	 * @param subsampleFactor Subsample factor of the image.
	 * @param usedScales Used scales to perform the registration.
	 * @param initialDeformation Initial deformation detail. (0 = Very coarse, 1 = Coarse, 2 = Fine, 3 = Very fine)
	 * @param finalDeformation Final deformation detail. (0 = Very coarse, 1 = Coarse, 2 = Fine, 3 = Very fine, 4 = Super fine)
	 * @param divWeight Divergence weight.
	 * @param curlWeight Curl weight.
	 * @param landmarkWeight Landmark weight.
	 * @param imageWeight Image weight.
	 * @param consistencyWeight Transformation's consistency weight.
	 * @param stopThreshold Threshold to stop registering.
	 * @param showProcess Flag to show intermediate results and processing log.
	 * @param mode Mode of registration. (0 = Fast, 1 = Accurate, 2 = Mono)
	 * @param plugin Reference to the BUnwarp plugin
	 */
	public BigBUnwarpper(String srcPath, String tgtPath, String transformedSrcPath, String transformedTgtPath,
	    String srcResultPath, String tgtResultPath, List<ROI2DPoint> srcLandmarks, List<ROI2DPoint> tgtLandmarks,
	    ROI2D srcMask, ROI2D tgtMask, int subsampleFactor, double[] usedScales, int initialDeformation, int finalDeformation, double divWeight,
	    double curlWeight, double landmarkWeight, double imageWeight, double consistencyWeight, double stopThreshold,
	    boolean showProcess, int mode, BUnwarp plugin) {
		this.srcPath = srcPath;
		this.tgtPath = tgtPath;
		this.transformedSrcPath = transformedSrcPath;
		this.transformedTgtPath = transformedTgtPath;
		this.srcResultPath = srcResultPath;
		this.tgtResultPath = tgtResultPath;
		this.srcLandmarks = srcLandmarks;
		this.tgtLandmarks = tgtLandmarks;
		this.srcMask = srcMask;
		this.tgtMask = tgtMask;
		this.subsampleFactor = subsampleFactor;
		this.usedScales = usedScales;
		this.initialDeformation = initialDeformation;
		this.finalDeformation = finalDeformation;
		this.divWeight = divWeight;
		this.curlWeight = curlWeight;
		this.landmarkWeight = landmarkWeight;
		this.imageWeight = imageWeight;
		this.consistencyWeight = consistencyWeight;
		this.stopThreshold = stopThreshold;
		this.showProcess = showProcess;
		this.mode = mode;
		this.plugin = plugin;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Thread#run()
	 */
	@Override
	public void run() {
		
	}
	
	
}
