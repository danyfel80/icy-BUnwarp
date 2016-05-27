package algorithms.danyfel80.registration.bunwarp;

import java.util.List;

import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import plugins.adufour.ezplug.EzPlug;
import plugins.danyfel80.registration.bunwarp.BUnwarp;
import plugins.kernel.roi.roi2d.ROI2DPoint;
import plugins.kernel.roi.roi2d.ROI2DPolygon;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class BUnwarpper extends Thread {
	// Images
	/** image representation for the source */
	private Sequence sourceSeq;
	/** image representation for the target */
	private Sequence targetSeq;
	/** source image model */
	private BSplineModel sourceModel;
	/** target image model */
	private BSplineModel targetModel;

	// Landmarks
	/** point handler for the landmarks in the source image */
	private List<ROI2DPoint> sourceLandmarks;
	/** point handler for the landmarks in the target image */
	private List<ROI2DPoint> targetLandmarks;

	// Masks for the images
	/** source image mask */
	private ROI2DPolygon sourceMask;
	/** target image mask */
	private ROI2DPolygon targetMask;

	// Transformation parameters
	/** maximum image subsampling factor */
	private int maxImageSubsamplingFactor;
	/** minimum scale deformation */
	private int minScaleDeformation;
	/** maximum scale deformation */
	private int maxScaleDeformation;
	/** minimum image scale */
	private int minScaleImage;
	/** flag to specify the level of resolution in the output */
	private int outputLevel;
	/** flag to show the optimizer */
	private boolean showMarquardtOptim;
	/** divergence weight */
	private double divWeight;
	/** curl weight */
	private double curlWeight;
	/** landmark weight */
	private double landmarkWeight;
	/** weight for image similarity */
	private double imageWeight;
	/** weight for the deformations consistency */
	private double consistencyWeight;
	/** stopping threshold */
	private double stopThreshold;
	/** level of accuracy */
	private int accurateMode;
	/** maximum depth for the image pyramid */
	private int imagePyramidDepth;

	// Ez Plugin
	/** Ez Plugin reference */
	private EzPlug plugin;
	/*
	 * List<ROI2DPoint> srcLandmarks, List<ROI2DPoint> tgtLandmarks, ROI2DPolygon
	 * srcMask, ROI2DPolygon tgtMask, RegistrationModeEnum mode, Integer
	 * maxSubsamplingFactor, MinimumScaleDeformationEnum minScaleDef,
	 * MaximumScaleDeformationEnum maxScaleDef, Double divWeight, Double
	 * curlWeight, Double landmarkWeight, Double consistencyWeight, Double
	 * imageWeight, Double stopThreshold, Boolean showProcess, EzPlug plugin
	 */
	public BUnwarpper(final Sequence sourceSequence, final Sequence targetSequence,
	    final List<ROI2DPoint> sourceLandmarks, final List<ROI2DPoint> targetLandmarks, final ROI2DPolygon sourceMask,
	    final ROI2DPolygon targetMask, final int maxImageSubsamplingFactor, final int minScaleDeformation,
	    final int maxScaleDeformation, final int minScaleImage, final double divWeight, final double curlWeight,
	    final double landmarkWeight, final double imageWeight, final double consistencyWeight, final double stopThreshold,
	    final int outputLevel, final boolean showMarquardtOptim, final int accurateMode, final EzPlug plugin) {
		this.sourceSeq = sourceSequence;
		this.targetSeq = targetSequence;
		this.sourceLandmarks = sourceLandmarks;
		this.targetLandmarks = targetLandmarks;
		this.sourceMask = sourceMask;
		this.targetMask = targetMask;
		this.maxImageSubsamplingFactor = maxImageSubsamplingFactor;
		this.minScaleDeformation = minScaleDeformation;
		this.maxScaleDeformation = maxScaleDeformation;
		this.minScaleImage = minScaleImage;
		this.divWeight = divWeight;
		this.curlWeight = curlWeight;
		this.landmarkWeight = landmarkWeight;
		this.imageWeight = imageWeight;
		this.consistencyWeight = consistencyWeight;
		this.stopThreshold = stopThreshold;
		this.outputLevel = outputLevel;
		this.showMarquardtOptim = showMarquardtOptim;
		this.accurateMode = accurateMode;
		this.plugin = plugin;

		createSourceImage(this.accurateMode < RegistrationModeEnum.MONO.getNumber());
		createTargetImage();
	}

	private void createSourceImage(boolean isReverse) {
		sourceModel = new BSplineModel(sourceSeq.getFirstImage(), isReverse, (int) Math.pow(2, maxImageSubsamplingFactor));
		computeImagePyramidDepth();
		sourceModel.setPyramidDepth(imagePyramidDepth + minScaleImage);
	}

	private void computeImagePyramidDepth() {
		imagePyramidDepth = maxScaleDeformation - minScaleDeformation + 1;
	}

	private void createTargetImage() {
		targetModel = new BSplineModel(targetSeq.getFirstImage(), true, (int) Math.pow(2, maxImageSubsamplingFactor));
		computeImagePyramidDepth();
		targetModel.setPyramidDepth(imagePyramidDepth + minScaleImage);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Thread#run()
	 */
	@Override
	public void run() {
		super.run();

		// Start pyramids
		System.out.println("Starting image pyramids...");
		if (targetModel.getWidth() > BSplineModel.MAX_OUTPUT_SIZE || targetModel.getHeight() > BSplineModel.MAX_OUTPUT_SIZE
		    || sourceModel.getWidth() > BSplineModel.MAX_OUTPUT_SIZE || sourceModel.getHeight() > BSplineModel.MAX_OUTPUT_SIZE)
			System.out.println("Starting image pyramids...");

		sourceModel.startPyramids();
		targetModel.startPyramids();

		try {
			sourceModel.join();
			targetModel.join();
		} catch(InterruptedException e) {
			System.err.println("Unhandled interruption: " + e);
		}
		
		// Create output image (source-target)
		final Sequence[] outputSeqs = initializeOutputSeqs();

		// If mono mode, reset consistency weight
		if (this.accurateMode == RegistrationModeEnum.MONO.getNumber())
			this.consistencyWeight = 0.0;

		// IJ.log("FinalAction: maxImageSubsamplingFactor = " +
		// maxImageSubsamplingFactor);

		// Prepare registration parameters
		final Transformation warp = new Transformation(sourceSeq, targetSeq, sourceModel, targetModel, sourceLandmarks, targetLandmarks, sourceMask,
		    targetMask, minScaleDeformation, maxScaleDeformation,
		    minScaleImage, divWeight, curlWeight, landmarkWeight, imageWeight, consistencyWeight, stopThreshold,
		    outputLevel, showMarquardtOptim, accurateMode, outputSeqs[0], outputSeqs[1], plugin);

		// Perform the registration
		System.out.println("Registering...");

		long start = System.currentTimeMillis(); // start timing

		if (accurateMode == RegistrationModeEnum.MONO.getNumber()) {
			// Do unidirectional registration
			warp.doUnidirectionalRegistration();
			// Show results.
			if (sourceModel.isSubOutput()) {
				System.out.println("Calculating final transformed source image");
			}
			warp.showDirectResults();

		} else {
			// Do bidirectional registration
			warp.doBidirectionalRegistration();
			if (sourceModel.isSubOutput()) {
				System.out.println("Calculating final transformed source image");
			}
			warp.showDirectResults();
			if (targetModel.isSubOutput()) {
				System.out.println("Calculating final transformed target image");
			}
			warp.showInverseResults();
		}

		long stop = System.currentTimeMillis(); // stop timing
		if (outputLevel == 2)
			System.out.println("\nRegistration time: " + (stop - start) + "ms"); // print
		
		((BUnwarp)plugin).restoreAll();
	}

	private Sequence[] initializeOutputSeqs() {
		int Ydimt = targetModel.getHeight();
    int Xdimt = targetModel.getWidth();
    int Xdims = sourceModel.getWidth();
    int Ydims = sourceModel.getHeight();
    double[] tImage = targetModel.isSubOutput() ? targetModel.getSubImage() : targetModel.getImage();
    double[] sImage = sourceModel.isSubOutput() ? sourceModel.getSubImage() : sourceModel.getImage(); 
    int sSubFactorX = 1;
    int sSubFactorY = 1;
    int tSubFactorX = 1;
    int tSubFactorY = 1;
    Sequence[] outputSeqs = new Sequence[2];
    	
    String extraTitleS = "";
    String extraTitleT = "";

    if(targetModel.isSubOutput() || sourceModel.isSubOutput())
    	System.out.println("Initializing output windows...");
    
    // If the output (difference) images are subsampled (because they were
    // larger than the maximum size), update variables.
    if(targetModel.isSubOutput())
    {        	
    	tSubFactorX = Xdimt / targetModel.getSubWidth();
    	tSubFactorY = Ydimt / targetModel.getSubHeight();
    	extraTitleT = " (Subsampled)";
    	Xdimt = targetModel.getSubWidth();
    	Ydimt = targetModel.getSubHeight();        	          	        	       			        	
    }
    
    if(sourceModel.isSubOutput())
	{
		sSubFactorX = Xdims / sourceModel.getSubWidth();
    	sSubFactorY = Ydims / sourceModel.getSubHeight();
    	extraTitleS = " (Subsampled)";
		Xdims = sourceModel.getSubWidth();
		Ydims = sourceModel.getSubHeight();
	} 
    
    // Float processor for the output source-target image.
    final IcyBufferedImage ibi = new IcyBufferedImage(Xdimt, Ydimt, 1, DataType.FLOAT);
    float[] ibiData = ibi.getDataXYAsFloat(0);               
                    
    for (int i=0; i<Ydimt; i++)
	{
		final int i_offset_t = i * Xdimt; 
		final int i_offset_s = i * Xdims; 
		final int i_s_sub = i * sSubFactorY;
		final int i_t_sub = i * tSubFactorY;
		
		for (int j=0; j<Xdimt; j++)
		{
			    				
			if (sourceMask.contains(j * sSubFactorX, i_s_sub) && targetMask.contains(j * tSubFactorX, i_t_sub)
					&& j < Xdims && i < Ydims)
				ibiData[j + i_offset_t] = (float) (tImage[i_offset_t + j] - sImage[i_offset_s + j]);
			else
			{
				ibiData[j + i_offset_t] = 0;
			}

		}
	}
	ibi.dataChanged();

    final Sequence seq1 = new Sequence("Output Source-Target" + extraTitleS, ibi);
    plugin.addSequence(seq1);
    
    outputSeqs[0] = seq1;

    // Create output image (target-source) if necessary                        
    
    if(this.accurateMode != RegistrationModeEnum.MONO.getNumber())
    {
    	final IcyBufferedImage ibi2 = new IcyBufferedImage(Xdims, Ydims, 1, DataType.FLOAT);
    	float[] ibi2Data = ibi2.getDataXYAsFloat(0);

    	for (int i=0; i<Ydims; i++)
    	{
    		int i_offset_t = i * Xdimt; 
    		int i_offset_s = i * Xdims; 
    		int i_s_sub = i * sSubFactorY;
    		int i_t_sub = i * tSubFactorY;
    		
    		for (int j=0; j<Xdims; j++)
    			if (targetMask.contains(j * tSubFactorX, i_t_sub) && sourceMask.contains(j * sSubFactorX, i_s_sub)
    					&& i < Ydimt && j < Xdimt)
    				ibi2Data[j + i_offset_s] = (float) (sImage[i_offset_s + j] - tImage[i_offset_t + j]);
    			else 
    				ibi2Data[j + i_offset_s] = 0;
    	}
    	ibi2.dataChanged();
    	
    	
    	final Sequence seq2 = new Sequence("Output Target-Source" + extraTitleT, ibi2);
    	plugin.addSequence(seq2);
    	
    	outputSeqs[1] = seq2;
    }
    else
    	outputSeqs[1] = null;
    
    return outputSeqs;
	}

}
