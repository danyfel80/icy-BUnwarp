package algorithms.danyfel80.registration.bunwarp;

import java.io.IOException;
import java.util.List;

import algorithms.danyfel80.bigimage.BigImageLoader;
import icy.common.exception.UnsupportedFormatException;
import icy.roi.ROI2D;
import icy.sequence.Sequence;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import plugins.danyfel80.registration.bunwarp.BUnwarp;
import plugins.kernel.roi.roi2d.ROI2DPoint;

/**
 * Class to perform BUnwarp on big images
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigBUnwarpper extends Thread {
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
	 * 
	 * @param srcPath
	 *          Source image path. Used for registering.
	 * @param tgtPath
	 *          Target image path. Used for registering.
	 * @param transformedSrcPath
	 *          Transformed source image path. Used to get final results.
	 * @param transformedTgtPath
	 *          Transformed target image path. Used to get final results.
	 * @param srcResultPath
	 *          Final result source image path.
	 * @param tgtResultPath
	 *          Final result target image path.
	 * @param srcLandmarks
	 *          Source image landmarks.
	 * @param tgtLandmarks
	 *          Source image landmarks.
	 * @param srcMask
	 *          Source mask. Indicates the area to register.
	 * @param tgtMask
	 *          Target mask. Indicates the area to register.
	 * @param subsampleFactor
	 *          Subsample factor of the image.
	 * @param usedScales
	 *          Used scales to perform the registration.
	 * @param initialDeformation
	 *          Initial deformation detail. (0 = Very coarse, 1 = Coarse, 2 =
	 *          Fine, 3 = Very fine)
	 * @param finalDeformation
	 *          Final deformation detail. (0 = Very coarse, 1 = Coarse, 2 = Fine,
	 *          3 = Very fine, 4 = Super fine)
	 * @param divWeight
	 *          Divergence weight.
	 * @param curlWeight
	 *          Curl weight.
	 * @param landmarkWeight
	 *          Landmark weight.
	 * @param imageWeight
	 *          Image weight.
	 * @param consistencyWeight
	 *          Transformation's consistency weight.
	 * @param stopThreshold
	 *          Threshold to stop registering.
	 * @param showProcess
	 *          Flag to show intermediate results and processing log.
	 * @param mode
	 *          Mode of registration. (0 = Fast, 1 = Accurate, 2 = Mono)
	 * @param plugin
	 *          Reference to the BUnwarp plugin
	 */
	public BigBUnwarpper(String srcPath, String tgtPath, String transformedSrcPath, String transformedTgtPath,
	    String srcResultPath, String tgtResultPath, List<ROI2DPoint> srcLandmarks, List<ROI2DPoint> tgtLandmarks,
	    ROI2D srcMask, ROI2D tgtMask, int subsampleFactor, double[] usedScales, int initialDeformation,
	    int finalDeformation, double divWeight, double curlWeight, double landmarkWeight, double imageWeight,
	    double consistencyWeight, double stopThreshold, boolean showProcess, int mode, BUnwarp plugin) {
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Thread#run()
	 */
	@Override
	public void run() {
		ProgressBar.setPlugin(plugin);
		ProgressBar.setProgressBarMessage("Loading source image...");

		Sequence srcSeq;
		Sequence tgtSeq;
//		Sequence srcTgtSeq;
//		Sequence tgtTgtSeq;
//		Dimension srcDim = BigImageTools.getSequenceSize(srcPath);
		// Dimension tgtDim = BigImageTools.getSequenceSize(tgtPath);

		// ---- First Scale Registration
		// Get sequences
		try {
			BigImageLoader.setPluginGUI(plugin.getUI());
			ProgressBar.setProgressBarMessage("Loading source image");
			srcSeq = BigImageLoader.loadDownsampledImage(srcPath, null, 1000, 1000);
			ProgressBar.setProgressBarMessage("Loading target image");
			tgtSeq = BigImageLoader.loadDownsampledImage(tgtPath, null, 1000, 1000);
			
//			srcTgtSeq = SequenceUtil.getCopy(srcSeq);
//			tgtTgtSeq = SequenceUtil.getCopy(tgtSeq);
		} catch (UnsupportedFormatException | IOException e1) {
			e1.printStackTrace();
			return;
		}

		// Register images
		BUnwarpper bu = new BUnwarpper(srcSeq, tgtSeq, srcLandmarks, tgtLandmarks, srcMask, tgtMask, subsampleFactor,
		    initialDeformation, finalDeformation, 0, divWeight, curlWeight, landmarkWeight, imageWeight, consistencyWeight,
		    stopThreshold, showProcess ? 2 : 1, showProcess, mode, plugin);
		bu.start();
		try {
			bu.join();
		} catch (InterruptedException e) {
			System.err.println("Thread interrupted: " + e.getMessage());
			return;
		}

		// Show results

//		bu.getRegisteredSource(srcTgtSeq);
//		Icy.getMainInterface().addSequence(srcTgtSeq);
		try {
			System.out.println("saving to " + srcResultPath + ", based on " + transformedSrcPath);
			bu.saveBigRegisteredSource(srcResultPath, transformedSrcPath, tgtPath, null);
		} catch (ServiceException | IOException | FormatException | InterruptedException e) {
			e.printStackTrace();
			return;
		}

		if (mode != RegistrationModeEnum.MONO.getNumber()) {
//			bu.getRegisteredTarget(tgtTgtSeq);
//			Icy.getMainInterface().addSequence(tgtTgtSeq);
			try {
				bu.saveBigRegister1edTarget(tgtResultPath, transformedTgtPath, srcPath, null);
			} catch (ServiceException | IOException | FormatException | InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		// Sequence result = bu.getRegisteredSource(srcResultPath,
		// transformedSrcPath, tgtPath);
		// Icy.getMainInterface().addSequence(result);
		// if (mode != RegistrationModeEnum.MONO.getNumber()) {
		// Sequence result1 = bu.getRegisteredTarget(tgtResultPath,
		// transformedTgtPath, srcPath);
		// Icy.getMainInterface().addSequence(result1);
		// }

		bu = null;
		System.gc();
		return;
//		// Next Scales Registration
//		for (int si = 0; si < usedScales.length; si++) {
//			double scale = usedScales[si];
//			int tileAmount = (int) Math.round(1.0 / scale);
//			Dimension tileDim = new Dimension(srcDim.width / tileAmount, srcDim.height / tileAmount);
//			Dimension tileSize = new Dimension(tileAmount + (srcDim.width % tileAmount > 0 ? 1 : 0),
//			    tileAmount + (srcDim.height % tileAmount > 0 ? 1 : 0));
//			int tileBorderSize = Math.max(tileDim.width, tileDim.height) / 8;
//			int nProc = Runtime.getRuntime().availableProcessors();
//			BUnwarpper[] bus = new BUnwarpper[nProc];
//			Rectangle[] rects = new Rectangle[nProc];
//
//			int usedThreads = 0;
//			int processedTiles = 0;
//			for (int i = 0; i < tileSize.width && !plugin.isPluginInterrumped(); i++) {
//				for (int j = 0; j < tileSize.height && !plugin.isPluginInterrumped(); j++) {
//					try {
//						rects[usedThreads] = new Rectangle(i * tileDim.width - tileBorderSize, j * tileDim.height - tileBorderSize,
//						    tileSize.width + tileBorderSize, tileSize.height + tileBorderSize);
//						srcSeq = BigImageLoader.loadDownsampledImage(srcResultPath, rects[usedThreads], 1023, 1023);
//						tgtSeq = BigImageLoader.loadDownsampledImage(tgtPath, rects[usedThreads], 1023, 1023);
//						srcTgtSeq = SequenceUtil.getCopy(srcSeq);
//						tgtTgtSeq = SequenceUtil.getCopy(tgtSeq);
//					} catch (UnsupportedFormatException | IOException e1) {
//						e1.printStackTrace();
//						return;
//					}
//
//					// Register images
//
//					bus[usedThreads] = new BUnwarpper(srcSeq, tgtSeq, srcLandmarks, tgtLandmarks, srcMask, tgtMask,
//					    subsampleFactor, initialDeformation, finalDeformation, 0, divWeight, curlWeight, landmarkWeight,
//					    imageWeight, consistencyWeight, stopThreshold, showProcess ? 2 : 1, showProcess,
//					    RegistrationModeEnum.MONO.getNumber(), plugin);
//					bus[usedThreads++].start();
//
//					if (usedThreads >= nProc || processedTiles == tileSize.width * tileSize.height) {
//						for (int t = 0; t < usedThreads; t++) {
//							try {
//								bus[t].join();
//							} catch (InterruptedException e) {
//								System.err.println("Thread interrupted: " + e.getMessage());
//								return;
//							}
//							// TODO Save registered tile
//							//bus[t].saveRegisteredSource(srcResultPath, transformedSrcPath, tgtPath, rects[t]);
//							bus[t].getRegisteredSource(srcTgtSeq);
//							Icy.getMainInterface().addSequence(srcTgtSeq);
//						}
//						usedThreads = 0;
//					}
//				}
//			}
//
//			// TODO inverse registration
//		}

	}

}
