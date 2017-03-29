package algorithms.danyfel80.registration.bunwarp;

import java.awt.Rectangle;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;

import algorithms.danyfel80.bigimage.BigImageLoader;
import icy.common.exception.UnsupportedFormatException;
import icy.main.Icy;
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
public class BigBUnwarpper implements Runnable {
	private String	srcPath;
	private String	tgtPath;
	private String	transformedSrcPath;
	private String	transformedTgtPath;

	private String	srcResultPath;
	private String	tgtResultPath;
	private String	transformedSrcResultPath;
	private String	transformedTgtResultPath;

	private List<ROI2DPoint>	srcLandmarks;
	private List<ROI2DPoint>	tgtLandmarks;

	private ROI2D	srcMask;
	private ROI2D	tgtMask;

	private int	subsampleFactor;
	private int	initialDeformation;
	private int	finalDeformation;

	private double	divWeight;
	private double	curlWeight;
	private double	landmarkWeight;
	private double	imageWeight;
	private double	consistencyWeight;

	private double stopThreshold;

	private boolean	showProcess;
	private int			mode;
	private BUnwarp	plugin;

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
	 * @param transformedSrcResultPath
	 *          Final result source image path.
	 * @param transformedTgtResultPath
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
			String srcResultPath, String tgtResultPath, String transformedSrcResultPath, String transformedTgtResultPath,
			List<ROI2DPoint> srcLandmarks, List<ROI2DPoint> tgtLandmarks, ROI2D srcMask, ROI2D tgtMask, int subsampleFactor,
			int initialDeformation, int finalDeformation, double divWeight, double curlWeight, double landmarkWeight,
			double imageWeight, double consistencyWeight, double stopThreshold, boolean showProcess, int mode,
			BUnwarp plugin) {
		this.srcPath = srcPath;
		this.tgtPath = tgtPath;
		this.transformedSrcPath = transformedSrcPath;
		this.transformedTgtPath = transformedTgtPath;
		this.srcResultPath = srcResultPath;
		this.tgtResultPath = tgtResultPath;
		this.transformedSrcResultPath = transformedSrcResultPath;
		this.transformedTgtResultPath = transformedTgtResultPath;
		this.srcLandmarks = srcLandmarks;
		this.tgtLandmarks = tgtLandmarks;
		this.srcMask = srcMask;
		this.tgtMask = tgtMask;
		this.subsampleFactor = subsampleFactor;
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
	 * @see java.lang.Thread#run()
	 */
	@Override
	public void run() {
		ProgressBar.setPlugin(plugin);

		Sequence srcSeq;
		Sequence tgtSeq;
		// Sequence srcTgtSeq;
		// Sequence tgtTgtSeq;
		// Dimension srcDim = BigImageTools.getSequenceSize(srcPath);
		// Dimension tgtDim = BigImageTools.getSequenceSize(tgtPath);

		// ---- First Scale Registration
		// Get sequences
		try {
			BigImageLoader loader = new BigImageLoader();
			loader.setPluginGUI(plugin.getUI());
			ProgressBar.setProgressBarMessage("Loading source image");
			srcSeq = loader.loadDownsampledImage(srcPath, null, 1000, 1000, true);
			ProgressBar.setProgressBarMessage("Loading source mask");
			srcMask = loader.loadDownsampledMask(srcSeq, srcPath, null, 1000, 1000, true);
			ProgressBar.setProgressBarMessage("Loading target image");
			tgtSeq = loader.loadDownsampledImage(tgtPath, null, 1000, 1000, true);
			ProgressBar.setProgressBarMessage("Loading target mask");
			tgtMask = loader.loadDownsampledMask(tgtSeq, tgtPath, null, 1000, 1000, true);

			// srcTgtSeq = SequenceUtil.getCopy(srcSeq);
			// tgtTgtSeq = SequenceUtil.getCopy(tgtSeq);
		} catch (UnsupportedFormatException | IOException e1) {
			e1.printStackTrace();
			return;
		}

		// Register images

		BUnwarpper bu = new BUnwarpper(srcSeq, tgtSeq, srcLandmarks, tgtLandmarks, srcMask, tgtMask, subsampleFactor,
				initialDeformation, finalDeformation, 0, divWeight, curlWeight, landmarkWeight, imageWeight, consistencyWeight,
				stopThreshold, showProcess ? 2 : 1, showProcess, mode, plugin);
		Thread but = new Thread(bu);
		but.start();
		try {
			but.join();
			but = null;
		} catch (InterruptedException e) {
			System.err.println("Thread interrupted: " + e.getMessage());
			return;
		}

		// Show results
		if (plugin.isPluginInterrumped()) return;
		// bu.getRegisteredSource(srcTgtSeq);
		// Icy.getMainInterface().addSequence(srcTgtSeq);
		try {
			System.out.println("saving to " + transformedSrcResultPath + ", based on " + transformedSrcPath);
			bu.saveBigRegisteredSource(srcResultPath, transformedSrcResultPath, srcPath, transformedSrcPath, tgtPath, null);
		} catch (ServiceException | IOException | FormatException | InterruptedException | ExecutionException e) {
			e.printStackTrace();
			return;
		}

		if (plugin.isPluginInterrumped()) return;

		if (mode != RegistrationModeEnum.MONO.getNumber()) {
			// bu.getRegisteredTarget(tgtTgtSeq);
			// Icy.getMainInterface().addSequence(tgtTgtSeq);
			try {
				bu.saveBigRegisteredTarget(tgtResultPath, transformedTgtResultPath, tgtPath, transformedTgtPath, srcPath, null);
			} catch (ServiceException | IOException | FormatException | InterruptedException | ExecutionException e) {
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

		// // Next Scales Registration
		// for (int si = 0; si < usedScales.length; si++) {
		// double scale = usedScales[si];
		// if (plugin != null && plugin.getUI() != null) {
		// plugin.getUI().setProgressBarMessage("Registering at scale " + scale);
		// plugin.getUI().setProgressBarValue(0);
		// }
		//
		// int tileAmount = (int) Math.ceil(1d / scale);
		// Dimension tileDim = new Dimension(srcDim.width / tileAmount,
		// srcDim.height / tileAmount);
		// Dimension residualTileDim = new Dimension(srcDim.width % tileDim.width,
		// srcDim.height % tileDim.height);
		// Dimension tileSize = new Dimension(tileAmount + (tileDim.width *
		// tileAmount < srcDim.width ? 1 : 0),
		// tileAmount + (tileDim.height * tileAmount < srcDim.height ? 1 : 0));
		// int tileBorderSize = Math.max(tileDim.width, tileDim.height) / 8;
		// // int nProc = Runtime.getRuntime().availableProcessors();
		//
		// ExecutorService threadPool = Executors.newFixedThreadPool(1/* nProc */);
		//
		// int tileNum = 0;
		// Dimension usedTileDim = tileDim;
		// for (int i = 0; i < tileSize.width && !plugin.isPluginInterrumped(); i++)
		// {
		// usedTileDim.width = (i + 1) * tileDim.width <= srcDim.width ?
		// tileDim.width : residualTileDim.width;
		// for (int j = 0; j < tileSize.height && !plugin.isPluginInterrumped();
		// j++) {
		// usedTileDim.height = (j + 1) * tileDim.height <= srcDim.height ?
		// tileDim.height : residualTileDim.height;
		//
		// Rectangle rect = new Rectangle(i * tileDim.width - tileBorderSize, j *
		// tileDim.height - tileBorderSize,
		// usedTileDim.width + tileBorderSize, usedTileDim.height + tileBorderSize);
		//
		// // Register images
		//
		// String sourceResultPath = FilenameUtils.getFullPath(srcResultPath);
		// sourceResultPath += FilenameUtils.getBaseName(srcResultPath) + "_" +
		// tileNum + "_BUnwarp.";
		// sourceResultPath += FilenameUtils.getExtension(srcResultPath);
		//
		// String transformedSourceResultPath =
		// FilenameUtils.getFullPath(transformedSrcResultPath);
		// transformedSourceResultPath +=
		// FilenameUtils.getBaseName(transformedSrcResultPath) + "_" + tileNum
		// + "_BUnwarp.";
		// transformedSourceResultPath +=
		// FilenameUtils.getExtension(transformedSrcResultPath);
		//
		// String sourcePath = srcResultPath;
		// String transformedSourcePath = transformedSrcResultPath;
		// String targetPath = tgtPath;
		//
		// threadPool.submit(new BUnwarpperTask(subsampleFactor, initialDeformation,
		// finalDeformation, 0, divWeight, curlWeight, landmarkWeight,
		// imageWeight, consistencyWeight, stopThreshold, showProcess ? 2 : 1,
		// showProcess,
		// RegistrationModeEnum.MONO.getNumber(), plugin, sourceResultPath,
		// transformedSourceResultPath, sourcePath,
		// transformedSourcePath, transformedSrcResultPath, targetPath, rect));
		// tileNum++;
		// }
		// }
		// threadPool.shutdown();
		// try {
		// threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.HOURS);
		// } catch (InterruptedException e) {
		// e.printStackTrace();
		// }
		// // TODO inverse registration
		// }

	}

	@SuppressWarnings("unused")
	private static class BUnwarpperTask implements Runnable {

		final int			maxImageSubsamplingFactor;
		final int			minScaleDeformation;
		final int			maxScaleDeformation;
		final int			minScaleImage;
		final double	divWeight;
		final double	curlWeight;
		final double	landmarkWeight;
		final double	imageWeight;
		final double	consistencyWeight;
		final double	stopThreshold;
		final int			outputLevel;
		final boolean	showMarquardtOptim;
		final int			accurateMode;
		final BUnwarp	plugin;

		BUnwarpper	unwarp;
		String			sourceResultPath;
		String			transformedSourceResultPath;
		String			sourcePath;
		String			transformedSourcePath;
		String			targetPath;
		String			srcPath;
		Rectangle		tile;

		public BUnwarpperTask(final int maxImageSubsamplingFactor, final int minScaleDeformation,
				final int maxScaleDeformation, final int minScaleImage, final double divWeight, final double curlWeight,
				final double landmarkWeight, final double imageWeight, final double consistencyWeight,
				final double stopThreshold, final int outputLevel, final boolean showMarquardtOptim, final int accurateMode,
				final BUnwarp plugin, String sourceResultPath, String transformedSourceResultPath, String sourcePath,
				String transformedSourcePath, String srcPath, String targetPath, Rectangle tile) {

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

			this.sourceResultPath = sourceResultPath;
			this.transformedSourceResultPath = transformedSourceResultPath;
			this.sourcePath = sourcePath;
			this.transformedSourcePath = transformedSourcePath;
			this.targetPath = targetPath;
			this.srcPath = srcPath;
			this.tile = tile;

		}

		@Override
		public void run() {

			BigImageLoader loader = new BigImageLoader();
			Sequence srcSeq = null;
			ROI2D srcMask = null;
			Sequence tgtSeq = null;
			ROI2D tgtMask = null;
			try {
				ProgressBar.setProgressBarMessage("Loading source tile image");
				srcSeq = loader.loadDownsampledImage(sourcePath, tile, 1023, 1023, true);
				ProgressBar.setProgressBarMessage("Loading source mask");
				srcMask = loader.loadDownsampledMask(srcSeq, srcPath, tile, 1023, 1023, true);
				ProgressBar.setProgressBarMessage("Loading target tile image");
				tgtSeq = loader.loadDownsampledImage(targetPath, tile, 1023, 1023, true);
				ProgressBar.setProgressBarMessage("Loading source mask");
				tgtMask = loader.loadDownsampledMask(tgtSeq, targetPath, tile, 1000, 1000, true);
			} catch (UnsupportedFormatException | IOException e1) {
				e1.printStackTrace();
			}

			if (srcSeq != null) {
				if (srcMask != null) srcSeq.addROI(srcMask);
				Icy.getMainInterface().addSequence(srcSeq);
			}
			if (tgtSeq != null) {
				if (tgtMask != null) tgtSeq.addROI(tgtMask);
				Icy.getMainInterface().addSequence(tgtSeq);
			}
			// srcTgtSeq = SequenceUtil.getCopy(srcSeq);
			// tgtTgtSeq = SequenceUtil.getCopy(tgtSeq);

			unwarp = new BUnwarpper(srcSeq, tgtSeq, new ArrayList<ROI2DPoint>(), new ArrayList<ROI2DPoint>(), srcMask,
					tgtMask, maxImageSubsamplingFactor, minScaleDeformation, maxScaleDeformation, minScaleImage, divWeight,
					curlWeight, landmarkWeight, imageWeight, consistencyWeight, stopThreshold, outputLevel, showMarquardtOptim,
					accurateMode, plugin);

			Thread thr = new Thread(unwarp);
			thr.start();
			try {
				thr.join();
				thr = null;
				System.out.println("saving " + tile);
				unwarp.saveBigRegisteredSource(sourceResultPath, transformedSourceResultPath, sourcePath, transformedSourcePath,
						targetPath, tile);
			} catch (InterruptedException | ServiceException | IOException | FormatException |ExecutionException e) {
				e.printStackTrace();
			}

		}

	}

}
