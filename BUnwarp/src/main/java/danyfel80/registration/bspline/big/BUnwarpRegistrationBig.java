package danyfel80.registration.bspline.big;

import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import algorithms.danyfel80.io.sequence.large.LargeSequenceExporter;
import algorithms.danyfel80.io.sequence.large.LargeSequenceHelper;
import algorithms.danyfel80.io.sequence.large.LargeSequenceImporter;
import danyfel80.registration.bspline.classic.BUnwarpRegistration;
import danyfel80.registration.bspline.classic.DeformationScale;
import danyfel80.registration.bspline.classic.RegistrationMode;
import danyfel80.registration.bspline.classic.Transformation;
import icy.common.exception.UnsupportedFormatException;
import icy.common.listener.DetailedProgressListener;
import icy.file.FileUtil;
import icy.roi.ROI;
import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;
import icy.util.XMLUtil;
import loci.formats.ome.OMEXMLMetadata;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DPoint;

public class BUnwarpRegistrationBig {
	private static final Dimension REGISTERED_MAX_DIMENSION = new Dimension(1500, 1000);
	/* Internal Variables */
	private File sourceFile;
	private File targetFile;
	private List<ROI2DPoint> sourceLandmarkROIs;
	private List<ROI2DPoint> targetLandmarkROIs;
	private RegistrationMode registrationMode;
	private int initialSubsampleFactor;
	private File transformedSourceFile;
	private File transformedTargetFile;
	private File transformedSourceOutputFile;
	private File transformedTargetOutputFile;
	private DeformationScale initialDeformationScale;
	private DeformationScale finalDeformationScale;
	private double divWeight;
	private double curlWeight;
	private double landmarkWeight;
	private double imageWeight;
	private double consistencyWeight;
	private double stopThreshold;
	private boolean showProcess;

	private Set<DetailedProgressListener> progressListeners;
	private Set<DetailedProgressListener> progressOutputListeners;

	private int subsampledSourceResolution;
	private int subsampledTargetResolution;
	private Sequence subsampledSourceSequence;
	private Sequence subsampledTargetSequence;
	private ROI2DArea sourceMask;
	private ROI2DArea targetMask;

	private BUnwarpRegistration registration;

	public BUnwarpRegistrationBig() {
		progressListeners = new HashSet<>();
		progressOutputListeners = new HashSet<>();
	}

	/**
	 * @return The sourceFile.
	 */
	public File getSourceFile() {
		return sourceFile;
	}

	/**
	 * @param sourceFile
	 *          The sourceFile to set.
	 */
	public void setSourceFile(File sourceFile) {
		this.sourceFile = sourceFile;
	}

	/**
	 * @return The targetFile.
	 */
	public File getTargetFile() {
		return targetFile;
	}

	/**
	 * @param targetFile
	 *          The targetFile to set.
	 */
	public void setTargetFile(File targetFile) {
		this.targetFile = targetFile;
	}

	public void setSourceLandmarkROIs(List<ROI2DPoint> sourceLandmarkROIs) {
		this.sourceLandmarkROIs = sourceLandmarkROIs;
	}

	public void setTargetLandmarkROIs(List<ROI2DPoint> targetLandmarkROIs) {
		this.targetLandmarkROIs = targetLandmarkROIs;
	}

	/**
	 * @return The registrationMode.
	 */
	public RegistrationMode getRegistrationMode() {
		return registrationMode;
	}

	/**
	 * @param registrationMode
	 *          The registrationMode to set.
	 */
	public void setRegistrationMode(RegistrationMode registrationMode) {
		this.registrationMode = registrationMode;
	}

	/**
	 * @return The initialSubsampleFactor.
	 */
	public int getInitialSubsampleFactor() {
		return initialSubsampleFactor;
	}

	/**
	 * @param initialSubsampleFactor
	 *          The initialSubsampleFactor to set.
	 */
	public void setInitialSubsampleFactor(int initialSubsampleFactor) {
		this.initialSubsampleFactor = initialSubsampleFactor;
	}

	/**
	 * @return The transformedSourceFile.
	 */
	public File getTransformedSourceFile() {
		return transformedSourceFile;
	}

	/**
	 * @param transformedSourceFile
	 *          The transformedSourceFile to set.
	 */
	public void setTransformedSourceFile(File transformedSourceFile) {
		this.transformedSourceFile = transformedSourceFile;
	}

	/**
	 * @return The transformedTargetFile.
	 */
	public File getTransformedTargetFile() {
		return transformedTargetFile;
	}

	/**
	 * @param transformedTargetFile
	 *          The transformedTargetFile to set.
	 */
	public void setTransformedTargetFile(File transformedTargetFile) {
		this.transformedTargetFile = transformedTargetFile;
	}

	/**
	 * @return The transformedSourceOutputFile.
	 */
	public File getTransformedSourceOutputFile() {
		return transformedSourceOutputFile;
	}

	/**
	 * @param transformedSourceOutputFile
	 *          The transformedSourceOutputFile to set.
	 */
	public void setTransformedSourceOutputFile(File transformedSourceOutputFile) {
		this.transformedSourceOutputFile = transformedSourceOutputFile;
	}

	/**
	 * @return The transformedTargetOutputFile.
	 */
	public File getTransformedTargetOutputFile() {
		return transformedTargetOutputFile;
	}

	/**
	 * @param transformedTargetOutputFile
	 *          The transformedTargetOutputFile to set.
	 */
	public void setTransformedTargetOutputFile(File transformedTargetOutputFile) {
		this.transformedTargetOutputFile = transformedTargetOutputFile;
	}

	/**
	 * @return The initialDeformationScale.
	 */
	public DeformationScale getInitialDeformationScale() {
		return initialDeformationScale;
	}

	/**
	 * @param initialDeformationScale
	 *          The initialDeformationScale to set.
	 */
	public void setInitialDeformationScale(DeformationScale initialDeformationScale) {
		this.initialDeformationScale = initialDeformationScale;
	}

	/**
	 * @return The finalDeformationScale.
	 */
	public DeformationScale getFinalDeformationScale() {
		return finalDeformationScale;
	}

	/**
	 * @param finalDeformationScale
	 *          The finalDeformationScale to set.
	 */
	public void setFinalDeformationScale(DeformationScale finalDeformationScale) {
		this.finalDeformationScale = finalDeformationScale;
	}

	/**
	 * @return The divWeight.
	 */
	public double getDivWeight() {
		return divWeight;
	}

	/**
	 * @param divWeight
	 *          The divWeight to set.
	 */
	public void setDivWeight(double divWeight) {
		this.divWeight = divWeight;
	}

	/**
	 * @return The curlWeight.
	 */
	public double getCurlWeight() {
		return curlWeight;
	}

	/**
	 * @param curlWeight
	 *          The curlWeight to set.
	 */
	public void setCurlWeight(double curlWeight) {
		this.curlWeight = curlWeight;
	}

	/**
	 * @return The landmarkWeight.
	 */
	public double getLandmarkWeight() {
		return landmarkWeight;
	}

	/**
	 * @param landmarkWeight
	 *          The landmarkWeight to set.
	 */
	public void setLandmarkWeight(double landmarkWeight) {
		this.landmarkWeight = landmarkWeight;
	}

	/**
	 * @return The imageWeight.
	 */
	public double getImageWeight() {
		return imageWeight;
	}

	/**
	 * @param imageWeight
	 *          The imageWeight to set.
	 */
	public void setImageWeight(double imageWeight) {
		this.imageWeight = imageWeight;
	}

	/**
	 * @return The consistencyWeight.
	 */
	public double getConsistencyWeight() {
		return consistencyWeight;
	}

	/**
	 * @param consistencyWeight
	 *          The consistencyWeight to set.
	 */
	public void setConsistencyWeight(double consistencyWeight) {
		this.consistencyWeight = consistencyWeight;
	}

	/**
	 * @return The stopThreshold.
	 */
	public double getStopThreshold() {
		return stopThreshold;
	}

	/**
	 * @param stopThreshold
	 *          The stopThreshold to set.
	 */
	public void setStopThreshold(double stopThreshold) {
		this.stopThreshold = stopThreshold;
	}

	/**
	 * @return The showProcess.
	 */
	public boolean isShowProcess() {
		return showProcess;
	}

	/**
	 * @param showProcess
	 *          The showProcess to set.
	 */
	public void setShowProcess(boolean showProcess) {
		this.showProcess = showProcess;
	}

	/**
	 * Adds a listener to notify when a progress on the registration is made.
	 * 
	 * @param listener
	 *          Listener to register.
	 */
	public void addProgressListener(DetailedProgressListener listener) {
		this.progressListeners.add(listener);
	}

	/**
	 * Adds a listener to notify when an output progress image is available.
	 * 
	 * @param listener
	 *          Listener to register.
	 */
	public void addProgressOutputListener(DetailedProgressListener listener) {
		this.progressOutputListeners.add(listener);
	}

	private void notifyProgress(double progress, String message) {
		for (DetailedProgressListener l: progressListeners) {
			l.notifyProgress(progress, message, null);
		}
	}

	private void notifyProgressOutput(Sequence sequence) {
		for (DetailedProgressListener l: progressOutputListeners) {
			l.notifyProgress(Double.NaN, "", sequence);
		}
	}

	public void compute() throws BUnwarpRegistrationBigException, InterruptedException {
		long startTime = System.currentTimeMillis();
		long endTime;
		try {
			importSubsampledSequences();
		} catch (Exception e) {
			throw new BUnwarpRegistrationBigException("Could not import input images.", e);
		}
		endTime = System.currentTimeMillis();
		System.out.println("Images loaded: " + (endTime - startTime) + "ms");

		if (Thread.interrupted())
			throw new InterruptedException();

		startTime = System.currentTimeMillis();
		try {
			performRegistration();
		} catch (Exception e) {
			throw new BUnwarpRegistrationBigException("Could not perform the registration.", e);
		}
		endTime = System.currentTimeMillis();
		System.out.println("Transformation computed: " + (endTime - startTime) + "ms");

		if (Thread.interrupted())
			throw new InterruptedException();

		startTime = System.currentTimeMillis();
		try {
			applyTransformationOnTransformedImages();
		} catch (InterruptedException e) {
			throw new BUnwarpRegistrationBigException("Could not apply transformation to large images: Interrupted", e);
		} catch (Exception e) {
			throw new BUnwarpRegistrationBigException("Could not apply transformation to large images:\n" + e.getMessage(),
					e);
		}
		endTime = System.currentTimeMillis();
		System.out.println("Registered images saved: " + (endTime - startTime) + "ms");
	}

	private void importSubsampledSequences() throws Exception {
		importSubsampledSourceSequence();
		importSubsampledTargetSequence();
	}

	private void importSubsampledSourceSequence() throws Exception {
		subsampledSourceResolution = LargeSequenceHelper.getResolutionLevel(getSourceFile(), REGISTERED_MAX_DIMENSION);
		subsampledSourceSequence = importSubsampledSequence(getSourceFile(), subsampledSourceResolution, "source");
		List<ROI> rois = subsampledSourceSequence.getROIs();
		subsampledSourceSequence = SequenceUtil.convertToType(subsampledSourceSequence, DataType.UBYTE, true, true);
		subsampledSourceSequence.addROIs(rois, false);
		notifyProgressOutput(subsampledSourceSequence);
	}

	private void importSubsampledTargetSequence() throws Exception {
		subsampledTargetResolution = LargeSequenceHelper.getResolutionLevel(getTargetFile(), REGISTERED_MAX_DIMENSION);
		subsampledTargetSequence = importSubsampledSequence(getTargetFile(), subsampledTargetResolution, "target");
		List<ROI> rois = subsampledTargetSequence.getROIs();
		subsampledTargetSequence = SequenceUtil.convertToType(subsampledTargetSequence, DataType.UBYTE, true, true);
		subsampledTargetSequence.addROIs(rois, false);
		notifyProgressOutput(subsampledTargetSequence);
	}

	private Sequence importSubsampledSequence(File file, int resolution, String label) throws Exception {
		LargeSequenceImporter importer = new LargeSequenceImporter();
		importer.setFilePath(file.toPath());
		importer.setTargetResolution(subsampledSourceResolution);
		importer.addProgressListener((progress, message, data) -> {
			notifyProgress(progress, String.format("Loading " + label + " image (%s)", message));
			return true;
		});
		return importer.call();
	}

	private void performRegistration() throws BUnwarpRegistrationBigException, InterruptedException {
		registration = new BUnwarpRegistration();
		setRegistrationParameters();
		for (DetailedProgressListener listener: progressListeners) {
			registration.addProgressListener(listener);
		}
		for (DetailedProgressListener listener: progressOutputListeners) {
			registration.addProgressOutputListener(listener);
		}
		registration.compute();
	}

	private void setRegistrationParameters() throws BUnwarpRegistrationBigException {
		this.registration.setSourceSequence(subsampledSourceSequence);
		this.registration.setTargetSequence(subsampledTargetSequence);
		this.registration.setRegistrationMode(getRegistrationMode());
		this.registration.setInitialSubsampleFactor(getInitialSubsampleFactor());
		this.registration.setTransformedSourceSequence(subsampledSourceSequence);
		this.registration.setTransformedTargetSequence(subsampledTargetSequence);
		this.registration.setInitialDeformationScale(getInitialDeformationScale());
		this.registration.setFinalDeformationScale(getFinalDeformationScale());
		this.registration.setDivWeight(getDivWeight());
		this.registration.setCurlWeight(getCurlWeight());
		this.registration.setLandmarkWeight(getLandmarkWeight());
		this.registration.setImageWeight(getImageWeight());
		this.registration.setConsistencyWeight(getConsistencyWeight());
		this.registration.setStopThreshold(getStopThreshold());
		this.registration.setShowProcess(isShowProcess());
		if (sourceLandmarkROIs == null || sourceLandmarkROIs.isEmpty() || targetLandmarkROIs == null
				|| targetLandmarkROIs.isEmpty()) {
			extractLandmarks();
		}
		this.registration.setSourceLandmarks(convertLandmarksToPoints(sourceLandmarkROIs));
		this.registration.setTargetLandmarks(convertLandmarksToPoints(targetLandmarkROIs));
		extractMasks();
		this.registration.setSourceMask(sourceMask);
		this.registration.setTargetMask(targetMask);
	}

	private void extractLandmarks() throws BUnwarpRegistrationBigException {
		sourceLandmarkROIs = subsampledSourceSequence.getROIs(ROI2DPoint.class, false).stream()
				.sorted(Comparator.comparing(ROI2DPoint::getName)).collect(Collectors.toList());
		targetLandmarkROIs = subsampledTargetSequence.getROIs(ROI2DPoint.class, false).stream()
				.sorted(Comparator.comparing(ROI2DPoint::getName)).collect(Collectors.toList());

		checkLandmarksAreCorrectInSizeAndName();
	}

	private void checkLandmarksAreCorrectInSizeAndName() throws BUnwarpRegistrationBigException {
		if (sourceLandmarkROIs.size() != targetLandmarkROIs.size())
			throw new BUnwarpRegistrationBigException(
					String.format("Source landmarks(%d) and target landmarks(%d) have different size.", sourceLandmarkROIs.size(),
							targetLandmarkROIs.size()));

		for (Iterator<ROI2DPoint> itSource = sourceLandmarkROIs.iterator(), itTarget = targetLandmarkROIs
				.iterator(); itSource.hasNext();) {
			ROI2DPoint sourceLandmark = itSource.next();
			ROI2DPoint targetLandmark = itTarget.next();
			if (!sourceLandmark.getName().equals(targetLandmark.getName()))
				throw new BUnwarpRegistrationBigException("No corresponding landmark for " + sourceLandmark.getName());
		}
	}

	private static List<Point2D> convertLandmarksToPoints(List<ROI2DPoint> landmarkROIs) {
		return landmarkROIs.stream().map(roi -> roi.getPoint()).collect(Collectors.toList());
	}

	private void extractMasks() {
		sourceMask = new ROI2DArea();
		sourceMask.addRect(0, 0, registration.getSourceSequence().getWidth(), registration.getSourceSequence().getHeight());
		targetMask = new ROI2DArea();
		targetMask.addRect(0, 0, registration.getTargetSequence().getWidth(), registration.getTargetSequence().getHeight());
	}

	private void applyTransformationOnTransformedImages() throws InterruptedException {
		ThreadPoolExecutor transformExportThreads = (ThreadPoolExecutor) Executors.newFixedThreadPool(1);
		Future<Void> transformationOnSourceFuture = transformExportThreads.submit(getTransformationOnSourceImageTask());
		transformExportThreads.shutdown();
		try {
			try {
				transformationOnSourceFuture.get();
			} catch (ExecutionException e) {
				throw new BUnwarpRegistrationBigException("Could not save transformed source image:\n" + e.getMessage(), e);
			}
		} finally {
			transformExportThreads.shutdownNow();
			if (!transformExportThreads.awaitTermination(10, TimeUnit.SECONDS)) {
				throw new BUnwarpRegistrationBigException(
						"Transformation writer thread still did not finished writing source image.");
			}
		}

		if (getRegistrationMode() != RegistrationMode.MONO) {
			transformExportThreads = (ThreadPoolExecutor) Executors.newFixedThreadPool(1);
			Future<Void> transformationOnTargetFuture = transformExportThreads.submit(getTransformationOnTargetImageTask());
			transformExportThreads.shutdown();
			try {
				try {
					transformationOnTargetFuture.get();
				} catch (ExecutionException e) {
					throw new BUnwarpRegistrationBigException("Could not save transformed target image:\n" + e.getMessage(), e);
				}

			} finally {
				transformExportThreads.shutdownNow();
				if (!transformExportThreads.awaitTermination(10, TimeUnit.SECONDS)) {
					throw new BUnwarpRegistrationBigException(
							"Transformation writer thread still did not finished writing target image.");
				}
			}
		}

	}

	private Callable<Void> getTransformationOnSourceImageTask() {
		return () -> {
			applyTransformationOnTransformedSourceImage();
			return null;
		};
	}

	private Callable<Void> getTransformationOnTargetImageTask() {
		return () -> {
			applyTransformationOnTransformedTargetImage();
			return null;
		};
	}

	private void applyTransformationOnTransformedSourceImage() throws InterruptedException {
		try (LargeSequenceExporter sourceExporter = new LargeSequenceExporter()) {
			sourceExporter.setOutputFilePath(getTransformedSourceOutputFile().toPath());
			OMEXMLMetadata outputImageMetadata = getTransformedSourceImageMetadata();
			sourceExporter.setOutputImageMetadata(outputImageMetadata);
			// TileProvider taking into account sourceToTargetTransformation and transformedSourceFile
			TransformationSourceTileProvider tileProvider = getTrasfromedSourceTileProvider(outputImageMetadata,
					sourceExporter.TILE_SIZE);
			sourceExporter.setTileProvider(tileProvider);
			sourceExporter.addProgressListener((progress, message, data) -> {
				notifyProgress(progress, String.format("Creating transformed source image (%s)", message));
				return true;
			});
			// start exporting
			try {
				sourceExporter.write();
			} finally {
				tileProvider.close();
			}
			applyTransformationOnSourceROIs();

		} catch (InterruptedException e) {
			throw e;
		} catch (Exception e) {
			throw new BUnwarpRegistrationBigException("Could not save transformed source image", e);
		}
	}

	private OMEXMLMetadata getTransformedSourceImageMetadata() throws IOException {
		String name;
		Dimension targetSize;
		int sourceSizeC;
		DataType sourceDataType;
		try {
			name = LargeSequenceHelper.getImageName(transformedSourceFile);
			targetSize = LargeSequenceHelper.getImageDimension(transformedTargetFile);
			sourceSizeC = LargeSequenceHelper.getImageChannelCount(transformedSourceFile);
			sourceDataType = LargeSequenceHelper.getImageDataType(transformedSourceFile);
		} catch (Exception e) {
			throw new IOException("Could not create transformed source image metadata.", e);
		}
		OMEXMLMetadata metaData = LargeSequenceExporter.createMetadata(targetSize.width, targetSize.height, sourceSizeC,
				sourceDataType);
		MetaDataUtil.setName(metaData, 0, name);
		return metaData;
	}

	private TransformationSourceTileProvider getTrasfromedSourceTileProvider(OMEXMLMetadata outputImageMetadata,
			Dimension tileSize) {
		TransformationSourceTileProvider sourceTileProvider = new TransformationSourceTileProvider();
		sourceTileProvider.setTransformedSourceFilePath(getTransformedSourceFile().toPath());
		sourceTileProvider.setOutputImageMetatada(outputImageMetadata);
		sourceTileProvider.setOutputTileSize(tileSize);
		sourceTileProvider.setTransformation(registration.getTransformation());
		sourceTileProvider.setSubsampledSourceSize(subsampledSourceSequence.getDimension2D());
		sourceTileProvider.setSubsampledTargetSize(subsampledTargetSequence.getDimension2D());
		return sourceTileProvider;
	}

	private void applyTransformationOnSourceROIs() throws InterruptedException, IOException {
		OMEXMLMetadata outputSourceImageMetadata;
		OMEXMLMetadata outputTargetImageMetadata;
		OMEXMLMetadata outputSourceTImageMetadata;
		try {
			outputSourceImageMetadata = (OMEXMLMetadata) LargeSequenceHelper.getImageMetadata(transformedSourceFile);
			outputTargetImageMetadata = (OMEXMLMetadata) LargeSequenceHelper.getImageMetadata(transformedTargetFile);
			outputSourceTImageMetadata = (OMEXMLMetadata) LargeSequenceHelper.getImageMetadata(transformedSourceOutputFile);
		} catch (UnsupportedFormatException | IOException e) {
			throw new IOException("Could not read source output image metadata", e);
		}
		Sequence tempTOutputS = new Sequence(outputSourceTImageMetadata,
				MetaDataUtil.getName(outputSourceTImageMetadata, 0));
		tempTOutputS.setFilename(transformedSourceOutputFile.getPath());

		List<ROI2DPoint> sourceRois = retrieveRoisFromFileXML(transformedSourceFile.toPath());
		Point2D outputSourceSize = new Point2D.Double(MetaDataUtil.getSizeX(outputSourceImageMetadata, 0),
				MetaDataUtil.getSizeY(outputSourceImageMetadata, 0));
		Point2D subsampledSourceSize = new Point2D.Double(subsampledSourceSequence.getWidth(),
				subsampledSourceSequence.getHeight());
		Point2D outputTargetSize = new Point2D.Double(MetaDataUtil.getSizeX(outputTargetImageMetadata, 0),
				MetaDataUtil.getSizeY(outputTargetImageMetadata, 0));
		Point2D subsampledTargetSize = new Point2D.Double(subsampledTargetSequence.getWidth(),
				subsampledTargetSequence.getHeight());
		sourceRois = transformSourceRois(sourceRois, registration.getTransformation(), outputSourceSize,
				subsampledSourceSize, outputTargetSize, subsampledTargetSize);
		tempTOutputS.addROIs(sourceRois, false);
		tempTOutputS.saveXMLData();
	}

	private List<ROI2DPoint> retrieveRoisFromFileXML(Path file) {
		String fileName = FileUtil.getFileName(file.toString(), false);
		Path xmlFile = file.resolveSibling(fileName + ".xml");
		if (Files.exists(xmlFile)) {
			Document xml = XMLUtil.loadDocument(xmlFile.toFile());
			Element rootElement = XMLUtil.getRootElement(xml);
			Element roisElement = XMLUtil.getElement(rootElement, "rois");
			List<ROI2DPoint> rois = ROI.loadROIsFromXML(roisElement).stream().filter(r -> r instanceof ROI2DPoint)
					.map(r -> (ROI2DPoint) r).collect(Collectors.toList());
			return rois;
		}
		return Collections.emptyList();
	}

	private List<ROI2DPoint> transformSourceRois(List<ROI2DPoint> sourceRois, Transformation transformation,
			Point2D outputSourceSize, Point2D subsampledSourceSize, Point2D outputTargetSize, Point2D subsampledTargetSize) {
		return sourceRois.stream().map(r -> {
			Point2D pos = r.getPosition2D();
			double[] tPos = new double[] {pos.getX() * subsampledTargetSize.getX() / outputTargetSize.getX(),
					pos.getY() * subsampledTargetSize.getY() / outputTargetSize.getY()};
			transformation.transform(tPos[0], tPos[1], tPos, true);
			ROI2DPoint nr = new ROI2DPoint((tPos[0] * outputSourceSize.getX()) / subsampledSourceSize.getX(),
					(tPos[1] * outputSourceSize.getY()) / subsampledSourceSize.getY());
			nr.setName(r.getName());
			return nr;
		}).collect(Collectors.toList());
	}

	private void applyTransformationOnTransformedTargetImage() throws InterruptedException {
		try (LargeSequenceExporter targetExporter = new LargeSequenceExporter()) {
			targetExporter.setOutputFilePath(getTransformedTargetOutputFile().toPath());
			OMEXMLMetadata outputImageMetadata = getTransformedTargetImageMetadata();
			targetExporter.setOutputImageMetadata(outputImageMetadata);
			// TileProvider taking into account sourceToTargetTransformation and transformedSourceFile
			TransformationTargetTileProvider tileProvider = getTrasfromedTargetTileProvider(outputImageMetadata,
					targetExporter.TILE_SIZE);
			targetExporter.setTileProvider(tileProvider);
			targetExporter.addProgressListener((progress, message, data) -> {
				notifyProgress(progress, String.format("Creating transformed target image (%s)", message));
				return true;
			});
			// start exporting
			try {
				targetExporter.write();
			} finally {
				tileProvider.close();
			}
			applyTransformationOnTargetROIs();

		} catch (InterruptedException e) {
			throw e;
		} catch (Exception e) {
			throw new BUnwarpRegistrationBigException("Could not save transformed target image", e);
		}
	}

	private OMEXMLMetadata getTransformedTargetImageMetadata() throws IOException {
		String name;
		Dimension sourceSize;
		int targetSizeC;
		DataType targetDataType;
		try {
			name = LargeSequenceHelper.getImageName(transformedTargetFile);
			sourceSize = LargeSequenceHelper.getImageDimension(transformedSourceFile);
			targetSizeC = LargeSequenceHelper.getImageChannelCount(transformedTargetFile);
			targetDataType = LargeSequenceHelper.getImageDataType(transformedTargetFile);
		} catch (Exception e) {
			throw new IOException("Could not create transformed target image metadata.", e);
		}
		OMEXMLMetadata metaData = LargeSequenceExporter.createMetadata(sourceSize.width, sourceSize.height, targetSizeC,
				targetDataType);
		MetaDataUtil.setName(metaData, 0, name);
		return metaData;
	}

	private TransformationTargetTileProvider getTrasfromedTargetTileProvider(OMEXMLMetadata outputImageMetadata,
			Dimension tileSize) {
		TransformationTargetTileProvider targetTileProvider = new TransformationTargetTileProvider();
		targetTileProvider.setTransformedSourceFilePath(getTransformedTargetFile().toPath());
		targetTileProvider.setOutputImageMetatada(outputImageMetadata);
		targetTileProvider.setOutputTileSize(tileSize);
		targetTileProvider.setTransformation(registration.getTransformation());
		targetTileProvider.setSubsampledSourceSize(subsampledTargetSequence.getDimension2D());
		targetTileProvider.setSubsampledTargetSize(subsampledSourceSequence.getDimension2D());
		return targetTileProvider;
	}

	private void applyTransformationOnTargetROIs() throws IOException {
		OMEXMLMetadata outputTargetImageMetadata;
		OMEXMLMetadata outputSourceImageMetadata;
		OMEXMLMetadata outputTargetTImageMetadata;
		try {
			outputTargetImageMetadata = (OMEXMLMetadata) LargeSequenceHelper.getImageMetadata(transformedTargetFile);
			outputSourceImageMetadata = (OMEXMLMetadata) LargeSequenceHelper.getImageMetadata(transformedSourceFile);
			outputTargetTImageMetadata = (OMEXMLMetadata) LargeSequenceHelper.getImageMetadata(transformedTargetOutputFile);
		} catch (UnsupportedFormatException | IOException e) {
			throw new IOException("Could not read target output image metadata", e);
		}
		Sequence tempTOutputT = new Sequence(outputTargetTImageMetadata,
				MetaDataUtil.getName(outputTargetTImageMetadata, 0));
		tempTOutputT.setFilename(transformedTargetOutputFile.getPath());

		List<ROI2DPoint> targetRois = retrieveRoisFromFileXML(transformedTargetFile.toPath());
		Point2D outputTargetSize = new Point2D.Double(MetaDataUtil.getSizeX(outputTargetImageMetadata, 0),
				MetaDataUtil.getSizeY(outputTargetImageMetadata, 0));
		Point2D subsampledTargetSize = new Point2D.Double(subsampledTargetSequence.getWidth(),
				subsampledTargetSequence.getHeight());
		Point2D outputSourceSize = new Point2D.Double(MetaDataUtil.getSizeX(outputSourceImageMetadata, 0),
				MetaDataUtil.getSizeY(outputSourceImageMetadata, 0));
		Point2D subsampledSourceSize = new Point2D.Double(subsampledSourceSequence.getWidth(),
				subsampledSourceSequence.getHeight());
		targetRois = transformTargetRois(targetRois, registration.getTransformation(), outputTargetSize,
				subsampledTargetSize, outputSourceSize, subsampledSourceSize);
		tempTOutputT.addROIs(targetRois, false);
		tempTOutputT.saveXMLData();
	}

	private List<ROI2DPoint> transformTargetRois(List<ROI2DPoint> targetRois, Transformation transformation,
			Point2D outputTargetSize, Point2D subsampledTargetSize, Point2D outputSourceSize, Point2D subsampledSourceSize) {
		return targetRois.stream().map(r -> {
			Point2D pos = r.getPosition2D();
			double[] tPos = new double[] {(pos.getX() * subsampledSourceSize.getX()) / outputSourceSize.getX(),
					(pos.getY() * subsampledSourceSize.getY()) / outputSourceSize.getY()};
			transformation.transform(tPos[0], tPos[1], tPos, false);
			return new ROI2DPoint((tPos[0] * outputTargetSize.getX()) / subsampledTargetSize.getX(),
					(tPos[1] * outputTargetSize.getY()) / subsampledTargetSize.getY());
		}).collect(Collectors.toList());
	}
}
