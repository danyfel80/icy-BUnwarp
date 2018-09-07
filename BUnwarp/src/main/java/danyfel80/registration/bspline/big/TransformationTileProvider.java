package danyfel80.registration.bspline.big;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.IOException;
import java.nio.file.Path;

import algorithms.danyfel80.io.sequence.cursor.IcyBufferedImageCursor;
import algorithms.danyfel80.io.sequence.tileprovider.CachedLargeSequenceTileProvider;
import algorithms.danyfel80.io.sequence.tileprovider.ITileProvider;
import danyfel80.registration.bspline.classic.BSplineModel;
import danyfel80.registration.bspline.classic.Transformation;
import icy.common.exception.UnsupportedFormatException;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.type.DataType;
import loci.formats.ome.OMEXMLMetadata;
import plugins.kernel.importer.LociImporterPlugin;

public abstract class TransformationTileProvider implements ITileProvider {

	private Path transformedSourceFilePath;
	private OMEXMLMetadata outputImageMetadata;
	private Dimension outputTileSize;
	private Transformation transformation;
	private Dimension subsampledSourceSize;
	private Dimension subsampledTargetSize;

	private boolean providerPrepared;
	private LociImporterPlugin transformedSourceSequenceImporter;
	private Dimension transformedSourceFullSize;
	private Rectangle transformedSourceFullRectangle;
	private int transformedSourceSizeC;
	private DataType transformedSourceDataType;
	private Dimension transformedSourceTileSize;
	private Rectangle transformedSourceFullTileGrid;
	private CachedLargeSequenceTileProvider transformedSourceTileProvider;
	private Dimension outputFullSize;
	private int outputTileSizeC;
	private DataType outputTileDataType;

	private Point tileIndex;
	private Point currentTilePosition;
	private Rectangle currentTileRectangle;
	private BSplineModel coefficientsX;
	private BSplineModel coefficientsY;
	private double sourceXScaleFactor;
	private double sourceYScaleFactor;
	private Rectangle neededSourceRectangle;
	private Sequence neededSourceSequence;
	private IcyBufferedImage outputTileImage;
	private double coefficientScaleFactorY;
	private double coefficientScaleFactorX;

	public Path getTransformedSourceFilePath() {
		return transformedSourceFilePath;
	}

	public void setTransformedSourceFilePath(Path filePath) {
		this.transformedSourceFilePath = filePath;
	}

	public OMEXMLMetadata getOutputImageMetadata() {
		return outputImageMetadata;
	}

	public void setOutputImageMetatada(OMEXMLMetadata metadata) {
		this.outputImageMetadata = metadata;
	}

	public Dimension getOutputTileSize() {
		return outputTileSize;
	}

	public void setOutputTileSize(Dimension tileSize) {
		this.outputTileSize = tileSize;
	}

	public Transformation getTransformation() {
		return transformation;
	}

	public void setTransformation(Transformation transformation) {
		this.transformation = transformation;
	}

	public Dimension getSubsampledSourceSize() {
		return subsampledSourceSize;
	}

	public void setSubsampledSourceSize(Dimension subsampledSourceSize) {
		this.subsampledSourceSize = subsampledSourceSize;
	}

	public Dimension getSubsampledTargetSize() {
		return subsampledTargetSize;
	}

	public void setSubsampledTargetSize(Dimension subsampledTargetSize) {
		this.subsampledTargetSize = subsampledTargetSize;
	}

	public IcyBufferedImage getOutputTileImage() {
		return outputTileImage;
	}

	protected void setOutputTileImage(IcyBufferedImage outputImageTile) {
		this.outputTileImage = outputImageTile;
	}

	@Override
	public IcyBufferedImage getTile(Point tileIndex) throws IOException {
		prepareProvider();
		setCurrentTileIndex(tileIndex);
		try {
			loadNeededInputForCurrentOutputTile();
		} catch (Exception e) {
			throw new IOException(String.format("Could not load needed input image for tile %s", tileIndex), e);
		}
		try {
			writeOutputTileImage();
		} catch (InterruptedException e) {
			throw new IOException(String.format("Interrupted while writing tile (%s) to output.", tileIndex), e);
		}

		return getOutputTileImage();
	}

	private void prepareProvider() throws IOException {
		if (!isProviderPrepared()) {
			setupTransformedSourceSequenceImporter();
			setupCacheTileProvider();
			setupOutputInformation();
			setupTransformationParameters();
			setProviderPrepared(true);
		}
	}

	protected boolean isProviderPrepared() {
		return providerPrepared;
	}

	protected void setProviderPrepared(boolean prepared) {
		providerPrepared = prepared;
	}

	private void setupTransformedSourceSequenceImporter() throws IOException {
		LociImporterPlugin transformedSourceSequenceImporter = new LociImporterPlugin();
		try {
			transformedSourceSequenceImporter.open(getTransformedSourceFilePath().toString(),
					LociImporterPlugin.FLAG_METADATA_ALL);
		} catch (UnsupportedFormatException | IOException e) {
			transformedSourceSequenceImporter.close();
			throw new IOException("Could not open importer", e);
		}
		setTransformedSourceSequenceImporter(transformedSourceSequenceImporter);
		try {
			computeTransformedSourceFullSize();
			computeTransformedSourceFullRectangle();
			computeTransformedSourceSizeC();
			computeTransformedSourceDataType();
			computeTransformedSourceTileSize();
			computeTransformedSourceFullTileGrid();
		} catch (UnsupportedFormatException | IOException e) {
			throw new IOException("Could not retreive the needed transformed source image metadata", e);
		}
	}

	protected void setTransformedSourceSequenceImporter(LociImporterPlugin largeSequenceImporter) {
		transformedSourceSequenceImporter = largeSequenceImporter;
	}

	protected LociImporterPlugin getTransformedSourceSequenceImporter() {
		return transformedSourceSequenceImporter;
	}

	protected void computeTransformedSourceFullSize() throws UnsupportedFormatException, IOException {
		setTransformedSourceFullSize(
				new Dimension(MetaDataUtil.getSizeX(getTransformedSourceSequenceImporter().getOMEXMLMetaData(), 0),
						MetaDataUtil.getSizeY(getTransformedSourceSequenceImporter().getOMEXMLMetaData(), 0)));
	}

	private Dimension getTransformedSourceFullSize() {
		return transformedSourceFullSize;
	}

	private void setTransformedSourceFullSize(Dimension dimension) {
		transformedSourceFullSize = dimension;
	}

	protected void computeTransformedSourceFullRectangle() {
		setTransformedSourceFullRectangle(new Rectangle(getTransformedSourceFullSize()));
	}

	private Rectangle getTransformedSourceFullRectangle() {
		return transformedSourceFullRectangle;
	}

	private void setTransformedSourceFullRectangle(Rectangle rectangle) {
		transformedSourceFullRectangle = rectangle;
	}

	protected void computeTransformedSourceSizeC() throws UnsupportedFormatException, IOException {
		setTransformedSourceSizeC(MetaDataUtil.getSizeC(getTransformedSourceSequenceImporter().getOMEXMLMetaData(), 0));
	}

	protected int getTransformedSourceSizeC() {
		return transformedSourceSizeC;
	}

	private void setTransformedSourceSizeC(int sizeC) {
		transformedSourceSizeC = sizeC;
	}

	protected void computeTransformedSourceDataType() throws UnsupportedFormatException, IOException {
		setTransformedSourceDataType(
				MetaDataUtil.getDataType(getTransformedSourceSequenceImporter().getOMEXMLMetaData(), 0));
	}

	private DataType getTransformedSourceDataType() {
		return transformedSourceDataType;
	}

	private void setTransformedSourceDataType(DataType dataType) {
		transformedSourceDataType = dataType;
	}

	protected void computeTransformedSourceTileSize() throws UnsupportedFormatException, IOException {
		setTransformedSourceTileSize(new Dimension(getTransformedSourceSequenceImporter().getTileWidth(0),
				getTransformedSourceSequenceImporter().getTileHeight(0)));
	}

	private Dimension getTransformedSourceTileSize() {
		return transformedSourceTileSize;
	}

	private void setTransformedSourceTileSize(Dimension dimension) {
		transformedSourceTileSize = dimension;
	}

	protected void computeTransformedSourceFullTileGrid() {
		int maxTileX = (getTransformedSourceFullSize().width / getTransformedSourceTileSize().width)
				+ ((getTransformedSourceFullSize().width % getTransformedSourceTileSize().width != 0)? 1: 0);
		int maxTileY = (getTransformedSourceFullSize().height / getTransformedSourceTileSize().height)
				+ ((getTransformedSourceFullSize().height % getTransformedSourceTileSize().height != 0)? 1: 0);
		setTransformedSourceFullTileGrid(new Rectangle(0, 0, maxTileX, maxTileY));

	}

	private Rectangle getTransformedSourceFullTileGrid() {
		return this.transformedSourceFullTileGrid;
	}

	private void setTransformedSourceFullTileGrid(Rectangle rectangle) {
		this.transformedSourceFullTileGrid = rectangle;
	}

	protected void setupCacheTileProvider() throws IOException {
		try {
			setTransformedSourceTileProvider(
					new CachedLargeSequenceTileProvider.Builder(getTransformedSourceSequenceImporter()).build());
		} catch (IllegalArgumentException | IOException e) {
			throw new IOException("Could set cached tile provider for transformed source image", e);
		}
	}

	private CachedLargeSequenceTileProvider getTransformedSourceTileProvider() {
		return transformedSourceTileProvider;
	}

	private void setTransformedSourceTileProvider(CachedLargeSequenceTileProvider tileProvider) {
		transformedSourceTileProvider = tileProvider;
	}

	protected void setupOutputInformation() {
		computeOutputFullSize();
		computeOutputTileSizeC();
		computeOutputTileDataType();
	}

	protected void computeOutputFullSize() {
		int outputSizeX = MetaDataUtil.getSizeX(getOutputImageMetadata(), 0);
		int outputSizeY = MetaDataUtil.getSizeY(getOutputImageMetadata(), 0);
		setOutputFullSize(new Dimension(outputSizeX, outputSizeY));
	}

	private Dimension getOutputFullSize() {
		return this.outputFullSize;
	}

	private void setOutputFullSize(Dimension dimension) {
		this.outputFullSize = dimension;
	}

	protected void computeOutputTileSizeC() {
		int sizeC = MetaDataUtil.getSizeC(getOutputImageMetadata(), 0);
		setOutputTileSizeC(sizeC);
	}

	private int getOutputTileSizeC() {
		return outputTileSizeC;
	}

	private void setOutputTileSizeC(int sizeC) {
		outputTileSizeC = sizeC;
	}

	protected void computeOutputTileDataType() {
		DataType dataType = MetaDataUtil.getDataType(getOutputImageMetadata(), 0);
		setOutputTileDataType(dataType);
	}

	private DataType getOutputTileDataType() {
		return outputTileDataType;
	}

	private void setOutputTileDataType(DataType dataType) {
		outputTileDataType = dataType;
	}

	protected void setupTransformationParameters() {
		computeTransformationCoefficientModels();
		computeTransformationCoefficientScaleFactors();
		computeTransformationInterpolationScaleFactors();
	}

	protected abstract void computeTransformationCoefficientModels();

	protected BSplineModel getCoefficientsX() {
		return coefficientsX;
	}

	protected void setCoefficientsX(BSplineModel bSplineModel) {
		coefficientsX = bSplineModel;
	}

	protected BSplineModel getCoefficientsY() {
		return coefficientsY;
	}

	protected void setCoefficientsY(BSplineModel bSplineModel) {
		coefficientsY = bSplineModel;
	}

	private void computeTransformationCoefficientScaleFactors() {
		setCoefficientScaleFactorY((getCoefficientsX().getHeight() - 3) / (double) (getOutputFullSize().height - 1));
		setCoefficientScaleFactorX((getCoefficientsX().getWidth() - 3) / (double) (getOutputFullSize().width - 1));
	}

	protected double getCoefficientScaleFactorY() {
		return coefficientScaleFactorY;
	}

	protected void setCoefficientScaleFactorY(double scaleFactor) {
		coefficientScaleFactorY = scaleFactor;
	}

	protected double getCoefficientScaleFactorX() {
		return coefficientScaleFactorX;
	}

	protected void setCoefficientScaleFactorX(double scaleFactor) {
		coefficientScaleFactorX = scaleFactor;
	}

	private void computeTransformationInterpolationScaleFactors() {
		setSourceXScaleFactor(getTransformedSourceFullSize().getWidth() / getSubsampledSourceSize().getWidth());
		setSourceYScaleFactor(getTransformedSourceFullSize().getHeight() / getSubsampledSourceSize().getHeight());
	}

	protected double getSourceXScaleFactor() {
		return sourceXScaleFactor;
	}

	protected void setSourceXScaleFactor(double factorX) {
		sourceXScaleFactor = factorX;
	}

	protected double getSourceYScaleFactor() {
		return sourceYScaleFactor;
	}

	protected void setSourceYScaleFactor(double factorY) {
		sourceYScaleFactor = factorY;
	}

	protected void setCurrentTileIndex(Point tileIndex) {
		this.tileIndex = tileIndex;
		computeCurrentTileRectangle();
	}

	public Point getCurrentTileIndex() {
		return tileIndex;
	}

	private void computeCurrentTileRectangle() {
		computeCurrentTilePosition();
		setCurrentTileRectangle(new Rectangle(currentTilePosition, getOutputTileSize()));
		clipCurrentTileRectangleToOutput();
	}

	protected Rectangle getCurrentTileRectangle() {
		return currentTileRectangle;
	}

	protected void setCurrentTileRectangle(Rectangle rectangle) {
		this.currentTileRectangle = rectangle;
	}

	private void computeCurrentTilePosition() {
		setCurrentTilePosition(new Point(getCurrentTileIndex().x * getOutputTileSize().width,
				getCurrentTileIndex().y * getOutputTileSize().height));
	}

	protected Point getCurrentTilePosition() {
		return currentTilePosition;
	}

	protected void setCurrentTilePosition(Point point) {
		this.currentTilePosition = point;
	}

	protected void clipCurrentTileRectangleToOutput() {
		if (getCurrentTileRectangle().getMaxX() >= getOutputFullSize().width) {
			getCurrentTileRectangle().width -= getCurrentTileRectangle().getMaxX() - getOutputFullSize().width;
		}
		if (getCurrentTileRectangle().getMaxY() >= getOutputFullSize().height) {
			getCurrentTileRectangle().height -= getCurrentTileRectangle().getMaxY() - getOutputFullSize().height;
		}
	}

	private void loadNeededInputForCurrentOutputTile() throws Exception {
		computeNeededInputTilesForCurrentOutputTile();
		if (!getNeededSourceRectangle().isEmpty()) {
			setNeededSourceSequence(new Sequence(buildNeededInputFromCache(getNeededSourceRectangle())));
		}
	}

	private void computeNeededInputTilesForCurrentOutputTile() throws IOException {
		double minX = Double.POSITIVE_INFINITY, minY = Double.POSITIVE_INFINITY, maxX = Double.NEGATIVE_INFINITY,
				maxY = Double.NEGATIVE_INFINITY;

		for (int j = 0; j < getCurrentTileRectangle().height; j++) {
			double scaledY = ((j + getCurrentTileRectangle().y) * getCoefficientScaleFactorY()) + 1;
			for (int i = 0; i < getCurrentTileRectangle().width; i++) {
				double scaledX = ((i + getCurrentTileRectangle().x) * getCoefficientScaleFactorX()) + 1;
				double interpolatedX = getCoefficientsX().prepareForInterpolationAndInterpolateI(scaledX, scaledY, false,
						false);
				double interpolatedY = getCoefficientsY().prepareForInterpolationAndInterpolateI(scaledX, scaledY, false,
						false);
				minX = Math.min(interpolatedX, minX);
				minY = Math.min(interpolatedY, minY);
				maxX = Math.max(interpolatedX, maxX);
				maxY = Math.max(interpolatedY, maxY);
			}
		}

		int scaledMinX = (int) Math.floor(minX * sourceXScaleFactor);
		int scaledMinY = (int) Math.floor(minY * sourceYScaleFactor);
		int scaledMaxX = (int) Math.ceil(maxX * sourceXScaleFactor);
		int scaledMaxY = (int) Math.ceil(maxY * sourceYScaleFactor);

		setNeededSourceRectangle(new Rectangle(scaledMinX, scaledMinY, scaledMaxX - scaledMinX, scaledMaxY - scaledMinY));
		Rectangle.intersect(getNeededSourceRectangle(), getTransformedSourceFullRectangle(), getNeededSourceRectangle());
	}

	protected Sequence getNeededSourceSequence() {
		return neededSourceSequence;
	}

	protected void setNeededSourceSequence(Sequence sequence) {
		neededSourceSequence = sequence;
	}

	private IcyBufferedImage buildNeededInputFromCache(Rectangle neededRectangle) throws IOException {
		IcyBufferedImage neededImage = new IcyBufferedImage(neededRectangle.width, neededRectangle.height,
				getTransformedSourceSizeC(), getTransformedSourceDataType());
		Rectangle neededTileGrid = computeNeededTileGrid(neededRectangle);

		Point currentNeededTileIndex = new Point();
		Point currentNeededTilePosition = new Point();
		neededImage.beginUpdate();
		for (int j = 0; j < neededTileGrid.height; j++) {
			currentNeededTileIndex.y = neededTileGrid.y + j;
			currentNeededTilePosition.y = currentNeededTileIndex.y * getTransformedSourceTileProvider().getTileSize().height
					- neededRectangle.y;

			if (0 < currentNeededTileIndex.y || currentNeededTileIndex.y <= getTransformedSourceFullTileGrid().height) {
				for (int i = 0; i < neededTileGrid.width; i++) {
					currentNeededTileIndex.x = neededTileGrid.x + i;
					currentNeededTilePosition.x = currentNeededTileIndex.x
							* getTransformedSourceTileProvider().getTileSize().width - neededRectangle.x;

					if (0 < currentNeededTileIndex.x || currentNeededTileIndex.x <= getTransformedSourceFullTileGrid().width) {
						IcyBufferedImage currentNeededTileImage;
						try {
							currentNeededTileImage = getTransformedSourceTileProvider().getTile(currentNeededTileIndex);
						} catch (IllegalArgumentException | IOException e) {
							throw new IOException(
									String.format("Could not retrieve transformed source tile (%s)", currentNeededTileIndex), e);
						}
						neededImage.copyData(currentNeededTileImage, null, currentNeededTilePosition);
					}
				}
			}
		}
		neededImage.endUpdate();

		return neededImage;
	}

	private Rectangle computeNeededTileGrid(Rectangle neededRectangle) {
		int startTileX = Math.floorDiv((int) neededRectangle.getMinX(), getTransformedSourceTileSize().width);
		int startTileY = Math.floorDiv((int) neededRectangle.getMinY(), getTransformedSourceTileSize().height);
		int endTileX = Math.floorDiv((int) neededRectangle.getMaxX(), getTransformedSourceTileSize().width)
				+ (Math.floorMod((int) neededRectangle.getMaxX(), getTransformedSourceTileSize().width) != 0? 1: 0);
		int endTileY = Math.floorDiv((int) neededRectangle.getMaxY(), getTransformedSourceTileSize().height)
				+ (Math.floorMod((int) neededRectangle.getMaxY(), getTransformedSourceTileSize().height) != 0? 1: 0);

		Rectangle neededTileGrid = new Rectangle(startTileX, startTileY, endTileX - startTileX, endTileY - startTileY);
		neededTileGrid = neededTileGrid.intersection(getTransformedSourceFullTileGrid());
		return neededTileGrid;
	}

	private Rectangle getNeededSourceRectangle() {
		return neededSourceRectangle;
	}

	private void setNeededSourceRectangle(Rectangle rectangle) {
		this.neededSourceRectangle = rectangle;
	}

	private void writeOutputTileImage() throws InterruptedException {
		initializeOutputTileImage();
		if (getNeededSourceSequence() != null) {
			BSplineModel[] neededImageModel = new BSplineModel[getNeededSourceSequence().getSizeC()];
			IcyBufferedImageCursor outputCursor = new IcyBufferedImageCursor(getOutputTileImage());
			for (int c = 0; c < neededImageModel.length; c++) {
				neededImageModel[c] = new BSplineModel(
						IcyBufferedImageUtil.extractChannel(getNeededSourceSequence().getFirstImage(), c), false, 0);
				neededImageModel[c].setPyramidDepth(0);
				neededImageModel[c].startPyramids();
			}
			try {
				for (int c = 0; c < neededImageModel.length; c++) {
					neededImageModel[c].getThread().join();
				}
			} catch (InterruptedException e) {
				System.out.println("Should never reach this");
				e.printStackTrace();
			}
			double targetYScaleFactor = (getCoefficientsX().getHeight() - 3) / (double) (getOutputFullSize().height - 1);
			double targetXScaleFactor = (getCoefficientsX().getWidth() - 3) / (double) (getOutputFullSize().width - 1);
			for (int j = 0; j < getCurrentTileRectangle().height; j++) {
				if (Thread.interrupted())
					throw new InterruptedException();

				double scaledY = ((j + getCurrentTileRectangle().y) * targetYScaleFactor) + 1;
				if (j < getOutputTileImage().getHeight()) {
					for (int i = 0; i < getCurrentTileRectangle().width; i++) {
						double scaledX = ((i + getCurrentTileRectangle().x) * targetXScaleFactor) + 1;
						double interpolatedX = getCoefficientsX().prepareForInterpolationAndInterpolateI(scaledX, scaledY, false,
								false) * getSourceXScaleFactor();
						double interpolatedY = getCoefficientsY().prepareForInterpolationAndInterpolateI(scaledX, scaledY, false,
								false) * getSourceYScaleFactor();
						if (0 <= interpolatedX && interpolatedX < getTransformedSourceFullSize().width && 0 <= interpolatedY
								&& interpolatedY < getTransformedSourceFullSize().height) {
							for (int c = 0; c < neededImageModel.length; c++) {
								double value = neededImageModel[c].prepareForInterpolationAndInterpolateI(
										interpolatedX - getNeededSourceRectangle().x, interpolatedY - getNeededSourceRectangle().y, false,
										false);
								try {
									outputCursor.setSafe(i, j, c, value);
								} catch (ArrayIndexOutOfBoundsException e) {
									throw new RuntimeException(String.format("Index (%d, %d) out of bounds (%d, %d)", i, j,
											getOutputTileImage().getWidth(), getOutputTileImage().getHeight()), e);
								}
							}
						}
					}
				}
			}
			outputCursor.commitChanges();
		}
	}

	private void initializeOutputTileImage() {
		setOutputTileImage(new IcyBufferedImage(getCurrentTileRectangle().width, getCurrentTileRectangle().height,
				getOutputTileSizeC(), getOutputTileDataType()));
	}

	public void close() throws IOException {
		if (getTransformedSourceSequenceImporter() != null) {
			getTransformedSourceSequenceImporter().close();
		}
		if (getTransformedSourceTileProvider() != null) {
			try {
				getTransformedSourceTileProvider().close();
			} catch (Exception e) {
				throw new IOException("Could not close transformed source tile provider", e);
			}
		}
	}

}
