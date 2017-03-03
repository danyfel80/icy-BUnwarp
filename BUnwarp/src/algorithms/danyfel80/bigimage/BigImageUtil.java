package algorithms.danyfel80.bigimage;

import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import icy.common.exception.UnsupportedFormatException;
import icy.painter.Anchor2D;
import icy.roi.ROI;
import icy.roi.ROIUtil;
import icy.type.DataType;
import loci.formats.ome.OMEXMLMetadataImpl;
import plugins.kernel.importer.LociImporterPlugin;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.kernel.roi.roi2d.ROI2DRectangle;
import plugins.kernel.roi.roi2d.ROI2DShape;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class BigImageUtil {

	/**
	 * Gets the size of the image specified by the given path.
	 * 
	 * @param path
	 *          Path of the image
	 * @return Image dimension
	 */
	public static Dimension getSequenceSize(String path) {
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadataImpl imgProps = importer.getMetaData();
			int imgSizeX = imgProps.getPixelsSizeX(0).getValue();
			int imgSizeY = imgProps.getPixelsSizeY(0).getValue();

			return new Dimension(imgSizeX, imgSizeY);

		} catch (UnsupportedFormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				importer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return null;
	}

	/**
	 * Gets the channel count of the image specified by the given path.
	 * 
	 * @param path
	 *          Path of the image
	 * @return Image dimension. 0 if image is not correctly read.
	 */
	public static int getSequenceChannelCount(String path) {
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadataImpl imgProps = importer.getMetaData();
			return imgProps.getPixelsSizeC(0).getValue();

		} catch (UnsupportedFormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				importer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return 0;
	}

	/**
	 * Gets the data type of the image specified by the given path.
	 * 
	 * @param path
	 *          Path of the image
	 * @return Data type of the image.
	 */
	public static DataType getSequenceDataType(String path) {
		LociImporterPlugin importer = new LociImporterPlugin();
		try {
			importer.open(path, 0);
			OMEXMLMetadataImpl imgProps = importer.getMetaData();
			return DataType.getDataTypeFromPixelType(imgProps.getPixelsType(0));

		} catch (UnsupportedFormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				importer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return null;
	}

	public static List<ROI> getROIsInTile(List<ROI> rois, Rectangle tile) {
		List<ROI> result = new ArrayList<>();
		ROI2DRectangle rectROI = new ROI2DRectangle(tile);
		for (ROI roi : rois) {
			List<ROI> ROIsToIntersect = new ArrayList<>();
			ROIsToIntersect.add(roi);
			ROIsToIntersect.add(rectROI);
			result.add(ROIUtil.getIntersection(ROIsToIntersect));
		}
		return result;
	}

	public static void scaleROI(ROI roi, double scale) {
		if (roi instanceof ROI2DShape) {
			scaleROI2DShape((ROI2DShape) roi, scale);
		} else if (roi instanceof ROI2DArea) {
			scaleROI2DArea((ROI2DArea) roi, scale);
		} else {
			throw new IllegalArgumentException("the specified element cannot be scaled");
		}
	}

	private static void scaleROI2DShape(ROI2DShape roi, double scale) {
		List<Anchor2D> pts = roi.getControlPoints();
		Point2D pos = roi.getPosition();
		for (Anchor2D pt : pts) {
			double posX = pt.getX() - pos.getX();
			double posY = pt.getY() - pos.getY();
			posX *= scale;
			posY *= scale;
			pt.moveTo(posX + pos.getX(), posY + pos.getY());
			roi.controlPointPositionChanged(pt);
		}
		pos.setLocation(pos.getX() * scale, pos.getY() * scale);
		roi.setPosition2D(pos);
	}

	private static void scaleROI2DArea(ROI2DArea roi, double scale) {
		Rectangle bounds = roi.getBounds();
		boolean[] mask = roi.getBooleanMask(bounds, true);

		Rectangle scaledBounds = new Rectangle((int) (bounds.x * scale), (int) (bounds.y * scale),
		    (int) (bounds.width * scale), (int) (bounds.height * scale));
		boolean[] scaledMask = new boolean[scaledBounds.width*scaledBounds.height];
		for (int y = 0; y < scaledBounds.height; y++) {
			int yOff = y*scaledBounds.width;
			int scaledYOff = (int)((y/scale)*bounds.width);
			for (int x = 0; x < scaledBounds.width; x++) {
				scaledMask[x+yOff] = mask[(int)(x/scale) + scaledYOff];
			}
		}
		
		roi.setAsBooleanMask(scaledBounds, scaledMask);
	}
}
