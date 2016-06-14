package algorithms.danyfel80.bigimage;

import java.awt.Dimension;
import java.io.IOException;

import icy.common.exception.UnsupportedFormatException;
import icy.type.DataType;
import loci.formats.ome.OMEXMLMetadataImpl;
import plugins.kernel.importer.LociImporterPlugin;

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
	
}
