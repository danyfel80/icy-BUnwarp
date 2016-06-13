package algorithms.danyfel80.bigimage;

import java.awt.Dimension;
import java.awt.Point;
import java.io.File;
import java.io.IOException;

import icy.file.FileUtil;
import icy.image.IcyBufferedImage;
import icy.sequence.MetaDataUtil;
import icy.type.DataType;
import icy.util.OMEUtil;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.IFormatWriter;
import loci.formats.ome.OMEXMLMetadataImpl;
import loci.formats.out.OMETiffWriter;
import loci.formats.out.TiffWriter;
import loci.formats.tiff.IFD;
import loci.formats.tiff.TiffCompression;

/**
 * This class allows to save big images by tiles.
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageSaver {
	private File outFile;

	private int sizeX;
	private int sizeY;
	private int sizeC;
	private DataType dataType;
	
	private int tileSizeX;
	private int tileSizeY;
	private IFD ifd;
	
	private boolean isSeparateChannels;
	private boolean isLittleEndian;

	private TiffWriter writer;

	// private final Dimension tgtDim;
	// private final IcyColorModel tgtColorModel;
	// private final OMEXMLMetadata tgtMetadata;

	public BigImageSaver(File outFile, Dimension tgtDim, int sizeC, DataType dataType, Dimension tileSize)
	    throws ServiceException, FormatException, IOException {
		this.outFile = outFile;
		this.sizeX = tgtDim.width;
		this.sizeY = tgtDim.height;
		this.sizeC = sizeC;
		this.dataType = dataType;
		this.tileSizeX = tileSize.width;
		this.tileSizeY = tileSize.height;
		initializeWriter();
		writer.setSeries(0);
	}

	private void initializeWriter() throws ServiceException, FormatException, IOException {
		writer = new TiffWriter();
		
		OMEXMLMetadataImpl mdi = OMEUtil.createOMEMetadata();
		MetaDataUtil.setMetaData(mdi, sizeX, sizeY, sizeC, 1, 1, -1, -1, dataType, isSeparateChannels);
		this.isLittleEndian = !mdi.getPixelsBinDataBigEndian(0, 0);
//		System.out.println("little endian = " + isLittleEndian);
		this.isSeparateChannels = getSeparateChannelFlag(writer, sizeC, dataType);
//		System.out.println("separate channels = " + isSeparateChannels);
		writer.setMetadataRetrieve(mdi);
		writer.setCompression(TiffCompression.LZW.getCodecName());
		writer.setWriteSequentially(true);
		writer.setInterleaved(false);
		writer.setBigTiff(true);
		if (FileUtil.exists(outFile.getAbsolutePath())) {
			FileUtil.delete(outFile, true);
		}
		writer.setId(outFile.getAbsolutePath());
		
		ifd = new IFD();
		long[] rowPerStrip = new long[1];
		rowPerStrip[0] = tileSizeY;
		ifd.put(IFD.TILE_WIDTH, tileSizeX);
		ifd.put(IFD.TILE_LENGTH, tileSizeY);
		ifd.put(IFD.ROWS_PER_STRIP, rowPerStrip);
	}
	
	/**
	 * Return the separate channel flag from specified writer and color space
	 */
	private static boolean getSeparateChannelFlag(IFormatWriter writer, int numChannel, DataType dataType) {
		if (writer instanceof OMETiffWriter)
			return (numChannel == 2) || (numChannel > 4) || (dataType.getSize() > 1);
		// return numChannel > 1;
		return false;
	}

	public void closeWriter() throws IOException {
		this.writer.close();
	}
	
	public void saveTile(IcyBufferedImage srcIBI, Point tgtPoint)
	    throws ServiceException, IOException, FormatException {

		byte[] data = null;

		// separated channel data
		if (isSeparateChannels) {
			for (int c = 0; c < sizeC; c++) {
				
				if (srcIBI != null) {
					data = srcIBI.getRawData(c, isLittleEndian);
					writer.saveBytes(c, data, ifd, tgtPoint.x, tgtPoint.y, srcIBI.getSizeX(), srcIBI.getSizeY());
				}
			}
		} else {
			if (srcIBI != null) {
				data = srcIBI.getRawData(isLittleEndian);
				writer.saveBytes(0, data, ifd, tgtPoint.x, tgtPoint.y, srcIBI.getSizeX(), srcIBI.getSizeY());
			}
		}
	}
}
