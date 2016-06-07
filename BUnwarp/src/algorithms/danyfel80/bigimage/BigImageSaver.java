package algorithms.danyfel80.bigimage;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;

import icy.file.FileUtil;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.image.colormodel.IcyColorModel;
import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.util.OMEUtil;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.IFormatWriter;
import loci.formats.ome.OMEXMLMetadata;
import loci.formats.ome.OMEXMLMetadataImpl;
import loci.formats.out.OMETiffWriter;
import loci.formats.tiff.IFD;

/**
 * This class allows to save big images by tiles.
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageSaver {

	private final OMETiffWriter writer;
	private final boolean littleEndian;
	private final int sizeC;
//	private final Dimension tgtDim;
//	private final IcyColorModel tgtColorModel;
//	private final OMEXMLMetadata tgtMetadata;
	private final boolean separateChannel;

	public BigImageSaver(File tgtFile, Dimension tgtDim, IcyColorModel tgtColorModel, OMEXMLMetadata tgtMetadata)
	    throws FormatException, IOException {
		writer = new OMETiffWriter();
//		this.tgtDim = tgtDim;
//		this.tgtColorModel = tgtColorModel;
//		this.tgtMetadata = tgtMetadata;
		try {
			writer.setCompression("LZW");
		} catch (FormatException e) {
			/* no compression */}

		// first delete the file else LOCI won't save it correctly
		if (tgtFile.exists())
			tgtFile.delete();
		// ensure parent directory exist
		FileUtil.ensureParentDirExist(tgtFile);

		sizeC = tgtColorModel.getNumComponents();
		this.separateChannel = getSeparateChannelFlag(writer, sizeC, tgtColorModel.getDataType_());

		// set settings
		writer.setFramesPerSecond(1);
		// generate metadata
		final OMEXMLMetadataImpl writerMetadata = OMEUtil.createOMEMetadata(tgtMetadata);

		try {
			MetaDataUtil.setMetaData(writerMetadata, (int) tgtDim.getWidth(), (int) tgtDim.getHeight(), sizeC, 1, 1,
			    tgtColorModel.getDataType_(), separateChannel);
		} catch (ServiceException e) {
			e.printStackTrace();
		}
		writer.setMetadataRetrieve(writerMetadata);
		// no interleave (XP default viewer want interleaved channel to correctly
		// read image)
		writer.setInterleaved(false);
		// save on big image
		writer.setBigTiff(true);
		// set id
		writer.setId(tgtFile.getAbsolutePath());
		// init
		writer.setSeries(0);
		// usually give better save performance
		writer.setWriteSequentially(true);

		// get endianess
		littleEndian = !writer.getMetadataRetrieve().getPixelsBinDataBigEndian(0, 0).booleanValue();
	}

	public void saveTile(Sequence seq, Rectangle srcRect, Point tgtPoint)
	    throws ServiceException, IOException, FormatException {

		byte[] data = null;

		int imageIndex = 0;
		if (srcRect == null)
			srcRect = new Rectangle(0, 0, seq.getSizeX(), seq.getSizeY());
		final IcyBufferedImage image = IcyBufferedImageUtil.getSubImage(seq.getFirstImage(), srcRect);
		
		System.out.println(srcRect);
		
		IFD ifd = new IFD();
		ifd.put(IFD.TILE_WIDTH, (long)srcRect.width);
		ifd.put(IFD.TILE_LENGTH, (long)srcRect.height);

		System.out.println(ifd.getTileWidth());
		// separated channel data
		if (separateChannel) {
			for (int c = 0; c < sizeC; c++) {
				if (image != null) {
					// avoid multiple allocation
					data = image.getRawData(c, data, 0, littleEndian);
					writer.saveBytes(imageIndex, data, ifd, tgtPoint.x, tgtPoint.y, image.getSizeX(), image.getSizeY());
				}

				imageIndex++;
			}
		} else {
			if (image != null) {
				// avoid multiple allocation
				data = image.getRawData(data, 0, littleEndian);
				writer.saveBytes(imageIndex, data, ifd, tgtPoint.x, tgtPoint.y, image.getSizeX(), image.getSizeY());
			}

			imageIndex++;
		}
	}

	public void close() throws IOException {
		this.writer.close();
	}

	/**
	 * Return the separate channel flag from specified writer and color space
	 */
	private static boolean getSeparateChannelFlag(IFormatWriter writer, int numChannel, DataType dataType) {
		if (writer instanceof OMETiffWriter)
			return (numChannel == 2) || (numChannel > 4) || (dataType.getSize() > 1);

		return false;
	}
}
