package algorithms.danyfel80.bigimage;

import java.awt.Dimension;
import java.awt.Point;
import java.io.File;
import java.io.IOException;

import icy.file.FileUtil;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import icy.sequence.MetaDataUtil;
import icy.type.DataType;
import icy.util.OMEUtil;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.IFormatWriter;
import loci.formats.meta.MetadataRetrieve;
import loci.formats.out.OMETiffWriter;
import loci.formats.out.TiffWriter;
import loci.formats.tiff.IFD;
import loci.formats.tiff.TiffCompression;
import ome.xml.meta.OMEXMLMetadata;

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

	private int tileSizeX = 256;
	private int tileSizeY = 256;
	private IFD ifd;

	private boolean isSeparateChannels;
	private boolean isLittleEndian;

	private TiffWriter writer;

	// private final Dimension tgtDim;
	// private final IcyColorModel tgtColorModel;
	// private final OMEXMLMetadata tgtMetadata;

	public BigImageSaver(File outFile, Dimension tgtDim, int sizeC, DataType dataType, Dimension tileSize)
			throws ServiceException, FormatException, IOException {

		int width = tileSize.width;
		int height = tileSize.height;
		if (width % 256 != 0 || height % 256 != 0) {
			throw new FormatException("tile size must be multiple of 256 on width and height.");
		}
		this.outFile = outFile;
		this.sizeX = tgtDim.width;
		this.sizeY = tgtDim.height;
		this.sizeC = sizeC;
		this.dataType = dataType;
		initializeWriter();
		writer.setSeries(0);
	}

	private void initializeWriter() throws ServiceException, FormatException, IOException {
		writer = new OMETiffWriter();
		writer.setCompression(TiffCompression.LZW.getCodecName());

		OMEXMLMetadata mdi = OMEUtil.createOMEXMLMetadata();
		this.isSeparateChannels = getSeparateChannelFlag(writer, sizeC, dataType);
		MetaDataUtil.setMetaData(mdi, sizeX, sizeY, sizeC, 1, 1, dataType, isSeparateChannels);

		writer.setMetadataRetrieve((MetadataRetrieve) mdi);
		writer.setWriteSequentially(true);
		writer.setInterleaved(false);
		writer.setBigTiff(true);
		if (FileUtil.exists(outFile.getAbsolutePath())) {
			FileUtil.delete(outFile, true);
		}
		writer.setId(outFile.getAbsolutePath());
		writer.setSeries(0);

		this.isLittleEndian = !writer.getMetadataRetrieve().getPixelsBinDataBigEndian(0, 0).booleanValue();

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
		// Only if we have more than 1 channel
		if (numChannel > 1) {
			// Only channel amount is different than three and data type size is more
			// than 1 byte
			return (numChannel != 3) || (dataType.getSize() > 1);
		}
		// Otherwise use fused channels
		return false;
	}

	public void closeWriter() throws IOException {
		this.writer.close();
	}

	public synchronized void saveTile(IcyBufferedImage srcIBI, Point tgtPoint)
			throws ServiceException, IOException, FormatException {

		byte[] data = null;

		int tx = 0, ty = 0;
		if (srcIBI != null) {
			tx = (srcIBI.getWidth() + tileSizeX - 1) / tileSizeX;
			ty = (srcIBI.getHeight() + tileSizeY - 1) / tileSizeY;
		}

		// separated channel data
		if (isSeparateChannels) {
			for (int c = 0; c < sizeC; c++) {
				if (srcIBI != null) {
					for (int ti = 0; ti < tx; ti++) {
						int currentTilePosX = ti * tileSizeX;
						int currTileSizeX = (currentTilePosX + tileSizeX < srcIBI.getWidth()) ? tileSizeX
								: srcIBI.getWidth() - currentTilePosX;
						for (int tj = 0; tj < ty; tj++) {
							int currentTilePosY = tj * tileSizeY;
							int currTileSizeY = (currentTilePosY + tileSizeY < srcIBI.getHeight()) ? tileSizeY
									: srcIBI.getHeight() - currentTilePosY;
							IcyBufferedImage srcTile = IcyBufferedImageUtil.getSubImage(srcIBI, currentTilePosX, currentTilePosY,
									currTileSizeX, currTileSizeY);
							data = srcIBI.getRawData(c, isLittleEndian);
							writer.saveBytes(c, data, ifd, tgtPoint.x + currentTilePosX, tgtPoint.y + currentTilePosY,
									srcTile.getSizeX(), srcTile.getSizeY());
						}
					}
				}
			}
		}
		// All Channels in the same data block
		else {
			if (srcIBI != null) {
				for (int ti = 0; ti < tx; ti++) {
					int currentTilePosX = ti * tileSizeX;
					int currTileSizeX = (currentTilePosX + tileSizeX < srcIBI.getWidth()) ? tileSizeX
							: srcIBI.getWidth() - currentTilePosX;
					for (int tj = 0; tj < ty; tj++) {
						int currentTilePosY = tj * tileSizeY;
						int currTileSizeY = (currentTilePosY + tileSizeY < srcIBI.getHeight()) ? tileSizeY
								: srcIBI.getHeight() - currentTilePosY;
						IcyBufferedImage srcTile = IcyBufferedImageUtil.getSubImage(srcIBI, currentTilePosX, currentTilePosY,
								currTileSizeX, currTileSizeY);
						data = srcTile.getRawData(isLittleEndian);
						try {
							writer.saveBytes(0, data, ifd, tgtPoint.x + currentTilePosX, tgtPoint.y + currentTilePosY,
									srcTile.getSizeX(), srcTile.getSizeY());
						} catch (Exception e) {
							System.out.println(data);
							System.out.println(ifd);
							System.out.println(tgtPoint);
							System.out.println(srcIBI.getBounds());
							throw e;
						}
					}
				}
			}
		}
	}
}
