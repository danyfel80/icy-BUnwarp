package plugins.danyfel80.bigimage;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;

import algorithms.danyfel80.bigimage.BigImageSaver;
import icy.image.IcyBufferedImage;
import icy.image.IcyBufferedImageUtil;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;

public class SaveBigImage extends EzPlug {

	EzVarSequence inSeq = new EzVarSequence("Sequence");
	EzVarFile inFile = new EzVarFile("File", "");
	EzVarInteger inTileSize = new EzVarInteger("Tile size");

	@Override
	protected void initialize() {
		addEzComponent(inSeq);
		addEzComponent(inFile);
		addEzComponent(inTileSize);
	}

	@Override
	protected void execute() {
		IcyBufferedImage ibi = inSeq.getValue().getFirstImage();
		
		File file = inFile.getValue();
		int tileSize = inTileSize.getValue();
		int tileSizeX = tileSize;
		int tileSizeY = tileSize;

		this.getUI().setProgressBarMessage("Saving");

		int x, y;
		int w, h;
		int sizeX, sizeY;
		int n, m;
		int diffWidth, diffHeight;

		sizeX = ibi.getWidth();
		sizeY = ibi.getHeight();
		if (tileSizeX <= 0)
			tileSizeX = sizeX;
		if (tileSizeY <= 0)
			tileSizeY = sizeY;
		n = sizeX / tileSizeX;
		m = sizeY / tileSizeY;
		if (n == 0) {
			tileSizeX = sizeX;
			n = 1;
		}
		if (m == 0) {
			tileSizeY = sizeY;
			m = 1;
		}
		diffWidth = sizeX - n * tileSizeX;
		diffHeight = sizeY - m * tileSizeY;
		if (diffWidth > 0)
			n++;
		if (diffHeight > 0)
			m++;

		BigImageSaver saver = null;
		try {
			saver = new BigImageSaver(file, new Dimension(ibi.getWidth(), ibi.getHeight()), ibi.getSizeC(),
					ibi.getDataType_(), new Dimension(tileSizeX, tileSizeY));
		} catch (ServiceException | FormatException | IOException e1) {
			e1.printStackTrace();
			return;
		}
		
		try {
			for (int i = 0; i < m; i++) {
				if (diffHeight > 0 && i == (m - 1)) {
					y = sizeY - diffHeight;
					h = diffHeight;
				} else {
					y = tileSizeY * i;
					h = tileSizeY;
				}
				for (int j = 0; j < n; j++) {
					if (diffWidth > 0 && j == (n - 1)) {
						x = sizeX - diffWidth;
						w = diffWidth;
					} else {
						x = tileSizeX * j;
						w = tileSizeX;
					}
					
					Rectangle rect = new Rectangle(x, y, w, h);
					Point point = new Point(x, y);
//					System.out.println("Saving tile " + rect + " at " + point);
					saver.saveTile(IcyBufferedImageUtil.getSubImage(ibi, rect), point);

				}
			}

		} catch (ServiceException | FormatException | IOException e) {
			e.printStackTrace();
		} finally {
			try {
				saver.closeWriter();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		this.getUI().setProgressBarMessage("Done");
	}

	@Override
	public void clean() {
	}
}
