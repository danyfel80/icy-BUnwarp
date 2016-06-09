package plugins.danyfel80.bigimage;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;

import algorithms.danyfel80.bigimage.BigImageSaver;
import icy.sequence.Sequence;
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
		Sequence seq = inSeq.getValue();
		File file = inFile.getValue();
		int tileSize = inTileSize.getValue();
		int tileSizeX = tileSize;
		int tileSizeY = tileSize;

		this.getUI().setProgressBarMessage("Saving");

		BigImageSaver saver;
		try {
			saver = new BigImageSaver(file, new Dimension(seq.getWidth(), seq.getHeight()), seq.getColorModel(),
			    seq.getMetadata());

			for (int i = 0; i < seq.getSizeX(); i += tileSize) {
				if (i + tileSizeX > seq.getSizeX())
					tileSizeX = seq.getSizeX() - i;
				for (int j = 0; j < seq.getSizeY(); j += tileSize) {
					if (j + tileSizeY > seq.getSizeY())
						tileSizeY = seq.getSizeY() - j;
					
					try {
						//addSequence(new Sequence(IcyBufferedImageUtil.getSubImage(seq.getFirstImage(), new Rectangle(i, j, tileSizeX, tileSizeY))));
						saver.saveTile(seq, new Rectangle(i, j, tileSizeX, tileSizeY), new Point(i, j));
					} catch (ServiceException | IOException | FormatException e) {
						e.printStackTrace();
					}
				}
			}
			
			saver.close();
		} catch (FormatException | IOException e1) {
			e1.printStackTrace();
		}
		this.getUI().setProgressBarMessage("Done");
	}

	@Override
	public void clean() {
	}
}
