/*
 * Copyright 2010-2016 Institut Pasteur.
 * 
 * This file is part of Icy.
 * 
 * Icy is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Icy is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Icy. If not, see <http://www.gnu.org/licenses/>.
 */
package plugins.danyfel80.bigimage;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.IntStream;

import org.apache.commons.io.FilenameUtils;

import algorithms.danyfel80.bigimage.BigImageToTiffConverter;
import algorithms.danyfel80.bigimage.BigImageUtil;
import icy.common.exception.UnsupportedFormatException;
import icy.common.listener.RichProgressListener;
import icy.gui.dialog.MessageDialog;
import icy.gui.frame.progress.AnnounceFrame;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarIntegerArrayNative;
import plugins.adufour.ezplug.EzVarListener;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class BigImageConvertToTiff extends EzPlug implements Block, EzStoppable {

	EzVarFile								varInFile;
	EzVarIntegerArrayNative	varInChannels;

	BigImageToTiffConverter converter;

	/*
	 * (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#initialize()
	 */
	@Override
	protected void initialize() {
		this.varInFile = new EzVarFile("File", null);
		this.varInChannels = new EzVarIntegerArrayNative("Using channels", new int[][] { { 1, 1, 1 } }, true);

		this.varInFile.addVarChangeListener(new EzVarListener<File>() {
			@Override
			public void variableChanged(EzVar<File> source, File newValue) {
				if (newValue != null && newValue.exists()) {
					int sizeC = 0;
					try {
						sizeC = BigImageUtil.getSequenceChannelCount(newValue.getPath());
					} catch (IOException | UnsupportedFormatException e) {
						MessageDialog.showDialog("Error", e.getMessage(), MessageDialog.ERROR_MESSAGE);
						return;
					}
					varInChannels.setValue(IntStream.range(0, sizeC).map(x -> 1).toArray());
				}
			}
		});

		addEzComponent(varInFile);
		addEzComponent(varInChannels);
	}

	/*
	 * (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#execute()
	 */
	@Override
	protected void execute() {
		try {
			validateInput();
		} catch (Exception e) {
			MessageDialog.showDialog("Error", e.getMessage(), MessageDialog.ERROR_MESSAGE);
			return;
		}

		File inFile = varInFile.getValue();
		File outFile = Paths.get(FilenameUtils.getFullPath(inFile.getAbsolutePath()),
				FilenameUtils.getBaseName(inFile.getAbsolutePath()) + "_ToTiff.tif").toFile();

		int[] intChannels = varInChannels.getValue();
		boolean[] channels = new boolean[intChannels.length];
		IntStream.range(0, intChannels.length).forEach(i -> channels[i] = intChannels[i] != 0);

		converter = new BigImageToTiffConverter(inFile, outFile, channels);
		converter.addProgressListener(new RichProgressListener() {
			@Override
			public boolean notifyProgress(double position, double length, String comment, Object data) {
				if (getUI() != null) {
					getUI().setProgressBarValue(position / length);
					getUI().setProgressBarMessage(comment);
				}
				return true;
			}
		});

		ExecutorService executor = Executors.newSingleThreadExecutor();
		Future<Void> futureRun = executor.submit(converter);

		try {
			futureRun.get();
		} catch (InterruptedException | ExecutionException e) {
			MessageDialog.showDialog("Error", e.getMessage(), MessageDialog.ERROR_MESSAGE);
			e.printStackTrace();
			return;
		} finally {
			executor.shutdownNow();
		}

		new AnnounceFrame("Conversion done", 3);
	}

	private void validateInput() throws IllegalArgumentException, UnsupportedFormatException, IOException {
		if (varInFile.getValue() == null || !varInFile.getValue().exists())
			throw new IllegalArgumentException("Invalid File: null file or doesn't exist.");
		int sizeC = BigImageUtil.getSequenceChannelCount(varInFile.getValue().getPath());
		if (varInChannels.getValue().length != sizeC)
			throw new IllegalArgumentException("Using channels has less channels than those in the image file");
	}

	/*
	 * (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#clean()
	 */
	@Override
	public void clean() {}

	/*
	 * (non-Javadoc)
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.
	 * VarList)
	 */
	@Override
	public void declareInput(VarList inputMap) {

	}

	/*
	 * (non-Javadoc)
	 * @see
	 * plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util
	 * .VarList)
	 */
	@Override
	public void declareOutput(VarList outputMap) {

	}

	@Override
	public void stopExecution() {
		if (converter != null) converter.interrupt();
	}
}
