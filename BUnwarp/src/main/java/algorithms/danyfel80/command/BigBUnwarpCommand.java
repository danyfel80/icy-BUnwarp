/*
 * Copyright 2010-2018 Institut Pasteur.
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
package algorithms.danyfel80.command;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class BigBUnwarpCommand implements CommandProcess<String> {

	private static final String COMMAND = "big-unwarp";
	private static final String NAME = "Big Unwarp";
	private static final String[] ARGS_DESCRIPTION = new String[] {"sourceImg", "targetImg", "transformedSourceImg",
			"transformedTargetImg", "subsampleFactor", "minimumDeformation(0-3)", "maximumDeformation(0-4)", "wDiv", "wCurl",
			"wImg", "wConsistency", "stopThreshold", "registrationMode(0: fast, 1:accurate, 2:mono)"};
	private static final String DESCRIPTION = "Performs large image registration using the BUnwarp method from Argandas";

	private String[] args;

	@Override
	public String getCommand() {
		return COMMAND;
	}

	@Override
	public String getName() {
		return NAME;
	}

	@Override
	public String[] getArgumentsDescription() {
		return ARGS_DESCRIPTION;
	}

	@Override
	public String getDescription() {
		return DESCRIPTION;
	}

	@Override
	public CommandProcess<String> setArguments(String[] args) {
		this.args = args;
		return this;
	}

	@Override
	public String call() throws Exception {
		return null;
		//		String srcPath = this.args[0];
		//		String tgtPath = this.args[1];
		//
		//		String transformedSrcPath = this.args[2];
		//		String transformedTgtPath = this.args[3];
		//
		//		String srcResultPath;
		//		String tgtResultPath;
		//		String transformedSrcResultPath;
		//		String transformedTgtResultPath;
		//
		//		int subsampleFactor = Integer.parseInt(this.args[4]);
		//		MinimumScaleDeformationEnum iniDef = MinimumScaleDeformationEnum.values()[Integer.parseInt(this.args[5])];
		//		MaximumScaleDeformationEnum finDef = MaximumScaleDeformationEnum.values()[Integer.parseInt(this.args[6])];
		//
		//		double divWeight = Double.parseDouble(this.args[7]);
		//		double curlWeight = Double.parseDouble(this.args[8]);
		//		double imageWeight = Double.parseDouble(this.args[9]);
		//		double consistencyWeight = Double.parseDouble(this.args[10]);
		//
		//		double stopThreshold = Double.parseDouble(this.args[11]);
		//
		//		RegistrationModeEnum mode = RegistrationModeEnum.values()[Integer.parseInt(this.args[12])];
		//
		//		srcResultPath = FilenameUtils.getFullPath(srcPath);
		//		srcResultPath += FilenameUtils.getBaseName(srcPath);
		//		srcResultPath += "_BUnwarp.tif";
		//		// srcResultPath += FilenameUtils.getExtension(srcPath);
		//
		//		tgtResultPath = FilenameUtils.getFullPath(tgtPath);
		//		tgtResultPath += FilenameUtils.getBaseName(tgtPath);
		//		tgtResultPath += "_BUnwarp.tif";
		//		// tgtResultPath += FilenameUtils.getExtension(tgtPath);
		//
		//		transformedSrcResultPath = FilenameUtils.getFullPath(transformedSrcPath);
		//		transformedSrcResultPath += FilenameUtils.getBaseName(transformedSrcPath);
		//		transformedSrcResultPath += "_BUnwarp.tif";
		//		// transformedSrcResultPath +=
		//		// FilenameUtils.getExtension(transformedSrcPath);
		//
		//		transformedTgtResultPath = FilenameUtils.getFullPath(transformedTgtPath);
		//		transformedTgtResultPath += FilenameUtils.getBaseName(transformedTgtPath);
		//		transformedTgtResultPath += "_BUnwarp.tif";
		//		// transformedTgtResultPath +=
		//		// FilenameUtils.getExtension(transformedTgtPath);
		//
		//		long startTime = System.nanoTime();
		//
		//		List<ROI2DPoint> srcLandmarks = null;
		//		List<ROI2DPoint> tgtLandmarks = null;
		//
		//		ROI2DPolygon srcMask = null;
		//		ROI2DPolygon tgtMask = null;
		//
		//		BigBUnwarpper bu = new BigBUnwarpper(srcPath, tgtPath, transformedSrcPath, transformedTgtPath, srcResultPath,
		//				tgtResultPath, transformedSrcResultPath, transformedTgtResultPath, srcLandmarks, tgtLandmarks, srcMask, tgtMask,
		//				subsampleFactor, iniDef.getNumber(), finDef.getNumber(), divWeight, curlWeight,
		//				0/* inLandmarkWeight.getValue() */, imageWeight, consistencyWeight, stopThreshold, false, mode.getNumber(),
		//				(progress, message, data) -> {
		//					return false;
		//				});
		//		Thread but = new Thread(bu);
		//		but.start();
		//		try {
		//			but.join();
		//			but = null;
		//		} catch (InterruptedException e) {
		//			System.err.println("Thread interrupted: " + e.getMessage());
		//		}
		//
		//		long endTime = System.nanoTime();
		//		long totalTime = (endTime - startTime);
		//		return String.format("Done (%d millisecs)", totalTime / 1000000);
	}

}
