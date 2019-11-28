package plugins.danyfel80.registration.elastic;

import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import algorithms.danyfel80.io.sequence.large.LargeSequenceHelper;
import danyfel80.registration.evaluation.PointCorrespondenceEvaluation;
import icy.common.exception.UnsupportedFormatException;
import icy.file.FileUtil;
import icy.roi.ROI;
import icy.util.XMLUtil;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarFile;
import plugins.kernel.roi.roi2d.ROI2DPoint;

public class PointRegistrationEvaluation extends EzPlug {

	private EzVarFile sourceImage;
	private EzVarFile targetImage;

	private List<? extends ROI2DPoint> sourceRois;
	private List<? extends ROI2DPoint> targetRois;
	private Dimension imageSize;

	@Override
	protected void initialize() {
		sourceImage = new EzVarFile("Source image", null);
		targetImage = new EzVarFile("Target image", null);

		addEzComponent(sourceImage);
		addEzComponent(targetImage);
	}

	@Override
	protected void execute() {
		sourceRois = getSequenceRois(sourceImage.getValue(true));
		targetRois = getSequenceRois(targetImage.getValue(true));
		try {
			imageSize = retrieveSequenceSize(targetImage.getValue(true));
		} catch (IOException | UnsupportedFormatException e) {
			throw new RuntimeException(e);
		}

		System.out.println(PointCorrespondenceEvaluation.evaluate(sourceRois, targetRois, imageSize));
	}

	private Dimension retrieveSequenceSize(File imageFile) throws IOException, UnsupportedFormatException {
		return LargeSequenceHelper.getImageDimension(imageFile);
	}

	private List<? extends ROI2DPoint> getSequenceRois(File imageFile) {
		String fileName = FileUtil.getFileName(imageFile.toString(), false);
		Path xmlFile = imageFile.toPath().resolveSibling(fileName + ".xml");
		if (Files.exists(xmlFile)) {
			Document xml = XMLUtil.loadDocument(xmlFile.toFile());
			Element rootElement = XMLUtil.getRootElement(xml);
			Element roisElement = XMLUtil.getElement(rootElement, "rois");
			List<? extends ROI> rois = ROI.loadROIsFromXML(roisElement);
			return rois.stream().filter(r -> (r instanceof ROI2DPoint)).map(r -> (ROI2DPoint) r).collect(Collectors.toList());
		} else {
			return Collections.emptyList();
		}
	}

	@Override
	public void clean() {}

}
