package danyfel80.registration.evaluation;

import java.awt.Dimension;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.util.Pair;

import plugins.kernel.roi.roi2d.ROI2DPoint;

public class PointCorrespondenceEvaluation {

	public static double evaluate(List<? extends ROI2DPoint> points1, List<? extends ROI2DPoint> points2,
			Dimension imageSize) {
		Objects.requireNonNull(points1, "Points 1 is null");
		Objects.requireNonNull(points2, "Points 2 is null");
		PointCorrespondenceEvaluation evaluator = new PointCorrespondenceEvaluation(points1, points2, imageSize);
		evaluator.compute();
		return evaluator.getScore();
	}

	private Map<String, ? extends ROI2DPoint> pts1, pts2;
	private Dimension imgSize;

	public PointCorrespondenceEvaluation(List<? extends ROI2DPoint> points1, List<? extends ROI2DPoint> points2,
			Dimension imageSize) {
		this.pts1 = convertToPointMap(points1);
		this.pts2 = convertToPointMap(points2);
		this.imgSize = imageSize;
	}

	private Map<String, ? extends ROI2DPoint> convertToPointMap(List<? extends ROI2DPoint> points) {
		return points.stream().collect(Collectors.toMap(p -> p.getName(), p -> p));
	}

	private void compute() {
		filterCommonPoints();
		evaluatePointDistances();
		computeMedianScore();
	}

	private Set<String> commonPtsNames;

	private void filterCommonPoints() {
		Set<String> pts1Names = pts1.keySet();
		Set<String> pts2Names = pts2.keySet();
		this.commonPtsNames = new HashSet<>(pts1Names);
		commonPtsNames.retainAll(pts2Names);
	}

	private Map<String, Double> distances;

	private void evaluatePointDistances() {
		double normalizer = Math.sqrt(imgSize.getWidth() * imgSize.getWidth() + imgSize.getHeight() * imgSize.getHeight());
		distances = commonPtsNames.stream()
				.map(k -> new Pair<String, Double>(k, euclideanDistance(pts1.get(k), pts2.get(k)) / normalizer))
				.collect(Collectors.toMap(pr -> pr.getKey(), pr -> pr.getValue()));
	}

	private double euclideanDistance(ROI2DPoint point1, ROI2DPoint point2) {
		return point1.getPoint().distance(point2.getPoint());
	}

	double score;

	private void computeMedianScore() {
		double[] sortedDists = distances.values().stream().mapToDouble(Double::doubleValue).sorted().toArray();
		System.out.println(sortedDists.length);
		if (sortedDists.length % 2 == 1)
			score = sortedDists[sortedDists.length / 2];
		else
			score = (sortedDists[(sortedDists.length / 2) - 1] + sortedDists[sortedDists.length / 2]) / 2d;
	}

	private double getScore() {
		return score;
	}
}
