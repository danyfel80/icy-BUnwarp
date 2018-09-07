package danyfel80.registration.bspline.big;

import danyfel80.registration.bspline.classic.BSplineModel;

public class TransformationTargetTileProvider extends TransformationTileProvider {

	@Override
	protected void computeTransformationCoefficientModels() {
		setCoefficientsX(new BSplineModel(getTransformation().getInverseDeformationCoefficientsX()));
		setCoefficientsY(new BSplineModel(getTransformation().getInverseDeformationCoefficientsY()));
	}
}
