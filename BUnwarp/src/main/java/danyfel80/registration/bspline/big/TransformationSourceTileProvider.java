package danyfel80.registration.bspline.big;

import danyfel80.registration.bspline.classic.BSplineModel;

public class TransformationSourceTileProvider extends TransformationTileProvider {

	@Override
	protected void computeTransformationCoefficientModels() {
		setCoefficientsX(new BSplineModel(getTransformation().getDirectDeformationCoefficientsX()));
		setCoefficientsY(new BSplineModel(getTransformation().getDirectDeformationCoefficientsY()));
	}
}
