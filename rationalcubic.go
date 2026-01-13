package jaeckel

import (
	"math"
)

const (
	// minimumControlParam = -(1 - √ε) where ε = 2^-52
	minimumControlParam = -0.99999998509883880615234375

	// maximumControlParam = 2/ε² where ε = 2^-52
	maximumControlParam = 4.0564819207303340847894502572032e+31
)

// interpolate performs rational cubic interpolation.
func interpolate(x, xL, xR, yL, yR, dL, dR, r float64) float64 {
	h := xR - xL
	if math.Abs(h) <= 0 {
		return 0.5 * (yL + yR)
	}

	t := (x - xL) / h
	if r < maximumControlParam {
		omt := 1.0 - t
		t2 := t * t
		omt2 := omt * omt
		numerator := yR*t2*t + (r*yR-h*dR)*t2*omt + (r*yL+h*dL)*t*omt2 + yL*omt2*omt
		denominator := 1.0 + (r-3.0)*t*omt
		return numerator / denominator
	}

	// linear interpolation for very large r
	return yR*t + yL*(1.0-t)
}

// controlParameterToFitSecondDerivativeAtLeftSide calculates r to match a target second derivative at the left boundary.
func controlParameterToFitSecondDerivativeAtLeftSide(xL, xR, yL, yR, dL, dR, secondDerivativeL float64) float64 {
	h := xR - xL
	numerator := 0.5*h*secondDerivativeL + (dR - dL)
	denominator := (yR-yL)/h - dL

	if isZero(denominator) {
		if numerator > 0 {
			return maximumControlParam
		}
		return minimumControlParam
	}
	return numerator / denominator
}

// controlParameterToFitSecondDerivativeAtRightSide calculates r to match a target second derivative at the right boundary.
func controlParameterToFitSecondDerivativeAtRightSide(xL, xR, yL, yR, dL, dR, secondDerivativeR float64) float64 {
	h := xR - xL
	numerator := 0.5*h*secondDerivativeR + (dR - dL)
	denominator := dR - (yR-yL)/h

	if isZero(denominator) {
		if numerator > 0 {
			return maximumControlParam
		}
		return minimumControlParam
	}
	return numerator / denominator
}

// minimumControlParameter ensures shape preservation (monotonicity/convexity).
func minimumControlParameter(dL, dR, s float64, preferShapePreservation bool) float64 {
	monotonic := dL*s >= 0 && dR*s >= 0
	convex := dL <= s && s <= dR
	concave := dL >= s && s >= dR

	if !monotonic && !convex && !concave {
		return minimumControlParam
	}

	dRMinusDl := dR - dL
	dRMinusS := dR - s
	sMinusDl := s - dL
	r1 := -math.MaxFloat64
	r2 := -math.MaxFloat64

	if monotonic {
		if !isZero(s) {
			r1 = (dR + dL) / s
		} else if preferShapePreservation {
			r1 = maximumControlParam
		}
	}

	if convex || concave {
		if !(isZero(sMinusDl) || isZero(dRMinusS)) {
			r2 = math.Max(math.Abs(dRMinusDl/dRMinusS), math.Abs(dRMinusDl/sMinusDl))
		} else if preferShapePreservation {
			r2 = maximumControlParam
		}
	} else if monotonic && preferShapePreservation {
		r2 = maximumControlParam
	}

	return math.Max(minimumControlParam, math.Max(r1, r2))
}

// convexControlParameterToFitSecondDerivativeAtLeftSide fits the second derivative while respecting the minimum r for shape preservation.
func convexControlParameterToFitSecondDerivativeAtLeftSide(xL, xR, yL, yR, dL, dR, secondDerivativeL float64, preferShapePreservation bool) float64 {
	r := controlParameterToFitSecondDerivativeAtLeftSide(xL, xR, yL, yR, dL, dR, secondDerivativeL)
	rMin := minimumControlParameter(dL, dR, (yR-yL)/(xR-xL), preferShapePreservation)
	return math.Max(r, rMin)
}

// convexControlParameterToFitSecondDerivativeAtRightSide fits the second derivative while respecting the minimum r for shape preservation.
func convexControlParameterToFitSecondDerivativeAtRightSide(xL, xR, yL, yR, dL, dR, secondDerivativeR float64, preferShapePreservation bool) float64 {
	r := controlParameterToFitSecondDerivativeAtRightSide(xL, xR, yL, yR, dL, dR, secondDerivativeR)
	rMin := minimumControlParameter(dL, dR, (yR-yL)/(xR-xL), preferShapePreservation)
	return math.Max(r, rMin)
}

func isZero(x float64) bool {
	return math.Abs(x) < math.SmallestNonzeroFloat64
}
