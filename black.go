package letsberational

import "math"

const (
	// eta and tau are the tuning parameters for the boundary of Region I and II.
	eta = -10.0
	tau = 1.15
)

func normalizedBlack(thetaX, s float64) float64 {
	if s <= 0 {
		return 0
	}

	// Exact ATM case
	if thetaX == 0 {
		return 2.0*normCDF(0.5*s) - 1.0
	}

	// Region I (Deep Tails)
	if isRegionI(thetaX, s) {
		return asymptoticExpansionOfScaledNormalisedBlack(thetaX/s, 0.5*s) * normalizedVega(thetaX, s)
	}

	// Region II (Small s)
	if isRegionII(thetaX, s) {
		return smallTExpansionOfScaledNormalizedBlack(thetaX/s, 0.5*s) * normalizedVega(thetaX, s)
	}

	// 6. Default (Region IV)
	return normalisedBlackWithOptimalUseOfCodysFunctions(thetaX, s)
}

// normalizedVega returns the derivative of the price with respect to s (volatility).
func normalizedVega(x, s float64) float64 {
	ax := math.Abs(x)

	// Case 1: x is zero (At-The-Money).
	if ax <= 0 {
		return invSqrtTwoPi * math.Exp(-0.125*s*s)
	}

	// Case 2: Safeguard against s being too small or zero.
	if s <= 0 || s <= ax*math.Sqrt(math.SmallestNonzeroFloat64) {
		return 0
	}

	// Case 3: Standard
	h := x / s
	t := 0.5 * s
	return invSqrtTwoPi * math.Exp(-0.5*(h*h+t*t))
}

// smallTExpansionOfScaledNormalizedBlack calculates (normalizedBlack / normalizedVega)
// for small values of t (s/2).
func smallTExpansionOfScaledNormalizedBlack(h, t float64) float64 {
	yH := sqrtPiOverTwo * erfcx(-h*invSqrtTwo)
	a := 1.0 + h*yH

	h2 := h * h
	t2 := t * t

	b0 := 2.0 * a
	b1 := (-1.0 + a*(3.0+h2)) / 3.0
	b2 := (-7.0 - h2 + a*(15.0+h2*(10.0+h2))) / 60.0
	b3 := (-57.0 + (-18.0-h2)*h2 + a*(105.0+h2*(105.0+h2*(21.0+h2)))) / 2520.0
	b4 := (-561.0 + h2*(-285.0+(-33.0-h2)*h2) + a*(945.0+h2*(1260.0+h2*(378.0+h2*(36.0+h2))))) / 181440.0
	b5 := (-6555.0 + h2*(-4680.0+h2*(-840.0+(-52.0-h2)*h2)) + a*(10395.0+h2*(17325.0+h2*(6930.0+h2*(990.0+h2*(55.0+h2)))))) / 19958400.0
	b6 := (-89055.0 + h2*(-82845.0+h2*(-20370.0+h2*(-1926.0+(-75.0-h2)*h2))) + a*(135135.0+h2*(270270.0+h2*(135135.0+h2*(25740.0+h2*(2145.0+h2*(78.0+h2))))))) / 3113510400.0

	return t * (b0 + t2*(b1+t2*(b2+t2*(b3+t2*(b4+t2*(b5+b6*t2))))))
}

func normalisedBlackWithOptimalUseOfCodysFunctions(thetaX, s float64) float64 {
	h := thetaX / s
	t := 0.5 * s

	q1 := -invSqrtTwo * (h + t)
	q2 := -invSqrtTwo * (h - t)

	expPos := math.Exp(0.5 * thetaX)
	expNeg := math.Exp(-0.5 * thetaX)
	expSq := math.Exp(-0.5 * (h*h + t*t))

	var twoB float64
	if q1 < codyThreshold {
		if q2 < codyThreshold {
			twoB = expPos*erfc(q1) - expNeg*erfc(q2)
		} else {
			twoB = expPos*erfc(q1) - expSq*erfcx(q2)
		}
	} else {
		if q2 < codyThreshold {
			twoB = expSq*erfcx(q1) - expNeg*erfc(q2)
		} else {
			twoB = expSq * (erfcx(q1) - erfcx(q2))
		}
	}

	return math.Max(0.5*twoB, 0.0)
}

func isRegionI(thetaX, s float64) bool {
	return thetaX < s*eta && s*(0.5*s-(tau+0.5+eta))+thetaX < 0
}

func isRegionII(thetaX, s float64) bool {
	return s*(s-(2*tau))-thetaX/eta < 0
}
