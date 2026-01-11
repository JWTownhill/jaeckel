package jaeckel

import (
	"math"
)

const (
	xNeg = -26.6287357137514
	xBig = 26.543
)

// Cody's rational approximations - hardcoded Horner form for performance.
// Coefficients ordered highest to lowest power for Horner's method: ((a*x+b)*x+c)*x+d...

// evalAB evaluates the rational function for |x| <= 0.46875.
func evalAB(z float64) float64 {
	// numAB: [0.185777706184603153, 3.1611237438705656, 113.864154151050156, 377.485237685302021, 3209.37758913846947]
	// denAB: [1.0, 23.6012909523441209, 244.024637934444173, 1282.61652607737228, 2844.23683343917062]
	return ((((0.185777706184603153*z+3.1611237438705656)*z+113.864154151050156)*z+377.485237685302021)*z + 3209.37758913846947) /
		((((1.0*z+23.6012909523441209)*z+244.024637934444173)*z+1282.61652607737228)*z + 2844.23683343917062)
}

// evalCD evaluates the rational function for 0.46875 < |x| <= 4.
func evalCD(y float64) float64 {
	// numCD: [2.15311535474403846e-8, .564188496988670089, 8.88314979438837594, 66.1191906371416295, 298.635138197400131, 881.95222124176909, 1712.04761263407058, 2051.07837782607147, 1230.33935479799725]
	// denCD: [1.0, 15.7449261107098347, 117.693950891312499, 537.181101862009858, 1621.38957456669019, 3290.79923573345963, 4362.61909014324716, 3439.36767414372164, 1230.33935480374942]
	return ((((((((2.15311535474403846e-8*y+0.564188496988670089)*y+8.88314979438837594)*y+66.1191906371416295)*y+298.635138197400131)*y+881.95222124176909)*y+1712.04761263407058)*y+2051.07837782607147)*y + 1230.33935479799725) /
		((((((((1.0*y+15.7449261107098347)*y+117.693950891312499)*y+537.181101862009858)*y+1621.38957456669019)*y+3290.79923573345963)*y+4362.61909014324716)*y+3439.36767414372164)*y + 1230.33935480374942)
}

// evalPQ evaluates the rational function for |x| > 4.
func evalPQ(z float64) float64 {
	// numPQ: [0.0163153871373020978, 0.305326634961232344, 0.360344899949804439, 0.125781726111229246, 0.0160837851487422766, 6.58749161529837803e-4]
	// denPQ: [1.0, 2.56852019228982242, 1.87295284992346047, 0.527905102951428412, 0.0605183413124413191, 0.00233520497626869185]
	return z * (((((0.0163153871373020978*z+0.305326634961232344)*z+0.360344899949804439)*z+0.125781726111229246)*z+0.0160837851487422766)*z + 6.58749161529837803e-4) /
		(((((1.0*z+2.56852019228982242)*z+1.87295284992346047)*z+0.527905102951428412)*z+0.0605183413124413191)*z + 0.00233520497626869185)
}

// evalOneMinusErfcx evaluates the rational function for oneMinusErfcx in [-0.2, 1/3].
func evalOneMinusErfcx(x float64) float64 {
	// numOneMinus: [1.4069285713634565e-2, 1.4069188744609651e-1, 5.7689001208873741e-1, 1.1514967181784756, 1.0000000000000002]
	// denOneMinus: [1.2463320728346347e-2, 1.358008134514386e-1, 6.2486081658640257e-1, 1.5089908593742723, 1.9037494962421563, 1.0]
	return ((((1.4069285713634565e-2*x+1.4069188744609651e-1)*x+5.7689001208873741e-1)*x+1.1514967181784756)*x + 1.0000000000000002) /
		(((((1.2463320728346347e-2*x+1.358008134514386e-1)*x+6.2486081658640257e-1)*x+1.5089908593742723)*x+1.9037494962421563)*x + 1.0)
}

func erf(x float64) float64 {
	y := math.Abs(x)
	if y <= codyThreshold {
		return x * evalAB(y*y)
	}
	var erfcAbsX float64
	if y >= xBig {
		erfcAbsX = 0
	} else {
		if y <= 4 {
			erfcAbsX = evalCD(y)
		} else {
			erfcAbsX = (oneOverSqrtPi - evalPQ(1/(y*y))) / y
		}
		erfcAbsX *= smoothenedExpNegSq(y)
	}
	if x < 0 {
		return erfcAbsX - 1
	}
	return 1 - erfcAbsX
}

func erfc(x float64) float64 {
	y := math.Abs(x)
	if y <= codyThreshold {
		return 1 - x*evalAB(y*y)
	}
	var erfcAbsX float64
	if y >= xBig {
		erfcAbsX = 0
	} else {
		if y <= 4 {
			erfcAbsX = evalCD(y)
		} else {
			erfcAbsX = (oneOverSqrtPi - evalPQ(1/(y*y))) / y
		}
		erfcAbsX *= smoothenedExpNegSq(y)
	}
	if x < 0 {
		return 2 - erfcAbsX
	}
	return erfcAbsX
}

// erfcx computes the scaled complementary error function exp(x^2) * erfc(x).
func erfcx(x float64) float64 {
	y := math.Abs(x)

	if y <= codyThreshold {
		z := y * y
		return math.Exp(z) * (1 - x*evalAB(z))
	}

	if x < xNeg {
		return math.MaxFloat64
	}

	result := erfcxAboveThreshold(y)

	if x < 0 {
		expx2 := smoothenedExpPosSq(x)
		return (expx2 + expx2) - result
	}

	return result
}

// f(x) := 1 - erfcx(x)  ≈  2/√π·x - x² + 4/(3√π)·x³ + ... for small x.
func oneMinusErfcx(x float64) float64 {
	if x < -0.2 || x > 1.0/3.0 {
		return 1.0 - erfcx(x)
	}
	// Remez-optimized minimax rational function of order (4,5) for g(x) := (2/√π-f(x)/x)/x.
	// The relative accuracy of f(x) ≈ x·(2/√π-x·g(x)) is better than 2.5E-17 (in perfect arithmetic) on x in [-0.2,1/3].
	ratio := evalOneMinusErfcx(x)
	return x * (twoOverSqrtPi - x*ratio)
}

func smoothenedExpPosSq(x float64) float64 {
	xTilde := math.Trunc(x*16) / 16.0
	return math.Exp(xTilde*xTilde) * math.Exp((x-xTilde)*(x+xTilde))
}

func smoothenedExpNegSq(y float64) float64 {
	yTilde := math.Trunc(y*16) / 16.0
	return math.Exp(-yTilde*yTilde) * math.Exp(-(y-yTilde)*(y+yTilde))
}

func erfcxAboveThreshold(y float64) float64 {
	if y <= 4.0 {
		return evalCD(y)
	}
	return (oneOverSqrtPi - evalPQ(1/(y*y))) / y
}

// Erf returns the error function.
// erf(x) = (2/√π) ∫₀ˣ exp(-t²) dt
func Erf(x float64) float64 {
	return erf(x)
}

// Erfc returns the complementary error function.
// erfc(x) = 1 - erf(x) = (2/√π) ∫ₓ^∞ exp(-t²) dt
func Erfc(x float64) float64 {
	return erfc(x)
}

// Erfcx returns the scaled complementary error function.
// erfcx(x) = exp(x²) · erfc(x)
//
// This function is numerically stable for large positive x where
// erfc(x) would underflow.
func Erfcx(x float64) float64 {
	return erfcx(x)
}
