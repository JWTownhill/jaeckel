package letsberational

import (
	"math"
)

const (
	xNeg = -26.6287357137514
	xBig = 26.543
)

var (
	// Constants from Cody's rational approximations for erfcx.
	// These are ordered so that index 0 is the highest power coefficient for Horner's method.
	numAB = [5]float64{0.185777706184603153, 3.1611237438705656, 113.864154151050156, 377.485237685302021, 3209.37758913846947}
	denAB = [5]float64{1.0, 23.6012909523441209, 244.024637934444173, 1282.61652607737228, 2844.23683343917062}

	numCD = [9]float64{2.15311535474403846e-8, .564188496988670089, 8.88314979438837594, 66.1191906371416295, 298.635138197400131, 881.95222124176909, 1712.04761263407058, 2051.07837782607147, 1230.33935479799725}
	denCD = [9]float64{1.0, 15.7449261107098347, 117.693950891312499, 537.181101862009858, 1621.38957456669019, 3290.79923573345963, 4362.61909014324716, 3439.36767414372164, 1230.33935480374942}

	numPQ = [6]float64{0.0163153871373020978, 0.305326634961232344, 0.360344899949804439, 0.125781726111229246, 0.0160837851487422766, 6.58749161529837803e-4}
	denPQ = [6]float64{1.0, 2.56852019228982242, 1.87295284992346047, 0.527905102951428412, 0.0605183413124413191, 0.00233520497626869185}

	// coefficients for oneMinusErfcx near zero
	numOneMinus = []float64{1.4069285713634565e-2, 1.4069188744609651e-1, 5.7689001208873741e-1, 1.1514967181784756, 1.0000000000000002}
	denOneMinus = []float64{1.2463320728346347e-2, 1.358008134514386e-1, 6.2486081658640257e-1, 1.5089908593742723, 1.9037494962421563, 1.0}
)

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
	ratio := evalRational(x, numOneMinus[:], denOneMinus[:])
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

func evalAB(z float64) float64 {
	return evalRational(z, numAB[:], denAB[:])
}

func evalCD(y float64) float64 {
	return evalRational(y, numCD[:], denCD[:])
}

func evalPQ(z float64) float64 {
	return z * evalRational(z, numPQ[:], denPQ[:])
}
