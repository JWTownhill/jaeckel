package letsberational

import (
	"math"
)

const (
	uMax           = 0.3413447460685429
	cdfAsympLimit1 = -10.0
	
)

var cdfAsympLimit2 = -1.0 / math.Sqrt(epsilon)

// Polynomial coefficients grouped by Branch (Rational Approximations)
// Each array is ordered [c_high ... c_low] for Horner's Method.
var (
	// Branch I: r < 2.05
	lowProbA1 = [6]float64{-13.054072340494093704, -83.383894003636969722, -74.594687726045926821, 65.451292110261454609, 47.170590600740689449, 3.691562302945566191}
	lowProbB1 = [6]float64{0.00018295174852053530579, 9.2216887978737432303, 59.270122556046077717, 71.813812182579255459, 20.837211328697753726, 1.0}

	// Branch II: 2.05 <= r < 3.41
	lowProbA2 = [6]float64{-1.2013147879435525574, -10.05916339568646151, -18.1254427791789183, 0.68397370256591532878, 14.49177828689122096, 3.2340179116317970288}
	lowProbB2 = [6]float64{0.000010957576098829595323, 0.84884892199149255469, 7.1369811056109768745, 14.656370665176799712, 8.8820931773304337525, 1.0}

	// Branch III: 3.41 <= r < 6.7
	lowProbA3 = [6]float64{-0.15414319494013597492, -2.8699061335882526744, -11.070534689309368061, -5.1633929115525534628, 9.9483724317036560676, 3.1252235780087584807}
	lowProbB3 = [6]float64{1.3565983564441297634e-7, 0.10897972234131828901, 2.0307076064309043613, 8.1086341122361532407, 7.076769154309171622, 1.0}

	// Branch IV: 6.7 <= r < 12.9
	lowProbA4 = [6]float64{-0.01612303318390145052, -0.47595169546783216436, -2.9644251353150605663, -3.688196041019692267, 2.250881388987032271, 2.6161264950897283681}
	lowProbB4 = [6]float64{3.0848093570966787291e-9, 0.011400087282177594359, 0.33663746405626400164, 2.1282030272153188194, 3.2517455169035921495, 1.0}

	// Branch V: r >= 12.9
	lowProbA5 = [6]float64{-0.0010566357727202585402, -0.065127593753781672404, -0.86385181219213758847, -2.5894451568465728432, -0.042799650734502094297, 2.3226849047872302955}
	lowProbB5 = [6]float64{2.3135343206304887818e-11, 0.0007471447992167225483, 0.046054974512474443189, 0.61320841329197493341, 1.9361316119254412206, 1.0}

	// Midrange coefficients
	midA = [7]float64{-7.58939881401259242, 134.233243502653864, 690.489242061408612, 749.97781456657924, 301.870541922933937, 50.260572167303103, 2.92958954698308805}
	midB = [6]float64{179.227008508102628, 479.123914509756757, 386.821208540417453, 129.404120448755281, 18.918538074574598, 1.0}
)

// normPDF returns the standard normal probability density function.
func normPDF(z float64) float64 {
	return invSqrtTwoPi * math.Exp(-0.5*z*z)
}

// normCDF implements Peter Jäckel's Normal Cumulative Distribution Function.
// It uses an asymptotic expansion for very negative values to maintain relative accuracy.
func normCDF(z float64) float64 {
	if z <= cdfAsympLimit1 {
		sum := 1.0
		if z >= cdfAsympLimit2 {
			zsqr := z * z
			i := 1.0
			g := 1.0
			a := math.MaxFloat64
			var lasta float64
			for {
				lasta = a
				x := (4*i - 3) / zsqr
				y := x * ((4*i - 1) / zsqr)
				a = g * (x - y)
				sum -= a
				g *= y
				i++
				a = math.Abs(a)
				if !(lasta > a && a >= math.Abs(sum*epsilon)) {
					break
				}
			}
		}
		return -normPDF(z) * sum / z
	}
	return 0.5 * math.Erfc(-z*invSqrtTwo)
}

// inverseNormCDF is the "PJ-2024-Inverse-Normal" algorithm.
func inverseNormCDF(p float64) float64 {
	u := p - 0.5
	if math.Abs(u) < uMax {
		return inverseNormCDFMidrange(u)
	}
	if u > 0 {
		return -inverseNormCDFLowProb(1.0 - p)
	}
	return inverseNormCDFLowProb(p)
}

// inverseNormCDFLowProb handles the tail regions.
func inverseNormCDFLowProb(p float64) float64 {
	r := math.Sqrt(-math.Log(p))

	var a, b []float64
	switch {
	case r < 2.05:
		a, b = lowProbA1[:], lowProbB1[:]
	case r < 3.41:
		a, b = lowProbA2[:], lowProbB2[:]
	case r < 6.7:
		a, b = lowProbA3[:], lowProbB3[:]
	case r < 12.9:
		a, b = lowProbA4[:], lowProbB4[:]
	default:
		a, b = lowProbA5[:], lowProbB5[:]
	}

	return evalRational(r, a, b)
}

// inverseNormCDFMidrange handles the center of the distribution.
func inverseNormCDFMidrange(u float64) float64 {
	return u * evalRational(uMax*uMax-u*u, midA[:], midB[:])
}
