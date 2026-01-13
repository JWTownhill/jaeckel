package jaeckel

import (
	"math"
)

const (
	// uMax is the boundary for midrange vs tail in InverseNormCDF.
	uMax = 0.3413447460685429

	// cdfAsympLimit1 is the threshold below which NormCDF uses asymptotic expansion.
	cdfAsympLimit1 = -10.0

	// cdfAsympLimit2 = -1/√ε where ε = 2^-52
	cdfAsympLimit2 = -67108864.0
)

// NormPDF returns the standard normal probability density function.
// This is φ(z) = exp(-z²/2) / √(2π).
func NormPDF(z float64) float64 {
	return invSqrtTwoPi * math.Exp(-0.5*z*z)
}

// NormCDF returns the standard normal cumulative distribution function.
// This is Φ(z) = P(Z ≤ z) where Z ~ N(0,1).
func NormCDF(z float64) float64 {
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
		return -NormPDF(z) * sum / z
	}
	return 0.5 * Erfc(-z*invSqrtTwo)
}

// InverseNormCDF returns the inverse of the standard normal CDF.
// This is the quantile function Φ⁻¹(p) such that Φ(Φ⁻¹(p)) = p.
func InverseNormCDF(p float64) float64 {
	u := p - 0.5
	if math.Abs(u) < uMax {
		return inverseNormCDFMidrange(u)
	}
	if u > 0 {
		return -inverseNormCDFLowProb(1.0 - p)
	}
	return inverseNormCDFLowProb(p)
}

// inverseNormCDFLowProb handles the tail regions with hardcoded Horner-form polynomials.
func inverseNormCDFLowProb(p float64) float64 {
	r := math.Sqrt(-math.Log(p))

	switch {
	case r < 2.05:
		// Branch I: r < 2.05
		return (((((-13.054072340494093704*r-83.383894003636969722)*r-74.594687726045926821)*r+65.451292110261454609)*r+47.170590600740689449)*r + 3.691562302945566191) /
			(((((0.00018295174852053530579*r+9.2216887978737432303)*r+59.270122556046077717)*r+71.813812182579255459)*r+20.837211328697753726)*r + 1.0)
	case r < 3.41:
		// Branch II: 2.05 <= r < 3.41
		return (((((-1.2013147879435525574*r-10.05916339568646151)*r-18.1254427791789183)*r+0.68397370256591532878)*r+14.49177828689122096)*r + 3.2340179116317970288) /
			(((((0.000010957576098829595323*r+0.84884892199149255469)*r+7.1369811056109768745)*r+14.656370665176799712)*r+8.8820931773304337525)*r + 1.0)
	case r < 6.7:
		// Branch III: 3.41 <= r < 6.7
		return (((((-0.15414319494013597492*r-2.8699061335882526744)*r-11.070534689309368061)*r-5.1633929115525534628)*r+9.9483724317036560676)*r + 3.1252235780087584807) /
			(((((1.3565983564441297634e-7*r+0.10897972234131828901)*r+2.0307076064309043613)*r+8.1086341122361532407)*r+7.076769154309171622)*r + 1.0)
	case r < 12.9:
		// Branch IV: 6.7 <= r < 12.9
		return (((((-0.01612303318390145052*r-0.47595169546783216436)*r-2.9644251353150605663)*r-3.688196041019692267)*r+2.250881388987032271)*r + 2.6161264950897283681) /
			(((((3.0848093570966787291e-9*r+0.011400087282177594359)*r+0.33663746405626400164)*r+2.1282030272153188194)*r+3.2517455169035921495)*r + 1.0)
	default:
		// Branch V: r >= 12.9
		return (((((-0.0010566357727202585402*r-0.065127593753781672404)*r-0.86385181219213758847)*r-2.5894451568465728432)*r-0.042799650734502094297)*r + 2.3226849047872302955) /
			(((((2.3135343206304887818e-11*r+0.0007471447992167225483)*r+0.046054974512474443189)*r+0.61320841329197493341)*r+1.9361316119254412206)*r + 1.0)
	}
}

func inverseNormCDFMidrange(u float64) float64 {
	r := uMax*uMax - u*u
	num := ((((((-7.58939881401259242*r+134.233243502653864)*r+690.489242061408612)*r+749.97781456657924)*r+301.870541922933937)*r+50.260572167303103)*r + 2.92958954698308805)
	den := ((((179.227008508102628*r+479.123914509756757)*r+386.821208540417453)*r+129.404120448755281)*r+18.918538074574598)*r + 1.0
	return u * num / den
}
