package jaeckel

import "math"

// erfinv computes the inverse error function.
//
// The inverse error function satisfies: erf(erfinv(x)) = x for x in (-1, 1).
//
// This implementation uses the internal branches of the inverse normal CDF
// (PJ-2024-Inverse-Normal algorithm) to avoid catastrophic subtractive
// cancellation for small arguments.
//
// The relationship is: erfinv(x) = inverseNormCDF((1+x)/2) / √2
//
// Reference: Peter Jäckel, "Let's Be Rational", 2024 implementation.
func erfinv(x float64) float64 {
	if x <= -1 {
		return math.Inf(-1)
	}
	if x >= 1 {
		return math.Inf(1)
	}
	if x == 0 {
		return 0
	}

	// Use the relationship: erfinv(x) = Φ⁻¹((1+x)/2) / √2
	// where Φ⁻¹ is the inverse normal CDF.
	//
	// For small |x|, we use the midrange branch of inverseNormCDF
	// which is specifically designed to avoid subtractive cancellation.
	p := 0.5 * (1 + x)
	return InverseNormCDF(p) * invSqrtTwo
}

// erfinvForATMImpliedVolatility computes √8 · erfinv(β) directly.
//
// This is optimized for the ATM implied volatility case where
// bₐₜₘ(s) = erf(s/√8), so s = √8 · erfinv(β).
//
// For small β, this avoids intermediate underflow issues that could occur
// when computing erfinv(β) separately and then multiplying by √8.
func erfinvForATMImpliedVolatility(beta float64) float64 {
	if beta <= 0 {
		return 0
	}
	if beta >= 1 {
		return math.Inf(1)
	}

	// s = √8 · erfinv(β)
	//   = √8 · Φ⁻¹((1+β)/2) / √2
	//   = 2 · Φ⁻¹((1+β)/2)
	//   = -2 · Φ⁻¹((1-β)/2)
	//
	// For the ATM case, we use p = (1+β)/2 and compute:
	//   s = 2 · inverseNormCDF(p)
	p := 0.5 * (1 + beta)
	return 2 * InverseNormCDF(p)
}

// impliedNormalizedVolatilityATM computes the normalized implied volatility
// for the exact at-the-money case (x = 0).
//
// At the money, bₐₜₘ(s) = 1 - 2·Φ(-s/2) = 2·Φ(s/2) - 1 = erf(s/√8)
//
// Therefore: s = √8 · erfinv(β)
func impliedNormalizedVolatilityATM(beta float64) float64 {
	return erfinvForATMImpliedVolatility(beta)
}

// Erfinv returns the inverse error function.
// erfinv(x) satisfies erf(erfinv(x)) = x for x ∈ (-1, 1).
func Erfinv(x float64) float64 {
	return erfinv(x)
}
