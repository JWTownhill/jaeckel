package jaeckel

import "math"

// Region boundary tuning parameters (Jaeckel 2024 values).
const (
	eta = -13.0
	tau = 0.21022410381342863 // 2 * math.Pow(epsilon, 1.0/16.0)
)

func normalizedBlack(thetaX, s float64) float64 {
	if s <= 0 {
		return 0
	}

	// Exact ATM case
	if thetaX == 0 {
		return 2.0*NormCDF(0.5*s) - 1.0
	}

	// Region I (Deep Tails)
	if isRegionI(thetaX, s) {
		return asymptoticExpansionOfScaledNormalizedBlack(thetaX/s, 0.5*s) * normalizedVega(thetaX, s)
	}

	// Region II (Small s)
	if isRegionII(thetaX, s) {
		return smallTExpansionOfScaledNormalizedBlack(thetaX/s, 0.5*s) * normalizedVega(thetaX, s)
	}

	// Regions III/IV (standard computation)
	return normalizedBlackStandard(thetaX, s)
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

func normalizedBlackStandard(thetaX, s float64) float64 {
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

// isRegionI checks if we should use the asymptotic expansion.
// Region I: h < η and t < (τ+½) + (|h|-|η|)
// where h = θx/s and t = s/2.
// Rewritten to avoid division by s:
//
//	s·(s/2-(τ+½+η)) + θx < 0
func isRegionI(thetaX, s float64) bool {
	return thetaX < s*eta && s*(0.5*s-(tau+0.5+eta))+thetaX < 0
}

// isRegionII checks if we should use the small-t expansion.
// Region II: t < τ + (½/|η|)·|h|
// where h = θx/s and t = s/2.
// Rewritten to avoid division by s:
//
//	s·(s-2·τ) - θx/η < 0
func isRegionII(thetaX, s float64) bool {
	return s*(s-2*tau)-thetaX/eta < 0
}

// Black computes the undiscounted Black option price.
//
// Parameters:
//   - F: Forward price
//   - K: Strike price
//   - sigma: Volatility (annualized)
//   - T: Time to expiry (in years)
//   - q: Option type indicator: +1 for call, -1 for put
//
// Returns the undiscounted option price.
func Black(F, K, sigma, T, q float64) float64 {
	s := sigma * math.Sqrt(T)

	// Specialisation for ATM: b(s) = 1 - 2·Φ(-s/2) = erf(s/√8)
	if K == F {
		return F * erf(s*0.5*invSqrtTwo)
	}

	// Map ITM to OTM using put-call parity
	x := math.Log(F / K)
	intrinsic := 0.0
	if q < 0 {
		if K > F {
			intrinsic = K - F
		}
	} else {
		if F > K {
			intrinsic = F - K
		}
	}

	if s <= 0 {
		return intrinsic
	}

	// Compute using geometric mean √(F·K) for numerical stability
	return intrinsic + math.Sqrt(F)*math.Sqrt(K)*normalizedBlack(-math.Abs(x), s)
}

// NormalizedBlack computes the normalized Black option price.
//
// The normalized price is defined as:
//
//	b(x,s,θ) = θ·[exp(θx/2)·Φ(θ·(x/s+s/2)) - exp(-θx/2)·Φ(θ·(x/s-s/2))]
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalized volatility
//   - q: Option type indicator: +1 for call, -1 for put
//
// Returns the normalized option price.
//
// To convert to actual price: price = √(F·K) · NormalizedBlack(log(F/K), σ√T, q)
func NormalizedBlack(x, s, q float64) float64 {
	// Specialisation for ATM: b(s) = erf(s/√8)
	if x == 0 {
		return erf(s * 0.5 * invSqrtTwo)
	}

	thetaX := x
	if q < 0 {
		thetaX = -x
	}

	return normalizedIntrinsic(thetaX) + normalizedBlack(-math.Abs(x), s)
}

// Vega computes the sensitivity of option price to volatility changes.
//
// Parameters:
//   - F: Forward price
//   - K: Strike price
//   - sigma: Volatility (annualized)
//   - T: Time to expiry (in years)
//
// Returns ∂Price/∂σ.
func Vega(F, K, sigma, T float64) float64 {
	x := math.Log(F / K)
	s := sigma * math.Sqrt(T)
	return math.Sqrt(F) * math.Sqrt(K) * NormalizedVega(x, s) * math.Sqrt(T)
}

// NormalizedVega computes the normalized vega.
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalized volatility
//
// Returns ∂b/∂s where b is the normalized Black price.
func NormalizedVega(x, s float64) float64 {
	ax := math.Abs(x)
	if ax <= 0 {
		return invSqrtTwoPi * math.Exp(-0.125*s*s)
	}
	if s <= 0 || s <= ax*math.Sqrt(math.SmallestNonzeroFloat64) {
		return 0
	}
	return normalizedVega(x, s)
}

// Volga computes the second derivative of option price with respect to volatility.
//
// Parameters:
//   - F: Forward price
//   - K: Strike price
//   - sigma: Volatility (annualized)
//   - T: Time to expiry (in years)
//
// Returns ∂²Price/∂σ².
func Volga(F, K, sigma, T float64) float64 {
	x := math.Log(F / K)
	s := sigma * math.Sqrt(T)
	return math.Sqrt(F) * math.Sqrt(K) * NormalizedVolga(x, s) * T
}

// NormalizedVolga computes the normalized volga (∂²b/∂s²).
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalized volatility
//
// Returns ∂²b/∂s² where b is the normalized Black price.
func NormalizedVolga(x, s float64) float64 {
	ax := math.Abs(x)
	if s <= 0 || (ax > 0 && s <= ax*math.Sqrt(math.SmallestNonzeroFloat64)) {
		return 0
	}
	h := x / s
	t := 0.5 * s
	h2 := h * h
	t2 := t * t
	return invSqrtTwoPi * math.Exp(-0.5*(h2+t2)) * (h2 - t2) / s
}
