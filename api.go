// Package letsberational provides high-precision Black-Scholes option pricing
// and implied volatility calculation using Peter Jäckel's "Let's Be Rational" algorithm.
//
// This implementation achieves full machine precision (approximately 16 significant digits)
// for implied volatility calculations across all input ranges, including extreme
// moneyness ratios and near-zero volatilities.
//
// # Key Features
//
//   - Implied volatility calculation with machine-precision accuracy
//   - Normalised and standard Black-Scholes pricing
//   - Greeks: Vega and Volga
//   - Numerically stable across extreme parameter ranges
//   - Based on the latest (2024) version of Peter Jäckel's algorithm
//
// # References
//
// Peter Jäckel, "Let's Be Rational", Wilmott Magazine, January 2015.
// Source code: www.jaeckel.org/LetsBeRational.7z
//
// # Example Usage
//
//	// Calculate implied volatility from option price
//	price := 10.0
//	F := 100.0  // Forward price
//	K := 105.0  // Strike price
//	T := 0.5    // Time to expiry in years
//	q := 1.0    // 1 for call, -1 for put
//	iv := letsberational.ImpliedBlackVolatility(price, F, K, T, q)
//
//	// Calculate option price from volatility
//	sigma := 0.20  // 20% volatility
//	optionPrice := letsberational.Black(F, K, sigma, T, q)
package letsberational

import "math"

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
//
// Note: To get the discounted price, multiply by exp(-r*T) where r is the
// risk-free rate.
func Black(F, K, sigma, T, q float64) float64 {
	s := sigma * math.Sqrt(T)

	// Specialisation for ATM: b(s) = 1 - 2·Φ(-s/2) = erf(s/√8)
	if K == F {
		return F * erf(s*0.5*invSqrtTwo)
	}

	// Map in-the-money to out-of-the-money using put-call parity
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

// NormalisedBlack computes the normalised Black option price.
//
// The normalised price is defined as:
//
//	b(x,s,θ) = θ·[exp(θx/2)·Φ(θ·(x/s+s/2)) - exp(-θx/2)·Φ(θ·(x/s-s/2))]
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalised volatility
//   - q: Option type indicator: +1 for call, -1 for put
//
// Returns the normalised option price.
//
// To convert to actual price: price = √(F·K) · NormalisedBlack(log(F/K), σ√T, q)
func NormalisedBlack(x, s, q float64) float64 {
	// Specialisation for ATM: b(s) = erf(s/√8)
	if x == 0 {
		return erf(s * 0.5 * invSqrtTwo)
	}

	thetaX := x
	if q < 0 {
		thetaX = -x
	}

	return normalisedIntrinsic(thetaX) + normalizedBlack(-math.Abs(x), s)
}

// ImpliedBlackVolatility computes the implied Black volatility from an option price.
//
// Parameters:
//   - price: The undiscounted option price
//   - F: Forward price
//   - K: Strike price
//   - T: Time to expiry (in years)
//   - q: Option type indicator: +1 for call, -1 for put
//
// Returns the implied volatility (annualized).
//
// Special return values:
//   - -Inf: price is below intrinsic value (arbitrage)
//   - +Inf: price is above maximum possible value
func ImpliedBlackVolatility(price, F, K, T, q float64) float64 {
	// Check if price is above maximum
	maxPrice := F
	if q < 0 {
		maxPrice = K
	}
	if price >= maxPrice {
		return volatilityValueToSignalPriceIsAboveMaximum
	}

	// Map in-the-money to out-of-the-money
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

	// Normalise the price
	sqrtFK := math.Sqrt(F) * math.Sqrt(K)
	beta := (price - intrinsic) / sqrtFK

	// Call the core algorithm
	normalisedVol := letsBerational(beta, -math.Abs(x), impliedVolatilityMaxIterations)

	// Convert from normalised to actual volatility
	return normalisedVol / math.Sqrt(T)
}

// NormalisedImpliedBlackVolatility computes the normalised implied volatility.
//
// Parameters:
//   - beta: The normalised option price
//   - x: log(F/K), the log-moneyness
//   - q: Option type indicator: +1 for call, -1 for put
//
// Returns the normalised implied volatility s = σ√T.
//
// Special return values:
//   - -Inf: price is below intrinsic value
//   - +Inf: price is above maximum possible value
func NormalisedImpliedBlackVolatility(beta, x, q float64) float64 {
	thetaX := x
	if q < 0 {
		thetaX = -x
	}
	return letsBerational(beta-normalisedIntrinsic(thetaX), -math.Abs(x), impliedVolatilityMaxIterations)
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
	return math.Sqrt(F) * math.Sqrt(K) * NormalisedVega(x, s) * math.Sqrt(T)
}

// NormalisedVega computes the normalised vega.
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalised volatility
//
// Returns ∂b/∂s where b is the normalised Black price.
func NormalisedVega(x, s float64) float64 {
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
	return math.Sqrt(F) * math.Sqrt(K) * NormalisedVolga(x, s) * T
}

// NormalisedVolga computes the normalised volga (∂²b/∂s²).
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalised volatility
//
// Returns ∂²b/∂s² where b is the normalised Black price.
func NormalisedVolga(x, s float64) float64 {
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

// ComplementaryNormalisedBlack computes bmax - b(x,s) without subtractive cancellation.
//
// This is useful when the option is deep in-the-money and the price is close
// to its maximum value bmax = exp(θx/2).
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalised volatility
//
// Returns exp(x/2) - b(x,s,+1) for calls, computed without subtractive cancellation.
func ComplementaryNormalisedBlack(x, s float64) float64 {
	return complementaryNormalisedBlack(x/s, s/2)
}

// NormCDF returns the standard normal cumulative distribution function.
// This is Φ(z) = P(Z ≤ z) where Z ~ N(0,1).
func NormCDF(z float64) float64 {
	return normCDF(z)
}

// NormPDF returns the standard normal probability density function.
// This is φ(z) = exp(-z²/2) / √(2π).
func NormPDF(z float64) float64 {
	return normPDF(z)
}

// InverseNormCDF returns the inverse of the standard normal CDF.
// This is the quantile function Φ⁻¹(p) such that Φ(Φ⁻¹(p)) = p.
func InverseNormCDF(p float64) float64 {
	return inverseNormCDF(p)
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

// Erfinv returns the inverse error function.
// erfinv(x) satisfies erf(erfinv(x)) = x for x ∈ (-1, 1).
func Erfinv(x float64) float64 {
	return erfinv(x)
}
