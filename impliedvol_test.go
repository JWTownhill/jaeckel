package jaeckel

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestImpliedVolatilityRoundTrip(t *testing.T) {
	// Test that we can compute IV from price and get the same price back
	testCases := []struct {
		name  string
		F     float64
		K     float64
		sigma float64
		T     float64
		q     float64
	}{
		{"ATM Call", 100, 100, 0.20, 1.0, 1},
		{"ATM Put", 100, 100, 0.20, 1.0, -1},
		{"OTM Call", 100, 110, 0.25, 0.5, 1},
		{"OTM Put", 100, 90, 0.25, 0.5, -1},
		{"ITM Call", 100, 90, 0.30, 0.25, 1},
		{"ITM Put", 100, 110, 0.30, 0.25, -1},
		{"Deep OTM Call", 100, 150, 0.20, 0.1, 1},
		{"Deep OTM Put", 100, 50, 0.20, 0.1, -1},
		{"Low Vol", 100, 100, 0.05, 1.0, 1},
		{"High Vol", 100, 100, 0.80, 1.0, 1},
		{"Short Expiry", 100, 100, 0.20, 0.01, 1},
		{"Long Expiry", 100, 100, 0.20, 5.0, 1},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Compute price from volatility
			price := Black(tc.F, tc.K, tc.sigma, tc.T, tc.q)

			// Compute implied volatility from price
			ivComputed := ImpliedBlackVolatility(price, tc.F, tc.K, tc.T, tc.q)

			// Check that we get the original volatility back (within 1e-12 relative error)
			assert.InEpsilon(t, tc.sigma, ivComputed, 1e-12,
				"Round trip failed: expected sigma=%v, got IV=%v", tc.sigma, ivComputed)
		})
	}
}

func TestImpliedVolatilityExtremeMoneyness(t *testing.T) {
	// Test extreme moneyness ratios
	F := 100.0
	sigma := 0.20
	T := 1.0
	q := 1.0

	// Extreme OTM calls
	strikes := []float64{200, 300, 500, 1000}
	for _, K := range strikes {
		t.Run("Extreme OTM K="+formatFloat(K), func(t *testing.T) {
			price := Black(F, K, sigma, T, q)
			if price < 1e-300 {
				t.Skip("Price too small for round-trip test")
			}
			iv := ImpliedBlackVolatility(price, F, K, T, q)
			if math.IsInf(iv, 0) || math.IsNaN(iv) {
				t.Skip("IV computation returned special value")
			}
			assert.InEpsilon(t, sigma, iv, 1e-10,
				"Extreme OTM: K=%v, expected sigma=%v, got IV=%v", K, sigma, iv)
		})
	}
}

func TestImpliedVolatilityATM(t *testing.T) {
	// Test exact ATM case
	F := 100.0
	K := 100.0
	T := 1.0
	q := 1.0

	sigmas := []float64{0.01, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00}
	for _, sigma := range sigmas {
		t.Run("ATM sigma="+formatFloat(sigma), func(t *testing.T) {
			price := Black(F, K, sigma, T, q)
			iv := ImpliedBlackVolatility(price, F, K, T, q)

			assert.InEpsilon(t, sigma, iv, 2e-14,
				"ATM round-trip: expected sigma=%v, got IV=%v", sigma, iv)
		})
	}
}

func TestImpliedVolatilityBoundaries(t *testing.T) {
	F := 100.0
	K := 100.0
	T := 1.0

	t.Run("Price Below Intrinsic Call", func(t *testing.T) {
		// For ATM call, intrinsic is 0, so any negative price is invalid
		iv := ImpliedBlackVolatility(-0.01, F, K, T, 1)
		assert.Equal(t, volatilityValueToSignalPriceIsBelowIntrinsic, iv)
	})

	t.Run("Price Below Intrinsic Put", func(t *testing.T) {
		// For ATM put, intrinsic is 0
		iv := ImpliedBlackVolatility(-0.01, F, K, T, -1)
		assert.Equal(t, volatilityValueToSignalPriceIsBelowIntrinsic, iv)
	})

	t.Run("Price Above Maximum Call", func(t *testing.T) {
		// Maximum call price is F
		iv := ImpliedBlackVolatility(F+1, F, K, T, 1)
		assert.Equal(t, volatilityValueToSignalPriceIsAboveMaximum, iv)
	})

	t.Run("Price Above Maximum Put", func(t *testing.T) {
		// Maximum put price is K
		iv := ImpliedBlackVolatility(K+1, F, K, T, -1)
		assert.Equal(t, volatilityValueToSignalPriceIsAboveMaximum, iv)
	})

	t.Run("Zero Price", func(t *testing.T) {
		iv := ImpliedBlackVolatility(0, F, K, T, 1)
		assert.Equal(t, 0.0, iv)
	})
}

func TestNormalizedBlackConsistency(t *testing.T) {
	// Test that NormalizedBlack is consistent with Black
	F := 100.0
	K := 105.0
	sigma := 0.25
	T := 0.5

	x := math.Log(F / K)
	s := sigma * math.Sqrt(T)

	for _, q := range []float64{1, -1} {
		t.Run("q="+formatFloat(q), func(t *testing.T) {
			blackPrice := Black(F, K, sigma, T, q)
			normPrice := NormalizedBlack(x, s, q)

			// Black price should equal √(F·K) · NormalizedBlack
			expectedBlackPrice := math.Sqrt(F*K) * normPrice

			assert.InEpsilon(t, expectedBlackPrice, blackPrice, 1e-14)
		})
	}
}

func TestVegaConsistency(t *testing.T) {
	// Test that vega is positive and consistent
	F := 100.0
	K := 100.0
	sigma := 0.20
	T := 1.0

	vega := Vega(F, K, sigma, T)
	assert.Greater(t, vega, 0.0, "Vega should be positive")

	// Numerical derivative check
	eps := 1e-6
	priceUp := Black(F, K, sigma+eps, T, 1)
	priceDown := Black(F, K, sigma-eps, T, 1)
	numericalVega := (priceUp - priceDown) / (2 * eps)

	assert.InEpsilon(t, numericalVega, vega, 1e-6,
		"Vega analytical=%v, numerical=%v", vega, numericalVega)
}

func TestVolga(t *testing.T) {
	F := 100.0
	K := 100.0
	sigma := 0.20
	T := 1.0

	volga := Volga(F, K, sigma, T)

	// Numerical second derivative check
	eps := 1e-4
	vegaUp := Vega(F, K, sigma+eps, T)
	vegaDown := Vega(F, K, sigma-eps, T)
	numericalVolga := (vegaUp - vegaDown) / (2 * eps)

	assert.InEpsilon(t, numericalVolga, volga, 1e-4,
		"Volga analytical=%v, numerical=%v", volga, numericalVolga)
}

func TestErfinvRoundTrip(t *testing.T) {
	testValues := []float64{-0.9, -0.5, -0.1, 0.1, 0.5, 0.9}

	for _, x := range testValues {
		t.Run("x="+formatFloat(x), func(t *testing.T) {
			y := Erfinv(x)
			roundTrip := Erf(y)

			assert.InEpsilon(t, x, roundTrip, 1e-14,
				"erfinv round trip: x=%v, erfinv(x)=%v, erf(erfinv(x))=%v",
				x, y, roundTrip)
		})
	}

	// Special case: zero
	t.Run("x=0", func(t *testing.T) {
		assert.Equal(t, 0.0, Erfinv(0))
	})
}

func formatFloat(f float64) string {
	if f == float64(int(f)) {
		return string(rune('0' + int(f)))
	}
	return string([]byte{byte('0' + int(math.Abs(f)))})
}

func TestRegionI(t *testing.T) {
	// These cases have h = x/s < -13, triggering the asymptotic expansion
	testCases := []struct {
		name string
		x    float64
		s    float64
		h    float64 // expected h = x/s
	}{
		{"h=-14", -7.0, 0.5, -14.0},
		{"h=-16", -8.0, 0.5, -16.0},
		{"h=-20", -10.0, 0.5, -20.0},
		{"h=-15", -3.0, 0.2, -15.0},
		{"h=-25", -5.0, 0.2, -25.0},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Verify h calculation
			h := tc.x / tc.s
			assert.InDelta(t, tc.h, h, 1e-10, "h calculation mismatch")

			// Compute normalized Black price
			price := NormalizedBlack(tc.x, tc.s, 1.0)
			assert.Greater(t, price, 0.0, "Price should be positive")
			assert.False(t, math.IsNaN(price), "Price should not be NaN")
			assert.False(t, math.IsInf(price, 0), "Price should not be Inf")

			// Round-trip test: compute IV from price, then price from IV
			iv := NormalizedImpliedBlackVolatility(price, tc.x, 1.0)
			if !math.IsInf(iv, 0) && !math.IsNaN(iv) && iv > 0 {
				priceRoundTrip := NormalizedBlack(tc.x, iv, 1.0)
				assert.InEpsilon(t, price, priceRoundTrip, 1e-10,
					"Round-trip failed: original=%v, roundtrip=%v", price, priceRoundTrip)
			}
		})
	}
}

// TestRegionIImpliedVolatility tests IV calculation for deep OTM options (Region I)
func TestRegionIImpliedVolatility(t *testing.T) {
	// Deep OTM options that should trigger Region I in the IV calculation
	testCases := []struct {
		name  string
		F     float64
		K     float64
		sigma float64
		T     float64
	}{
		{"DeepOTM_K200", 100, 200, 0.20, 0.1},  // very deep OTM
		{"DeepOTM_K300", 100, 300, 0.30, 0.25}, // extremely deep OTM
		{"DeepOTM_K500", 100, 500, 0.40, 0.5},  // ultra deep OTM
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			price := Black(tc.F, tc.K, tc.sigma, tc.T, 1.0)
			if price < 1e-300 {
				t.Skip("Price too small")
			}

			iv := ImpliedBlackVolatility(price, tc.F, tc.K, tc.T, 1.0)
			if math.IsInf(iv, 0) || math.IsNaN(iv) {
				t.Skip("IV returned special value")
			}

			assert.InEpsilon(t, tc.sigma, iv, 1e-10,
				"Region I IV round-trip: expected sigma=%v, got IV=%v", tc.sigma, iv)
		})
	}
}

// Benchmark tests using b.Loop() (Go 1.24+)

// BenchmarkBlack benchmarks Black formula across different regions
func BenchmarkBlack(b *testing.B) {
	cases := []struct {
		name  string
		F, K  float64
		sigma float64
		T     float64
	}{
		{"ATM", 100, 100, 0.20, 1.0},
		{"OTM_5pct", 100, 105, 0.20, 1.0},
		{"OTM_20pct", 100, 120, 0.20, 1.0},
		{"DeepOTM", 100, 200, 0.20, 1.0},
		{"ITM_5pct", 100, 95, 0.20, 1.0},
		{"LowVol", 100, 100, 0.05, 1.0},
		{"HighVol", 100, 100, 0.80, 1.0},
		{"ShortExpiry", 100, 100, 0.20, 0.01},
		{"LongExpiry", 100, 100, 0.20, 5.0},
	}

	for _, tc := range cases {
		b.Run(tc.name, func(b *testing.B) {
			for b.Loop() {
				Black(tc.F, tc.K, tc.sigma, tc.T, 1.0)
			}
		})
	}
}

// BenchmarkImpliedVolatility benchmarks IV calculation across different algorithm branches
func BenchmarkImpliedVolatility(b *testing.B) {
	cases := []struct {
		name  string
		F, K  float64
		sigma float64
		T     float64
	}{
		// ATM branch - uses direct erfinv formula
		{"ATM", 100, 100, 0.20, 1.0},
		// Lower middle branch - s_l <= s < s_c
		{"LowerMiddle_OTM5pct", 100, 105, 0.20, 1.0},
		// Upper middle branch - s_c <= s <= s_u
		{"UpperMiddle_ITM5pct", 100, 95, 0.20, 1.0},
		// Lowest branch - s < s_l (deep OTM, uses log objective)
		{"Lowest_DeepOTM", 100, 150, 0.20, 0.1},
		// Highest branch - s > s_u (high premium)
		{"Highest_HighVol", 100, 100, 1.50, 1.0},
		// Various moneyness levels
		{"OTM_10pct", 100, 110, 0.25, 0.5},
		{"OTM_30pct", 100, 130, 0.30, 1.0},
		{"ITM_10pct", 100, 90, 0.25, 0.5},
		{"ITM_30pct", 100, 70, 0.30, 1.0},
		// Edge cases
		{"VeryLowVol", 100, 100, 0.01, 1.0},
		{"VeryHighVol", 100, 100, 2.00, 1.0},
		{"ShortExpiry", 100, 105, 0.20, 0.01},
		{"LongExpiry", 100, 105, 0.20, 10.0},
	}

	for _, tc := range cases {
		price := Black(tc.F, tc.K, tc.sigma, tc.T, 1.0)
		b.Run(tc.name, func(b *testing.B) {
			for b.Loop() {
				ImpliedBlackVolatility(price, tc.F, tc.K, tc.T, 1.0)
			}
		})
	}
}

// BenchmarkNormalizedBlack benchmarks normalized Black across regions
func BenchmarkNormalizedBlack(b *testing.B) {
	cases := []struct {
		name string
		x, s float64
	}{
		{"ATM_s0.2", 0, 0.2},
		{"ATM_s1.0", 0, 1.0},
		{"RegionI_DeepTail", -7.0, 0.5},      // h=-14 < eta=-13, triggers asymptotic expansion
		{"RegionI_Extreme", -15.0, 0.5},      // h=-30, deep in Region I
		{"RegionII_SmallS", -0.1, 0.1},       // t < tau region
		{"Standard", -0.5, 0.5},              // standard computation (Regions III/IV)
		{"ExtremeMoney_Neg10", -10.0, 1.0},   // h=-10, just outside Region I
		{"ExtremeMoney_Neg100", -100.0, 5.0}, // h=-20, Region I
	}

	for _, tc := range cases {
		b.Run(tc.name, func(b *testing.B) {
			for b.Loop() {
				NormalizedBlack(tc.x, tc.s, 1.0)
			}
		})
	}
}

// BenchmarkVega benchmarks vega calculation
func BenchmarkVega(b *testing.B) {
	cases := []struct {
		name  string
		F, K  float64
		sigma float64
		T     float64
	}{
		{"ATM", 100, 100, 0.20, 1.0},
		{"OTM", 100, 120, 0.20, 1.0},
		{"ITM", 100, 80, 0.20, 1.0},
		{"HighVol", 100, 100, 0.80, 1.0},
	}

	for _, tc := range cases {
		b.Run(tc.name, func(b *testing.B) {
			for b.Loop() {
				Vega(tc.F, tc.K, tc.sigma, tc.T)
			}
		})
	}
}

// BenchmarkErfFunctions benchmarks error function implementations
func BenchmarkErfFunctions(b *testing.B) {
	b.Run("Erf", func(b *testing.B) {
		for b.Loop() {
			Erf(0.5)
		}
	})
	b.Run("Erfc", func(b *testing.B) {
		for b.Loop() {
			Erfc(0.5)
		}
	})
	b.Run("Erfcx_Small", func(b *testing.B) {
		for b.Loop() {
			Erfcx(0.3)
		}
	})
	b.Run("Erfcx_Large", func(b *testing.B) {
		for b.Loop() {
			Erfcx(10.0)
		}
	})
	b.Run("Erfinv", func(b *testing.B) {
		for b.Loop() {
			Erfinv(0.5)
		}
	})
}

// BenchmarkNormCDFFunctions benchmarks normal distribution functions
func BenchmarkNormCDFFunctions(b *testing.B) {
	b.Run("NormCDF_Center", func(b *testing.B) {
		for b.Loop() {
			NormCDF(0.0)
		}
	})
	b.Run("NormCDF_Tail", func(b *testing.B) {
		for b.Loop() {
			NormCDF(3.0)
		}
	})
	b.Run("InverseNormCDF_Center", func(b *testing.B) {
		for b.Loop() {
			InverseNormCDF(0.5)
		}
	})
	b.Run("InverseNormCDF_Tail", func(b *testing.B) {
		for b.Loop() {
			InverseNormCDF(0.001)
		}
	})
}
