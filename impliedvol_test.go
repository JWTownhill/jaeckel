package jaeckel

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

// TestPutRoundTrip verifies put options work correctly in round-trip pricing.
// Internally, puts are converted to calls via put-call parity, then the call
// implied vol is computed. This test ensures that transformation is accurate.
// Call options are comprehensively tested in TestComprehensiveRoundTrip.
func TestPutRoundTrip(t *testing.T) {
	testCases := []struct {
		name  string
		F     float64
		K     float64
		sigma float64
		T     float64
	}{
		{"ATM", 100, 100, 0.20, 1.0},
		{"OTM Put", 100, 90, 0.25, 0.5},
		{"ITM Put", 100, 110, 0.30, 0.25},
		{"Deep OTM Put", 100, 50, 0.20, 0.1},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			price := Black(tc.F, tc.K, tc.sigma, tc.T, Put)
			iv := ImpliedBlackVolatility(price, tc.F, tc.K, tc.T, Put)
			assert.InEpsilon(t, tc.sigma, iv, 1e-12)
		})
	}
}

func TestImpliedVolatilityBoundaries(t *testing.T) {
	F := 100.0
	K := 100.0
	T := 1.0

	t.Run("Price Below Intrinsic Call", func(t *testing.T) {
		iv := ImpliedBlackVolatility(-0.01, F, K, T, Call)
		assert.Equal(t, volatilityValueToSignalPriceIsBelowIntrinsic, iv)
	})

	t.Run("Price Below Intrinsic Put", func(t *testing.T) {
		iv := ImpliedBlackVolatility(-0.01, F, K, T, Put)
		assert.Equal(t, volatilityValueToSignalPriceIsBelowIntrinsic, iv)
	})

	t.Run("Price Above Maximum Call", func(t *testing.T) {
		iv := ImpliedBlackVolatility(F+1, F, K, T, Call)
		assert.Equal(t, volatilityValueToSignalPriceIsAboveMaximum, iv)
	})

	t.Run("Price Above Maximum Put", func(t *testing.T) {
		iv := ImpliedBlackVolatility(K+1, F, K, T, Put)
		assert.Equal(t, volatilityValueToSignalPriceIsAboveMaximum, iv)
	})

	t.Run("Zero Price", func(t *testing.T) {
		iv := ImpliedBlackVolatility(0, F, K, T, Call)
		assert.Equal(t, 0.0, iv)
	})
}

func TestImpliedVolatilityITMParity(t *testing.T) {
	// Test ITM options - they should map to OTM internally
	F := 100.0
	sigma := 0.25
	T := 1.0

	t.Run("ITM Call", func(t *testing.T) {
		K := 90.0
		price := Black(F, K, sigma, T, Call)
		iv := ImpliedBlackVolatility(price, F, K, T, Call)
		assert.InEpsilon(t, sigma, iv, 1e-12)
	})

	t.Run("ITM Put", func(t *testing.T) {
		K := 110.0
		price := Black(F, K, sigma, T, Put)
		iv := ImpliedBlackVolatility(price, F, K, T, Put)
		assert.InEpsilon(t, sigma, iv, 1e-12)
	})
}

func TestNormalizedImpliedVolatilityBoundaries(t *testing.T) {
	t.Run("zero beta returns zero", func(t *testing.T) {
		iv := NormalizedImpliedBlackVolatility(0, -0.5, Call)
		assert.Equal(t, 0.0, iv, "Zero price should return zero IV")
	})

	t.Run("negative beta returns below-intrinsic signal", func(t *testing.T) {
		// For OTM call (x<0), intrinsic=0, so beta goes directly to letsBeRational
		iv := NormalizedImpliedBlackVolatility(-0.1, -0.5, Call)
		assert.Equal(t, volatilityValueToSignalPriceIsBelowIntrinsic, iv, "Negative price should signal below intrinsic")
	})

	t.Run("beta >= bMax returns above-maximum signal", func(t *testing.T) {
		x := -0.5
		// For OTM call, thetaX = x, and letsBeRational gets -|x| = x
		// bMax = exp(0.5 * (-|x|)) = exp(-0.25)
		bMax := math.Exp(0.5 * x)
		iv := NormalizedImpliedBlackVolatility(bMax+0.01, x, Call)
		assert.Equal(t, volatilityValueToSignalPriceIsAboveMaximum, iv, "Price above max should signal above maximum")
	})

	t.Run("Put option parity", func(t *testing.T) {
		x := 0.5 // Put is OTM when x > 0
		s := 0.3
		price := NormalizedBlack(x, s, Put)
		iv := NormalizedImpliedBlackVolatility(price, x, Put)
		assert.InEpsilon(t, s, iv, 1e-10)
	})
}

func TestComprehensiveRoundTrip(t *testing.T) {
	// Log-moneyness values from Jäckel's timing tests
	xValues := []float64{
		0, 1e-14, 1e-8, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5,
		1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
	}

	// For each |x|, test multiple volatility levels
	for _, absX := range xValues {
		x := -absX

		t.Run(formatTestName(absX), func(t *testing.T) {
			sMid := math.Sqrt(math.Max(2*absX, 1e-16))
			sValues := generateSValues(sMid, absX)

			for _, s := range sValues {
				if s <= 0 {
					continue
				}

				t.Run(formatSName(s), func(t *testing.T) {
					price := NormalizedBlack(x, s, Call)

					// Skip degenerate cases where no unique IV exists:
					// - Price underflowed to zero (vol would be 0)
					// - Price at/above maximum bMax = exp(0.5*x) (vol would be +inf)
					bMax := math.Exp(0.5 * x)
					if price <= 0 || price >= bMax {
						t.Skipf("Degenerate case: price=%v, bMax=%v for x=%v, s=%v", price, bMax, x, s)
						return
					}

					iv := NormalizedImpliedBlackVolatility(price, x, Call)

					// Adaptive tolerance based on regime difficulty
					tolerance := 1e-10
					if s < 1e-6 {
						tolerance = 1e-5
					} else if s < 1e-4 {
						tolerance = 1e-8
					} else if price > 0.99*bMax {
						tolerance = 1e-6 // Near-maximum prices have reduced precision
					} else if s > 10 {
						tolerance = 2e-9
					}

					assert.InEpsilon(t, s, iv, tolerance,
						"Volatility round-trip: expected s=%v, got iv=%v", s, iv)
				})
			}
		})
	}
}

// generateSValues creates a range of volatility values to test for a given x
func generateSValues(sMid, absX float64) []float64 {
	seen := make(map[float64]bool)
	add := func(vals ...float64) []float64 {
		var result []float64
		for _, v := range vals {
			if v > 0 && !seen[v] {
				seen[v] = true
				result = append(result, v)
			}
		}
		return result
	}

	var sValues []float64

	// Very small s (approaching zero but still valid)
	if absX > 0 {
		sMin := math.Max(absX*1e-6, 1e-15)
		sValues = append(sValues, add(sMin, sMin*10, sMin*100, sMin*1000)...)
	}

	// Around the critical point s_mid
	if sMid > 0 {
		sValues = append(sValues, add(
			sMid*0.01, sMid*0.1, sMid*0.5,
			sMid*0.9, sMid, sMid*1.1,
			sMid*2, sMid*5, sMid*10,
		)...)
	}

	// Fixed reasonable values
	sValues = append(sValues, add(
		0.01, 0.05, 0.1, 0.2, 0.3, 0.5,
		1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0,
	)...)

	// For large |x|, also test large s
	if absX > 10 {
		sValues = append(sValues, add(absX*0.5, absX*0.8, absX*1.0, absX*1.2)...)
	}

	return sValues
}

func formatTestName(x float64) string {
	if x == 0 {
		return "x=0"
	}
	if x < 0.001 {
		return fmt.Sprintf("x=%.0e", x)
	}
	if x == float64(int(x)) {
		return fmt.Sprintf("x=%d", int(x))
	}
	return fmt.Sprintf("x=%.3f", x)
}

func formatSName(s float64) string {
	if s < 0.001 {
		return fmt.Sprintf("s=%.1e", s)
	}
	if s >= 100 {
		return fmt.Sprintf("s=%.0f", s)
	}
	return fmt.Sprintf("s=%.3f", s)
}

func BenchmarkImpliedVolatility(b *testing.B) {
	cases := []struct {
		name  string
		F, K  float64
		sigma float64
		T     float64
	}{
		// ATM branch - uses direct erfinv formula
		{"ATM", 100, 100, 0.20, 1.0},
		// Lower middle branch
		{"LowerMiddle_OTM5pct", 100, 105, 0.20, 1.0},
		// Upper middle branch
		{"UpperMiddle_ITM5pct", 100, 95, 0.20, 1.0},
		// Lowest branch - deep OTM, uses log objective
		{"Lowest_DeepOTM", 100, 150, 0.20, 0.1},
		// Highest branch - high premium
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
		price := Black(tc.F, tc.K, tc.sigma, tc.T, Call)
		b.Run(tc.name, func(b *testing.B) {
			for b.Loop() {
				ImpliedBlackVolatility(price, tc.F, tc.K, tc.T, Call)
			}
		})
	}
}
