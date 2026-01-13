package jaeckel

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestBlack(t *testing.T) {
	t.Run("ATM Call", func(t *testing.T) {
		F, K, sigma, T := 100.0, 100.0, 0.20, 1.0
		price := Black(F, K, sigma, T, Call)
		assert.Greater(t, price, 0.0)
		assert.Less(t, price, F)
	})

	t.Run("Zero Vol Returns Intrinsic", func(t *testing.T) {
		F, K, T := 100.0, 95.0, 1.0
		callPrice := Black(F, K, 0, T, Call)
		putPrice := Black(F, K, 0, T, Put)
		assert.Equal(t, F-K, callPrice)
		assert.Equal(t, 0.0, putPrice)
	})
}

func TestNormalizedBlackConsistency(t *testing.T) {
	F := 100.0
	K := 105.0
	sigma := 0.25
	T := 0.5

	x := math.Log(F / K)
	s := sigma * math.Sqrt(T)

	for _, optType := range []OptionType{Call, Put} {
		name := "Call"
		if optType == Put {
			name = "Put"
		}
		t.Run(name, func(t *testing.T) {
			blackPrice := Black(F, K, sigma, T, optType)
			normPrice := NormalizedBlack(x, s, optType)

			expectedBlackPrice := math.Sqrt(F*K) * normPrice

			assert.InEpsilon(t, expectedBlackPrice, blackPrice, 1e-14)
		})
	}
}

func TestVegaConsistency(t *testing.T) {
	F := 100.0
	K := 100.0
	sigma := 0.20
	T := 1.0

	vega := Vega(F, K, sigma, T)
	assert.Greater(t, vega, 0.0, "Vega should be positive")

	eps := 1e-6
	priceUp := Black(F, K, sigma+eps, T, Call)
	priceDown := Black(F, K, sigma-eps, T, Call)
	numericalVega := (priceUp - priceDown) / (2 * eps)

	assert.InEpsilon(t, numericalVega, vega, 1e-6,
		"Vega analytical=%v, numerical=%v", vega, numericalVega)
}

func TestNormalizedVega(t *testing.T) {
	t.Run("ATM (x=0)", func(t *testing.T) {
		s := 0.2
		expected := (1 / math.Sqrt(2*math.Pi)) * math.Exp(-0.125*s*s)
		actual := NormalizedVega(0, s)
		assert.InDelta(t, expected, actual, 1e-15)
	})

	t.Run("Zero vol returns zero", func(t *testing.T) {
		assert.Equal(t, 0.0, NormalizedVega(0.5, 0))
		assert.Equal(t, 0.0, NormalizedVega(-0.5, 0))
	})

	t.Run("Tiny vol relative to |x| returns zero", func(t *testing.T) {
		// s <= ax * sqrt(SmallestNonzeroFloat64) should return 0
		ax := 1.0
		tinyS := ax * math.Sqrt(math.SmallestNonzeroFloat64) * 0.5
		assert.Equal(t, 0.0, NormalizedVega(ax, tinyS))
		assert.Equal(t, 0.0, NormalizedVega(-ax, tinyS))
	})

	t.Run("Standard case", func(t *testing.T) {
		x, s := 0.2, 0.3
		expected := math.Exp(-0.5*(math.Pow(x/s, 2)+math.Pow(s/2, 2))) / math.Sqrt(2*math.Pi)
		actual := NormalizedVega(x, s)
		assert.InDelta(t, expected, actual, 1e-15)
	})

	t.Run("Extreme OTM underflow", func(t *testing.T) {
		assert.Equal(t, 0.0, NormalizedVega(600.0, 0.5))
		assert.Equal(t, 0.0, NormalizedVega(-600.0, 0.5))
	})
}

func TestNormalizedVolga(t *testing.T) {
	t.Run("Zero vol returns zero", func(t *testing.T) {
		assert.Equal(t, 0.0, NormalizedVolga(0.5, 0))
		assert.Equal(t, 0.0, NormalizedVolga(0, 0))
	})

	t.Run("Tiny vol relative to |x| returns zero", func(t *testing.T) {
		ax := 1.0
		tinyS := ax * math.Sqrt(math.SmallestNonzeroFloat64) * 0.5
		assert.Equal(t, 0.0, NormalizedVolga(ax, tinyS))
	})

	t.Run("ATM (x=0)", func(t *testing.T) {
		s := 0.5
		t2 := 0.25 * s * s
		expected := (1 / math.Sqrt(2*math.Pi)) * math.Exp(-0.5*t2) * -t2 / s
		actual := NormalizedVolga(0, s)
		assert.InDelta(t, expected, actual, 1e-15)
	})

	t.Run("Standard case", func(t *testing.T) {
		x, s := 0.2, 0.3
		actual := NormalizedVolga(x, s)
		assert.False(t, math.IsNaN(actual))
		assert.False(t, math.IsInf(actual, 0))
	})
}

func TestVolga(t *testing.T) {
	F := 100.0
	K := 100.0
	sigma := 0.20
	T := 1.0

	volga := Volga(F, K, sigma, T)

	eps := 1e-4
	vegaUp := Vega(F, K, sigma+eps, T)
	vegaDown := Vega(F, K, sigma-eps, T)
	numericalVolga := (vegaUp - vegaDown) / (2 * eps)

	assert.InEpsilon(t, numericalVolga, volga, 1e-4,
		"Volga analytical=%v, numerical=%v", volga, numericalVolga)
}

func TestNormalizedBlackStandard(t *testing.T) {
	t.Run("At the Money (ATM)", func(t *testing.T) {
		thetaX := 0.0
		s := 0.2
		expected := math.Erf(s / (2 * math.Sqrt(2)))
		actual := normalizedBlackStandard(thetaX, s)
		assert.InDelta(t, expected, actual, 1e-15)
	})

	t.Run("Deep In The Money", func(t *testing.T) {
		// As thetaX becomes large, call price approaches 1 - exp(-thetaX)
		thetaX := 10.0
		s := 0.1
		actual := normalizedBlackStandard(thetaX, s)
		assert.Greater(t, actual, 0.99)
	})

	t.Run("Deep Out Of The Money stability", func(t *testing.T) {
		thetaX, s := -20.0, 0.1
		actual := normalizedBlackStandard(thetaX, s)

		assert.False(t, math.IsNaN(actual), "Result should not be NaN")
		assert.False(t, math.IsInf(actual, 0), "Result should not be Inf")
		assert.GreaterOrEqual(t, actual, 0.0, "Option price cannot be negative")
		assert.Less(t, actual, 1e-10, "Deep OTM price should be near zero but stable")
	})

	t.Run("Intrinsic Value Floor", func(t *testing.T) {
		thetaX := -500.0
		s := 0.0001
		actual := normalizedBlackStandard(thetaX, s)
		assert.Equal(t, 0.0, actual)
	})

	t.Run("q1 < threshold, q2 >= threshold", func(t *testing.T) {
		// Need: q1 < codyThreshold and q2 >= codyThreshold
		// q1 = -invSqrtTwo * (h + t), q2 = -invSqrtTwo * (h - t)
		// where h = thetaX/s and t = s/2
		thetaX := -2.0
		s := 1.5
		actual := normalizedBlackStandard(thetaX, s)
		assert.GreaterOrEqual(t, actual, 0.0)
		assert.False(t, math.IsNaN(actual))
	})

	t.Run("q1 >= threshold, q2 < threshold", func(t *testing.T) {
		// Need: q1 >= codyThreshold and q2 < codyThreshold
		thetaX := 3.0
		s := 0.8
		actual := normalizedBlackStandard(thetaX, s)
		assert.GreaterOrEqual(t, actual, 0.0)
		assert.False(t, math.IsNaN(actual))
	})
}

func TestBlackConsistency(t *testing.T) {
	t.Run("Monotonicity", func(t *testing.T) {
		s := 0.2
		price1 := normalizedBlackStandard(-0.1, s)
		price2 := normalizedBlackStandard(0.0, s)
		price3 := normalizedBlackStandard(0.1, s)

		assert.Less(t, price1, price2)
		assert.Less(t, price2, price3)
	})
}

func TestNormalizedBlackEdgeCases(t *testing.T) {
	t.Run("s <= 0 returns zero", func(t *testing.T) {
		testCases := []struct {
			name   string
			thetaX float64
			s      float64
		}{
			{"s=0, thetaX=0", 0.0, 0.0},
			{"s=0, thetaX positive", 1.0, 0.0},
			{"s=0, thetaX negative", -1.0, 0.0},
			{"s negative", 0.5, -0.1},
			{"s very negative", 0.5, -10.0},
		}

		for _, tc := range testCases {
			t.Run(tc.name, func(t *testing.T) {
				result := normalizedBlack(tc.thetaX, tc.s)
				assert.Equal(t, 0.0, result, "normalizedBlack should return 0 when s <= 0")
			})
		}
	})

	t.Run("ATM formula across volatilities", func(t *testing.T) {
		// At the money: b(0, s) = 2*Phi(s/2) - 1 = erf(s/(2*sqrt(2)))
		// This tests both the mathematical identity and call/put symmetry at ATM
		testCases := []float64{0.01, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0}

		for _, s := range testCases {
			t.Run(fmt.Sprintf("s=%.2f", s), func(t *testing.T) {
				result := normalizedBlack(0.0, s)
				expected := 2.0*NormCDF(0.5*s) - 1.0

				assert.InDelta(t, expected, result, 1e-15,
					"ATM case should match 2*Phi(s/2) - 1 for s=%v", s)

				erfExpected := math.Erf(s / (2 * math.Sqrt(2)))
				assert.InDelta(t, erfExpected, result, 1e-15,
					"ATM case should match erf(s/(2*sqrt(2))) for s=%v", s)
			})
		}
	})
}

func TestRegionI(t *testing.T) {
	testCases := []struct {
		name string
		x    float64
		s    float64
	}{
		{"h=-14", -7.0, 0.5},
		{"h=-16", -8.0, 0.5},
		{"h=-20", -10.0, 0.5},
		{"h=-15", -3.0, 0.2},
		{"h=-25", -5.0, 0.2},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			price := NormalizedBlack(tc.x, tc.s, Call)
			assert.Greater(t, price, 0.0, "Price should be positive")
			assert.False(t, math.IsNaN(price), "Price should not be NaN")
			assert.False(t, math.IsInf(price, 0), "Price should not be Inf")

			iv := NormalizedImpliedBlackVolatility(price, tc.x, Call)
			if !math.IsInf(iv, 0) && !math.IsNaN(iv) && iv > 0 {
				priceRoundTrip := NormalizedBlack(tc.x, iv, Call)
				assert.InEpsilon(t, price, priceRoundTrip, 1e-10,
					"Round-trip failed: original=%v, roundtrip=%v", price, priceRoundTrip)
			}
		})
	}
}

func TestRegionII(t *testing.T) {
	testCases := []struct {
		name string
		x    float64
		s    float64
	}{
		// Small s values that trigger Region II
		{"SmallS_ATM", 0.0, 0.1},
		{"SmallS_SlightlyOTM", -0.05, 0.08},
		{"SmallS_OTM", -0.1, 0.1},
		{"VerySmallS", -0.02, 0.03},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			price := NormalizedBlack(tc.x, tc.s, Call)
			assert.Greater(t, price, 0.0, "Price should be positive")
			assert.False(t, math.IsNaN(price), "Price should not be NaN")
			assert.False(t, math.IsInf(price, 0), "Price should not be Inf")

			iv := NormalizedImpliedBlackVolatility(price, tc.x, Call)
			if !math.IsInf(iv, 0) && !math.IsNaN(iv) && iv > 0 {
				priceRoundTrip := NormalizedBlack(tc.x, iv, Call)
				assert.InEpsilon(t, price, priceRoundTrip, 1e-10,
					"Round-trip failed: original=%v, roundtrip=%v", price, priceRoundTrip)
			}
		})
	}
}

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
				Black(tc.F, tc.K, tc.sigma, tc.T, Call)
			}
		})
	}
}

func BenchmarkNormalizedBlack(b *testing.B) {
	cases := []struct {
		name string
		x, s float64
	}{
		{"ATM_s0.2", 0, 0.2},
		{"ATM_s1.0", 0, 1.0},
		{"RegionI_DeepTail", -7.0, 0.5},      // h=-14, triggers asymptotic expansion
		{"RegionI_Extreme", -15.0, 0.5},      // h=-30, deep in Region I
		{"RegionII_SmallS", -0.05, 0.08},     // small-t expansion
		{"Standard", -0.5, 0.5},              // Regions III/IV (standard computation)
		{"ExtremeMoney_Neg10", -10.0, 1.0},   // h=-10, just outside Region I
		{"ExtremeMoney_Neg100", -100.0, 5.0}, // h=-20, Region I
	}

	for _, tc := range cases {
		b.Run(tc.name, func(b *testing.B) {
			for b.Loop() {
				NormalizedBlack(tc.x, tc.s, Call)
			}
		})
	}
}

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
