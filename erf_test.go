package jaeckel

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestErf(t *testing.T) {
	tests := []struct {
		name     string
		x        float64
		expected float64
		delta    float64
	}{
		{"zero", 0, 0, 1e-15},
		{"small positive", 0.1, math.Erf(0.1), 1e-15},
		{"at threshold", 0.46875, math.Erf(0.46875), 1e-15},
		{"one", 1.0, math.Erf(1.0), 1e-15},
		{"medium", 2.5, math.Erf(2.5), 1e-15},
		{"larger", 5.0, math.Erf(5.0), 1e-15},
		{"small negative", -0.1, math.Erf(-0.1), 1e-15},
		{"negative at threshold", -0.46875, math.Erf(-0.46875), 1e-15},
		{"negative one", -1.0, math.Erf(-1.0), 1e-15},
		{"negative medium", -2.5, math.Erf(-2.5), 1e-15},
		{"negative larger", -5.0, math.Erf(-5.0), 1e-15},
		{"large positive saturates", 30.0, 1.0, 0},
		{"large negative saturates", -30.0, -1.0, 0},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.InDelta(t, tt.expected, Erf(tt.x), tt.delta, "Erf(%v)", tt.x)
		})
	}
}

func TestErfSpecialValues(t *testing.T) {
	t.Run("NaN input", func(t *testing.T) {
		assert.True(t, math.IsNaN(Erf(math.NaN())), "Erf(NaN) should be NaN")
	})

	t.Run("+Inf input", func(t *testing.T) {
		assert.Equal(t, 1.0, Erf(math.Inf(1)), "Erf(+Inf) should be 1")
	})

	t.Run("-Inf input", func(t *testing.T) {
		assert.Equal(t, -1.0, Erf(math.Inf(-1)), "Erf(-Inf) should be -1")
	})
}

func TestErfc(t *testing.T) {
	tests := []struct {
		name        string
		x           float64
		expected    float64
		useRelative bool // use relative error for very small values
		tolerance   float64
	}{
		{"zero", 0, math.Erfc(0), false, 1e-15},
		{"small positive", 0.1, math.Erfc(0.1), false, 1e-15},
		{"at threshold", 0.46875, math.Erfc(0.46875), false, 1e-15},
		{"one", 1.0, math.Erfc(1.0), false, 1e-15},
		{"four", 4.0, math.Erfc(4.0), true, 1e-13},
		{"ten", 10.0, math.Erfc(10.0), true, 1e-13},
		// Negative values - erfc(-x) = 2 - erfc(x)
		{"small negative", -0.1, math.Erfc(-0.1), false, 1e-15},
		{"negative at threshold", -0.46875, math.Erfc(-0.46875), false, 1e-15},
		{"negative one", -1.0, math.Erfc(-1.0), false, 1e-15},
		{"negative two", -2.0, math.Erfc(-2.0), false, 1e-15},
		{"negative four", -4.0, math.Erfc(-4.0), false, 1e-15},
		{"negative ten", -10.0, math.Erfc(-10.0), false, 1e-15},
		// Boundary cases
		{"large positive underflows", 30.0, 0.0, false, 0},
		{"large negative saturates", -30.0, 2.0, false, 0},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := Erfc(tt.x)
			if tt.useRelative && tt.expected > 0 {
				assert.InEpsilon(t, tt.expected, actual, tt.tolerance, "Erfc(%v)", tt.x)
			} else {
				assert.InDelta(t, tt.expected, actual, tt.tolerance, "Erfc(%v)", tt.x)
			}
		})
	}
}

func TestErfcSpecialValues(t *testing.T) {
	t.Run("NaN input", func(t *testing.T) {
		assert.True(t, math.IsNaN(Erfc(math.NaN())), "Erfc(NaN) should be NaN")
	})

	t.Run("+Inf input", func(t *testing.T) {
		assert.Equal(t, 0.0, Erfc(math.Inf(1)), "Erfc(+Inf) should be 0")
	})

	t.Run("-Inf input", func(t *testing.T) {
		assert.Equal(t, 2.0, Erfc(math.Inf(-1)), "Erfc(-Inf) should be 2")
	})
}

// naiveErfcx is the unstable implementation used for comparison in small ranges.
func naiveErfcx(x float64) float64 {
	return math.Exp(x*x) * math.Erfc(x)
}

func TestErfcx(t *testing.T) {
	t.Run("Stable Range Comparison", func(t *testing.T) {
		tests := []struct {
			name string
			x    float64
		}{
			{"zero", 0.0},
			{"small positive", 0.1},
			{"below threshold", 0.4},
			{"negative below threshold", -0.4},
			{"medium", 2.0},
			{"larger", 5.0},
		}
		for _, tt := range tests {
			t.Run(tt.name, func(t *testing.T) {
				assert.InEpsilon(t, naiveErfcx(tt.x), Erfcx(tt.x), 1e-13, "Relative error too high at x=%f", tt.x)
			})
		}
	})

	t.Run("Overflow Avoidance", func(t *testing.T) {
		x := 30.0
		actual := Erfcx(x)

		require.False(t, math.IsInf(actual, 0), "Large positive x=%f should not be Infinity", x)
		require.False(t, math.IsNaN(actual), "Large positive x=%f should not be NaN", x)

		assert.Positive(t, actual, "erfcx should always be positive")

		// Expected value from Wolfram Alpha: erfcx(30) ≈ 0.018795888861...
		assert.InDelta(t, 0.01879588886, actual, 1e-10, "Large positive x=%f: value out of expected bounds", x)
	})

	t.Run("Negative Large Value", func(t *testing.T) {
		x := -25.0
		actual := Erfcx(x)
		assert.Greater(t, actual, 1e100, "Large negative x=%f should result in a massive value", x)
		assert.Less(t, actual, math.MaxFloat64, "x=-25 should still be under the MaxFloat64 cap")
	})

	t.Run("XNEG Boundary", func(t *testing.T) {
		// xNeg is -26.628... any x smaller should cap at MaxFloat64
		assert.Equal(t, math.MaxFloat64, Erfcx(xNeg-1e-14), "Should return MaxFloat64 for values below xNeg")
	})
}

func TestErfcxSpecialValues(t *testing.T) {
	t.Run("NaN input", func(t *testing.T) {
		assert.True(t, math.IsNaN(Erfcx(math.NaN())), "Erfcx(NaN) should be NaN")
	})

	t.Run("+Inf input", func(t *testing.T) {
		result := Erfcx(math.Inf(1))
		assert.Equal(t, 0.0, result, "Erfcx(+Inf) should be 0")
	})

	t.Run("-Inf input", func(t *testing.T) {
		result := Erfcx(math.Inf(-1))
		assert.Equal(t, math.MaxFloat64, result, "Erfcx(-Inf) should be MaxFloat64")
	})
}

func TestOneMinusErfcx(t *testing.T) {
	t.Run("Zero and Near-Zero", func(t *testing.T) {
		assert.InDelta(t, 0.0, oneMinusErfcx(0.0), 1e-16)

		smallX := 1e-7
		// f(x) ≈ 2/√π·x - x²
		expected := (twoOverSqrtPi * smallX) - (smallX * smallX)
		assert.InEpsilon(t, expected, oneMinusErfcx(smallX), 1e-14)
	})

	t.Run("Negative Range", func(t *testing.T) {
		x := -0.1
		actual := oneMinusErfcx(x)
		expected := 1.0 - naiveErfcx(x)

		assert.InEpsilon(t, expected, actual, 1e-14, "Should match naive calculation at x=-0.1")
	})

	t.Run("Boundary Consistency", func(t *testing.T) {
		upper := 1.0 / 3.0
		valUpperInside := oneMinusErfcx(upper - 1e-16)
		valUpperOutside := oneMinusErfcx(upper + 1e-16)
		assert.InDelta(t, valUpperInside, valUpperOutside, 1e-15, "Discontinuity at upper boundary 1/3")

		lower := -0.2
		valLowerInside := oneMinusErfcx(lower + 1e-16)
		valLowerOutside := oneMinusErfcx(lower - 1e-16)
		assert.InDelta(t, valLowerInside, valLowerOutside, 1e-15, "Discontinuity at lower boundary -0.2")
	})
}

func TestErfinvRoundTrip(t *testing.T) {
	tests := []struct {
		name string
		x    float64
	}{
		{"x=-0.9", -0.9},
		{"x=-0.5", -0.5},
		{"x=-0.1", -0.1},
		{"x=0", 0},
		{"x=0.1", 0.1},
		{"x=0.5", 0.5},
		{"x=0.9", 0.9},
		{"x=0.99", 0.99},
		{"x=-0.99", -0.99},
		{"x=0.999", 0.999},
		{"x=-0.999", -0.999},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			y := Erfinv(tt.x)
			roundTrip := Erf(y)

			if tt.x == 0 {
				assert.Equal(t, 0.0, y)
				assert.Equal(t, 0.0, roundTrip)
			} else {
				assert.InEpsilon(t, tt.x, roundTrip, 1e-14,
					"erfinv round trip: x=%v, erfinv(x)=%v, erf(erfinv(x))=%v",
					tt.x, y, roundTrip)
			}
		})
	}
}

func TestErfinvBoundaries(t *testing.T) {
	t.Run("x <= -1 returns -Inf", func(t *testing.T) {
		assert.True(t, math.IsInf(Erfinv(-1), -1))
		assert.True(t, math.IsInf(Erfinv(-1.5), -1))
	})

	t.Run("x >= 1 returns +Inf", func(t *testing.T) {
		assert.True(t, math.IsInf(Erfinv(1), 1))
		assert.True(t, math.IsInf(Erfinv(1.5), 1))
	})

	t.Run("near boundaries", func(t *testing.T) {
		// Values very close to boundaries should still work
		assert.False(t, math.IsInf(Erfinv(0.9999999), 0))
		assert.False(t, math.IsInf(Erfinv(-0.9999999), 0))
	})
}

func TestErfinvSpecialValues(t *testing.T) {
	t.Run("NaN input", func(t *testing.T) {
		assert.True(t, math.IsNaN(Erfinv(math.NaN())), "Erfinv(NaN) should be NaN")
	})

	t.Run("+Inf input", func(t *testing.T) {
		// +Inf is outside domain (-1, 1), should return +Inf
		assert.True(t, math.IsInf(Erfinv(math.Inf(1)), 1), "Erfinv(+Inf) should be +Inf")
	})

	t.Run("-Inf input", func(t *testing.T) {
		// -Inf is outside domain (-1, 1), should return -Inf
		assert.True(t, math.IsInf(Erfinv(math.Inf(-1)), -1), "Erfinv(-Inf) should be -Inf")
	})
}

func TestErfinvForATMImpliedVolatility(t *testing.T) {
	t.Run("beta <= 0 returns 0", func(t *testing.T) {
		assert.Equal(t, 0.0, erfinvForATMImpliedVolatility(0))
		assert.Equal(t, 0.0, erfinvForATMImpliedVolatility(-0.5))
	})

	t.Run("beta >= 1 returns +Inf", func(t *testing.T) {
		assert.True(t, math.IsInf(erfinvForATMImpliedVolatility(1), 1))
		assert.True(t, math.IsInf(erfinvForATMImpliedVolatility(1.5), 1))
	})

	t.Run("round-trip with ATM Black formula", func(t *testing.T) {
		// At ATM: beta = erf(s / sqrt(8)), so s = sqrt(8) * erfinv(beta)
		testSigmas := []float64{0.1, 0.2, 0.5, 1.0, 2.0}
		for _, s := range testSigmas {
			beta := math.Erf(s / math.Sqrt(8))
			sComputed := erfinvForATMImpliedVolatility(beta)
			assert.InEpsilon(t, s, sComputed, 1e-12,
				"s=%v, beta=%v, sComputed=%v", s, beta, sComputed)
		}
	})
}

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
