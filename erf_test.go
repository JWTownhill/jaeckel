package letsberational

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestErf(t *testing.T) {
	for _, x := range []float64{0, 0.1, 0.46875, 1.0, 2.5, 5.0} {
		assert.InDelta(t, math.Erf(x), erf(x), 1e-15, "Erf failed at x=%f", x)
		assert.InDelta(t, math.Erf(-x), erf(-x), 1e-15, "Erf failed at x=%f", -x)
	}
	assert.Equal(t, 1.0, erf(30.0))
	assert.Equal(t, -1.0, erf(-30.0))
}

func TestErfc(t *testing.T) {
	for _, x := range []float64{0, 0.1, 0.46875, 1.0, 4.0, 10.0} {
		expected := math.Erfc(x)
		actual := erfc(x)

		if expected > 1e-10 {
			assert.InDelta(t, expected, actual, 1e-15, "Erfc absolute error at x=%f", x)
		} else {
			// For very small values, relative error is more important
			assert.InEpsilon(t, expected, actual, 1e-13, "Erfc relative error at x=%f", x)
		}
	}

	assert.Equal(t, 0.0, erfc(30.0), "Erfc should underflow to 0")
	assert.Equal(t, 2.0, erfc(-30.0), "Erfc(-large) should be 2.0")
}

// naiveErfcx is the unstable implementation used for comparison in small ranges.
func naiveErfcx(x float64) float64 {
	return math.Exp(x*x) * math.Erfc(x)
}
func TestErfcx(t *testing.T) {
	t.Run("Stable Range Comparison", func(t *testing.T) {
		smallValues := []float64{0.0, 0.1, 0.4, -0.4, 2.0, 5.0}
		for _, x := range smallValues {
			assert.InEpsilon(t, naiveErfcx(x), erfcx(x), 1e-13, "Relative error too high at x=%f", x)
		}
	})

	t.Run("Overflow Avoidance", func(t *testing.T) {
		x := 30.0
		actual := erfcx(x)

		require.False(t, math.IsInf(actual, 0), "Large positive x=%f should not be Infinity", x)
		require.False(t, math.IsNaN(actual), "Large positive x=%f should not be NaN", x)

		assert.Positive(t, actual, "erfcx should always be positive")

		assert.InDelta(t, 0.0188, actual, 1e-5, "Large positive x=%f: value out of expected bounds", x)
	})

	t.Run("Negative Large Value", func(t *testing.T) {
		x := -25.0
		actual := erfcx(x)
		assert.Greater(t, actual, 1e100, "Large negative x=%f should result in a massive value", x)
		assert.Less(t, actual, math.MaxFloat64, "x=-25 should still be under the MaxFloat64 cap")
	})

	t.Run("XNEG Boundary", func(t *testing.T) {
		// xNeg is -26.628... any x smaller should cap at MaxFloat64
		assert.Equal(t, math.MaxFloat64, erfcx(xNeg-1e-14), "Should return MaxFloat64 for values below xNeg")
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
