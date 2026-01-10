package letsberational

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestNormalizedVega(t *testing.T) {
	tests := []struct {
		name     string
		x        float64
		s        float64
		expected float64
		delta    float64
	}{
		{
			"Standard Case",
			0.2, 0.3,
			math.Exp(-0.5*(math.Pow(0.2/0.3, 2)+math.Pow(0.3/2, 2))) / (math.Sqrt(2 * math.Pi)),
			1e-15,
		},
		{
			"ATM Case (x=0)",
			0.0, 0.2,
			(1 / math.Sqrt(2*math.Pi)) * math.Exp(-0.125*0.2*0.2),
			1e-15,
		},
		{
			"Exact Zero Vol",
			0.5, 0.0,
			0.0,
			0.0,
		},
		{
			"Tiny Vol Safeguard (s <= ax * sqrt(DBL_MIN))",
			1.0, 1e-160,
			0.0,
			0.0,
		},
		{
			"Extreme OTM Underflow (Large x)",
			600.0, 0.5,
			0.0,
			0.0,
		},
		{
			"Deep ITM Underflow (Negative x)",
			-600.0, 0.5,
			0.0,
			0.0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := normalizedVega(tt.x, tt.s)
			if tt.delta == 0 {
				assert.Equal(t, tt.expected, actual)
			} else {
				assert.InDelta(t, tt.expected, actual, tt.delta)
			}
		})
	}
}

func TestNormalisedBlackWithOptimalUseOfCodysFunctions(t *testing.T) {
	t.Run("At the Money (ATM)", func(t *testing.T) {
		thetaX := 0.0
		s := 0.2
		expected := math.Erf(s / (2 * math.Sqrt(2)))
		actual := normalisedBlackWithOptimalUseOfCodysFunctions(thetaX, s)
		assert.InDelta(t, expected, actual, 1e-15)
	})

	t.Run("Deep In The Money", func(t *testing.T) {
		// As thetaX becomes large, call price approaches 1 - exp(-thetaX)
		thetaX := 10.0
		s := 0.1
		actual := normalisedBlackWithOptimalUseOfCodysFunctions(thetaX, s)
		assert.Greater(t, actual, 0.99)
	})

	t.Run("Deep Out Of The Money stability", func(t *testing.T) {
		thetaX, s := -20.0, 0.1
		actual := normalisedBlackWithOptimalUseOfCodysFunctions(thetaX, s)

		assert.False(t, math.IsNaN(actual), "Result should not be NaN")
		assert.False(t, math.IsInf(actual, 0), "Result should not be Inf")
		assert.GreaterOrEqual(t, actual, 0.0, "Option price cannot be negative")
		assert.Less(t, actual, 1e-10, "Deep OTM price should be near zero but stable")
	})

	t.Run("Intrinsic Value Floor", func(t *testing.T) {
		thetaX := -500.0
		s := 0.0001
		actual := normalisedBlackWithOptimalUseOfCodysFunctions(thetaX, s)
		assert.Equal(t, 0.0, actual)
	})
}

func TestBlackConsistency(t *testing.T) {
	t.Run("Monotonicity", func(t *testing.T) {
		s := 0.2
		price1 := normalisedBlackWithOptimalUseOfCodysFunctions(-0.1, s)
		price2 := normalisedBlackWithOptimalUseOfCodysFunctions(0.0, s)
		price3 := normalisedBlackWithOptimalUseOfCodysFunctions(0.1, s)

		assert.Less(t, price1, price2)
		assert.Less(t, price2, price3)
	})
}
