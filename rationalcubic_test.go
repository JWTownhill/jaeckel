package letsberational

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestInterpolate(t *testing.T) {
	tests := []struct {
		name     string
		x        float64
		xL, xR   float64
		yL, yR   float64
		dL, dR   float64
		r        float64
		expected float64
		delta    float64
	}{
		{"Midpoint Standard Cubic", 0.5, 0, 1, 0, 1, 1, 1, 3.0, 0.5, 1e-15},
		{"Left Boundary", 0.0, 0, 1, 10, 20, 1, 1, 3.0, 10.0, 1e-15},
		{"Right Boundary", 1.0, 0, 1, 10, 20, 1, 1, 3.0, 20.0, 1e-15},
		{"High R (Linear)", 0.5, 0, 1, 0, 10, 50, 50, 1e12, 5.0, 1e-7},
		{"Zero Interval", 0.5, 1, 1, 10, 20, 1, 1, 3.0, 15.0, 1e-15},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := interpolate(tt.x, tt.xL, tt.xR, tt.yL, tt.yR, tt.dL, tt.dR, tt.r)
			assert.InDelta(t, tt.expected, actual, tt.delta)
		})
	}
}

func TestMinimumControlParameter(t *testing.T) {
	tests := []struct {
		name     string
		dL       float64
		dR       float64
		s        float64
		prefer   bool
		expected float64
	}{
		{
			"Perfectly Linear Monotonic",
			1.0, 1.0, 1.0,
			false,
			2.0,
		},
		{
			"Increasing Monotonic",
			1.0, 3.0, 2.0,
			false,
			2.0,
		},
		{
			"Convex Shape (Non-Monotonic)",
			-1.0, 1.0, 0.5,
			false,
			4.0,
		},
		{
			"Flat Segment with Shape Preservation",
			0.0, 0.0, 0.0,
			true,
			maximumControlParam,
		},
		{
			"Oscillatory Fallback (Standard Cubic)",
			10.0, 10.0, 0.0,
			false,
			minimumControlParam,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := minimumControlParameter(tt.dL, tt.dR, tt.s, tt.prefer)
			assert.InDelta(t, tt.expected, actual, 1e-15)
		})
	}
}

func TestSecondDerivativeFit(t *testing.T) {
	tests := []struct {
		name     string
		xL, xR   float64
		yL, yR   float64
		dL, dR   float64
		deriv    float64
		expected float64
	}{
		{"Left Side Flat", 0, 1, 0, 1, 1, 1, 0.0, minimumControlParam},
		{"Left Side Positive Curve", 0, 1, 0, 1, 0, 2, 2.0, 3.0},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := controlParameterToFitSecondDerivativeAtLeftSide(
				tt.xL, tt.xR, tt.yL, tt.yR, tt.dL, tt.dR, tt.deriv,
			)
			assert.InDelta(t, tt.expected, actual, 1e-15)
		})
	}
}
