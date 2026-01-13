package jaeckel

import (
	"math"
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
		{"MaxR Linear Fallback", 0.5, 0, 1, 0, 10, 1, 1, maximumControlParam, 5.0, 1e-15},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := interpolate(tt.x, tt.xL, tt.xR, tt.yL, tt.yR, tt.dL, tt.dR, tt.r)
			assert.InDelta(t, tt.expected, actual, tt.delta)
		})
	}
}

func TestControlParameterToFitSecondDerivativeAtLeftSide(t *testing.T) {
	tests := []struct {
		name                string
		xL, xR              float64
		yL, yR              float64
		dL, dR              float64
		secondDerivativeL   float64
		expected            float64
		useInDeltaAssertion bool
	}{
		{
			name: "Flat case returns minimumControlParam",
			// denominator = (yR-yL)/h - dL = (1-0)/1 - 1 = 0
			// numerator = 0.5*h*0 + (dR - dL) = 0 + (1-1) = 0
			xL: 0, xR: 1, yL: 0, yR: 1, dL: 1, dR: 1,
			secondDerivativeL: 0.0,
			expected:          minimumControlParam,
		},
		{
			name: "Positive curvature",
			xL:   0, xR: 1, yL: 0, yR: 1, dL: 0, dR: 2,
			secondDerivativeL:   2.0,
			expected:            3.0,
			useInDeltaAssertion: true,
		},
		{
			name: "Denominator zero with positive numerator returns maximumControlParam",
			// h = 1, denominator = (yR-yL)/h - dL = 1 - 1 = 0
			// numerator = 0.5*1*2 + (2-0) = 1 + 2 = 3 > 0
			xL: 0, xR: 1, yL: 0, yR: 1, dL: 1, dR: 2,
			secondDerivativeL: 2.0,
			expected:          maximumControlParam,
		},
		{
			name: "Denominator zero with negative numerator returns minimumControlParam",
			// h = 1, denominator = (yR-yL)/h - dL = 1 - 1 = 0
			// numerator = 0.5*1*(-4) + (0-0) = -2 + 0 = -2 < 0
			xL: 0, xR: 1, yL: 0, yR: 1, dL: 1, dR: 0,
			secondDerivativeL: -4.0,
			expected:          minimumControlParam,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := controlParameterToFitSecondDerivativeAtLeftSide(
				tt.xL, tt.xR, tt.yL, tt.yR, tt.dL, tt.dR, tt.secondDerivativeL,
			)
			if tt.useInDeltaAssertion {
				assert.InDelta(t, tt.expected, actual, 1e-15)
			} else {
				assert.Equal(t, tt.expected, actual)
			}
		})
	}
}

func TestControlParameterToFitSecondDerivativeAtRightSide(t *testing.T) {
	tests := []struct {
		name                string
		xL, xR              float64
		yL, yR              float64
		dL, dR              float64
		secondDerivativeR   float64
		expected            float64
		useInDeltaAssertion bool
	}{
		{
			name: "Flat case returns minimumControlParam",
			// denominator = dR - (yR-yL)/h = 1 - (1-0)/1 = 0
			// numerator = 0.5*1*0 + (1-1) = 0
			xL: 0, xR: 1, yL: 0, yR: 1, dL: 1, dR: 1,
			secondDerivativeR: 0.0,
			expected:          minimumControlParam,
		},
		{
			name: "Positive curvature",
			xL:   0, xR: 1, yL: 0, yR: 1, dL: 0, dR: 2,
			secondDerivativeR:   2.0,
			expected:            3.0,
			useInDeltaAssertion: true,
		},
		{
			name: "Denominator zero with positive numerator returns maximumControlParam",
			// h = 1, denominator = dR - (yR-yL)/h = 1 - 1 = 0
			// numerator = 0.5*1*4 + (1-0) = 2 + 1 = 3 > 0
			xL: 0, xR: 1, yL: 0, yR: 1, dL: 0, dR: 1,
			secondDerivativeR: 4.0,
			expected:          maximumControlParam,
		},
		{
			name: "Denominator zero with negative numerator returns minimumControlParam",
			// h = 1, denominator = dR - (yR-yL)/h = 1 - 1 = 0
			// numerator = 0.5*1*(-8) + (1-0) = -4 + 1 = -3 < 0
			xL: 0, xR: 1, yL: 0, yR: 1, dL: 0, dR: 1,
			secondDerivativeR: -8.0,
			expected:          minimumControlParam,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := controlParameterToFitSecondDerivativeAtRightSide(
				tt.xL, tt.xR, tt.yL, tt.yR, tt.dL, tt.dR, tt.secondDerivativeR,
			)
			if tt.useInDeltaAssertion {
				assert.InDelta(t, tt.expected, actual, 1e-15)
			} else {
				assert.Equal(t, tt.expected, actual)
			}
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
			"Monotonic with Zero Slope",
			// dL*s = 0 >= 0 ✓, dR*s = 0 >= 0 ✓, so monotonic
			// Also convex and concave (dL <= s <= dR and dL >= s >= dR both true when all zero)
			10.0, 10.0, 0.0,
			false,
			minimumControlParam,
		},
		{
			"Concave Only (Non-Monotonic)",
			// dL*s = 0.5 >= 0, dR*s = -0.5 < 0, so not monotonic
			// convex: 1 <= 0.5 is false
			// concave: 1 >= 0.5 && 0.5 >= -1 is true
			1.0, -1.0, 0.5,
			false,
			4.0,
		},
		{
			"Not Monotonic Not Convex Not Concave",
			// dL*s = -2 < 0, so not monotonic
			// convex: 2 <= -1 is false
			// concave: 2 >= -1 && -1 >= 3 is false
			// Hits early return branch
			2.0, 3.0, -1.0,
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

func TestConvexControlParameterToFitSecondDerivativeAtLeftSide(t *testing.T) {
	tests := []struct {
		name                     string
		xL, xR                   float64
		yL, yR                   float64
		dL, dR                   float64
		secondDerivativeL        float64
		preferShapePreservation  bool
		expectedGreaterOrEqualTo float64
	}{
		{
			name: "Basic convex fit",
			xL:   0, xR: 1, yL: 0, yR: 1, dL: 0, dR: 2,
			secondDerivativeL:        2.0,
			preferShapePreservation:  false,
			expectedGreaterOrEqualTo: minimumControlParam,
		},
		{
			name: "With shape preservation",
			xL:   0, xR: 1, yL: 0, yR: 1, dL: 1, dR: 1,
			secondDerivativeL:        0.0,
			preferShapePreservation:  true,
			expectedGreaterOrEqualTo: minimumControlParam,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := convexControlParameterToFitSecondDerivativeAtLeftSide(
				tt.xL, tt.xR, tt.yL, tt.yR, tt.dL, tt.dR, tt.secondDerivativeL, tt.preferShapePreservation,
			)
			assert.GreaterOrEqual(t, actual, tt.expectedGreaterOrEqualTo)
		})
	}
}

func TestConvexControlParameterToFitSecondDerivativeAtRightSide(t *testing.T) {
	tests := []struct {
		name                     string
		xL, xR                   float64
		yL, yR                   float64
		dL, dR                   float64
		secondDerivativeR        float64
		preferShapePreservation  bool
		expectedGreaterOrEqualTo float64
	}{
		{
			name: "Basic convex fit",
			xL:   0, xR: 1, yL: 0, yR: 1, dL: 0, dR: 2,
			secondDerivativeR:        2.0,
			preferShapePreservation:  false,
			expectedGreaterOrEqualTo: minimumControlParam,
		},
		{
			name: "With shape preservation",
			xL:   0, xR: 1, yL: 0, yR: 1, dL: 1, dR: 1,
			secondDerivativeR:        0.0,
			preferShapePreservation:  true,
			expectedGreaterOrEqualTo: minimumControlParam,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := convexControlParameterToFitSecondDerivativeAtRightSide(
				tt.xL, tt.xR, tt.yL, tt.yR, tt.dL, tt.dR, tt.secondDerivativeR, tt.preferShapePreservation,
			)
			assert.GreaterOrEqual(t, actual, tt.expectedGreaterOrEqualTo)
		})
	}
}

func TestIsZero(t *testing.T) {
	tests := []struct {
		name     string
		x        float64
		expected bool
	}{
		{"Exact zero", 0.0, true},
		{"Smallest positive", math.SmallestNonzeroFloat64, false},
		{"Negative smallest", -math.SmallestNonzeroFloat64, false},
		{"Half smallest", math.SmallestNonzeroFloat64 / 2, true},
		{"Normal value", 1.0, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.Equal(t, tt.expected, isZero(tt.x))
		})
	}
}
