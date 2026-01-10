package letsberational

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestEvalRational(t *testing.T) {
	tests := []struct {
		name     string
		x        float64
		num      []float64
		den      []float64
		expected float64
	}{
		{
			"Simple ratio 10/2",
			1.0,
			[]float64{10},
			[]float64{2},
			5.0,
		},
		{
			"Linear: (3x + 1) / (x + 1) at x=2",
			2.0,
			[]float64{3, 1}, // 3(2)+1 = 7
			[]float64{1, 1}, // 1(2)+1 = 3
			7.0 / 3.0,
		},
		{
			"Quadratic: (x^2 + 2x + 1) / (x + 1) at x=3",
			3.0,
			[]float64{1, 2, 1}, // (x+1)^2 = 16
			[]float64{1, 1},    // (x+1) = 4
			4.0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.InDelta(t, tt.expected, evalRational(tt.x, tt.num, tt.den), 1e-9, "Evaluation should match within delta")
		})
	}
}
