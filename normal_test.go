package letsberational

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"gonum.org/v1/gonum/stat/distuv"
)

func TestNormPDF(t *testing.T) {
	tests := []struct {
		name     string
		z        float64
		expected float64
	}{
		{"Zero", 0.0, 0.3989422804014327},
		{"One Sigma", 1.0, 0.24197072451914337},
		{"Negative Two", -2.0, 0.05399096651318806},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			assert.InDelta(t, tt.expected, normPDF(tt.z), 1e-15)
		})
	}
}

func TestNormCDF(t *testing.T) {
	tests := []struct {
		name      string
		z         float64
		expected  float64
		isEpsilon bool
	}{
		{"Median", 0.0, 0.5, false},
		{"1.96 Sigma", 1.96, 0.9750021048517795, false},
		{"-1.96 Sigma", -1.96, 0.024997895148220435, false},
		// calculated using mpmath Python library
		// 3.6709661993127508e-51 = 0.5 * mpmath.erfc(15 / mpmath.sqrt(2))
		{"Deep Tail", -15.0, 3.6709661993127508e-51, true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			actual := normCDF(tt.z)
			if tt.isEpsilon {
				assert.InEpsilon(t, tt.expected, actual, 1e-14)
			} else {
				assert.InDelta(t, tt.expected, actual, 1e-15)
			}
		})
	}
}

func TestInverseNormCDF(t *testing.T) {
	t.Run("RoundTripIdentity", func(t *testing.T) {
		for _, p := range []float64{0.000001, 0.01, 0.5, 0.99, 0.999999} {
			t.Run(fmt.Sprintf("p=%v", p), func(t *testing.T) {
				z := inverseNormCDF(p)
				pRoundTrip := normCDF(z)
				assert.InEpsilon(t, p, pRoundTrip, 1e-13)
			})
		}
	})
}

var result float64

func BenchmarkCDFComparison(b *testing.B) {
	z := -12.0
	b.Run("Custom", func(b *testing.B) {
		for b.Loop() {
			result = normCDF(z)
		}
	})
	b.Run("MathErfc", func(b *testing.B) {
		for b.Loop() {
			result = 0.5 * math.Erfc(-z*invSqrtTwo)
		}
	})
}

func BenchmarkInverseCDFComparison(b *testing.B) {
	p := 0.00001
	b.Run("Custom", func(b *testing.B) {
		for b.Loop() {
			result = inverseNormCDF(p)
		}
	})
	b.Run("Gonum", func(b *testing.B) {
		dist := distuv.UnitNormal
		for b.Loop() {
			result = dist.Quantile(p)
		}
	})
}
