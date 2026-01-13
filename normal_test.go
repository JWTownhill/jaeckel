package jaeckel

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
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
			assert.InDelta(t, tt.expected, NormPDF(tt.z), 1e-15)
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
			actual := NormCDF(tt.z)
			if tt.isEpsilon {
				assert.InEpsilon(t, tt.expected, actual, 1e-14)
			} else {
				assert.InDelta(t, tt.expected, actual, 1e-15)
			}
		})
	}
}

func TestInverseNormCDF(t *testing.T) {
	// Midrange: p=0.5 should give exactly 0
	t.Run("Midrange", func(t *testing.T) {
		assert.InDelta(t, 0.0, InverseNormCDF(0.5), 1e-15)
	})

	// Round-trip tests covering lower tail, midrange, and upper tail
	for _, p := range []float64{0.000001, 0.01, 0.5, 0.99, 0.999999} {
		t.Run(fmt.Sprintf("RoundTrip_p=%v", p), func(t *testing.T) {
			z := InverseNormCDF(p)
			assert.InEpsilon(t, p, NormCDF(z), 1e-13)
		})
	}
}

// TestInverseNormCDFLowProbBranches tests all branches of inverseNormCDFLowProb
// via round-trip verification. Branches are determined by r = sqrt(-log(p)):
//   - Branch I:   r < 2.05   (p > ~0.015)
//   - Branch II:  r < 3.41   (p > ~9e-6)
//   - Branch III: r < 6.7    (p > ~2e-20)
//   - Branch IV:  r < 12.9   (p > ~1e-73)
//   - Branch V:   r >= 12.9  (p <= ~1e-73)
func TestInverseNormCDFLowProbBranches(t *testing.T) {
	testCases := []struct {
		name      string
		p         float64
		tolerance float64
	}{
		{"Branch I (r < 2.05)", 0.05, 1e-10},
		{"Branch II (r < 3.41)", 1e-4, 1e-10},
		{"Branch III (r < 6.7)", 1e-10, 1e-10},
		{"Branch IV (r < 12.9)", 1e-50, 1e-5},
		{"Branch V (r >= 12.9)", 1e-100, 1e-5},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			z := InverseNormCDF(tc.p)
			pRoundTrip := NormCDF(z)
			assert.InEpsilon(t, tc.p, pRoundTrip, tc.tolerance)
		})
	}
}

var result float64

func BenchmarkCDFComparison(b *testing.B) {
	z := -12.0
	b.Run("Custom", func(b *testing.B) {
		for b.Loop() {
			result = NormCDF(z)
		}
	})
	b.Run("MathErfc", func(b *testing.B) {
		for b.Loop() {
			result = 0.5 * math.Erfc(-z*invSqrtTwo)
		}
	})
}

func BenchmarkInverseNormCDF(b *testing.B) {
	p := 0.00001
	for b.Loop() {
		result = InverseNormCDF(p)
	}
}

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
