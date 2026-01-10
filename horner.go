package letsberational

// evalRational evaluates a rational function P(x)/Q(x) using Horner's Method.
// num and den are the coefficients of P and Q ordered highest to lowest degree coefficient
func evalRational(x float64, num, den []float64) float64 {
	p := num[0]
	for i := 1; i < len(num); i++ {
		p = p*x + num[i]
	}
	q := den[0]
	for i := 1; i < len(den); i++ {
		q = q*x + den[i]
	}
	return p / q
}
