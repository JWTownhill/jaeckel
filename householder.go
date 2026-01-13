package jaeckel

// householder3Factor computes the third-order Householder correction factor.
//
// Given Newton step ν = -g/g', and the ratios h₂ = g”/g' and h₃ = g”'/g',
// the Householder(3) iteration is:
//
//	s_{n+1} = s_n + ν · (1 + ν·h₂/2) / (1 + ν·(h₂ + ν·h₃/6))
//
// This function returns the factor: (1 + ν·h₂/2) / (1 + ν·(h₂ + ν·h₃/6))
func householder3Factor(nu, h2, h3 float64) float64 {
	return (1 + 0.5*h2*nu) / (1 + nu*(h2+h3*nu*(1.0/6.0)))
}

// householder4Factor computes the fourth-order Householder correction factor.
//
// Given Newton step ν = -g/g', and the ratios h₂ = g"/g', h₃ = g"'/g', h₄ = g""/g',
// the Householder(4) iteration uses a more complex factor that provides
// fifth-order convergence accuracy.
func householder4Factor(nu, h2, h3, h4 float64) float64 {
	return (1 + nu*(h2+nu*h3*(1.0/6.0))) / (1 + nu*(1.5*h2+nu*(h2*h2*0.25+h3*(1.0/3.0)+nu*h4*(1.0/24.0))))
}
