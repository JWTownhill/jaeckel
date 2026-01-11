package jaeckel

import "math"

// Sentinel values for implied volatility results
const (
	volatilityValueToSignalPriceIsBelowIntrinsic = -math.MaxFloat64
	volatilityValueToSignalPriceIsAboveMaximum   = math.MaxFloat64
)

// Default maximum iterations for implied volatility calculation
const impliedVolatilityMaxIterations = 2

// letsBerational computes the normalized implied Black volatility using
// Peter Jäckel's "Let's Be Rational" algorithm.
//
// Parameters:
//   - beta: The normalized option price (must be positive)
//   - thetaX: θ·x where θ=±1 (call/put) and x=ln(F/K). Must satisfy θx ≤ 0.
//   - n: Maximum number of iterations (typically 2 is sufficient)
//
// Returns the normalized implied volatility s = σ√T.
//
// Special return values:
//   - volatilityValueToSignalPriceIsBelowIntrinsic: if price is below intrinsic
//   - volatilityValueToSignalPriceIsAboveMaximum: if price is above maximum
//
// Reference: Peter Jäckel, "Let's Be Rational", Wilmott Magazine, January 2015.
func letsBerational(beta, thetaX float64, n int) float64 {
	// Input validation
	if beta <= 0 {
		if beta == 0 {
			return 0
		}
		return volatilityValueToSignalPriceIsBelowIntrinsic
	}

	bMax := math.Exp(0.5 * thetaX)
	if beta >= bMax {
		return volatilityValueToSignalPriceIsAboveMaximum
	}

	// Specialise for exact ATM case
	if thetaX == 0 {
		return impliedNormalizedVolatilityATM(beta)
	}

	// Critical point calculation
	sqrtAx := math.Sqrt(-thetaX)
	sC := math.Sqrt2 * sqrtAx
	ome := oneMinusErfcx(sqrtAx)
	bC := 0.5 * bMax * ome

	var s float64
	ds := -math.MaxFloat64 // Initialize to -MaxFloat64 so loop condition is true initially

	// Four branches based on price level
	if beta < bC {
		// LOWER HALF: s < s_c
		sL := sC - sqrtPiOverTwo*ome

		bL := bLOverBMax(sC) * bMax

		if beta < bL {
			// LOWEST BRANCH: s < s_l
			s = lowestBranchInitialGuess(thetaX, sL, bL, beta)

			// Iterate using logarithmic objective function
			lnBeta := math.Log(beta)
			for i := 0; i < n && math.Abs(ds) > epsilon*s; i++ {
				bx, lnVega := scaledNormalizedBlackAndLnVega(thetaX, s)
				lnB := math.Log(bx) + lnVega
				bpob := 1 / bx // b'/b in scaled form

				h := thetaX / s
				x2OverS3 := h * h / s
				bH2 := x2OverS3 - s/4
				lambda := 1 / lnB
				otLambda := 1 + 2*lambda
				h2 := bH2 - bpob*otLambda
				c := 3 * (x2OverS3 / s) // (h/s)²

				bH3 := bH2*bH2 - c - 0.25
				sqBpob := bpob * bpob
				bppob := bH2 * bpob
				mu := 6 * lambda * (1 + lambda)
				h3 := bH3 + sqBpob*(2+mu) - bppob*3*otLambda

				nu := (lnBeta - lnB) * lnB / lnBeta / bpob

				if thetaX < -190 {
					// Use Householder(4) for extreme moneyness
					bH4 := bH2*(bH3-0.5) - (bH2-2/s)*2*c
					h4 := bH4 - bpob*(sqBpob*(6+lambda*(22+lambda*(36+lambda*24)))-bppob*(12+6*mu)) -
						bppob*bH2*3*otLambda - bH3*bpob*4*otLambda
					ds = nu * householder4Factor(nu, h2, h3, h4)
				} else {
					ds = nu * householder3Factor(nu, h2, h3)
				}
				s += ds
			}
			return s

		} else {
			// LOWER MIDDLE: s_l ≤ s < s_c
			invVC := sqrtTwoPi / bMax
			invVL := invNormalizedVega(thetaX, sL)
			rLM := convexControlParameterToFitSecondDerivativeAtRightSide(bL, bC, sL, sC, invVL, invVC, 0.0, false)
			s = interpolate(beta, bL, bC, sL, sC, invVL, invVC, rLM)
		}
	} else {
		// UPPER HALF: s_c ≤ s
		sU := sC + sqrtPiOverTwo*(2-ome)
		bU := bUOverBMax(sC) * bMax

		if beta <= bU {
			// UPPER MIDDLE: s_c ≤ s ≤ s_u
			invVC := sqrtTwoPi / bMax
			invVU := invNormalizedVega(thetaX, sU)
			rUM := convexControlParameterToFitSecondDerivativeAtLeftSide(bC, bU, sC, sU, invVC, invVU, 0.0, false)
			s = interpolate(beta, bC, bU, sC, sU, invVC, invVU, rUM)
		} else {
			// HIGHEST BRANCH: s_u < s and β > b_max/2
			s = highestBranchInitialGuess(thetaX, sU, bU, bMax, beta)

			if beta > 0.5*bMax {
				// Use complementary Black function for better accuracy
				betaBar := bMax - beta
				for i := 0; i < n && math.Abs(ds) > epsilon*s; i++ {
					h := thetaX / s
					t := s / 2
					// gp = b'/b̄ computed directly without subtractive cancellation
					gp := (2 / sqrtTwoPi) / (erfcx((t+h)*invSqrtTwo) + erfcx((t-h)*invSqrtTwo))
					bBar := normalizedVega(thetaX, s) / gp

					g := math.Log(betaBar / bBar)
					x2OverS3 := (h * h) / s
					bH2 := x2OverS3 - s/4
					c := 3 * (x2OverS3 / s)
					bH3 := bH2*bH2 - c - 0.25

					nu := -g / gp
					h2 := bH2 + gp
					h3 := bH3 + gp*(2*gp+3*bH2)

					if thetaX < -580 {
						// Use Householder(4) for extreme moneyness
						bH4 := bH2*(bH3-0.5) - (bH2-2/s)*2*c
						h4 := bH4 + gp*(6*gp*(gp+2*bH2)+3*bH2*bH2+4*bH3)
						ds = nu * householder4Factor(nu, h2, h3, h4)
					} else {
						ds = nu * householder3Factor(nu, h2, h3)
					}
					s += ds
				}
				return s
			}
		}
	}

	// MIDDLE BRANCHES (ITERATION): s_l ≤ s and (s < s_u or β ≤ b_max/2)
	// Simple objective function: g(s) = b(θx,s) - β
	for i := 0; i < n && math.Abs(ds) > epsilon*s; i++ {
		b := normalizedBlack(thetaX, s)
		invBp := invNormalizedVega(thetaX, s)
		nu := (beta - b) * invBp
		h := thetaX / s
		x2OverS3 := (h * h) / s
		h2 := x2OverS3 - s*0.25
		h3 := h2*h2 - 3*(x2OverS3/s) - 0.25
		ds = nu * householder3Factor(nu, h2, h3)
		s += ds
	}
	return s
}

// lowestBranchInitialGuess computes the initial guess for the lowest branch
// using rational cubic interpolation.
func lowestBranchInitialGuess(thetaX, sL, bL, beta float64) float64 {
	fLowerMapL, dFLowerMapLDBeta, d2FLowerMapLDBeta2 := computeFLowerMapAndDerivatives(thetaX, sL)

	rLL := convexControlParameterToFitSecondDerivativeAtRightSide(0, bL, 0, fLowerMapL, 1, dFLowerMapLDBeta, d2FLowerMapLDBeta2, true)
	f := interpolate(beta, 0, bL, 0, fLowerMapL, 1, dFLowerMapLDBeta, rLL)

	if f <= 0 {
		// Switch to quadratic interpolation for extreme values
		t := beta / bL
		f = (fLowerMapL*t + bL*(1-t)) * t
	}

	return inverseFLowerMap(thetaX, f)
}

// highestBranchInitialGuess computes the initial guess for the highest branch
// using rational cubic interpolation.
func highestBranchInitialGuess(thetaX, sU, bU, bMax, beta float64) float64 {
	fUpperMapH, dFUpperMapHDBeta, d2FUpperMapHDBeta2 := computeFUpperMapAndDerivatives(thetaX, sU)

	sqrtDBLMax := math.Sqrt(math.MaxFloat64)
	if d2FUpperMapHDBeta2 > -sqrtDBLMax && d2FUpperMapHDBeta2 < sqrtDBLMax {
		rUU := convexControlParameterToFitSecondDerivativeAtLeftSide(bU, bMax, fUpperMapH, 0, dFUpperMapHDBeta, -0.5, d2FUpperMapHDBeta2, true)
		f := interpolate(beta, bU, bMax, fUpperMapH, 0, dFUpperMapHDBeta, -0.5, rUU)
		if f > 0 {
			return inverseFUpperMap(f)
		}
	}

	// Fallback to quadratic interpolation
	h := bMax - bU
	t := (beta - bU) / h
	f := (fUpperMapH*(1-t) + 0.5*h*t) * (1 - t)
	return inverseFUpperMap(f)
}

// scaledNormalizedBlackAndLnVega returns both the scaled normalized Black price
// and the natural log of vega for efficiency in the iteration.
func scaledNormalizedBlackAndLnVega(thetaX, s float64) (bx, lnVega float64) {
	h := thetaX / s
	t := 0.5 * s
	lnVega = -0.5*math.Log(twoPi) - 0.5*(h*h+t*t)

	if isRegionI(thetaX, s) {
		return asymptoticExpansionOfScaledNormalizedBlack(h, t), lnVega
	}
	if isRegionII(thetaX, s) {
		return smallTExpansionOfScaledNormalizedBlack(h, t), lnVega
	}

	// Region III/IV: standard computation
	b := normalizedBlackStandard(thetaX, s)
	return b * math.Exp(-lnVega), lnVega
}

// invNormalizedVega returns 1/vega for numerical stability.
func invNormalizedVega(x, s float64) float64 {
	ax := math.Abs(x)
	if ax <= 0 {
		return sqrtTwoPi * math.Exp(0.125*s*s)
	}
	if s <= 0 || s <= ax*math.Sqrt(math.SmallestNonzeroFloat64) {
		return math.MaxFloat64
	}
	h := x / s
	t := 0.5 * s
	return sqrtTwoPi * math.Exp(0.5*(h*h+t*t))
}

// computeFLowerMapAndDerivatives computes f_lower_map and its first two derivatives.
// This is used for rational cubic interpolation in the lowest branch.
func computeFLowerMapAndDerivatives(x, s float64) (f, fp, fpp float64) {
	ax := math.Abs(x)
	z := sqrtOneOverThree * ax / s
	y := z * z
	s2 := s * s

	phi := NormCDF(-z)
	phiFunc := NormPDF(z)

	fpp = piOverSix * y / (s2 * s) * phi * (8*sqrtThree*s*ax + (3*s2*(s2-8)-8*x*x)*phi/phiFunc) * math.Exp(2*y+0.25*s2)

	phi2 := phi * phi
	fp = twoPi * y * phi2 * math.Exp(y+0.125*s*s)
	f = twoPiOverSqrtTwentySeven * ax * (phi2 * phi)

	return f, fp, fpp
}

// inverseFLowerMap inverts the lower mapping function.
func inverseFLowerMap(x, f float64) float64 {
	if f <= 0 {
		return math.SmallestNonzeroFloat64
	}
	ax := math.Abs(x)
	// Use math.Abs on result as inverseNormCDF can return negative for small probabilities
	return math.Abs(x / (sqrtThree * InverseNormCDF(sqrtThreeOverThirdRootTwoPi*math.Cbrt(f)/math.Cbrt(ax))))
}

// computeFUpperMapAndDerivatives computes f_upper_map and its first two derivatives.
// This is used for rational cubic interpolation in the highest branch.
func computeFUpperMapAndDerivatives(x, s float64) (f, fp, fpp float64) {
	f = NormCDF(-0.5 * s)
	w := (x / s) * (x / s)
	fp = -0.5 * math.Exp(0.5*w)
	fpp = sqrtPiOverTwo * math.Exp(w+0.125*s*s) * w / s
	return f, fp, fpp
}

// inverseFUpperMap inverts the upper mapping function.
func inverseFUpperMap(f float64) float64 {
	return -2 * InverseNormCDF(f)
}

// bLOverBMax computes b_l(x)/b_max(x) using Remez-optimized rational approximations.
// This is a univariate function of s_c = √(2|x|).
func bLOverBMax(sC float64) float64 {
	// Four branches based on s_c value
	if sC < 2.6267851073127395 {
		if sC < 0.7099295739719539 {
			// Branch I: smallest |x|
			g := (8.0741072372882856924e-2 + sC*(9.8078911786358897272e-2+sC*(3.9760631445677058375e-2+sC*(5.9716928459589189876e-3+sC*(-6.4036399341479799981e-6+4.5425102093616062245e-7*sC))))) /
				(1 + sC*(1.8594977672287664353+sC*(1.3658801475711790419+sC*(4.6132707108655653215e-1+6.1254597049831720643e-2*sC))))
			return (sC * sC) * (0.07560996640296361767172 + sC*(sC*g-0.09672719281339436290858))
		}
		// Branch II
		return (1.9795737927598581235e-9 + sC*(-2.7081288564685588037e-8+sC*(7.5610142272549044609e-2+sC*(6.917130174466834016e-2+sC*(2.9537058950963019803e-2+sC*(6.5849252702302307774e-3+6.9711400639834715731e-4*sC)))))) /
			(1 + sC*(2.1941448525586579756+sC*(2.1297103549995181357+sC*(1.1571483187179784072+sC*(3.7831622253060456794e-1+sC*(7.1714862448829349869e-2+6.6361975827861200167e-3*sC))))))
	}
	if sC < 7.348469228349534 {
		// Branch III
		return (-9.3325115354837883291e-5 + sC*(5.3118033972794648837e-4+sC*(7.4114855448345002595e-2+sC*(7.4039658186822817454e-2+sC*(3.9225177407687604785e-2+sC*(1.0022913378254090083e-2+1.7012579407246055469e-3*sC)))))) /
			(1 + sC*(2.2217238132228132256+sC*(2.3441816707087403282+sC*(1.3912323646271141826+sC*(5.3231258443501838354e-1+sC*(1.1744005919716101572e-1+1.6195405895930935811e-2*sC))))))
	}
	// Branch IV: largest |x|
	return (1.4500072297240603183e-3 + sC*(-1.5116692485011195757e-3+sC*(7.1682178310936334831e-2+sC*(3.921610857820463493e-2+sC*(2.9342405658628443931e-2+sC*(5.1832526171631521426e-3+1.6930208078421474854e-3*sC)))))) /
		(1 + sC*(1.6176313502305414664+sC*(1.6823159175281531664+sC*(8.4878307567372222113e-1+sC*(3.7543742137375791321e-1+sC*(7.126137099644302999e-2+1.6116992546788676159e-2*sC))))))
}

// bUOverBMax computes b_u(x)/b_max(x) using Remez-optimized rational approximations.
// This is a univariate function of s_c = √(2|x|).
func bUOverBMax(sC float64) float64 {
	// Four branches based on s_c value
	if sC < 1.7888543819998317 {
		if sC < 0.7745966692414833 {
			// Branch I: smallest |x|
			g := (-6.063099881233561706e-2 + sC*(-8.1011946637120604985e-2+sC*(-4.2505564862438753828e-2+sC*(-8.9880000946868691788e-3+sC*(-7.5603072110443268356e-6+4.3879556621540147458e-7*sC))))) /
				(1 + sC*(1.8400371530721828756+sC*(1.5709283443886143691+sC*(6.8913245453611400484e-1+1.4703173061720980923e-1*sC))))
			return 0.7899085945560627246288 + (sC*sC)*(0.0614616805805147403487+sC*g)
		}
		// Branch II
		return (7.8990944435755287611e-1 + sC*(-1.2655410534988972886+sC*(-2.8803040699221003256+sC*(-2.6936198689113258727+sC*(-1.1213067281643205754+sC*(-2.1277793801691629892e-1+5.1486445905299802703e-6*sC)))))) /
			(1 + sC*(-1.6021222722060444448+sC*(-3.7242680976480704555+sC*(-3.2083117718907365085+sC*(-1.2922333835930958583-2.3762328334050001161e-1*sC)))))
	}
	if sC < 6.164414002968976 {
		// Branch III
		return (7.8990640048967596475e-1 + sC*(1.5993699253596663678+sC*(1.6481729039140370242+sC*(9.8227188109869200166e-1+sC*(3.6313557966186936883e-1+sC*(7.8277036261179606301e-2+9.3404307364538726214e-3*sC)))))) /
			(1 + sC*(2.0247407005640401446+sC*(2.0087454279103740489+sC*(1.1627561803056961973+sC*(4.2004672123723823581e-1+sC*(8.9130862793887234546e-2+1.0436767768858021717e-2*sC))))))
	}
	// Branch IV: largest |x|
	return (7.91133825948419359e-1 + sC*(1.24653733210880042+sC*(1.32747426980537386+sC*(6.95009705717846778e-1+sC*(3.05965944268228457e-1+sC*(6.02200363391352887e-2+1.29050244454344842e-2*sC)))))) /
		(1 + sC*(1.58117486714634672+sC*(1.60144713247629644+sC*(8.30040185836882436e-1+sC*(3.53071863813401531e-1+sC*(6.95901684131758475e-2+1.44197580643890011e-2*sC))))))
}

// ImpliedBlackVolatility computes the implied Black volatility from an option price.
//
// Parameters:
//   - price: The undiscounted option price
//   - F: Forward price
//   - K: Strike price
//   - T: Time to expiry (in years)
//   - q: Option type indicator: +1 for call, -1 for put
//
// Returns the implied volatility (annualized).
//
// Special return values:
//   - -Inf: price is below intrinsic value (arbitrage)
//   - +Inf: price is above maximum possible value
func ImpliedBlackVolatility(price, F, K, T, q float64) float64 {
	maxPrice := F
	if q < 0 {
		maxPrice = K
	}
	if price >= maxPrice {
		return volatilityValueToSignalPriceIsAboveMaximum
	}

	// Map ITM to OTM
	x := math.Log(F / K)
	intrinsic := 0.0
	if q < 0 {
		if K > F {
			intrinsic = K - F
		}
	} else {
		if F > K {
			intrinsic = F - K
		}
	}

	sqrtFK := math.Sqrt(F) * math.Sqrt(K)
	beta := (price - intrinsic) / sqrtFK

	normalizedVol := letsBerational(beta, -math.Abs(x), impliedVolatilityMaxIterations)
	return normalizedVol / math.Sqrt(T) // annualize
}

// NormalizedImpliedBlackVolatility computes the normalized implied volatility.
//
// Parameters:
//   - beta: The normalized option price
//   - x: log(F/K), the log-moneyness
//   - q: Option type indicator: +1 for call, -1 for put
//
// Returns the normalized implied volatility s = σ√T.
//
// Special return values:
//   - -Inf: price is below intrinsic value
//   - +Inf: price is above maximum possible value
func NormalizedImpliedBlackVolatility(beta, x, q float64) float64 {
	thetaX := x
	if q < 0 {
		thetaX = -x
	}
	return letsBerational(beta-normalizedIntrinsic(thetaX), -math.Abs(x), impliedVolatilityMaxIterations)
}
