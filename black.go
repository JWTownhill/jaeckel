package jaeckel

import "math"

// Region boundary tuning parameters (Jaeckel 2024 values).
const (
	eta = -13.0
	tau = 0.21022410381342863 // 2 * math.Pow(epsilon, 1.0/16.0)
)

// Black computes the undiscounted Black 1976 option price.
//
// Parameters:
//   - F: Forward price
//   - K: Strike price
//   - sigma: Volatility (annualized)
//   - T: Time to expiry (in years)
//   - optionType: Call or Put
//
// Returns the undiscounted option price.
func Black(F, K, sigma, T float64, optionType OptionType) float64 {
	s := sigma * math.Sqrt(T)

	// Specialisation for ATM: b(s) = 1 - 2·Φ(-s/2) = erf(s/√8)
	if K == F {
		return F * Erf(s*0.5*invSqrtTwo)
	}

	// Map ITM to OTM using put-call parity
	x := math.Log(F / K)
	intrinsic := 0.0
	if optionType == Put {
		if K > F {
			intrinsic = K - F
		}
	} else {
		if F > K {
			intrinsic = F - K
		}
	}

	if s <= 0 {
		return intrinsic
	}

	// Compute using geometric mean √(F·K) for numerical stability
	return intrinsic + math.Sqrt(F)*math.Sqrt(K)*normalizedBlack(-math.Abs(x), s)
}

// NormalizedBlack computes the normalized Black option price.
//
// The normalized price is defined as:
//
//	b(x,s,θ) = θ·[exp(θx/2)·Φ(θ·(x/s+s/2)) - exp(-θx/2)·Φ(θ·(x/s-s/2))]
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalized volatility
//   - optionType: Call or Put
//
// Returns the normalized option price.
func NormalizedBlack(x, s float64, optionType OptionType) float64 {
	// Specialisation for ATM: b(s) = erf(s/√8)
	if x == 0 {
		return Erf(s * 0.5 * invSqrtTwo)
	}

	thetaX := x
	if optionType == Put {
		thetaX = -x
	}

	return normalizedIntrinsic(thetaX) + normalizedBlack(-math.Abs(x), s)
}

// Vega computes the sensitivity of option price to volatility changes.
//
// Parameters:
//   - F: Forward price
//   - K: Strike price
//   - sigma: Volatility (annualized)
//   - T: Time to expiry (in years)
//
// Returns ∂Price/∂σ.
func Vega(F, K, sigma, T float64) float64 {
	x := math.Log(F / K)
	s := sigma * math.Sqrt(T)
	return math.Sqrt(F) * math.Sqrt(K) * NormalizedVega(x, s) * math.Sqrt(T)
}

// NormalizedVega computes the normalized vega (∂b/∂s).
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalized volatility
//
// Returns ∂b/∂s where b is the normalized Black price.
func NormalizedVega(x, s float64) float64 {
	ax := math.Abs(x)
	// ATM case (x = 0)
	if ax <= 0 {
		return invSqrtTwoPi * math.Exp(-0.125*s*s)
	}
	// Safeguard against s being too small or zero
	if s <= 0 || s <= ax*math.Sqrt(math.SmallestNonzeroFloat64) {
		return 0
	}
	h := x / s
	t := 0.5 * s
	return invSqrtTwoPi * math.Exp(-0.5*(h*h+t*t))
}

// Volga computes the second derivative of option price with respect to volatility.
//
// Parameters:
//   - F: Forward price
//   - K: Strike price
//   - sigma: Volatility (annualized)
//   - T: Time to expiry (in years)
//
// Returns ∂²Price/∂σ².
func Volga(F, K, sigma, T float64) float64 {
	x := math.Log(F / K)
	s := sigma * math.Sqrt(T)
	return math.Sqrt(F) * math.Sqrt(K) * NormalizedVolga(x, s) * T
}

// NormalizedVolga computes the normalized volga (∂²b/∂s²).
//
// Parameters:
//   - x: log(F/K), the log-moneyness
//   - s: σ√T, the normalized volatility
//
// Returns ∂²b/∂s² where b is the normalized Black price.
func NormalizedVolga(x, s float64) float64 {
	ax := math.Abs(x)
	if s <= 0 || (ax > 0 && s <= ax*math.Sqrt(math.SmallestNonzeroFloat64)) {
		return 0
	}
	h := x / s
	t := 0.5 * s
	h2 := h * h
	t2 := t * t
	return invSqrtTwoPi * math.Exp(-0.5*(h2+t2)) * (h2 - t2) / s
}

func normalizedBlack(thetaX, s float64) float64 {
	if s <= 0 {
		return 0
	}

	// Exact ATM case
	if thetaX == 0 {
		return 2.0*NormCDF(0.5*s) - 1.0
	}

	// Region I (Deep Tails)
	if isRegionI(thetaX, s) {
		return asymptoticExpansionOfScaledNormalizedBlack(thetaX/s, 0.5*s) * NormalizedVega(thetaX, s)
	}

	// Region II (Small s)
	if isRegionII(thetaX, s) {
		return smallTExpansionOfScaledNormalizedBlack(thetaX/s, 0.5*s) * NormalizedVega(thetaX, s)
	}

	// Regions III/IV (standard computation)
	return normalizedBlackStandard(thetaX, s)
}

// smallTExpansionOfScaledNormalizedBlack calculates (normalizedBlack / normalizedVega)
// for small values of t (s/2).
func smallTExpansionOfScaledNormalizedBlack(h, t float64) float64 {
	yH := sqrtPiOverTwo * Erfcx(-h*invSqrtTwo)
	a := 1.0 + h*yH

	h2 := h * h
	t2 := t * t

	b0 := 2.0 * a
	b1 := (-1.0 + a*(3.0+h2)) / 3.0
	b2 := (-7.0 - h2 + a*(15.0+h2*(10.0+h2))) / 60.0
	b3 := (-57.0 + (-18.0-h2)*h2 + a*(105.0+h2*(105.0+h2*(21.0+h2)))) / 2520.0
	b4 := (-561.0 + h2*(-285.0+(-33.0-h2)*h2) + a*(945.0+h2*(1260.0+h2*(378.0+h2*(36.0+h2))))) / 181440.0
	b5 := (-6555.0 + h2*(-4680.0+h2*(-840.0+(-52.0-h2)*h2)) + a*(10395.0+h2*(17325.0+h2*(6930.0+h2*(990.0+h2*(55.0+h2)))))) / 19958400.0
	b6 := (-89055.0 + h2*(-82845.0+h2*(-20370.0+h2*(-1926.0+(-75.0-h2)*h2))) + a*(135135.0+h2*(270270.0+h2*(135135.0+h2*(25740.0+h2*(2145.0+h2*(78.0+h2))))))) / 3113510400.0

	return t * (b0 + t2*(b1+t2*(b2+t2*(b3+t2*(b4+t2*(b5+b6*t2))))))
}

// normalizedBlackStandard computes the normalized Black price using Cody's erfc/erfcx
// functions with optimal branch selection to minimize exponential evaluations.
// The four branches select between erfc and erfcx based on q1 and q2 vs codyThreshold.
func normalizedBlackStandard(thetaX, s float64) float64 {
	h := thetaX / s
	t := 0.5 * s

	q1 := -invSqrtTwo * (h + t)
	q2 := -invSqrtTwo * (h - t)

	expPos := math.Exp(0.5 * thetaX)
	expNeg := math.Exp(-0.5 * thetaX)
	expSq := math.Exp(-0.5 * (h*h + t*t))

	var twoB float64
	if q1 < codyThreshold {
		if q2 < codyThreshold {
			twoB = expPos*Erfc(q1) - expNeg*Erfc(q2)
		} else {
			twoB = expPos*Erfc(q1) - expSq*Erfcx(q2)
		}
	} else {
		// NOTE: The branch "q1 >= codyThreshold && q2 < codyThreshold" appears
		// unreachable for s > 0.
		// Kept for consistency with Jäckel's original C++ (which includes
		// all four branches), and as a defensive safeguard.
		if q2 < codyThreshold {
			twoB = expSq*Erfcx(q1) - expNeg*Erfc(q2)
		} else {
			twoB = expSq * (Erfcx(q1) - Erfcx(q2))
		}
	}

	return math.Max(0.5*twoB, 0.0)
}

// isRegionI checks if we should use the asymptotic expansion.
// Region I: h < η and t < (τ+½) + (|h|-|η|)
// where h = θx/s and t = s/2.
// Rewritten to avoid division by s:
//
//	s·(s/2-(τ+½+η)) + θx < 0
func isRegionI(thetaX, s float64) bool {
	return thetaX < s*eta && s*(0.5*s-(tau+0.5+eta))+thetaX < 0
}

// isRegionII checks if we should use the small-t expansion.
// Region II: t < τ + (½/|η|)·|h|
// where h = θx/s and t = s/2.
// Rewritten to avoid division by s:
//
//	s·(s-2·τ) - θx/η < 0
func isRegionII(thetaX, s float64) bool {
	return s*(s-2*tau)-thetaX/eta < 0
}

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

// asymptoticExpansionOfScaledNormalizedBlack computes the scaled normalized Black function
// using an asymptotic expansion for Region I (deep tails).
//
// The scaled normalized Black function is defined as:
//
//	bx := b / (∂b/∂s)
//
// where b is the normalized Black price and ∂b/∂s is the vega.
//
// For h < η (typically η = -13) and t < (τ+½) + (|h|-|η|), we use an asymptotic
// expansion based on Abramowitz & Stegun (26.2.12).
//
// The expansion is in terms of:
//
//	h := θx/s (normalized log-moneyness)
//	t := s/2 (half of normalized volatility)
func asymptoticExpansionOfScaledNormalizedBlack(h, t float64) float64 {
	// Note that e := (t/h)² ∈ (0,1).
	e := (t / h) * (t / h)
	r := (h + t) * (h - t)
	q := (h / r) * (h / r)

	// Asymptotic expansion coefficients A0-A16.
	a := [17]float64{
		2,                              // a0
		-6 - 2*e,                       // a1
		30 + e*(60+6*e),                // a2
		-210 + e*(-1050+e*(-630-30*e)), // a3
		1890 + e*(17640+e*(26460+e*(7560+210*e))),                                                                                                                                                                                                              // a4
		-20790 + e*(-311850+e*(-873180+e*(-623700+e*(-103950-1890*e)))),                                                                                                                                                                                        // a5
		270270 + e*(5945940+e*(26756730+e*(35675640+e*(14864850+e*(1621620+20790*e))))),                                                                                                                                                                        // a6
		-4054050 + e*(-122972850+e*(-811620810+e*(-1739187450+e*(-1352701350+e*(-368918550+e*(-28378350-270270*e)))))),                                                                                                                                         // a7
		68918850 + e*(2756754000+e*(25086461400+e*(78843164400+e*(98553955500+e*(50172922800+e*(9648639000+e*(551350800+4054050*e))))))),                                                                                                                       // a8
		-1309458150 + e*(-66782365650+e*(-801388387800+e*(-3472683013800+e*(-6366585525300+e*(-5209024520700+e*(-1869906238200+e*(-267129462600+e*(-11785123350-68918850*e)))))))),                                                                             // a9
		27498621150 + e*(1741579339500+e*(26646163894350+e*(152263793682000+e*(384889034029500+e*(461866840835400+e*(266461638943500+e*(71056437051600+e*(7837107027750+e*(274986211500+1309458150*e))))))))),                                                  // a10
		-632468286450 + e*(-48700058056650+e*(-925301103076350+e*(-6741479465270550+e*(-22471598217568500+e*(-37180280687249700+e*(-31460237504595900+e*(-13482958930541100+e*(-2775903309229050+e*(-243500290283250+e*(-6957151150950-27498621150*e)))))))))), // a11
		15811707161250 + e*(1454677058835000+e*(33603040059088500+e*(304027505296515000+e*(1292116897510188750+e*(2819164140022230000+e*(3289024830025935000+e*(2067387036016302000+e*(684061886917158750+e*(112010133530295000+e*(80007238235925000+e*(189740485935000+632468286450*e))))))))))),                                                                                                                                                  // a12
		-426916093353750 + e*(-46249243446656250+e*(-1276479119127712500+e*(-14041270310404837500+e*(-74106704416025531200+e*(-206151377739125569000+e*(-317155965752500875000+e*(-274868503652167425000+e*(-133392067948845956000+e*(-35103175776012093800+e*(-4680423436801612500+e*(-277495460679937500+e*(-554990921359875000-15811707161250*e)))))))))))),                                                                                     // a13
		12380566707258750 + e*(1559951405114602500+e*(50698420666224581200+e*(66632210018466592500+e*(427556680951827302000+e*(1477013988742676130000+e*(2897219747149095490000+e*(3311108282456109140000+e*(2215520983114014200000+e*(855113361903654604000+e*(183238577550783129000+e*(20279368266489832500+e*(10139684133244916200+e*(173327933901622500+426916093353750*e))))))))))))),                                                         // a14
		-383797567925021250 + e*(-55650647349128081200+e*(-210359446979704147000+e*(-325556286992399275000+e*(-2495931533608394440000+e*(-10482912441155256700000+e*(-25535299536147420100000+e*(-37208579324100526400000+e*(-32831099403618111500000+e*(-17471520735258761100000+e*(-5491049373938467780000+e*(-976668860977197826000+e*(-91155760357871797100+e*(-38955453144389656900+e*(-5756963518875318750-12380566707258750*e)))))))))))))), // a15
		12665319741525701200 + e*(2093999530598915940000+e*(91088979581052843400+e*(16396016324589511800+e*(14801959181921087100+e*(74278922440185818700+e*(21997988568824261700+e*(39805884076920092600+e*(44781619586535104100+e*(31425697955463231000+e*(13617802447367400100+e*(35524702036610608900+e*(53287053054915913400+e*(42508190471157993600+e*(15704996479491869600+e*(20264511586441122000+383797567925021250*e))))))))))))))),       // a16
	}

	// Use thresholds to determine how many terms to include.
	thresholds := [12]float64{12.347, 12.958, 13.729, 14.718, 16.016, 17.769, 20.221, 23.816, 29.419, 38.93, 57.171, 99.347}

	// We always include terms a0-a4, and add more based on distance from boundary.
	distance := -h - t + tau + 0.5
	numExtraTerms := 0
	for i := range len(thresholds) {
		if distance < thresholds[i] {
			numExtraTerms = 12 - i
			break
		}
	}
	omega := 0.0
	startIdx := min(4+numExtraTerms, 16)

	for i := startIdx; i >= 2; i-- {
		omega = q * (a[i] + omega)
	}
	omega = a[0] + q*(a[1]+omega)
	return (t / r) * omega
}

// normalizedIntrinsic computes the intrinsic value of a normalized option.
// For θx > 0: intrinsic = exp(θx/2) - exp(-θx/2) = 2·sinh(θx/2)
// For θx ≤ 0: intrinsic = 0
func normalizedIntrinsic(thetaX float64) float64 {
	if thetaX <= 0 {
		return 0
	}
	return 2 * math.Sinh(0.5*thetaX) // Use sinh for numerical stability
}
