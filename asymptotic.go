package letsberational

import "math"

// asymptoticExpansionOfScaledNormalisedBlack computes the scaled normalised Black function
// using an asymptotic expansion for Region I (deep tails).
//
// The scaled normalised Black function is defined as:
//   bx := b / (∂b/∂s)
// where b is the normalised Black price and ∂b/∂s is the vega.
//
// For h < η (typically η = -13) and t < (τ+½) + (|h|-|η|), we use an asymptotic
// expansion based on Abramowitz & Stegun (26.2.12).
//
// The expansion is in terms of:
//   h := θx/s (normalized log-moneyness)
//   t := s/2 (half of normalized volatility)
//   e := (t/h)²
//   r := (h+t)·(h-t)
//   q := (h/r)²
//
// Reference: Peter Jäckel, "Let's Be Rational", Wilmott Magazine, January 2015.
func asymptoticExpansionOfScaledNormalisedBlack(h, t float64) float64 {
	// Note that e := (t/h)² ∈ (0,1).
	e := (t / h) * (t / h)
	r := (h + t) * (h - t)
	q := (h / r) * (h / r)

	// Asymptotic expansion coefficients A0-A16.
	// These are computed as polynomials in e (= (t/h)²).
	// Coefficients follow Jaeckel's optimized formulation.
	a := [17]float64{
		2,                                                                                                                                                                                                                                                                                             // a0
		-6 - 2*e,                                                                                                                                                                                                                                                                                      // a1
		30 + e*(60+6*e),                                                                                                                                                                                                                                                                               // a2
		-210 + e*(-1050+e*(-630-30*e)),                                                                                                                                                                                                                                                                // a3
		1890 + e*(17640+e*(26460+e*(7560+210*e))),                                                                                                                                                                                                                                                     // a4
		-20790 + e*(-311850+e*(-873180+e*(-623700+e*(-103950-1890*e)))),                                                                                                                                                                                                                               // a5
		270270 + e*(5945940+e*(26756730+e*(35675640+e*(14864850+e*(1621620+20790*e))))),                                                                                                                                                                                                                // a6
		-4054050 + e*(-122972850+e*(-811620810+e*(-1739187450+e*(-1352701350+e*(-368918550+e*(-28378350-270270*e)))))),                                                                                                                                                                                 // a7
		68918850 + e*(2756754000+e*(25086461400+e*(78843164400+e*(98553955500+e*(50172922800+e*(9648639000+e*(551350800+4054050*e))))))),                                                                                                                                                                // a8
		-1309458150 + e*(-66782365650+e*(-801388387800+e*(-3472683013800+e*(-6366585525300+e*(-5209024520700+e*(-1869906238200+e*(-267129462600+e*(-11785123350-68918850*e)))))))),                                                                                                                      // a9
		27498621150 + e*(1741579339500+e*(26646163894350+e*(152263793682000+e*(384889034029500+e*(461866840835400+e*(266461638943500+e*(71056437051600+e*(7837107027750+e*(274986211500+1309458150*e))))))))),                                                                                           // a10
		-632468286450 + e*(-48700058056650+e*(-925301103076350+e*(-6741479465270550+e*(-22471598217568500+e*(-37180280687249700+e*(-31460237504595900+e*(-13482958930541100+e*(-2775903309229050+e*(-243500290283250+e*(-6957151150950-27498621150*e)))))))))),                                          // a11
		15811707161250 + e*(1454677058835000+e*(33603040059088500+e*(304027505296515000+e*(1292116897510188750+e*(2819164140022230000+e*(3289024830025935000+e*(2067387036016302000+e*(684061886917158750+e*(112010133530295000+e*(80007238235925000+e*(189740485935000+632468286450*e))))))))))),       // a12
		-426916093353750 + e*(-46249243446656250+e*(-1276479119127712500+e*(-14041270310404837500+e*(-74106704416025531200+e*(-206151377739125569000+e*(-317155965752500875000+e*(-274868503652167425000+e*(-133392067948845956000+e*(-35103175776012093800+e*(-4680423436801612500+e*(-277495460679937500+e*(-554990921359875000-15811707161250*e)))))))))))), // a13
		12380566707258750 + e*(1559951405114602500+e*(50698420666224581200+e*(66632210018466592500+e*(427556680951827302000+e*(1477013988742676130000+e*(2897219747149095490000+e*(3311108282456109140000+e*(2215520983114014200000+e*(855113361903654604000+e*(183238577550783129000+e*(20279368266489832500+e*(10139684133244916200+e*(173327933901622500+426916093353750*e))))))))))))), // a14
		-383797567925021250 + e*(-55650647349128081200+e*(-210359446979704147000+e*(-325556286992399275000+e*(-2495931533608394440000+e*(-10482912441155256700000+e*(-25535299536147420100000+e*(-37208579324100526400000+e*(-32831099403618111500000+e*(-17471520735258761100000+e*(-5491049373938467780000+e*(-976668860977197826000+e*(-91155760357871797100+e*(-38955453144389656900+e*(-5756963518875318750-12380566707258750*e)))))))))))))), // a15
		12665319741525701200 + e*(2093999530598915940000+e*(91088979581052843400+e*(16396016324589511800+e*(14801959181921087100+e*(74278922440185818700+e*(21997988568824261700+e*(39805884076920092600+e*(44781619586535104100+e*(31425697955463231000+e*(13617802447367400100+e*(35524702036610608900+e*(53287053054915913400+e*(42508190471157993600+e*(15704996479491869600+e*(20264511586441122000+383797567925021250*e))))))))))))))), // a16
	}

	// Use thresholds to determine how many terms to include.
	// This follows Jaeckel's optimized approach for accuracy vs. performance.
	thresholds := [12]float64{12.347, 12.958, 13.729, 14.718, 16.016, 17.769, 20.221, 23.816, 29.419, 38.93, 57.171, 99.347}
	distance := -h - t + tau + 0.5

	// Find which threshold bracket we're in (determines number of terms needed)
	numExtraTerms := 0
	for i := 0; i < len(thresholds); i++ {
		if distance < thresholds[i] {
			numExtraTerms = 12 - i // How many extra terms beyond base 5 (a0-a4)
			break
		}
	}

	// Build omega using Horner's method.
	// We always include terms a0-a4 (5 base terms), and add more based on distance from boundary.
	// Start from highest needed term and work down.
	omega := 0.0
	startIdx := 4 + numExtraTerms // Start from this index and go down to 0

	// Clamp startIdx to valid range
	if startIdx > 16 {
		startIdx = 16
	}

	// Build polynomial from high to low using Horner's method
	for i := startIdx; i >= 2; i-- {
		omega = q * (a[i] + omega)
	}
	// Final terms: a0 + q*(a1 + omega)
	omega = a[0] + q*(a[1]+omega)

	// bx := (t/r) * omega
	return (t / r) * omega
}

// yprimeTailExpansionRationalFunctionPart computes the rational function part
// of the Y'(h) tail expansion for h < -4.
func yprimeTailExpansionRationalFunctionPart(w float64) float64 {
	num := w * (-2.9999999999994663866 + w*(-1.7556263323542206288e2+w*(-3.4735035445495633334e3+w*(-2.7805745693864308643e4+w*(-8.3836021460741980839e4-6.6818249032616849037e4*w)))))
	den := 1 + w*(6.3520877744831739102e1+w*(1.4404389037604337538e3+w*(1.4562545638507033944e4+w*(6.6886794165651675684e4+w*(1.2569970380923908488e5+6.9286518679803751694e4*w)))))
	return num / den
}

// yprime computes Y'(h) = 1 + h·Y(h) avoiding subtractive cancellation.
// Y(h) := Φ(h)/φ(h) where Φ is the standard normal CDF and φ is the PDF.
func yprime(h float64) float64 {
	// Thresholds copied from Cody's implementation
	if h < -4 {
		// Nonlinear-Remez optimized minimax rational function of order (5,6)
		// for g(w) := (Y'(h)/h²-1)/h² with w:=1/h².
		// Relative accuracy better than 9.8E-17 on h in [-∞,-4].
		w := 1 / (h * h)
		return w * (1 + yprimeTailExpansionRationalFunctionPart(w))
	}
	if h <= -codyThreshold {
		// Remez-optimized minimax rational function of order (7,7)
		// Relative accuracy better than 1.6E-16 on h in [-4,-0.46875].
		num := 1.0000000000594317229 - h*(6.1911449879694112749e-1-h*(2.2180844736576013957e-1-h*(4.5650900351352987865e-2-h*(5.545521007735379052e-3-h*(3.0717392274913902347e-4-h*(4.2766597835908713583e-8+8.4592436406580605619e-10*h))))))
		den := 1 - h*(1.8724286369589162071-h*(1.5685497236077651429-h*(7.6576489836589035112e-1-h*(2.3677701403094640361e-1-h*(4.6762548903194957675e-2-h*(5.5290453576936595892e-3-3.0822020417927147113e-4*h))))))
		return num / den
	}
	return 1 + h*sqrtPiOverTwo*erfcx(-invSqrtTwo*h)
}

// smallTExpansionOfScaledNormalizedBlackWithYprime is an alternative implementation
// using the optimized Yprime function for better accuracy in Region II.
func smallTExpansionOfScaledNormalizedBlackWithYprime(h, t float64) float64 {
	a := yprime(h)
	h2 := h * h
	t2 := t * t

	b0 := 2 * a
	b1 := (-1 + a*(3+h2)) / 3
	b2 := (-7 - h2 + a*(15+h2*(10+h2))) / 60
	b3 := (-57 + (-18-h2)*h2 + a*(105+h2*(105+h2*(21+h2)))) / 2520
	b4 := (-561 + h2*(-285+(-33-h2)*h2) + a*(945+h2*(1260+h2*(378+h2*(36+h2))))) / 181440
	b5 := (-6555 + h2*(-4680+h2*(-840+(-52-h2)*h2)) + a*(10395+h2*(17325+h2*(6930+h2*(990+h2*(55+h2)))))) / 19958400
	b6 := (-89055 + h2*(-82845+h2*(-20370+h2*(-1926+(-75-h2)*h2))) + a*(135135+h2*(270270+h2*(135135+h2*(25740+h2*(2145+h2*(78+h2))))))) / 3113510400

	return t * (b0 + t2*(b1+t2*(b2+t2*(b3+t2*(b4+t2*(b5+b6*t2))))))
}

// complementaryNormalisedBlack computes bₘₐₓ - b(x,s) without subtractive cancellation.
// This is crucial for numerical stability in the upper branch of implied volatility.
//
// The complementary normalised Black function is:
//   b̄(x,s) = exp(θx/2)·Φ(-x/s-s/2) + exp(-θx/2)·Φ(x/s-s/2)
//          = ½ · (erfcx((t+h)/√2) + erfcx((t-h)/√2)) · exp(-½(t²+h²))
//
// where h = x/s and t = s/2.
func complementaryNormalisedBlack(h, t float64) float64 {
	sum := erfcx((t+h)*invSqrtTwo) + erfcx((t-h)*invSqrtTwo)
	return 0.5 * sum * math.Exp(-0.5*(t*t+h*h))
}

// normalisedIntrinsic computes the intrinsic value of a normalised option.
// For θx > 0: intrinsic = exp(θx/2) - exp(-θx/2) = 2·sinh(θx/2)
// For θx ≤ 0: intrinsic = 0
func normalisedIntrinsic(thetaX float64) float64 {
	if thetaX <= 0 {
		return 0
	}
	// Use sinh for numerical stability
	return 2 * math.Sinh(0.5*thetaX)
}
