package jaeckel

import "math"

// We use the full precision numbers from Jaeckel
const (
	invSqrtTwo                  = 0.70710678118654752440084436210485
	codyThreshold               = 0.46875
	twoPi                       = 6.283185307179586476925286766559005768394338798750
	sqrtPiOverTwo               = 1.253314137315500251207882642405522626503493370305
	sqrtThree                   = 1.732050807568877293527446341505872366942805253810
	sqrtOneOverThree            = 0.577350269189625764509148780501957455647601751270
	twoPiOverSqrtTwentySeven    = 1.209199576156145233729385505094770488189377498728
	sqrtThreeOverThirdRootTwoPi = 0.938643487427383566075051356115075878414688769574
	piOverSix                   = 0.523598775598298873077107230546583814032861566563
	epsilon                     = 1.0 / (1 << 52)
	oneOverSqrtPi               = 1.0 / math.SqrtPi
	twoOverSqrtPi               = 2.0 / math.SqrtPi
	invSqrtTwoPi                = 0.398942280401432677939946
)

var (
	sqrtTwoPi = math.Sqrt(2.0 * math.Pi)
)
