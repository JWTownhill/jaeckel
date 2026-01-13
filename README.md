# jaeckel

A Go implementation of Peter Jäckel's "Let's Be Rational" algorithm for Black implied volatility calculation.

## Overview

This package provides machine-precision (approximately 16 significant digits) implied volatility calculations across all input ranges, including extreme moneyness ratios and near-zero volatilities. The implementation is based on the 2024 version of Peter Jäckel's algorithm.

## Features

- **Implied volatility calculation** with full machine-precision accuracy
- **Black pricing** (normalized and standard forms)
- **Greeks**: Vega and Volga
- **Numerically stable** across extreme parameter ranges
- **Zero allocations**: All operations are allocation-free

## Installation

```bash
go get github.com/JWTownhill/jaeckel
```

## Usage

```go
package main

import (
    "fmt"

    "github.com/JWTownhill/jaeckel"
)

func main() {
    F := 100.0              // Forward price
    K := 105.0              // Strike price
    sigma := 0.20           // 20% volatility
    T := 1.0                // Time to expiry in years
    optionType := jaeckel.Call  // Call or Put

    price := jaeckel.Black(F, K, sigma, T, optionType)
    fmt.Printf("Option price: %.6f\n", price)

    // Calculate implied volatility from option price
    iv := jaeckel.ImpliedBlackVolatility(price, F, K, T, optionType)
    fmt.Printf("Implied volatility: %.6f\n", iv)

    // Calculate vega
    vega := jaeckel.Vega(F, K, sigma, T)
    fmt.Printf("Vega: %.6f\n", vega)
}
```

## API Reference

### Option Pricing

- `Black(F, K, sigma, T float64, optionType OptionType) float64` - Black 1976 option price
- `NormalizedBlack(x, s float64, optionType OptionType) float64` - Normalized Black price where x=ln(F/K) and s=σ√T

### Implied Volatility

- `ImpliedBlackVolatility(price, F, K, T float64, optionType OptionType) float64` - Compute implied volatility from option price
- `NormalizedImpliedBlackVolatility(beta, x float64, optionType OptionType) float64` - Normalized implied volatility

### Greeks

- `Vega(F, K, sigma, T float64) float64` - ∂Price/∂σ
- `NormalizedVega(x, s float64) float64` - Normalized vega (∂b/∂s)
- `Volga(F, K, sigma, T float64) float64` - ∂²Price/∂σ²
- `NormalizedVolga(x, s float64) float64` - Normalized volga (∂²b/∂s²)

### Utility Functions

- `NormCDF(z float64) float64` - Standard normal CDF
- `NormPDF(z float64) float64` - Standard normal PDF
- `InverseNormCDF(p float64) float64` - Inverse normal CDF
- `Erf(x float64) float64` - Error function
- `Erfc(x float64) float64` - Complementary error function
- `Erfcx(x float64) float64` - Scaled complementary error function
- `Erfinv(x float64) float64` - Inverse error function

## Algorithm Details

The algorithm uses a multi-branch structure based on the normalized price level:

1. **ATM branch**: Direct formula using inverse error function
2. **Lower middle branch**: Rational cubic interpolation with Householder refinement
3. **Upper middle branch**: Rational cubic interpolation with Householder refinement
4. **Lowest branch**: Logarithmic objective with Householder(3/4) iterations
5. **Highest branch**: Complementary objective with Householder(3/4) iterations

Key optimizations include:

- Remez-optimized minimax rational approximations for boundary functions
- Region-specific expansions (asymptotic for deep tails, small-t for near-zero volatility)
- Householder iterations of order 3-4 for rapid convergence
- Careful numerical treatment to avoid subtractive cancellation

## Benchmarks

Performance on Apple M1 (arm64):

| Operation                    | Time       | Allocations |
| ---------------------------- | ---------- | ----------- |
| Black (ATM)                  | 3.5 ns     | 0           |
| Black (OTM)                  | 45-55 ns   | 0           |
| NormalizedBlack (ATM)        | 3.5 ns     | 0           |
| NormalizedBlack (Standard)   | 34 ns      | 0           |
| NormalizedBlack (Region I)   | 59-64 ns   | 0           |
| ImpliedVol (ATM)             | 18 ns      | 0           |
| ImpliedVol (Middle branches) | 252-282 ns | 0           |
| ImpliedVol (Lowest branch)   | 421 ns     | 0           |
| Vega                         | 15-23 ns   | 0           |
| Erf/Erfc                     | 16 ns      | 0           |
| Erfinv                       | 4.6 ns     | 0           |
| NormCDF                      | 3.5-21 ns  | 0           |
| InverseNormCDF               | 3.2-11 ns  | 0           |

Run benchmarks with:

```bash
go test -bench=. -benchmem ./...
```

## Credits and Attribution

This implementation is based on the work of **Peter Jäckel**, whose algorithm achieves full machine precision for implied volatility calculations.

### References

- Source code and paper available at: [www.jaeckel.org/LetsBeRational.7z](http://www.jaeckel.org/LetsBeRational.7z)

## License

MIT License - See [LICENSE](LICENSE) file for details.

Note: The original C++ implementation by Peter Jäckel is available under its own license terms. This Go implementation is an independent port that follows the published algorithm.
