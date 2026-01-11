# jaeckel

A high-precision Go implementation of Peter Jäckel's "Let's Be Rational" algorithm for Black implied volatility calculation.

## Overview

This package provides machine-precision (approximately 16 significant digits) implied volatility calculations across all input ranges, including extreme moneyness ratios and near-zero volatilities. The implementation is based on the latest (2024) version of Peter Jäckel's algorithm.

## Features

- **Implied volatility calculation** with full machine-precision accuracy
- **Black pricing** (normalized and standard forms)
- **Greeks**: Vega and Volga
- **Numerically stable** across extreme parameter ranges
- **High performance**: ~18ns for ATM, ~250-400ns for OTM/ITM cases

## Installation

```bash
go get jaeckel
```

## Usage

```go
package main

import (
    "fmt"
    "jaeckel"
)

func main() {
    // Calculate option price from volatility
    F := 100.0     // Forward price
    K := 105.0     // Strike price
    sigma := 0.20  // 20% volatility
    T := 1.0       // Time to expiry in years
    q := 1.0       // 1 for call, -1 for put

    price := jaeckel.Black(F, K, sigma, T, q)
    fmt.Printf("Option price: %.6f\n", price)

    // Calculate implied volatility from option price
    iv := jaeckel.ImpliedBlackVolatility(price, F, K, T, q)
    fmt.Printf("Implied volatility: %.6f\n", iv)

    // Calculate vega
    vega := jaeckel.Vega(F, K, sigma, T)
    fmt.Printf("Vega: %.6f\n", vega)
}
```

## API Reference

### Option Pricing

- `Black(F, K, sigma, T, q float64) float64` - Undiscounted Black option price
- `NormalizedBlack(x, s, q float64) float64` - Normalized Black price where x=ln(F/K) and s=σ√T

### Implied Volatility

- `ImpliedBlackVolatility(price, F, K, T, q float64) float64` - Compute implied volatility from option price
- `NormalizedImpliedBlackVolatility(beta, x, q float64) float64` - Normalized implied volatility

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

The algorithm uses a four-branch structure based on the normalized price level:

1. **ATM branch**: Direct formula using inverse error function
2. **Lower middle branch**: Rational cubic interpolation with Householder refinement
3. **Upper middle branch**: Rational cubic interpolation with Householder refinement
4. **Lowest/Highest branches**: Logarithmic objective function with higher-order Householder iterations

Key optimizations include:

- Remez-optimized minimax rational approximations for boundary functions
- Region-specific expansions (asymptotic for deep tails, small-t for near-zero volatility)
- Householder iterations of order 3-4 for rapid convergence
- Careful numerical treatment to avoid subtractive cancellation

## Benchmarks

Typical performance on Apple M1:

| Operation                     | Time        |
| ----------------------------- | ----------- |
| Black (ATM)                   | ~3.5 ns     |
| Black (OTM)                   | ~45-55 ns   |
| ImpliedVol (ATM)              | ~18 ns      |
| ImpliedVol (Middle branches)  | ~250-270 ns |
| ImpliedVol (Extreme branches) | ~370-390 ns |
| Vega                          | ~15-23 ns   |

## Credits and Attribution

This implementation is based on the work of **Peter Jäckel**, whose algorithm achieves full machine precision for implied volatility calculations.

### References

- Source code and paper available at: [www.jaeckel.org/LetsBeRational.7z](http://www.jaeckel.org/LetsBeRational.7z)

## License

MIT License - See [LICENSE](LICENSE) file for details.

Note: The original C++ implementation by Peter Jäckel is available under its own license terms. This Go implementation is an independent port that follows the published algorithm.

## Contributing

Contributions are welcome! Please ensure that any changes maintain the numerical accuracy guarantees of the original algorithm.
