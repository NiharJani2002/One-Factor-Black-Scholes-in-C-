# One-Factor-Black-Scholes-in-C-

# Black-Scholes Option Pricing Model - Complete Implementation in C++

## Table of Contents
1. [Overview](#overview)
2. [Mathematical Foundation](#mathematical-foundation)
3. [Black-Scholes Formula Derivation](#black-scholes-formula-derivation)
4. [Greeks and Risk Management](#greeks-and-risk-management)
5. [Code Architecture](#code-architecture)
6. [Implementation Details](#implementation-details)
7. [Mathematical Proof Verification](#mathematical-proof-verification)
8. [Optimization Techniques](#optimization-techniques)
9. [Usage Examples](#usage-examples)
10. [Limitations and Assumptions](#limitations-and-assumptions)
11. [Compilation and Execution](#compilation-and-execution)
12. [Testing and Validation](#testing-and-validation)

## Overview

This C++ implementation provides a complete Black-Scholes option pricing model for European options. The Black-Scholes model, developed by Fischer Black, Myron Scholes, and Robert Merton in the early 1970s, revolutionized financial derivatives pricing and earned Scholes and Merton the 1997 Nobel Prize in Economics.

### Key Features
- **Complete Black-Scholes Formula**: Both call and put option pricing
- **All Greeks Calculations**: Delta, Gamma, Theta, Vega, and Rho
- **Numerical Stability**: Optimized mathematical functions
- **Interactive Interface**: User-friendly input/output system
- **Scenario Analysis**: Multiple pricing scenarios
- **Input Validation**: Robust error handling

## Mathematical Foundation

### The Black-Scholes Differential Equation

The Black-Scholes model is based on the following stochastic differential equation for stock price movement:

```
dS = μS dt + σS dW
```

Where:
- `S` = Stock price
- `μ` = Expected return (drift)
- `σ` = Volatility
- `dW` = Wiener process (Brownian motion)
- `dt` = Infinitesimal time increment

### Key Assumptions

1. **Constant Parameters**: Risk-free rate `r` and volatility `σ` are constant
2. **Geometric Brownian Motion**: Stock prices follow a lognormal distribution
3. **No Dividends**: The stock pays no dividends during the option's life
4. **European Exercise**: Options can only be exercised at expiration
5. **No Transaction Costs**: Perfect market with no bid-ask spreads
6. **Continuous Trading**: Markets operate continuously
7. **Risk-Free Borrowing**: Unlimited borrowing/lending at risk-free rate

## Black-Scholes Formula Derivation

### Step 1: Risk-Neutral Valuation

Under risk-neutral measure, the stock price follows:
```
dS = rS dt + σS dW*
```

Where `W*` is a Wiener process under the risk-neutral measure.

### Step 2: Solution to the SDE

The stock price at time `T` is:
```
S(T) = S(0) * exp((r - σ²/2)T + σ√T * Z)
```

Where `Z ~ N(0,1)` is a standard normal random variable.

### Step 3: Option Payoff

For a European call option:
```
Payoff = max(S(T) - K, 0)
```

For a European put option:
```
Payoff = max(K - S(T), 0)
```

### Step 4: Risk-Neutral Expectation

The option price is the discounted expected payoff:
```
V = e^(-rT) * E[Payoff]
```

### Step 5: Final Black-Scholes Formulas

**Call Option Price:**
```
C = S₀ * N(d₁) - K * e^(-rT) * N(d₂)
```

**Put Option Price:**
```
P = K * e^(-rT) * N(-d₂) - S₀ * N(-d₁)
```

**Where:**
```
d₁ = [ln(S₀/K) + (r + σ²/2)T] / (σ√T)
d₂ = d₁ - σ√T
```

**And:**
- `N(x)` = Cumulative standard normal distribution function
- `S₀` = Current stock price
- `K` = Strike price
- `T` = Time to expiration
- `r` = Risk-free rate
- `σ` = Volatility

## Greeks and Risk Management

### Delta (Δ) - Price Sensitivity
Measures the rate of change of option price with respect to stock price.

**Call Delta:**
```cpp
Δ_call = N(d₁)
```

**Put Delta:**
```cpp
Δ_put = N(d₁) - 1 = -N(-d₁)
```

**Mathematical Interpretation:**
- Range: Call [0, 1], Put [-1, 0]
- Delta-neutral portfolios have Δ = 0
- Approximates hedge ratio for dynamic hedging

### Gamma (Γ) - Delta Sensitivity
Measures the rate of change of delta with respect to stock price.

```cpp
Γ = φ(d₁) / (S * σ * √T)
```

Where `φ(x)` is the standard normal probability density function.

**Properties:**
- Always positive for both calls and puts
- Highest for at-the-money options
- Approaches zero for deep ITM/OTM options

### Theta (Θ) - Time Decay
Measures the rate of change of option price with respect to time.

**Call Theta:**
```cpp
Θ_call = -[S * φ(d₁) * σ / (2√T) + r * K * e^(-rT) * N(d₂)]
```

**Put Theta:**
```cpp
Θ_put = -[S * φ(d₁) * σ / (2√T) - r * K * e^(-rT) * N(-d₂)]
```

### Vega (ν) - Volatility Sensitivity
Measures the rate of change of option price with respect to volatility.

```cpp
ν = S * φ(d₁) * √T
```

**Properties:**
- Always positive for both calls and puts
- Maximum for at-the-money options
- Higher for longer-term options

### Rho (ρ) - Interest Rate Sensitivity
Measures the rate of change of option price with respect to risk-free rate.

**Call Rho:**
```cpp
ρ_call = K * T * e^(-rT) * N(d₂)
```

**Put Rho:**
```cpp
ρ_put = -K * T * e^(-rT) * N(-d₂)
```

## Code Architecture

### Class Design

```cpp
class BlackScholes {
private:
    double S, K, T, r, sigma;  // Model parameters
    
    // Core mathematical functions
    double normalCDF(double x) const;
    double normalPDF(double x) const;
    double calculateD1() const;
    double calculateD2() const;
    
public:
    // Constructor and pricing methods
    BlackScholes(double, double, double, double, double);
    double callPrice() const;
    double putPrice() const;
    
    // Greeks calculations
    double delta(bool isCall = true) const;
    double gamma() const;
    double theta(bool isCall = true) const;
    double vega() const;
    double rho(bool isCall = true) const;
    
    // Utility functions
    void displayResults() const;
};
```

### Design Principles

1. **Encapsulation**: All parameters are private with controlled access
2. **Const Correctness**: Mathematical functions don't modify state
3. **Single Responsibility**: Each method has a specific purpose
4. **Reusability**: d₁ and d₂ calculations are centralized
5. **Error Handling**: Input validation and edge case management

## Implementation Details

### 1. Normal Distribution Functions

#### Cumulative Distribution Function (CDF)
```cpp
double normalCDF(double x) const {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}
```

**Mathematical Basis:**
The relationship between the error function and normal CDF:
```
N(x) = ½[1 + erf(x/√2)]
```

Where:
```
erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
```

#### Probability Density Function (PDF)
```cpp
double normalPDF(double x) const {
    return (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
}
```

**Mathematical Formula:**
```
φ(x) = (1/√(2π)) * e^(-x²/2)
```

### 2. Parameter Calculations

#### d₁ Calculation
```cpp
double calculateD1() const {
    return (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
}
```

**Derivation:**
Starting from the risk-neutral stock price evolution:
```
S(T) = S₀ * exp((r - σ²/2)T + σ√T * Z)
```

For the option to be in-the-money:
```
S(T) > K
S₀ * exp((r - σ²/2)T + σ√T * Z) > K
(r - σ²/2)T + σ√T * Z > ln(K/S₀)
Z > [ln(K/S₀) - (r - σ²/2)T] / (σ√T)
Z > -d₁
```

Therefore:
```
d₁ = [ln(S₀/K) + (r + σ²/2)T] / (σ√T)
```

#### d₂ Calculation
```cpp
double calculateD2() const {
    return calculateD1() - sigma * sqrt(T);
}
```

**Mathematical Relationship:**
```
d₂ = d₁ - σ√T
```

This relationship comes from the drift adjustment in the risk-neutral measure.

### 3. Edge Case Handling

#### Zero Time to Expiration
```cpp
if (T <= 0) return std::max(S - K, 0.0);  // Call
if (T <= 0) return std::max(K - S, 0.0);  // Put
```

When `T = 0`, options become worth their intrinsic value.

#### Greeks at Expiration
```cpp
if (T <= 0) {
    if (isCall) return (S > K) ? 1.0 : 0.0;  // Delta
    else return (S < K) ? -1.0 : 0.0;
}
```

At expiration, delta becomes a step function.

## Mathematical Proof Verification

### 1. Put-Call Parity Verification

The Black-Scholes formulas satisfy put-call parity:
```
C - P = S₀ - K * e^(-rT)
```

**Proof:**
```
C - P = S₀ * N(d₁) - K * e^(-rT) * N(d₂) - [K * e^(-rT) * N(-d₂) - S₀ * N(-d₁)]
      = S₀ * N(d₁) - K * e^(-rT) * N(d₂) - K * e^(-rT) * N(-d₂) + S₀ * N(-d₁)
      = S₀ * [N(d₁) + N(-d₁)] - K * e^(-rT) * [N(d₂) + N(-d₂)]
      = S₀ * 1 - K * e^(-rT) * 1
      = S₀ - K * e^(-rT)
```

This verification is implemented in the test suite.

### 2. Delta Relationship Verification

For calls and puts on the same underlying:
```
Δ_call - Δ_put = 1
```

**Proof:**
```
Δ_call - Δ_put = N(d₁) - [N(d₁) - 1] = N(d₁) - N(d₁) + 1 = 1
```

### 3. Gamma Equivalence

Both calls and puts have the same gamma:
```
Γ_call = Γ_put = φ(d₁) / (S * σ * √T)
```

This is verified through numerical testing in our implementation.

## Optimization Techniques

### 1. Mathematical Optimizations

#### Efficient Error Function Implementation
```cpp
// Using standard library erf() function for accuracy and performance
return 0.5 * (1.0 + erf(x / sqrt(2.0)));
```

The standard library `erf()` function uses optimized polynomial approximations.

#### Cached Calculations
```cpp
double d1 = calculateD1();
double d2 = calculateD2();  // Uses d1 - σ√T relationship
```

d₂ is calculated from d₁ to avoid redundant logarithm calculations.

### 2. Numerical Stability

#### Preventing Overflow/Underflow
- Using `exp(-r * T)` instead of `pow(e, -r * T)`
- Proper handling of extreme values in normal distribution functions
- Input validation to prevent mathematical errors

#### Precision Considerations
```cpp
std::cout << std::fixed << std::setprecision(4);
```

Using appropriate precision for financial calculations (typically 4 decimal places).

### 3. Memory Efficiency

#### Stack Allocation
All calculations use stack-allocated variables, avoiding heap fragmentation.

#### Minimal Memory Footprint
```cpp
class BlackScholes {
private:
    double S, K, T, r, sigma;  // Only 5 double values (40 bytes on 64-bit)
    // ...
};
```

### 4. Computational Complexity

- **Time Complexity**: O(1) for all calculations
- **Space Complexity**: O(1) constant space usage
- **Function Calls**: Minimal function call overhead with inline optimizations

## Usage Examples

### Basic Usage
```cpp
// Create option with parameters
BlackScholes option(100.0, 105.0, 0.25, 0.05, 0.20);

// Get prices
double callPrice = option.callPrice();
double putPrice = option.putPrice();

// Get Greeks
double delta = option.delta(true);  // Call delta
double gamma = option.gamma();
double theta = option.theta(true);  // Call theta
```

### Scenario Analysis
```cpp
// At-the-money option
BlackScholes atmOption(100.0, 100.0, 0.25, 0.05, 0.20);

// In-the-money call
BlackScholes itmCall(100.0, 90.0, 0.25, 0.05, 0.20);

// Out-of-the-money call
BlackScholes otmCall(100.0, 110.0, 0.25, 0.05, 0.20);
```

### Interactive Mode
The program provides an interactive interface for real-time calculations:

```
=== Black-Scholes Option Pricing Calculator ===
Enter the following parameters:
Current Stock Price: $100
Strike Price: $105
Time to Expiration (years): 0.25
Risk-free Rate (as decimal, e.g., 0.05 for 5%): 0.05
Volatility (as decimal, e.g., 0.20 for 20%): 0.20
```

## Limitations and Assumptions

### Model Limitations

1. **Constant Volatility**: Real volatility changes over time
2. **Constant Interest Rates**: Interest rates fluctuate
3. **No Dividends**: Many stocks pay dividends
4. **European Exercise**: American options can be exercised early
5. **Perfect Liquidity**: Real markets have transaction costs
6. **Continuous Trading**: Markets have gaps and closures

### Numerical Limitations

1. **Floating Point Precision**: Limited to double precision (15-17 significant digits)
2. **Extreme Parameters**: Very large or small values may cause numerical issues
3. **Time to Expiration**: T = 0 requires special handling

### Practical Considerations

1. **Market Data Quality**: Garbage in, garbage out
2. **Parameter Estimation**: Volatility and risk-free rate estimation challenges
3. **Model Risk**: Black-Scholes may not capture all market dynamics

## Compilation and Execution

### Compilation Options

#### Basic Compilation
```bash
g++ -o blackscholes blackscholes.cpp -lm
```

#### Optimized Compilation
```bash
g++ -O2 -o blackscholes blackscholes.cpp -lm
```

#### Debug Compilation
```bash
g++ -g -Wall -Wextra -o blackscholes_debug blackscholes.cpp -lm
```

#### Production Compilation
```bash
g++ -O3 -DNDEBUG -march=native -o blackscholes_prod blackscholes.cpp -lm
```

### Compiler Requirements

- **C++ Standard**: C++11 or later
- **Math Library**: `-lm` flag required for mathematical functions
- **Platform**: Cross-platform (Linux, macOS, Windows with MinGW)

### Execution
```bash
./blackscholes
```

## Testing and Validation

### Unit Tests

#### Put-Call Parity Test
```cpp
void testPutCallParity() {
    BlackScholes bs(100, 105, 0.25, 0.05, 0.20);
    double call = bs.callPrice();
    double put = bs.putPrice();
    double parity = call - put;
    double expected = 100 - 105 * exp(-0.05 * 0.25);
    assert(abs(parity - expected) < 1e-10);
}
```

#### Delta Sum Test
```cpp
void testDeltaSum() {
    BlackScholes bs(100, 105, 0.25, 0.05, 0.20);
    double callDelta = bs.delta(true);
    double putDelta = bs.delta(false);
    assert(abs((callDelta - putDelta) - 1.0) < 1e-10);
}
```

#### Greeks Consistency Test
```cpp
void testGamma() {
    BlackScholes bs(100, 100, 0.25, 0.05, 0.20);
    double gamma = bs.gamma();
    assert(gamma > 0);  // Gamma should always be positive
}
```

### Benchmark Tests

#### Known Values Verification
Test against published Black-Scholes values from financial literature.

#### Numerical Stability Tests
Test with extreme parameter values to ensure numerical stability.

#### Performance Benchmarks
Measure calculation time for large numbers of option calculations.

### Integration Tests

#### Market Data Validation
Compare results with real market option prices (considering bid-ask spreads).

#### Cross-Platform Testing
Verify consistent results across different operating systems and compilers.

## Advanced Extensions

### Potential Enhancements

1. **American Options**: Binomial tree or Monte Carlo methods
2. **Dividend Adjustments**: Modified Black-Scholes for dividend-paying stocks
3. **Implied Volatility**: Newton-Raphson solver for market-implied volatility
4. **Interest Rate Models**: Stochastic interest rate incorporation
5. **Jump Diffusion**: Merton jump-diffusion model
6. **Exotic Options**: Asian, barrier, and binary options

### Performance Optimizations

1. **SIMD Instructions**: Vectorized calculations for multiple options
2. **GPU Computing**: CUDA implementation for massive parallel calculations
3. **Lookup Tables**: Pre-computed normal distribution values
4. **Approximation Methods**: Faster approximations for real-time systems

## Mathematical References

1. Black, F., & Scholes, M. (1973). "The Pricing of Options and Corporate Liabilities"
2. Merton, R. C. (1973). "Theory of Rational Option Pricing"
3. Hull, J. C. "Options, Futures, and Other Derivatives" (10th Edition)
4. Wilmott, P. "Paul Wilmott Introduces Quantitative Finance"

## Conclusion

This implementation provides a complete, mathematically rigorous, and computationally efficient Black-Scholes option pricing system. The code demonstrates deep understanding of both the underlying mathematical theory and practical software engineering principles. While the Black-Scholes model has limitations in real-world applications, it remains the foundation for modern derivatives pricing and risk management.

The implementation serves as both an educational tool for understanding option pricing theory and a practical foundation for more advanced financial modeling systems.

---

*This implementation is for educational and research purposes. Financial decisions should not be made solely based on this model without considering its limitations and consulting with qualified financial professionals.*
