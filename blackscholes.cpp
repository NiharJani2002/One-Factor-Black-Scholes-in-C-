#include <iostream>
#include <cmath>
#include <iomanip>

class BlackScholes {
private:
    double S;  // Current stock price
    double K;  // Strike price
    double T;  // Time to expiration (in years)
    double r;  // Risk-free rate
    double sigma; // Volatility
    
    // Cumulative standard normal distribution function
    double normalCDF(double x) const {
        return 0.5 * (1.0 + erf(x / sqrt(2.0)));
    }
    
    // Standard normal probability density function
    double normalPDF(double x) const {
        return (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
    }
    
    // Calculate d1 parameter
    double calculateD1() const {
        return (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    }
    
    // Calculate d2 parameter
    double calculateD2() const {
        return calculateD1() - sigma * sqrt(T);
    }

public:
    // Constructor
    BlackScholes(double stockPrice, double strikePrice, double timeToExpiry, 
                 double riskFreeRate, double volatility) 
        : S(stockPrice), K(strikePrice), T(timeToExpiry), r(riskFreeRate), sigma(volatility) {}
    
    // Calculate call option price
    double callPrice() const {
        if (T <= 0) return std::max(S - K, 0.0);
        
        double d1 = calculateD1();
        double d2 = calculateD2();
        
        return S * normalCDF(d1) - K * exp(-r * T) * normalCDF(d2);
    }
    
    // Calculate put option price
    double putPrice() const {
        if (T <= 0) return std::max(K - S, 0.0);
        
        double d1 = calculateD1();
        double d2 = calculateD2();
        
        return K * exp(-r * T) * normalCDF(-d2) - S * normalCDF(-d1);
    }
    
    // Calculate Greeks
    double delta(bool isCall = true) const {
        if (T <= 0) {
            if (isCall) return (S > K) ? 1.0 : 0.0;
            else return (S < K) ? -1.0 : 0.0;
        }
        
        double d1 = calculateD1();
        if (isCall) {
            return normalCDF(d1);
        } else {
            return normalCDF(d1) - 1.0;
        }
    }
    
    double gamma() const {
        if (T <= 0) return 0.0;
        
        double d1 = calculateD1();
        return normalPDF(d1) / (S * sigma * sqrt(T));
    }
    
    double theta(bool isCall = true) const {
        if (T <= 0) return 0.0;
        
        double d1 = calculateD1();
        double d2 = calculateD2();
        
        double term1 = -(S * normalPDF(d1) * sigma) / (2.0 * sqrt(T));
        
        if (isCall) {
            double term2 = r * K * exp(-r * T) * normalCDF(d2);
            return (term1 - term2) / 365.0; // Convert to per day
        } else {
            double term2 = r * K * exp(-r * T) * normalCDF(-d2);
            return (term1 + term2) / 365.0; // Convert to per day
        }
    }
    
    double vega() const {
        if (T <= 0) return 0.0;
        
        double d1 = calculateD1();
        return S * normalPDF(d1) * sqrt(T) / 100.0; // Convert to per 1% change
    }
    
    double rho(bool isCall = true) const {
        if (T <= 0) return 0.0;
        
        double d2 = calculateD2();
        
        if (isCall) {
            return K * T * exp(-r * T) * normalCDF(d2) / 100.0; // Convert to per 1% change
        } else {
            return -K * T * exp(-r * T) * normalCDF(-d2) / 100.0; // Convert to per 1% change
        }
    }
    
    // Display all option information
    void displayResults() const {
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "\n=== Black-Scholes Option Pricing Results ===" << std::endl;
        std::cout << "Parameters:" << std::endl;
        std::cout << "  Stock Price (S): $" << S << std::endl;
        std::cout << "  Strike Price (K): $" << K << std::endl;
        std::cout << "  Time to Expiry (T): " << T << " years" << std::endl;
        std::cout << "  Risk-free Rate (r): " << r * 100 << "%" << std::endl;
        std::cout << "  Volatility (σ): " << sigma * 100 << "%" << std::endl;
        
        std::cout << "\nOption Prices:" << std::endl;
        std::cout << "  Call Price: $" << callPrice() << std::endl;
        std::cout << "  Put Price: $" << putPrice() << std::endl;
        
        std::cout << "\nGreeks:" << std::endl;
        std::cout << "  Call Delta: " << delta(true) << std::endl;
        std::cout << "  Put Delta: " << delta(false) << std::endl;
        std::cout << "  Gamma: " << gamma() << std::endl;
        std::cout << "  Call Theta: " << theta(true) << " (per day)" << std::endl;
        std::cout << "  Put Theta: " << theta(false) << " (per day)" << std::endl;
        std::cout << "  Vega: " << vega() << " (per 1% vol change)" << std::endl;
        std::cout << "  Call Rho: " << rho(true) << " (per 1% rate change)" << std::endl;
        std::cout << "  Put Rho: " << rho(false) << " (per 1% rate change)" << std::endl;
    }
};

// Function to get user input
void getUserInput(double& S, double& K, double& T, double& r, double& sigma) {
    std::cout << "=== Black-Scholes Option Pricing Calculator ===" << std::endl;
    std::cout << "Enter the following parameters:" << std::endl;
    
    std::cout << "Current Stock Price: $";
    std::cin >> S;
    
    std::cout << "Strike Price: $";
    std::cin >> K;
    
    std::cout << "Time to Expiration (years): ";
    std::cin >> T;
    
    std::cout << "Risk-free Rate (as decimal, e.g., 0.05 for 5%): ";
    std::cin >> r;
    
    std::cout << "Volatility (as decimal, e.g., 0.20 for 20%): ";
    std::cin >> sigma;
}

int main() {
    double S, K, T, r, sigma;
    char choice;
    
    do {
        // Get user input
        getUserInput(S, K, T, r, sigma);
        
        // Validate inputs
        if (S <= 0 || K <= 0 || T < 0 || sigma <= 0) {
            std::cout << "Error: Invalid input parameters. Please ensure all values are positive (T can be zero)." << std::endl;
            continue;
        }
        
        // Create BlackScholes object and display results
        BlackScholes bs(S, K, T, r, sigma);
        bs.displayResults();
        
        // Example calculations with different scenarios
        std::cout << "\n=== Scenario Analysis ===" << std::endl;
        
        // ATM option
        BlackScholes atmOption(S, S, T, r, sigma);
        std::cout << "At-the-Money (K = S = $" << S << "):" << std::endl;
        std::cout << "  Call Price: $" << std::fixed << std::setprecision(4) << atmOption.callPrice() << std::endl;
        std::cout << "  Put Price: $" << atmOption.putPrice() << std::endl;
        
        // ITM Call / OTM Put
        BlackScholes itmCall(S, S * 0.9, T, r, sigma);
        std::cout << "In-the-Money Call (K = $" << S * 0.9 << "):" << std::endl;
        std::cout << "  Call Price: $" << itmCall.callPrice() << std::endl;
        
        // OTM Call / ITM Put
        BlackScholes otmCall(S, S * 1.1, T, r, sigma);
        std::cout << "Out-of-the-Money Call (K = $" << S * 1.1 << "):" << std::endl;
        std::cout << "  Call Price: $" << otmCall.callPrice() << std::endl;
        
        std::cout << "\nDo you want to calculate another option? (y/n): ";
        std::cin >> choice;
        
    } while (choice == 'y' || choice == 'Y');
    
    std::cout << "Thank you for using the Black-Scholes Calculator!" << std::endl;
    return 0;
}