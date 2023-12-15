#ifndef MONTE_CARLO_H 
#define MONTE_CARLO_H

#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <random>
#include <functional>
#include <limits>

// Class for generating quasi-random numbers from Sobol sequence
class SobolGenerator{
private:
    static const int MAX_DIMENSIONS = 5;
    static const int MAX_BITS = 30;
    std::vector<std::vector<unsigned>> directionNumbers;
    std::vector<unsigned> index;

public:
    SobolGenerator(){
        directionNumbers = {
            {1 << (MAX_BITS -1)},
            {1 << (MAX_BITS -2), 3 << (MAX_BITS -3)},
            {1 << (MAX_BITS -3), 3 << (MAX_BITS -4)},
            {1 << (MAX_BITS -4), 3 << (MAX_BITS -5)},
            {1 << (MAX_BITS -5), 3 << (MAX_BITS -6)}
        };
    
        for (auto &vec : directionNumbers){
            vec.resize(MAX_BITS, 0);
        }

        index.resize(MAX_DIMENSIONS, 0);

        std::cout << "Construction for SobolGenerator complete" << std::endl;
    }

    ~SobolGenerator(){
        //std::cout << "Destruction for SobolGenerator complete" << std::endl;
    }

    std::vector<double> next() {
        std::vector<double> point(MAX_DIMENSIONS, 0.0);

        for (int dim = 0; dim < MAX_DIMENSIONS; ++dim){
            int i = 0;
            int m = 1;

            while (index[dim] & m) {
                index[dim] ^= directionNumbers[dim][i];
                i++;
                m <<= 1;
            }

            index[dim] ^= directionNumbers[dim][i];

            point[dim] = static_cast<double>(index[dim]) / (1U << MAX_BITS);
        }

        return point;
    }
};

// Function to apply importance sampling to a path of an asset
void ImportanceSampling(std::vector<double>& assetPrices, const std::vector<double>& shifts) {
    if (assetPrices.size() != shifts.size()) {
        throw std::invalid_argument("assetPrices and shifts must be of the same size.");
    }

    for (size_t i = 0; i < assetPrices.size(); ++i) {
        assetPrices[i] += shifts[i];
    }
}

// Function to compute the inverse of the standard normal cumulative distribution function (CDF)
double normInv(double p) {
    if (p < 0.0 || p > 1.0) {
        std::cerr << "Invalid input for normInv: " << p << std::endl;
        return NAN; // Not a Number
    }

    std::normal_distribution<double> dist(0.0, 1.0);
    std::default_random_engine generator;

    // Lower and upper bounds for the search
    double lower = -10.0; // Lower bound of the distribution
    double upper = 10.0;  // Upper bound of the distribution
    double x;
    double precision = 1e-6; // Precision of the search

    // Binary search for finding the inverse
    while (upper - lower > precision) {
        x = (lower + upper) / 2;
        double p_x = dist(generator);
        if (p_x < p) {
            lower = x;
        } else {
            upper = x;
        }
    }

    return (lower + upper) / 2;
}

// Function to simulate a path of an asset using the Geometric Brownian Motion (GBM) model
void SimulateAssetPath(std::vector<double>& assetPath,
                       SobolGenerator sobolGenerator,
                       double initialPrice,
                       double drift,
                       double volatility,
                       double timeToMaturity,
                       int numTimeSteps) {
    // Check if assetPath is the correct size, resize if necessary
    if (assetPath.size() != static_cast<size_t>(numTimeSteps)) {
        assetPath.resize(numTimeSteps);
    }

    double dt = timeToMaturity / numTimeSteps; // Time step size
    double currentPrice = initialPrice;

    // Starting the path with the initial price
    assetPath[0] = initialPrice;

    for (int i = 1; i < numTimeSteps; ++i) {
        // Generate a quasi-random number from Sobol sequence
        std::vector<double> sobolPoint = sobolGenerator.next();
        double z = normInv(sobolPoint[0]); // Assuming normInv is a function that computes the inverse of the standard normal CDF

        // Apply Geometric Brownian Motion (GBM) formula
        currentPrice *= exp((drift - 0.5 * volatility * volatility) * dt + volatility * sqrt(dt) * z);

        assetPath[i] = currentPrice;
    }

    //std::cout << "SimulateAssetPath complete" << std::endl;
}

// Function to calculate the payoff of an autocallable option
double AutocallableOption(const std::vector<std::vector<double>>& assetPaths,
                            double strikePrice,
                            const std::vector<int>& autocallDates,
                            double autocallBarrier,
                            double couponRate) {
    double payoff = 0.0;
    bool isAutocalled = false;

    for (int date : autocallDates) {
        bool aboveBarrier = true;

        for (const auto& path : assetPaths) {
            if (path[date] < autocallBarrier * strikePrice) {
                aboveBarrier = false;
                break;
            }
        }

        if (aboveBarrier) {
            payoff = couponRate * strikePrice; 
            isAutocalled = true;
            break;
        }
    }

    if (!isAutocalled) {
        double minAssetPrice = std::numeric_limits<double>::max();
        for (const auto& path : assetPaths) {
            minAssetPrice = std::min(minAssetPrice, path.back());
        }
        payoff = std::max(0.0, minAssetPrice - strikePrice); 
    }

    //std::cout << "AutocallableOption complete" << std::endl;

    return payoff;
}

// Structure for passing arguments to the thread function
struct ThreadArgs {
    int startPath;
    int numPaths;
    int numTimeSteps;
    int numAssets;
    std::vector<double> initialPrices;
    double drift;
    double volatility;
    double timeToMaturity;
    double strikePrice;
    std::vector<int> autocallDates;
    double autocallBarrier;
    double couponRate;
    double totalPayoff;  
    SobolGenerator sobolGenerator;  

    ThreadArgs(int start, int num, int timeSteps, int assets, const std::vector<double>& prices, double dr, double vol, double maturity, double strike, const std::vector<int>& dates, double barrier, double coupon)
        : startPath(start), numPaths(num), numTimeSteps(timeSteps), numAssets(assets), initialPrices(prices), drift(dr), volatility(vol), timeToMaturity(maturity), strikePrice(strike), autocallDates(dates), autocallBarrier(barrier), couponRate(coupon), totalPayoff(0.0), sobolGenerator() {
            std::cout << "ThreadArgs constructor called" << std::endl;
        }
};

// Thread function for simulating paths
void* simulatePaths(void* args) {
    auto* threadArgs = static_cast<ThreadArgs*>(args);

    double threadPayoff = 0.0;  // Store the payoff calculated by this thread

    for (int i = threadArgs->startPath; i < threadArgs->startPath + threadArgs->numPaths; ++i) {
        std::vector<std::vector<double>> assetPaths(threadArgs->numAssets, std::vector<double>(threadArgs->numTimeSteps));

        for (int asset = 0; asset < threadArgs->numAssets; ++asset) {
            SimulateAssetPath(assetPaths[asset], threadArgs->sobolGenerator, threadArgs->initialPrices[asset], threadArgs->drift, threadArgs->volatility, threadArgs->timeToMaturity, threadArgs->numTimeSteps);

            // Apply importance sampling to the path of each asset
            // Assuming shifts vector is part of threadArgs or calculated here
            std::vector<double> shifts(threadArgs->numTimeSteps, 0.0); // Example: No shift applied
            ImportanceSampling(assetPaths[asset], shifts);
        }

        // Calculate payoff for the path and accumulate in threadPayoff
        threadPayoff += AutocallableOption(assetPaths, threadArgs->strikePrice, threadArgs->autocallDates, threadArgs->autocallBarrier, threadArgs->couponRate);
    }

    // Store the result in the ThreadArgs structure
    threadArgs->totalPayoff = threadPayoff;

    std::cout << "Thread " << threadArgs->startPath << " complete" << std::endl;

    return nullptr;
}

#endif // MONTE_CARLO_H