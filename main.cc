#include <iostream>
#include <cmath>
#include <vector>
#include <pthread.h>
#include <chrono>
#include "monte_carlo.h"

int main(){

    // Parameters for the simulation
    int numPaths = 10000;                 // Number of Monte Carlo paths
    int numTimeSteps = 252;             // Number of time steps (e.g., days in a year)
    int numAssets = 5;                  // Number of underlying assets
    std::vector<double> initialPrices(numAssets, 100);  // Initial prices of assets
    double drift = 0.05;                // Drift rate
    double volatility = 0.2;            // Volatility
    double timeToMaturity = 1.0;        // Time to maturity in years
    double strikePrice = 180;           // Strike price of the option
    double autocallBarrier = 1.1;       // Autocall barrier (110% of strike price)
    double couponRate = 0.03;           // Coupon rate
    std::vector<int> autocallDates = {30, 90, 180, 270, 360}; // Autocall dates (in days)

    const int numThreads = 4;           // Number of threads
    //std::vector<std::thread> threads;
    pthread_t threads[numThreads];
    std::vector<ThreadArgs> args;
    int pathsPerThread = numPaths / numThreads;

    // Start timer 
    auto startTime = std::chrono::high_resolution_clock::now();

    // Create threads
    for (int i = 0; i < numThreads; ++i) {
        args.emplace_back(i * pathsPerThread, pathsPerThread, numTimeSteps, numAssets, initialPrices, drift, volatility, timeToMaturity, strikePrice, autocallDates, autocallBarrier, couponRate);
        pthread_create(&threads[i], nullptr, simulatePaths, &args.back());
        //threads.emplace_back(simulatePaths, &args.back());

        std::cout << "Thread " << i << " created" << std::endl;
    }

    // Join threads and accumulate results
    double totalPayoff = 0.0;
    for (int i = 0; i < numThreads; ++i) {
        pthread_join(threads[i], nullptr);
        //threads[i].join();
        totalPayoff += args[i].totalPayoff;

        std::cout << "Thread " << i << " joined" << std::endl;
    }

    double optionPrice = totalPayoff / numPaths;

    // Stop timer
    auto endTime = std::chrono::high_resolution_clock::now();

    // Print results
    std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms" << std::endl;
    std::cout << "Estimated price of the autocallable option: " << optionPrice << std::endl;


    return 0;
}