# Monte Carlo Simulation for Autocallable Options

## Overview
This project implements a Monte Carlo simulation for pricing autocallable options using C++. It utilizes multithreading to enhance performance and is designed to simulate various scenarios in financial markets.

## Features
- Sobol sequence for quasi-random number generation.
- Geometric Brownian Motion model for asset price simulation.
- Multithreaded simulations using C++ standard library's `std::thread`.
- Calculation of payoff for autocallable options.
- Performance measurement using `std::chrono`.

## Requirements
- C++ compiler (C++17 or later)
- POSIX-compliant operating system (for pthreads, if not using `std::thread`)

## Installation
Clone the repository to your local machine:

## Compiling and Running the Code

### Using Makefile
The project includes a `Makefile` for easy compilation. Follow these steps:

1. To clean up the compiled files, you can use:
    ```
    make clean
    ```

2. Run the following command to compile the code:
    ```
    make
    ```
    This will create an executable named `main`.

3. To run the simulation, use:
    ```
    ./main
    ```
