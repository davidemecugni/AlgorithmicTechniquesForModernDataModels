# AlgorithmicTechniquesForModernDataModels
Implementation of the algos seen in the DTU course 02289 Algorithmic techniques for modern data models in Fall 2025

## Project Structure

This project contains implementations for various algorithmic techniques organized in the following folders:

- **BloomFilters** - Bloom filter implementations
- **DistanceOracles** - Distance oracle algorithms
- **ApproxNNSearch** - Approximate nearest neighbor search
- **DistributedComputing** - Distributed computing algorithms
  - Intro - Introduction to distributed computing
  - SSSPandAPSP - Single-Source Shortest Path and All-Pairs Shortest Path
  - RandomizedColoring - Randomized coloring algorithms
- **Streaming** - Streaming algorithms
  - majority - Majority element finding
  - Mistra-Gries - Misra-Gries algorithm
  - Approx - Approximation algorithms
  - ApproxCounting - Approximate counting
  - FrequencyEstimation - Frequency estimation
  - CountMin - Count-Min sketch

## Building the Project

This project uses CMake as its build system. To build all executables:

```bash
# Create build directory
mkdir build
cd build

# Configure the project
cmake ..

# Build all executables
cmake --build .
```

## Running Executables

After building, executables will be generated in their respective subfolders within the build directory:

```bash
# Example: Run BloomFilters
./BloomFilters/BloomFilters

# Example: Run a streaming algorithm
./Streaming/majority/majority

# Example: Run a distributed computing algorithm
./DistributedComputing/Intro/Intro
```

## Requirements

- CMake 3.10 or higher
- C++17 compatible compiler (GCC, Clang, or MSVC)
