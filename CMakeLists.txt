cmake_minimum_required(VERSION 3.10)
project(AlgorithmicTechniquesForModernDataModels)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add subdirectories
add_subdirectory(BloomFilters)
add_subdirectory(DistanceOracles)
add_subdirectory(ApproxNNSearch)
add_subdirectory(DistributedComputing)
add_subdirectory(Streaming)
