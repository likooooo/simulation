cmake_minimum_required(VERSION 3.14)
project(simulation 
    DESCRIPTION "simulation"
    HOMEPAGE_URL "https://github.com/likooooo/simulation"
    LANGUAGES CXX CUDA
)
find_package(cpp_cuda REQUIRED)
find_package(py_visualizer REQUIRED)
include_directories(include)
add_subdirectory(src/optical)