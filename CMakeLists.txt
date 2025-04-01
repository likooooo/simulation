cmake_minimum_required(VERSION 3.14)
project(simulation 
    DESCRIPTION "simulation"
    HOMEPAGE_URL "https://github.com/likooooo/simulation"
    LANGUAGES CXX
)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# find_package(cpp_cuda REQUIRED)
find_package(py_visualizer REQUIRED)
find_package(kernels REQUIRED)
include_directories(include)
add_subdirectory(src/optical)
include(cmake/auto_ctest.cmake)