cmake_minimum_required(VERSION 3.14)
project(simulation 
    DESCRIPTION "cuda cpp wrapper"
    HOMEPAGE_URL "https://github.com/likooooo/cuda"
    LANGUAGES CXX CUDA
)
set(cpp_cuda_DIR /usr/local/cmake/cpp_cuda)
find_package(cpp_cuda REQUIRED)

include_directories(include)
add_subdirectory(src/optical)