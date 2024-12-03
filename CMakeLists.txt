cmake_minimum_required(VERSION 3.14)
project(simulation 
    DESCRIPTION "cuda cpp wrapper"
    HOMEPAGE_URL "https://github.com/likooooo/cuda"
    LANGUAGES CXX CUDA
)
set(cpp_cuda_DIR /usr/local/cmake/cpp_cuda)
find_package(cpp_cuda REQUIRED)

# enable_language(CUDA CXX)
# # 查找 CUDA 包
# find_package(CUDAToolkit REQUIRED)

# # 设置 CUDA 编译器
# set(CMAKE_CUDA_STANDARD 17)
# set(CMAKE_CUDA_STANDARD_REQUIRED ON)
# # 查找 CUDA 包
# find_package(CUDAToolkit REQUIRED)


# find_package(Threads REQUIRED)

# # 查找 cpp_cuda 包
# find_package(cpp_cuda REQUIRED)
add_executable(diffusion_simulate fluid_simulate.cu)
target_include_directories(diffusion_simulate PUBLIC cpp_cuda::cpp_cuda)
target_link_libraries(diffusion_simulate PUBLIC cpp_cuda::cpp_cuda)
