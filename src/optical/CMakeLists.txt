
# add_executable(diffusion_simulate fluid_simulate.cu)
# target_include_directories(diffusion_simulate PUBLIC cpp_cuda::cpp_cuda)
# target_link_libraries(diffusion_simulate PUBLIC cpp_cuda::cpp_cuda)

add_executable(test_source test_source.cpp)
target_link_libraries(test_source PUBLIC cpp_cuda::cpp_cuda)