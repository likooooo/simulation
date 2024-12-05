#include <source.hpp>
#include <cuda_test_common.hpp>
#include <cuda_vector.hpp>
using namespace cuda;
using TDisplay = pycallback_update_frame;
void test_traditional_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<OLreal> source(xsize * ysize);
    std::vector<OLreal> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<OLreal>(source.size());
        get_traditional_source(source.data(), xsize, ysize, {sigma, 0, 0});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_annular_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<OLreal> source(xsize * ysize);
    std::vector<OLreal> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<OLreal>(source.size());
        get_annular_source(source.data(), xsize, ysize, {1.0, sigma, 0, 0});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_dipole_fan_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<OLreal> source(xsize * ysize);
    std::vector<OLreal> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<OLreal>(source.size());
        get_dipole_fan_source(source.data(), xsize, ysize, {{0.5f + 0.5f * sigma, 0.5f, 0, 0}, M_PIf * 2 * sigma, M_PI_4f});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_quadratic_fan_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<OLreal> source(xsize * ysize);
    std::vector<OLreal> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<OLreal>(source.size());
        get_quadratic_fan_source(source.data(), xsize, ysize, {{0.5f + 0.5f * sigma, 0.5f, 0, 0}, M_PIf * sigma, M_PI_4f});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_dipole_leaf_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<OLreal> source(xsize * ysize);
    std::vector<OLreal> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<OLreal>(source.size());
        get_dipole_leaf_source(source.data(), xsize, ysize, {0.2f + 0.8f * sigma, 0.3f + 0.7f * sigma, 0, 0, 0});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_quadratic_leaf_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<OLreal> source(xsize * ysize);
    std::vector<OLreal> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<OLreal>(source.size());
        get_quadratic_leaf_source(source.data(), xsize, ysize, {{sigma, 0.7f, 0, 0, 0}, 0.2f});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<OLreal>(source.size());
        get_quadratic_leaf_source(source.data(), xsize, ysize, {{0.5f, sigma, 0, 0, 0}, 0.2f});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<OLreal>(source.size());
        get_quadratic_leaf_source(source.data(), xsize, ysize, {{0.5f, 0.7f, 0, 0, 0}, sigma});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
int main()
{
    py_loader::init();
    py_plot::get_default_visualizer_dir() = "/usr/local/bin";
    auto callback = py_plot::create_callback_simulation_fram_done(py::object(overload_click));
    int xsize,ysize,frame_count; 
    xsize = ysize= 470;
    frame_count = 20;
    test_traditional_source(xsize, ysize, frame_count, callback);
    test_annular_source(xsize, ysize, frame_count, callback);
    test_dipole_fan_source(xsize, ysize, frame_count, callback);
    test_quadratic_fan_source(xsize, ysize, frame_count, callback);
    test_dipole_leaf_source(xsize, ysize, frame_count, callback);
    test_quadratic_leaf_source(xsize, ysize, frame_count, callback);
}   