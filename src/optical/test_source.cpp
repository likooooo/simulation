#include <optical/source/source.hpp>
#include <cuda_test_common.hpp>
#include <cuda_vector.hpp>
using namespace cuda;
using TDisplay = pycallback_update_frame;
using real = float;
using ps = parametric_source<real>;
void test_traditional_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<real> source(xsize * ysize);
    std::vector<real> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<real>(source.size());
        ps::get_traditional_source(source.data(), xsize, ysize, {sigma, 0, 0});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_annular_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<real> source(xsize * ysize);
    std::vector<real> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<real>(source.size());
        ps::get_annular_source(source.data(), xsize, ysize, {1.0, sigma, 0, 0});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_dipole_fan_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<real> source(xsize * ysize);
    std::vector<real> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<real>(source.size());
        ps::get_dipole_fan_source(source.data(), xsize, ysize, {{0.5f + 0.5f * sigma, 0.5f, 0, 0}, M_PIf * 2 * sigma, M_PI_4f});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_quadratic_fan_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<real> source(xsize * ysize);
    std::vector<real> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<real>(source.size());
        ps::get_quadratic_fan_source(source.data(), xsize, ysize, {{0.5f + 0.5f * sigma, 0.5f, 0, 0}, M_PIf * sigma, M_PI_4f});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_dipole_leaf_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<real> source(xsize * ysize);
    std::vector<real> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<real>(source.size());
        ps::get_dipole_leaf_source(source.data(), xsize, ysize, {0.2f + 0.8f * sigma, 0.3f + 0.7f * sigma, 0, 0, 0});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_quadratic_leaf_source(int xsize, int ysize, int test_count, TDisplay& display)
{
    pageable_vector<real> source(xsize * ysize);
    std::vector<real> test(test_count);
    std::iota(test.begin(), test.end(), 0);
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<real>(source.size());
        ps::get_quadratic_leaf_source(source.data(), xsize, ysize, {{sigma, 0.7f, 0, 0, 0}, 0.2f});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<real>(source.size());
        ps::get_quadratic_leaf_source(source.data(), xsize, ysize, {{0.5f, sigma, 0, 0, 0}, 0.2f});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
    for(auto sigma : test){
        sigma /= test_count;
        source = pageable_vector<real>(source.size());
        ps::get_quadratic_leaf_source(source.data(), xsize, ysize, {{0.5f, 0.7f, 0, 0, 0}, sigma});
        display(create_ndarray_from_vector(source, {xsize, ysize}));
    }
}
void test_source_point_shift()
{
    using src = source<float>;
    using row1 = std::tuple<float, float, matrix3x3<float>>;
    std::vector<row1> rotate_matrixs;
    for(float crao : {0.0, 0.25, 0.5})
    for(float azimuth : {0.0, 0.25, 0.5, 0.75, 1.0})
        rotate_matrixs.push_back(row1(crao, azimuth, src::rotate_matrix(crao*M_PI, azimuth * M_PI)));
    print_table(rotate_matrixs, {"crao(PI)", "azimuth(PI)", "rotate-matrix"}, 1024);
    vec2<size_t> shape{2, 2};
    std::vector<float> srcimage(shape[0] * shape[1], 1);
    for(float crao : {0.0, 0.25})
    for(float azimuth : {0.0, 0.25})
    {
        auto sp = src::get_source_points(srcimage.data(), shape, vec2<float>{1, 1} / shape, 1e-6, crao*M_PI, azimuth * M_PI);
        std::cout << "crao(PI), azimuth(PI) = " << std::make_tuple(crao, azimuth) << std::endl;
        print_table(reinterpret_cast<std::vector<src::source_point::print_type>&>(sp), {"kr", "sigma", "intensity"}, 1024);
    }
    exit(0);       
}
int main()
{
    // test_source_point_shift();
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