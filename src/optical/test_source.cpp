#include <source.hpp>
#include <cuda_test_common.hpp>
#include <cuda_vector.hpp>

int main()
{
    py_loader::init();
    py_plot::get_default_visualizer_dir() = "/usr/local/bin";

    using namespace cuda;
    int xsize,ysize; 
    xsize = ysize= 470;
    pageable_vector<OLreal> source(xsize * ysize);
    // get_traditional_source(source.data(), xsize, ysize, traditional_source_params{1, 0, 0});
    // get_annular_source(source.data(), xsize, ysize, {1.0, 0.5, 0, 0});
    // get_annular_source(source.data(), xsize, ysize, {1.0, 0.5, 0, 0});
    // get_dipole_fan_source(source.data(), xsize, ysize, {{1.0, 0.5, 0, 0}, M_PI_4, M_PI_4});
    // get_quadratic_fan_source(source.data(), xsize, ysize, {{1.0, 0.5, 0, 0}, M_PI_4, M_PI_4});
    // get_dipole_leaf_source(source.data(), xsize, ysize, {0.2, 0.3, 0, 0, 0});
    get_quadratic_leaf_source(source.data(), xsize, ysize, {{0.5, 0.7, 0, 0, 0}, 0.2});
    imshow(source, {xsize, ysize});
}   