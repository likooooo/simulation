
#include <optical/simulation_grid_info.hpp>
#include <py_helper.hpp>

using rT = float;
using grid_info_t = grid_info<rT>;
using point_dbu_t = typename grid_info_t::point_dbu_t;
using point_physical_t = typename grid_info_t::point_physical_t;

grid_info_t create_grid_info(int create_mode, vec<size_t, grid_info_t::dim> shape, rT lambda, rT sigma, rT NA, vec2<point_physical_t> roi, rT dbu)
{
    using func = grid_info_t (*)(vec<size_t, grid_info_t::dim> shape, rT lambda, rT sigma, rT NA, vec2<point_physical_t> roi, rT dbu);
    constexpr vec3<func> mode_call = {grid_info_t::create_grid_info, grid_info_t::create_grid_info_bloch_mode, grid_info_t::create_grid_info_opc_mode};
    constexpr vec3<const char*> mode_name = {"default mode", "bloch mode", "OPC mode"};
    std::cout << "# create grid info in " << mode_name.at(create_mode) << ".\n";
    grid_info_t info = mode_call.at(create_mode)(shape, lambda, sigma, NA, roi, dbu);
    std::cout << info;
    return info;
}
//== ./test_optical 0 256 256 193 0 1.2  0 -200 400 200
int main(int argc, char** argv)
{
    size_t i = 1;
    int create_mode = std::stoi(argv[i++]);
    vec<size_t, grid_info_t::dim> shape{
        std::stoul(argv[i]),
        std::stoul(argv[i + 1])
    };
    i += 2;
    rT lambda = std::stof(argv[i++]);
    rT sigma = std::stof(argv[i++]);
    rT NA = std::stof(argv[i++]);
    vec2<point_physical_t> roi{
        std::stof(argv[i]),
        std::stof(argv[i+1]),
        std::stof(argv[i+2]),
        std::stof(argv[i+3])
    };
    i += 4;
    assert(argc == i);
    create_grid_info(create_mode, shape, lambda, sigma, NA, roi, 1e-6);
}

BOOST_PYTHON_MODULE(grid_info) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("create_grid_info", &create_grid_info);
}


