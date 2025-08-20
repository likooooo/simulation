#include <optical/simulation_grid_info.hpp>
#include <py_helper.hpp>

template<class T, size_t dim = 2> void regist_grid_info()
{
    using grid_info_t = grid_info<T, dim>;
    std::string suffix =  to_string(dim) + "d_" + get_numerical_type_str_suffix<T>();
    py::class_<grid_info_t>(("grid_info_" + suffix).c_str()).def(py::init<>())       
        .def_readwrite("dbu", &grid_info_t::dbu)
        .def_readwrite("spatial",&grid_info_t::spatial)
        .def_readwrite("fourier",&grid_info_t::fourier)
        .def_readwrite("tilesize",&grid_info_t::tilesize)
        .def("__repr__", (std::string (*)(const grid_info_t&))&to_string<grid_info_t>)
        .def("create_grid_info_bloch_mode", grid_info_t::create_grid_info_bloch_mode)
        .def("create_grid_info", grid_info_t::create_grid_info)
        .def("nyquist_sampling_rate", grid_info_t::nyquist_sampling_rate)
        .def("wafer_pov_k_space_boundary", grid_info_t::wafer_pov_k_space_boundary)
        .def("mask_pov_k_space_boundary", grid_info_t::mask_pov_k_space_boundary)
        .def("display", &grid_info_t::template display<real_t<T>>)
        .def("display", &grid_info_t::template display<complex_t<T>>)
    ;
    using FT = decltype(grid_info_t::fourier);
    py::class_<FT>(("fourier_info_" + suffix).c_str()).def(py::init<>())  
        .def_readwrite("start",&FT::start)
        .def_readwrite("step",&FT::step)
    ;
    using ST = decltype(grid_info_t::spatial);
    py::class_<ST>(("spatial_info_" + suffix).c_str()).def(py::init<>())  
        .def_readwrite("start",&ST::start)
        .def_readwrite("step",&ST::step)
    ;
}
regist_py(
    using std::array;
    init_stl_converters<       
        array<array<float, 1>, 2>,  array<size_t, 1>,
        array<array<float, 3>, 2>,
        array<array<std::complex<float>, 3>, 2>
    >();
    regist_grid_info<float, 1>();
    regist_grid_info<float, 2>();
    regist_grid_info<float, 3>();
);