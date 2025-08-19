#include <optical/source/source.hpp>
#include <py_helper.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

template<class T>
void regist_source_class()
{
    using source_point_t = source_point<T>;
    std::string suffix = get_numerical_type_str_suffix<T>();
    py::class_<source_point_t>(("source_point_" + suffix).c_str()).def(py::init<>())       
        .def_readwrite("intensity", &source_point_t::intensity)
        .def_readwrite("e_field_direction",&source_point_t::e_field_direction)
        .def_readwrite("DOP",&source_point_t::DOP)
        .def_readwrite("ellipticity",&source_point_t::ellipticity)
        .def_readwrite("sigmaxy",&source_point_t::sigmaxy)
        .def("k_vector", &source_point_t::k_vector)
        .def("polarization_state", &source_point_t::polarization_state)
        .def("sp_polarization_state", &source_point_t::sp_polarization_state)
        .def("get_crao_azimuth", &source_point_t::get_crao_azimuth)
        .def("get_current_polarized_dir", &source_point_t::get_current_polarized_dir)
        .def("__repr__", (std::string (*)(const source_point_t&))&to_string<source_point_t>)
    ;
    //== parametric source params
    using parametric_source_t = parametric_source<T>;
    using traditional_source = typename parametric_source_t::traditional_source_params;
    py::class_<traditional_source>(("traditional_source_" + suffix).c_str()).def(py::init<>())
        .def_readwrite("sigma", &traditional_source::sigma)
        .def_readwrite("center_x", &traditional_source::centerX)
        .def_readwrite("center_y", &traditional_source::centerY)
    ;
    using annular_source = typename parametric_source_t::annular_source_params;
    py::class_<annular_source>(("annular_source_" + suffix).c_str()).def(py::init<>())
        .def_readwrite("sigma_out", &annular_source::sigmaOut)
        .def_readwrite("sigma_in", &annular_source::sigmaIn)
        .def_readwrite("sigma_inner_shift_x", &annular_source::sigmaInnerShiftX)
        .def_readwrite("sigma_inner_shift_y", &annular_source::sigmaInnerShiftY)
        .def_readwrite("sigma_shift_x", &annular_source::sigmaShiftX)
        .def_readwrite("sigma_shift_y", &annular_source::sigmaShiftY)
    ;
    using dipole_fan_source = typename parametric_source_t::dipole_fan_source_params;
    py::class_<dipole_fan_source>(("dipole_fan_source_" + suffix).c_str()).def(py::init<>())
        .def_readwrite("annular", &dipole_fan_source::args)
        .def_readwrite("rotate_rad", &dipole_fan_source::rotAngle)
        .def_readwrite("span_rad", &dipole_fan_source::spanAngle)
    ;
    using quadratic_fan_source = typename parametric_source_t::quadratic_fan_source_params;
    py::class_<quadratic_fan_source>(("quadratic_fan_source_" + suffix).c_str()).def(py::init<>())
        .def_readwrite("annular", &quadratic_fan_source::args)
        .def_readwrite("rotate_rad", &quadratic_fan_source::rotAngle)
        .def_readwrite("span_rad", &quadratic_fan_source::spanAngle)
    ;
    using dipole_leaf_source = typename parametric_source_t::dipole_leaf_source_params;
    py::class_<dipole_leaf_source>(("dipole_leaf_source_" + suffix).c_str()).def(py::init<>())
        .def_readwrite("sigma_D", &dipole_leaf_source::sigma_D)
        .def_readwrite("sigma_d", &dipole_leaf_source::sigma_d)
        .def_readwrite("rotate_rad", &dipole_leaf_source::rotAngle)
        .def_readwrite("sigma_shift_x", &dipole_leaf_source::sigmaShiftX)
        .def_readwrite("sigma_shift_y", &dipole_leaf_source::sigmaShiftY)
    ;
    using quadratic_leaf_source = typename parametric_source_t::quadratic_leaf_source_params;
    py::class_<quadratic_leaf_source>(("quadratic_leaf_source_" + suffix).c_str()).def(py::init<>())
        .def_readwrite("dipole_leaf", &quadratic_leaf_source::params)
        .def_readwrite("aspect_ratio", &quadratic_leaf_source::aspectRatio)
    ;
    using source_grid_t = source_grid<T>;
    py::class_<source_grid_t>(("source_grid_" + suffix).c_str())
        .def(py::init<>())
        .def(py::init<size_t>())
        .def("create_traditional_source", source_grid_t::template create<traditional_source>,
            (py::arg("size"), py::arg("source_params"), py::arg("e_field_direction") = 0, py::arg("ellipticity") = 0, py::arg("polarization") = 0))
        .def("create_annular_source", source_grid_t::template create<annular_source>,
            (py::arg("size"), py::arg("source_params"), py::arg("e_field_direction") = 0, py::arg("ellipticity") = 0, py::arg("polarization") = 0))
        .def("create_dipole_fan_source", source_grid_t::template create<dipole_fan_source>,
            (py::arg("size"), py::arg("source_params"), py::arg("e_field_direction") = 0, py::arg("ellipticity") = 0, py::arg("polarization") = 0))
        .def("create_quadratic_fan_source", source_grid_t::template create<quadratic_fan_source>,
            (py::arg("size"), py::arg("source_params"), py::arg("e_field_direction") = 0, py::arg("ellipticity") = 0, py::arg("polarization") = 0))
        .def("create_dipole_leaf_source", source_grid_t::template create<dipole_leaf_source>,
            (py::arg("size"), py::arg("source_params"), py::arg("e_field_direction") = 0, py::arg("ellipticity") = 0, py::arg("polarization") = 0))
        .def("create_quadratic_leaf_source", source_grid_t::template create<quadratic_leaf_source>,
            (py::arg("size"), py::arg("source_params"), py::arg("e_field_direction") = 0, py::arg("ellipticity") = 0, py::arg("polarization") = 0))
        .def("__repr__", (std::string (*)(const source_grid_t&))&to_string<source_grid_t>)  
        .def("plot_wafer_pov", &source_grid_t::plot_wafer_pov, py::arg("grid_info") = grid_info<T>())
        .def("plot_mask_pov", &source_grid_t::plot_mask_pov)
        .def_readwrite("source_points", &source_grid_t::source_points)
        .def("shift_dc", (void (source_grid_t::*)(T, T, T))&source_grid_t::shift_dc)
        .def("clear_invalid_source_points", &source_grid_t::clear_invalid_source_points)
        .def("get_dc_from_chief_ray", source_grid_t::get_dc_from_chief_ray)
    ;
    init_stl_converters<std::vector<source_point_t>>();
}
regist_py
(
    regist_source_class<float>();
);