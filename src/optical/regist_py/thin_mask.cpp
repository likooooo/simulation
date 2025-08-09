#include <optical/near_field/thin_mask/binary_mask.hpp>
#include <optical/near_field/thin_mask/thin_mask.hpp>
#include <optical/near_field/diffraction.hpp>
#include <py_helper.hpp>
#include "backend.hpp"
template<class T, size_t dim = 2> void regist_thin_mask()
{
    using bm_t = binary_mask<T>;
    using cT = typename bm_t::cT;
    using rT = typename bm_t::rT;
    py::class_<bm_t>("binary_mask")
        .def("create", (std::vector<cT> (*)(cT , cT , const grid_info<rT, dim>&, const polys_vertex_dbu& ))&bm_t::create)
    ;
    py::class_<thin_mask<T, dim>>("thin_mask")
        .def(py::init<>())       
        .def("create", ( std::vector<cT> (*)(cT , cT , const grid_info<rT, dim>&, const lines_dbu& ))&thin_mask<T, dim>::create)
    ;

    if constexpr(dim == 2){
        py::class_<diffraction<T>>("diffraction")
            .def(py::init<const grid_info<rT, 2>&>())       
            .def("update_diffraction_source_points", &diffraction<T>::update_diffraction_source_points,  py::return_value_policy<py::reference_existing_object>())
            .def("get_imaging_pupil_intensity", &diffraction<T>::get_imaging_pupil_intensity)
            .def_readonly("grid_info", &diffraction<T>::gi)
            .def_readonly("source_points", &diffraction<T>::source_points)
        ;
    }
}
regist_py(
    regist_thin_mask<float, 2>();
    py::def("get_diffraction_order", get_diffraction_order<float>, (py::arg("grid_info"), py::arg("sigmaxy") = vec2<float>{0, 0}));
);