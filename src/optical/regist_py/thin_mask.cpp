#include <optical/near_field/thin_mask/binary_mask.hpp>
#include <optical/near_field/thin_mask/thin_mask.hpp>
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
    py::class_<bm_t>("thin_mask")
        .def(py::init<>())       
        .def("create", ( std::vector<cT> (*)(cT , cT , const grid_info<rT, dim>&, const lines_dbu& ))&thin_mask<T, dim>::create)
    ;
}
regist_py(
    regist_thin_mask<float, 2>();
);