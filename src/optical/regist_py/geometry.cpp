#include <optical/geometry/geo_manager.hpp>
#include <py_helper.hpp>

regist_py(
    
    using namespace std;
    init_stl_converters<vec2<array<int64_t, 1>>, vec2<array<float, 1>>, vec2<array<double, 1>>>();
    py::class_<geo_manager>("geo_manager")
        .def(py::init<>())    
        .def(py::init<const polys_vertex_dbu&>())      
        .def(py::init<const std::string&>())    
        .def(py::init<const std::string&, int>())    
        .def(py::init<const std::string&, const std::string&, int>())    
        .def("get_vertex", &geo_manager::get_vertex,  py::return_value_policy<py::reference_existing_object>())
        .def("get_all_edges", &geo_manager::get_all_edges)
        .def("get_poly_edges", &geo_manager::get_poly_edges)
        .def_readwrite("polys",&geo_manager::polys)
    ;
);