
#include <optical/optical_numerics.hpp>
#include <py_helper.hpp>

template<class T>void optical_numerics_test(rectangle<T> roi={0, -200, 400, 200},
    vec2<T> ambit={0, 0}, vec2<size_t> tilesize={160, 160}, 
    T maxNA = 0.5, T sigma = 0,  T lambda = 13.5
)
{
    std::cout  << std::endl << "* OPC mode" << std::endl;
    optical_numerics(roi, ambit, tilesize, maxNA, lambda).print();
    std::cout << "* Bloch mode" << std::endl;
    bloch_optical_numerics(roi, tilesize, maxNA, lambda).print();
}
int main()
{
    optical_numerics_test<float>();
}

BOOST_PYTHON_MODULE(lib_test_optical) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("optical_numerics_test", &optical_numerics_test<float>);
    py::def("optical_numerics_testV1", &main);
}