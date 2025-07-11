#include <optical/near_field/thin_mask/binary_mask.hpp>
#include <py_plugin.h>

 void test_thin_mask(const std::string& path = "/home/like/doc/Documents/YuWei/gds/gds/case11.gds")
{
    py_plugin::call<void>("gds_io", "plot_gds", path);
}

int main()
{
    py_engine::init();
    // thin_mask();
}
BOOST_PYTHON_MODULE(lib_test_nearfield) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    // py::def("thin_mask", &thin_mask);
}