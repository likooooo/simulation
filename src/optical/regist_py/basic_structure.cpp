#include <optical/geometry/geo_manager.hpp>
#include <py_helper.hpp>

template<class T>
void regist_geometry_structure()
{
    using point        = vec<T, 2>;
    using poly_vertex  = std::vector<point>;
    using polys_vertex = std::vector<poly_vertex>;

    using line         = vec2<point>;
    using lines        = std::vector<line>;
    using polys_lines  = std::vector<lines>;
    using std::array;
    init_stl_converters<
        polys_vertex, polys_lines
    >();
}
regist_py(
    using namespace std;
    init_stl_converters<array<int64_t, 1>, array<float, 1>, array<double, 1>>();
    regist_geometry_structure<dbu_t>();
    regist_geometry_structure<float>();
    regist_geometry_structure<double>();
);