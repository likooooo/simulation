#pragma once
#include <optical/geometry.hpp>

template<class TCallback> inline void foreach_poly_point(np::array2df poly, TCallback&& callback_input_two_points)
{
    auto [pVertex, vertex_size] = ndarray_ref_no_padding<point_dbu>(poly);
    for(size_t i = 0; i < vertex_size; i++) callback_input_two_points(pVertex[i]);
}
template<class TCallback> inline void foreach_poly_lines(np::array2df poly, TCallback&& callback_input_two_points)
{
    auto [pVertex, vertex_size] = ndarray_ref_no_padding<point_dbu>(poly);
    for(size_t i = 0; i < vertex_size - 1; i++) callback_input_two_points(pVertex[i], pVertex[i + 1]);
    callback_input_two_points(pVertex[vertex_size - 1], pVertex[0]);
}
template<class TCallback> inline void foreach_shapes(np::list_array2d shapes, TCallback&& callback_input_array2di){
    for(size_t i = 0; i < len(shapes); i++){
        np::array2di poly = py::extract<np::array2di>(shapes[i]);
        callback_input_array2di(poly);
    }
}

struct geo_manager
{
    polys_vertex_dbu polys;
    geo_manager() = default;

    explicit geo_manager(const polys_vertex_dbu& poly_list) : polys(poly_list)
    {
    }
private:
    explicit geo_manager(const py::tuple& poly_and_holes)
    {
        auto poly_list = convert_to<np::list_array2d>(poly_and_holes[0]);
        size_t n = len(poly_list);
        polys.reserve(n);
        foreach_shapes(poly_list, [&](np::array2di poly){
            polys.emplace_back().reserve(poly.shape(0));
            foreach_poly_point(poly, [&](const point_dbu& point){
                polys.back().push_back(point);
            });
        });
        // TODO : apply dbu
    }
public:
    geo_manager(const std::string& path):
        geo_manager(py_plugin::call<py::tuple>("klayout_op", "load_oas_vertexs", path, py::object(), py::object()))
    {
    }
    geo_manager(const std::string& path, int layer_id):
        geo_manager(py_plugin::call<py::tuple>("klayout_op", "load_oas_vertexs", path, py::object(), layer_id))
    {
    }
    geo_manager(const std::string& path, const std::string& cell_name, int layer_id) :
        geo_manager(py_plugin::call<py::tuple>("klayout_op", "load_oas_vertexs", path, cell_name, layer_id))
    {
    }

    const polys_vertex_dbu& get_vertex() const
    {
        return polys;
    }
    lines_dbu get_all_edges() const
    {
        lines_dbu lines;
        lines.reserve(polys.size() * 5);
        for(const auto& poly : polys){
            size_t n = poly.size();
            for(size_t i = 0; i < n; i++){
                lines.push_back(line_dbu{poly.at(i), poly.at((i + 1) % n)});
            }
        }
        return lines;
    }
    polys_lines_dbu get_poly_edges() const
    {
        polys_lines_dbu shapes; 
        shapes.reserve(polys.size());
        for(const auto& poly : polys){
            size_t n = poly.size();
            shapes.emplace_back().reserve(n);
            for(size_t i = 0; i < n; i++){
                shapes.back().push_back(line_dbu{poly.at(i), poly.at((i + 1) % n)});
            }
        }
        return shapes;
    }
};