#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
template<class Image, class TMeta = grid_start_step<real_t<typename Image::value_type>>> 
struct init_image
{
    using MetaData = TMeta;
    using rT = real_t<typename Image::value_type>;
    constexpr static std::pair<Image, MetaData> operator()(MetaData info, const size_t USF = 1)
    {
        info.tilesize *= USF;
        info.spatial.step /= rT(USF);
        info.fourier.step /= rT(USF);
        auto prod = std::accumulate(
            info.tilesize.begin(), info.tilesize.end(), 1, 
            [](auto a, auto b){return a * b;}
        );
        return Image(prod);
    }
};

namespace near_filed
{
    shapes_dbu load_shapes_from_file(const std::string& path, const std::string& cell_name, int layer_id)
    {
        py::tuple poly_and_holes = py_plugin::call<py::tuple>("klayout_op", "load_oas_vertexs", path, cell_name, layer_id);
        auto polys = convert_to<np::list_array2d>(poly_and_holes[0]);
        auto holes = convert_to<np::list_array2d>(poly_and_holes[1]);
        polys += holes;
        shapes_dbu shapes; shapes.reserve(len(polys));
        foreach_shapes(polys, [&](np::array2di poly){
            shapes.emplace_back().reserve(poly.shape(0));
            foreach_poly_lines(poly, [&](point_dbu from, point_dbu to){
                shapes.back().push_back(cutline_dbu{from, to});
            });
        });
        return shapes;
    }
}

template<class T, class Image = std::vector<T>>struct thin_mask
{
    using cT = complex_t<T>;
    using rT = real_t<T>;
    
};