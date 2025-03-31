#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>
extern bool verbose;
template<class Image, class MetaData> 
struct init_image
{
    constexpr static std::pair<Image, MetaData> operator()(MetaData info, const size_t USF = 1)
    {
        info.tilesize *= USF;
        info.spatial.step /= typename MetaData::value_type(USF);
        info.fourier.step /= typename MetaData::value_type(USF);
        auto prod = std::accumulate(
            info.tilesize.begin(), info.tilesize.end(), 1, 
            [](auto a, auto b){return a * b;}
        );
        return {Image(prod), info};
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

template<class T, class TMeta, class Image = std::vector<T>> struct thin_mask
{
    using cT = complex_t<T>;
    using rT = real_t<T>;
    static std::tuple<Image, Image, TMeta> intergral_image(const TMeta& info, const std::vector<poly_dbu>& polys)
    {
        auto [x, mask_info] = init_image<Image, TMeta>{}(info, 8);
        print_grid_start_step<TMeta, debug_print<thin_mask>>(mask_info, "intergral image");
        auto y = x;
        const auto& start = mask_info.spatial.start;
        const auto& step = mask_info.spatial.step;
        auto roi = cutline_dbu{point_dbu{0,0}, step * mask_info.tilesize};
        auto interpolate_to = [&](Image& im, cutline_dbu edge, point_dbu norm_dir){
            int64_t sign = norm_dir[0] + norm_dir[1];
            dissect_loop<point_dbu::value_type, 2>(edge, step, 
                [&](point_dbu current){
                    auto index = convert_to<point_dbu>(floor(current / mask_info.spatial.step));
                    if(numerics_logic::operator<(index, (mask_info.tilesize - 1)))
                    {
                        auto delta = convert_to<vec2<rT>>(current - index * mask_info.spatial.step);
                        delta /= step;
                        vec2<rT> coefx = linear_interpolate<rT>::get_coef(delta[0]);
                        vec2<rT> coefy = linear_interpolate<rT>::get_coef(delta[1]);
                        auto [ix, iy] = index;
                        im.at(iy * mask_info.tilesize[0] + ix)           += coefx[0] * coefy[0] * sign;
                        im.at(iy * mask_info.tilesize[0] + ix + 1)       += coefx[1] * coefy[0] * sign;
                        im.at((iy + 1) * mask_info.tilesize[0] + ix)     += coefx[0] * coefy[1] * sign;
                        im.at((iy + 1) * mask_info.tilesize[0] + ix + 1) += coefx[1] * coefy[1] * sign;
                    }
                }
            );
        };
        debug_print<thin_mask>::out("step", mask_info.spatial.step);
        debug_print<thin_mask>::out("roi", roi);
        for(const auto& poly : polys){
            for(auto edge : poly){
                edge -= start;
                if(!is_edge_inside_domain(edge, roi)) continue;
                auto [dx, dy] = unit_vector(edge);
                if(0 != dx && 0 == dy) interpolate_to(x, edge, {-dy, dx});
                else if(0 == dx && 0 != dy) interpolate_to(y, edge, {-dy, dx});
                else throw std::runtime_error("intergral_image failed at " + to_string(edge));
            }
        }
        return {x, y, mask_info};
    }

};