#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>

template<class Image, class MetaData> 
struct init_image
{
    constexpr std::pair<Image, MetaData> operator()(MetaData info, const size_t USF = 1)
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
    static std::tuple<Image, Image, TMeta> intergral_image(const TMeta& info, const std::vector<poly_dbu>& polys, rT dissect_coef = 0.5)
    {
        auto [x, mask_info] = init_image<Image, TMeta>{}(info, 8);
        print_grid_start_step<TMeta, debug_print<thin_mask>>(mask_info, "intergral image");
        auto y = x;
        const auto& start = mask_info.spatial.start;
        const auto& step = mask_info.spatial.step;
        auto roi = cutline_dbu{point_dbu{0,0}, step * mask_info.tilesize};
        auto interpolate_to = [&](Image& im, cutline_dbu edge, point_dbu norm_dir){
            rT sign = (norm_dir[0] + norm_dir[1]);
            dissect_loop<point_dbu::value_type, 2>(edge, step * dissect_coef, 
                [&](point_dbu current){
                    point_dbu index;
                    if((norm_dir == point_dbu{-1, 0}) || (norm_dir == point_dbu{0, -1}))
                        index = convert_to<point_dbu>(current / mask_info.spatial.step);
                    else if((norm_dir == point_dbu{1, 0}) || (norm_dir == point_dbu{0, 1}))
                        index = convert_to<point_dbu>((current + mask_info.spatial.step) / mask_info.spatial.step);
                    else throw std::runtime_error("dissect_loop failed at " + to_string(current));

                    // auto index = convert_to<point_dbu>(floor(current / mask_info.spatial.step));
                    if(full_compare<point_dbu, vec2<size_t>>::less(index, (mask_info.tilesize - 1)))
                    {
                        auto delta = convert_to<vec2<rT>>(current - index * mask_info.spatial.step);
                        delta /= step;
                        vec2<rT> coefx = linear_interpolate<rT>::get_coef(delta[0]);
                        vec2<rT> coefy = linear_interpolate<rT>::get_coef(delta[1]);
                        auto [ix, iy] = index;
                        im.at(iy * mask_info.tilesize[0] + ix)           += coefx[0] * coefy[0] * dissect_coef * 0.5 * sign;
                        im.at(iy * mask_info.tilesize[0] + ix + 1)       += coefx[1] * coefy[0] * dissect_coef * 0.5 * sign;
                        im.at((iy + 1) * mask_info.tilesize[0] + ix)     += coefx[0] * coefy[1] * dissect_coef * 0.5 * sign;
                        im.at((iy + 1) * mask_info.tilesize[0] + ix + 1) += coefx[1] * coefy[1] * dissect_coef * 0.5 * sign;
                    }
                    // else
                    //     error_unclassified::out("out of tilebox ", to_string(vec2<rT>{0.00025, 0.00025} * (start + current)));
                }
            );
        };
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