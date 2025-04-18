#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>
#include <fft/spectrum_analysis.hpp>
template<class Image, class MetaData> 
struct init_image
{
    constexpr std::pair<Image, MetaData> operator()(MetaData info, const size_t USF = 1) const
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
    template<class TShape>constexpr Image operator()(TShape shape) const
    {
        auto prod = std::accumulate(
            shape.begin(), shape.end(), 1, 
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
extern const uca::backend<double>& backend;
template<class T, class Image = std::vector<T>> struct thin_mask
{
    using cT = complex_t<T>;
    using rT = real_t<T>;
    constexpr static init_image<Image, grid_info_in_dbu> gen_image{};
    static std::tuple<Image, Image, grid_info_in_dbu> edge_pixelization(const grid_info_in_dbu& info, const std::vector<poly_dbu>& polys, size_t USF = 8, rT dissect_coef = 0.5)
    {
        auto [x, mask_info] = gen_image(info, USF);
        auto y = x;
        const auto& start = mask_info.spatial.start;
        const auto& step = mask_info.spatial.step;
        auto roi = cutline_dbu{point_dbu{0,0}, step * mask_info.tilesize};
        auto interpolate_to = [&](Image& im, cutline_dbu edge, point_dbu norm_dir){
            rT sign = rT(-1) * (norm_dir[0] + norm_dir[1]);
            int dirx = 1;
            int diry = 1;
            dissect_loop<point_dbu::value_type, 2>(edge, step * dissect_coef, 
                [&](point_dbu current){
                    point_dbu index;
                    index = convert_to<point_dbu>(current / mask_info.spatial.step);
                    if(!full_compare<point_dbu, vec2<size_t>>::less(index, (mask_info.tilesize - point_dbu{dirx, diry}))) return;
                    auto delta = convert_to<vec2<rT>>(current - index * mask_info.spatial.step) / step;
                    auto [ix, iy] = index;
                    auto coefs = linear_interpolate<rT>:: template get_coefs<2>(delta) * (dissect_coef * 0.5 * sign);
                    im.at(iy * mask_info.tilesize[0] + ix)                 += coefs[0][0];
                    im.at(iy * mask_info.tilesize[0] + ix + dirx)          += coefs[1][0];
                    im.at((iy + diry) * mask_info.tilesize[0] + ix)        += coefs[0][1];
                    im.at((iy + diry) * mask_info.tilesize[0] + ix + dirx) += coefs[1][1];
                }
            );
        };
        for(const auto& poly : polys){
            for(auto edge : poly){
                edge -= start;
                auto [dx, dy] = unit_vector(edge);
                if(0 != dx && 0 == dy) interpolate_to(x, edge, {-dy, dx});
                else if(0 == dx && 0 != dy) interpolate_to(y, edge, {-dy, dx});
                else throw std::runtime_error("edge_pixelization failed at " + to_string(edge));
            }
        }
        return {x, y, mask_info};
    }
    static std::tuple<Image, grid_info_in_dbu> get_edge_from_rasterization(const grid_info_in_dbu& info, const Image& image, const cutline_dbu& cutline, size_t USF = 1)
    {
        const auto& [start, step] = info.spatial;
        const auto& [from, to] = (cutline - start);
        grid_info_in_dbu cutline_meta = info;
        cutline_meta.spatial.start = cutline[0];
        cutline_meta.spatial.step = step / USF;
        
        cutline_meta.tilesize = convert_to<vec2<size_t>>(((to - from + step)/step) * USF); 
        Image line = gen_image(cutline_meta.tilesize);
        dissect_loop<point_dbu::value_type, 2>(cutline, step * (rT(1.0)/USF), 
            [&](point_dbu current){
                point_dbu idx_cutline = convert_to<point_dbu>((current - cutline[0]) / cutline_meta.spatial.step);
                if(0 != idx_cutline[0] && 0 != idx_cutline[1]){
                    error_unclassified::out("invalid cutline :", cutline, " index : ", idx_cutline);
                    return;
                }
                point_dbu idx_image = convert_to<point_dbu>((current - start) / info.spatial.step);
                if(!full_compare<point_dbu, vec2<size_t>>::less(idx_image, info.tilesize)){
                    error_unclassified::out("cutline is too close to image border. ", std::make_tuple(cutline, current, start, step, info.tilesize));
                    return;
                }
                auto delta = convert_to<vec2<rT>>(current - start - idx_image * info.spatial.step);
                delta /= step;
                const auto coefx = linear_interpolate<rT>::get_coef(delta[0]);
                const auto coefy = linear_interpolate<rT>::get_coef(delta[1]);
                auto [ix, iy] = idx_image;
                line.at(idx_cutline[0] + idx_cutline[1]) = 
                image.at(iy * info.tilesize[0] + ix)           * coefx[0] * coefy[0] +
                image.at(iy * info.tilesize[0] + ix + 1)       * coefx[1] * coefy[0] +
                image.at((iy + 1) * info.tilesize[0] + ix)     * coefx[0] * coefy[1] +
                image.at((iy + 1) * info.tilesize[0] + ix + 1) * coefx[1] * coefy[1];
            }
        );
        auto image_center = cutline_meta.spatial.step * (cutline_meta.tilesize - 1) /2 ;
        auto cutline_center = (cutline[1] - cutline[0] - 1)/2 ;
        rT offset =  rT(image_center[0] -  cutline_center[0]) / cutline_meta.spatial.step[0];
        if(0 != offset)
        {
            //== shift 反而会导致点的偏移 
            // line.reserve(line.size() + 2);
            // shift<T, complex_t<T>, rT>(line.data(), line.size(), 1, offset, 0);
            // error_unclassified::out("    TODO shift ", dbu_to_um(double(cutline_center[0] - image_center[0]), 0.25), "(nm)", " image_center=", image_center, "(dbu) cutline_center=", cutline_center, "(dbu)");
        }
        return {line, cutline_meta};
    }

    static std::tuple<Image, grid_info_in_dbu> mask_image(const grid_info_in_dbu& info, const std::vector<poly_dbu>& polys, size_t USF, rT dissect_coef)
    {
        auto [x, y, mask_info] = edge_pixelization(info, polys, USF, dissect_coef);
        print_grid_start_step<grid_info_in_dbu, debug_print<thin_mask>>(mask_info, "    intergral image");
        backend.integral_y(mask_info.tilesize, x.data());
        backend.integral_x(mask_info.tilesize, y.data());
        backend.VtAdd(x.size(), x.data(), y.data());
        return {y, mask_info};
    }
};