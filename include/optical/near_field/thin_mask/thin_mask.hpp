#pragma once
#include <optical/simulation_grid_info.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>

template<class T, size_t DIM>
struct thin_mask
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    //== shapes 必须要是封闭图形
    static std::vector<cT> create(cT clear_transmition, cT absorber_transmition, const grid_info<rT, DIM>& info, const lines_dbu& shapes)
    {
        std::vector<cT> x(info.total_size());
        create(clear_transmition, absorber_transmition, x.data(), info, shapes);
        return x;
    }
    static void create(cT clear_transmition, cT absorber_transmition, cT* p, const grid_info<rT, DIM>& info, const lines_dbu& shapes)
    {
        std::vector<cT> y(info.total_size());
        thin_mask<cT, DIM>::edge_pixelization(p, y.data(), info, shapes);
        vec2<size_t> shape;
        if constexpr (1 == DIM)
            shape = {info.tilesize[0], 1};
        else if constexpr(2 == DIM)
            shape = convert_to<vec2<size_t>>(info.tilesize);
        else 
            unreachable_constexpr_if<T>();
        // info.display(std::vector<cT>(p, p + info.total_size()));
        // info.display(y);
        get_math_backend<cT>().integral_y(shape, p);
        get_math_backend<cT>().integral_x(shape, y.data());
        absorber_transmition -= clear_transmition;
        for(size_t i = 0; i < info.total_size(); i++)
            p[i] = (p[i] + y[i]) * absorber_transmition + clear_transmition;
    }

    static void edge_pixelization(T* x, T* y, const grid_info<rT, DIM>& mask_info, const lines_dbu& polys)
    {
        constexpr rT dissect_coef = 1;
        const auto [start, step] = mask_info.spatial;
        auto interpolate_to = [&](T* im, line_dbu edge, rT sign){
            constexpr int dirx = 1;
            constexpr int diry = 1;
            edge -= start;
            dissect_loop<point_dbu::value_type, DIM>(edge, step * dissect_coef, 
                [&](point_dbu current){
                    point_dbu index;
                    index = convert_to<point_dbu>(current / step);
                    auto [ix, iy] = index;
                    auto delta = convert_to<vec2<rT>>(current - index * step) / step;
                    if(ix < mask_info.tilesize[0] - 1 && iy < mask_info.tilesize[1] - 1)
                    {
                        auto coefs = linear_interpolate<rT>:: template get_coefs<2>(delta) * (dissect_coef * (rT(1) / DIM) * sign);
                        im[iy * mask_info.tilesize[0] + ix]                 += coefs[0][0];
                        im[iy * mask_info.tilesize[0] + ix + dirx]          += coefs[1][0];
                        im[(iy + diry) * mask_info.tilesize[0] + ix]        += coefs[0][1];
                        im[(iy + diry) * mask_info.tilesize[0] + ix + dirx] += coefs[1][1];
                    }
                    else if(ix < mask_info.tilesize[0] - 1 && iy == mask_info.tilesize[1] - 1)
                    {
                        auto coefs = linear_interpolate<rT>::get_coef(delta[0]) * (dissect_coef * (rT(1) / DIM) * sign);
                        im[iy * mask_info.tilesize[0] + ix]                 += coefs[0];
                        im[iy * mask_info.tilesize[0] + ix + dirx]          += coefs[1];
                    }
                    else if (ix == mask_info.tilesize[0] - 1 && iy < mask_info.tilesize[1] - 1)
                    {
                        auto coefs = linear_interpolate<rT>::get_coef(delta[1]) * (dissect_coef * (rT(1) / DIM) * sign);
                        im[iy * mask_info.tilesize[0] + ix]          += coefs[0];
                        im[(iy + diry) * mask_info.tilesize[0] + ix] += coefs[1];
                    }
                }
            );
        };
        for(line_dbu edge : polys){
            const auto [dx, dy] = (edge[1] - edge[0]);
            //==  interpolate to x   interpolate to y
            //       -----                 +    -
            //     x:                  y:  +    -
            //       +++++                 +    -
            int sign = is_horizon(edge[0], edge[1]);
            if(0 != sign) {
                interpolate_to(x, edge, -sign);
                continue;
            }
            sign = is_vertical(edge[0], edge[1]);
            if(0 != sign) {
                interpolate_to(y, edge, sign);
                continue;
            }
            throw std::invalid_argument("support manhattan only. invalid edge " + to_string(edge));
        }
    }
};