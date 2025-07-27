#pragma once
#include <optical/simulation_grid_info.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>
// #include <kernels/kernel_loop.hpp>
#include <optical/kernel_loop_v2.hpp>

template<class T>
struct binary_mask
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    static std::vector<cT> create(cT clear_transmition, cT absorber_transmition, const grid_info<rT, 2>& info, const polys_vertex_dbu& shapes)
    {
        std::vector<cT> x(info.total_size(), clear_transmition);
        create(clear_transmition, absorber_transmition, x.data(), info, shapes);
        return x;
    }
    static void create(cT clear_transmition, cT absorber_transmition, cT* p, const grid_info<rT, 2>& info, const polys_vertex_dbu& shapes)
    {
        polys_vertex_dbu polys = shapes - info.spatial.start;
        for(poly_vertex_dbu& poly : polys){
            size_t n = poly.size();
            for (size_t i = 0; i < n; ++i) {
                point_dbu& p1 = poly.at(i);
                point_dbu& p2 = poly.at((i + 1) % n);
                if(1 == is_horizon(p1, p2)){
                    p1[1] -= info.spatial.step[1];
                    p2[1] -= info.spatial.step[1];
                }
                else if(-1 == is_vertical(p1, p2)){
                    p1[0] -= info.spatial.step[0];
                    p2[0] -= info.spatial.step[0];
                }
            }
        }
        kernels_v2::kernel_loop<rT, 2>(info.tilesize, 
            [&](const vec<rT, 2>& unused_center, const vec<size_t, 2>& indices){
                point_dbu current = info.spatial.step * indices;
                size_t i = info(indices);
                for(const poly_vertex_dbu& poly : polys){
                    if(is_point_in_clockwise_convex_polygon(poly, current)){
                        p[i] = absorber_transmition;
                        return;
                    }
                }
                p[i] = clear_transmition;
            }
        ); 
    }
    static bool is_point_in_clockwise_convex_polygon(const poly_vertex_dbu& polygon, const point_dbu p) 
    {
        size_t n = polygon.size();
        if(n == 0) return false;
        if(n < 3) {
            return is_point_on_segment(p, polygon.front(), polygon.back());
        }
        
        for (size_t i = 0; i < n; ++i) {
            const point_dbu& p1 = polygon.at(i);
            const point_dbu& p2 = polygon.at((i + 1) % n);
            if (cross_product(p2 - p1, p - p1) > 0)  return false;

        }
        return true;
    }

};