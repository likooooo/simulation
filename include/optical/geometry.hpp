#pragma once
#include "simulation_grid_info.hpp"
#include <py_helper.hpp>
#include <type_traist_notebook/type_traist.hpp>

using point_dbu        = point_nd_dbu<2>;
using poly_vertex_dbu  = std::vector<point_dbu>;
using polys_vertex_dbu = std::vector<poly_vertex_dbu>;

using line_dbu         = vec2<point_dbu>;
using lines_dbu        = std::vector<line_dbu>;
using polys_lines_dbu  = std::vector<lines_dbu>;

template<class T, size_t N> vec<T, N> unit_vector(const vec<T, N>& from, const vec<T, N>& to)
{
    auto vec = to - from;
    vec /= std::sqrt(vector_norm(vec));
    return vec;
}
template<class T, size_t N> vec<T, N> unit_vector(const vec2<vec<T, N>>& line)
{
    const auto& [from, to] = line;
    return unit_vector(from, to);
}
template<class T, size_t N> vec<T, N> norm_vector(vec<T, N> from, vec<T, N> to)
{
    if constexpr(N == 2){
        //== rotate matrix (theta = PI/2)
        //[0, -1]
        //[1,  0]
        from = {-from[1], from[0]};
        to = {-to[1], to[0]};
        auto vec = to - from;
        return unit_vector(from, to);
    }else{
        unreachable_constexpr_if<>();
    }
}
template<class T, size_t N> vec<T, N> norm_vector(const vec2<vec<T, N>>& line)
{
    const auto& [from, to] = line;
    return norm_vector(from, to);
}
inline typename point_dbu::value_type cross_product(point_dbu v1, point_dbu v2) 
{
    return v1[0] * v2[1] - v1[1] * v2[0];
}
//== p ∈[p1, p2)
inline bool is_point_on_segment(point_dbu p, point_dbu p1, point_dbu p2) 
{
    point_dbu n = (p1 - p) * (p2 - p);
    return (p1 == p) || (n[0] < 0 && n[1] == 0) || (n[1] < 0 && n[0] == 0);
}
inline int is_horizon(point_dbu from, point_dbu to) 
{
    from = to - from;
    return std::clamp<int>(int(0 == from[1]) * from[0], -1, 1);
}
inline int is_vertical(point_dbu from, point_dbu to) 
{
    from = to - from;
    return std::clamp<int>(int(0 == from[0]) * from[1], -1, 1);
}

template<class T, size_t N, class PixelFunc, size_t Dim = 0>
inline void __dissect_loop_impl(const vec2<vec<T, N>>& p, const vec<T, N>& step, PixelFunc&& func, std::array<T, N>& indices) 
{
    const auto&[from, to] = p;
    if constexpr (Dim < N) 
    {
        assert(
            [&]()->bool
            { 
                if(step[Dim] < 0) return false; 
                if(step[Dim] == 0 && from[Dim] != to[Dim]) return false;
                return true;
            }()
        );
        T lb = std::min(from[Dim], to[Dim]);
        T ub = std::max(from[Dim], to[Dim]);
        indices[Dim] = lb;
        if(ub == lb) ub = lb + step[Dim];
        for (indices[Dim] = lb; indices[Dim] < ub; indices[Dim] += step[Dim]) 
        {
            __dissect_loop_impl<T, N, PixelFunc, Dim + 1>(p, step, std::forward<PixelFunc>(func), indices);
        }
    } 
    else 
    {
        func(indices);
    }
}

template<class T, size_t N, class TCallback> void dissect_loop(const vec2<vec<T, N>>& p, const vec<T, N>& step, TCallback&& callback_in_point)
{
    std::array<T, N> indices{};
    __dissect_loop_impl<T, N>(p, step, std::forward<TCallback>(callback_in_point), indices);
}   

template<class T, size_t N> bool point_in_domain(const vec<T, N>& p, const vec<T, N>& from, const vec<T, N>& to){
    return from <= p && p <= to;
}
template<class T, size_t N> vec<T, N> domain_clamp(vec<T, N> p, const vec<T, N>& from, const vec<T, N>& to){
    for(size_t i = 0; i < N; i++) {
        if constexpr(is_real_or_complex_v<T>){
            p.at(i) = std::clamp(p.at(i), from.at(i), to.at(i));
        }
        else{
            p.at(i) = domain_clamp(p.at(i), from.at(i), to.at(i));
        }
    }
    return p;
}
template<class T, size_t N> bool is_point_inside_domain(const vec<T, N>& p, const vec2<vec<T, N>>& domain)
{
    const auto& [from, to] = domain;
    return full_compare<vec<T, N>, vec<T, N>>::less_equal(from, p) && full_compare<vec<T, N>,vec<T, N>>::greater(p, from);
}
template<class T, size_t N> bool is_edge_inside_domain(const vec2<vec<T, N>>& p, const vec2<vec<T, N>>& domain)
{
    const auto& [p1, p2] = p;
    return is_point_inside_domain(p1, domain) && is_point_inside_domain(p2, domain);
}