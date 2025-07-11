#pragma once
#include <py_helper.hpp>
#include <type_traist_notebook/type_traist.hpp>

using point_dbu         = point_nd_dbu<2>;
using poly_vertex_dbu   = std::vector<point_dbu>;
using shapes_vertex_dbu = std::vector<poly_vertex_dbu>;

using cutline_dbu = vec2<point_dbu>;
using poly_dbu    = std::vector<cutline_dbu>;
using shapes_dbu  = std::vector<poly_dbu>;


template<class TCallback> inline void foreach_poly_points(np::array2df poly, TCallback&& callback_input_point)
{
    auto [pVertex, vertex_size] = ndarray_ref_no_padding<point_dbu>(poly);
    for(size_t i = 0; i < vertex_size; i++, pVertex++) callback_input_point(*pVertex);
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
template<class TCallback> inline void foreach_shape_points(np::list_array2d shapes, TCallback&& callback_input_point){
    foreach_shapes(shapes, [&](np::array2di poly){
        foreach_poly_points(poly, std::forward<TCallback>(callback_input_point));
    });
}
template<class TCallback> inline void foreach_shape_lines(np::list_array2d shapes, TCallback&& callback_input_two_points){
    foreach_shapes(shapes, [&](np::array2di poly){
        foreach_poly_lines(poly, std::forward<TCallback>(callback_input_two_points));
    });
}
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


template<class T, size_t N, class PixelFunc, size_t Dim = 0>
inline void __dissect_loop_impl(const vec2<vec<T, N>>& p, const vec<T, N>& step, PixelFunc&& func, std::array<T, N>& indices) 
{
    assert(step > 0);
    const auto&[from, to] = p;
    if constexpr (Dim < N) 
    {
        T lb = std::min(from[Dim], to[Dim]);
        T ub = std::max(from[Dim], to[Dim]);
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