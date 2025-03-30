#pragma once
#include <py_helper.hpp>
#include <type_traist_notebook/type_traist.hpp>

using point_dbu = vec2<int64_t>;
using cutlinei = rectangle<int64_t>;
template<class T, class T1> void dbu_to_um(T& t, T1 dbu){t *= dbu;}
template<class T, class T1> void um_to_dbu(T& t, T1 dbu){t /= dbu;}
template<class T, class T1> T dbu_to_um(const T& t, T1 dbu){return t * dbu;}
template<class T, class T1> T um_to_dbu(const T& t, T1 dbu){return t / dbu;}

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
template<class T, size_t N> bool point_in_domain(const vec<T, N>& p, const vec<T, N>& from, const vec<T, N>& to){
    return from <= p && p <= to;
}
template<class T, size_t N> vec<T, N> point_in_domain_clamp(vec<T, N> p, const vec<T, N>& from, const vec<T, N>& to){
    for(size_t i = 0; i < N; i++) p.at(i) = std::clamp(p.at(i), from.at(i), to.at(i));
    return p;
}