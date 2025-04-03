#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
template<class T, class TMeta> 
inline vec<T, 5> get_feature_intensity_from_cutline(const cutline_dbu& cutline, const std::vector<T>& y, const TMeta& )
{

}


template<class T>
struct resist_cnn
{
    using XN = std::vector<T>;
    using MatrixP = std::vector<XN>;
    using VectorQ = std::vector<T>;
};

template<class T>
struct resist_least_squares
{

};