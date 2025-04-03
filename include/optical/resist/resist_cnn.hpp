#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>

template<class T> 
inline std::vector<vec<T, 5>> get_feature_intensity_from_cutline(const cutline_dbu& cutline, int64_t design_cd_in_dbu,const std::vector<std::vector<T>>& yArray, const point_dbu start, point_dbu step)
{
    int axis = -1;
    if(cutline[0][1] == cutline[1][1]) axis = 0;
    else if(cutline[0][0] == cutline[1][0]) axis = 1;
    assert(-1 != axis);
    auto center = (cutline[0][axis] + cutline[0][axis]) / 2 - start[axis]; 
    vec<int64_t, 5> pos {
        center,
        center - design_cd_in_dbu/2,
        center + design_cd_in_dbu/2,
        cutline[0][axis] - start[axis],
        cutline[1][axis] - start[axis],
    };
    std::vector<vec<T, 5>> features;
    features.reserve(yArray.size());
    for(const auto& y : yArray){
        vec<T, 5> feature;
        std::transform(pos.begin(), pos.end(), feature.begin(), [&](int64_t n){
            T current = T(n) / T(step[axis]);
            return cubic_interpolate<T>::eval(current, y);
        });
        features.push_back(feature);
    }
    return features;
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