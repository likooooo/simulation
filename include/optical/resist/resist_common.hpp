#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>

inline std::tuple<vec<int64_t, 5>, int> get_feature_pos_from_cutline(const cutline_dbu& cutline, int64_t target_cd_in_dbu, const point_dbu start, point_dbu step)
{
    int axis = -1;
    if(cutline[0][1] == cutline[1][1]) axis = 0;
    else if(cutline[0][0] == cutline[1][0]) axis = 1;
    assert(-1 != axis);
    auto center = (cutline[0][axis] + cutline[1][axis]) / 2 - start[axis]; 
    vec<int64_t, 5> pos {
        center,
        center - target_cd_in_dbu/2,
        center + target_cd_in_dbu/2,
        cutline[0][axis] - start[axis],
        cutline[1][axis] - start[axis],
    };
    return {pos , axis};
}
template<class T> inline std::vector<vec<T, 5>> get_feature_intensity_from_cutline(const cutline_dbu& cutline, int64_t measured_cd_in_dbu,const std::vector<std::vector<T>>& yArray, const point_dbu start, point_dbu step)
{
    auto [pos, axis] = get_feature_pos_from_cutline(cutline, measured_cd_in_dbu, start, step);
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
