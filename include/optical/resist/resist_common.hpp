#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>
using cutline_feature_pos_dbu = vec<double, 5>;
template<class T> using terms_cutline = std::vector<std::vector<T>>;
template<class T> using terms_features_intensity = std::vector<vec<T, 5>>;
inline std::tuple<cutline_feature_pos_dbu, int> get_feature_pos_from_cutline(const cutline_dbu& cutline, double target_cd_in_dbu, const point_dbu start, point_dbu step)
{
    int axis = -1;
    if(cutline[0][1] == cutline[1][1]) axis = 0;
    else if(cutline[0][0] == cutline[1][0]) axis = 1;
    assert(-1 != axis);
    auto center = double(cutline[0][axis] + cutline[1][axis]) / 2 - start[axis]; 
    vec<double, 5> pos {
        center,
        center - target_cd_in_dbu/2,
        center + target_cd_in_dbu/2,
        double(cutline[0][axis] - start[axis]),
        double(cutline[1][axis] - start[axis]),
    };
    return {pos , axis};
}
template<class T> inline terms_features_intensity<T> get_feature_intensity_from_cutline(const cutline_dbu& cutline, double measured_cd_in_dbu,const terms_cutline<T>& yArray, const point_dbu start, point_dbu step)
{
    auto [pos, axis] = get_feature_pos_from_cutline(cutline, measured_cd_in_dbu, start, step);
    terms_features_intensity<T> features;
    features.reserve(yArray.size());
    for(const auto& y : yArray){
        vec<T, 5> feature;
        std::transform(pos.begin(), pos.end(), feature.begin(), [&](T n){
            T current = T(n) / T(step[axis]);
            return cubic_interpolate<T>::eval(current, y);
        });
        features.push_back(feature);
    }
    return features;
}

void display_cutline(const cutline_data& data, const std::vector<double>& cutline_image,const point_dbu& start_dbu,  const point_dbu& step_dbu, float dbu, float threshold = 0.5)
{
    auto cutline = data.cutline;
    auto start = dbu_to_um(convert_to<vec2<float>>(start_dbu), dbu);
    auto step = dbu_to_um(convert_to<vec2<float>>(step_dbu), dbu);
    auto center = dbu_to_um((cutline[0] + cutline[1]) / 2, dbu);
    auto [features_in_dbu, dir] = get_feature_pos_from_cutline(cutline, um_to_dbu(data.measured_cd, dbu), start_dbu, step_dbu);
    auto features = dbu_to_um(convert_to<vec<float, 5>>(features_in_dbu), dbu);
    features += start[dir];
    plot_curves(std::vector<std::vector<double>>{
        cutline_image, std::vector<double>{1}, std::vector<double>{threshold, threshold}, std::vector<double>{0.25, 0.25}}, 
        {start[dir], features[0], features[1], features[3]}, 
        {step[dir], step[dir] , features[2] - features[1], features[4] - features[3]}, 
        {"cutline (um)", "center", "on", "out"}, {"b--", "r-x", "g--x", "r--x"});  
}