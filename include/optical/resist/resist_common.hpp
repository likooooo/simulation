#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/geometry.hpp>
#include <optical/polynomials.hpp>
#include <optical/clip.hpp>

using cutline_feature_pos_dbu = vec<double, 5>;
template<class T> using terms_cutline = std::vector<std::vector<T>>;
template<class T> using terms_features_intensity = std::vector<vec<T, 5>>;
template<class T> using aerial_slice = std::vector<T>;

inline int get_cutline_dir(const cutline_dbu& cutline)
{
    int axis = -1;
    if(cutline[0][1] == cutline[1][1]) axis = 0;
    else if(cutline[0][0] == cutline[1][0]) axis = 1;
    assert(-1 != axis);
    return axis;
}

inline std::tuple<cutline_feature_pos_dbu, int> get_feature_pos_from_cutline(const cutline_dbu& cutline, double target_cd_in_dbu, const point_dbu start, point_dbu step, double dbu)
{
    int axis = get_cutline_dir(cutline);
    // auto center = double(cutline[0][axis] + cutline[1][axis] -1) / 2 - start[axis]; 
    auto center = double(cutline[1][axis] - cutline[0][axis]) / 2 -(dbu_to_um<double>((const double)step[axis], dbu) * 1e3); 

    cutline_feature_pos_dbu pos {
        center,
        center - target_cd_in_dbu/2,
        center + target_cd_in_dbu/2,
        double(cutline[0][axis] - start[axis]),
        double(cutline[1][axis] - start[axis]),
    };
    return {pos , axis};
}
template<class T> inline std::tuple<vec2<T>, int> get_feature_on_measured_from_cutline(const cutline_dbu& cutline, double target_cd_in_dbu, const point_dbu start, point_dbu step, double dbu)
{
    int axis = get_cutline_dir(cutline);
    T center = T(cutline[1][axis] - cutline[0][axis]) / 2 -(dbu_to_um<T>((const T)step[axis], dbu) * 1e3); 
    T lhs = T(center - target_cd_in_dbu/2 - start[axis]) / T(step[axis]);
    T rhs = T(center + target_cd_in_dbu/2 - start[axis]) / T(step[axis]);
    return {vec2<T>{lhs, rhs}, axis};
}
template<class T> inline terms_features_intensity<T> get_feature_intensity_from_cutline(const cutline_dbu& cutline, double measured_cd_in_dbu, const terms_cutline<T>& yArray, const point_dbu start, point_dbu step, double dbu)
{
    auto [pos, axis] = get_feature_pos_from_cutline(cutline, measured_cd_in_dbu, start, step, dbu);
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
template<class T> inline vec3<std::vector<T>> get_dense_feature_intensity_from_cutline(const cutline_dbu& cutline, double measured_cd_in_dbu, const terms_cutline<T>& yArray, const point_dbu start, point_dbu step, double dbu)
{
    auto [pos, axis] = get_feature_on_measured_from_cutline(cutline, measured_cd_in_dbu, start, step, dbu);
    std::vector<T> in, on, out;

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

template<class T> inline aerial_slice<T> get_resist_slice_image(const terms_cutline<T>& terms, const std::vector<T>& coefs)
{
    assert(terms.size() == coefs.size());
    aerial_slice<T> sum = terms.front() * coefs.front();
    for(size_t i = 1; i < terms.size(); i++)
        sum += (terms.at(i) * coefs.at(i));
    return sum;
}
template<class T> inline std::vector<T> find_root_on_slice(aerial_slice<T> slice, T threshold, T start = 0, T step = 1)
{
    slice -= threshold;
    aerial_slice<T> shifted; 
    shifted.reserve(slice.size());
    std::copy(slice.begin() + 1, slice.end(), std::back_insert_iterator(shifted));
    //== period-border
    // shifted.push_back(slice.back());
    //== absorber-border
    shifted.push_back(slice.back());
    shifted *= slice;

    std::vector<T> roots; roots.reserve(4);
    for(size_t i = 0; i < shifted.size() - 1; i++){
        if(shifted.at(i) > 0) continue;
        T x = start + (step * T(i) + step * slice.at(i)/(slice.at(i) - slice.at(i + 1)));
        roots.push_back(x);
    }
    return roots;
}
template<class T> std::tuple<vec2<T>, bool> get_on_measured_cd_positive_slice(aerial_slice<T> slice, T threshold, const cutline_data& data, T start = 0, T step = 1)
{
    auto find_actual_symmetry_point =[&](const std::vector<T>& roots){
        T symmetry_point_guess = start + (step * (slice.size() - 1) / 2);
        for(int i = 0; i < roots.size() - 1; i++){
            T lhs = roots.at(i);
            T rhs = roots.at(i + 1);
            if((lhs - symmetry_point_guess) * (rhs - symmetry_point_guess) < 0) return i; 
        }
        return -1;
    };

    auto [it_lb, it_ub] = std::minmax_element(slice.begin(), slice.end());
    T lb = *it_lb, ub = *it_ub;
    bool is_ghost_cd = (lb >= threshold && ub <= threshold);
    std::vector<T> roots = find_root_on_slice(slice,  is_ghost_cd ?  (lb + ub) / 2 : threshold, start, step);
    // if(roots.size() == 1 && !is_ghost_cd) 
    //     roots = find_root_on_slice(slice, (lb + ub) / 2, start, step);
    if(roots.size() == 0) {
        error_unclassified::out("find root failed. ", " threshold is ",threshold, " mean is ",(lb + ub) / 2," aerial slice is ",slice);
        return {vec2<T>(), false};
    }
    if(roots.size() == 1){
        // error_unclassified::out("find singal-root. ", " threshold is ",threshold, " mean is ",(lb + ub) / 2," aerial slice is ",slice);
        T symmetry_point = start + step * std::distance(slice.begin(), it_ub);
        T half_len = symmetry_point - roots.at(0);
        auto result = half_len > 0 ? vec2<T>{roots.at(0), symmetry_point + half_len} : vec2<T>{symmetry_point + half_len, roots.at(0)};
        return {result, is_ghost_cd};
    }
    int center_index = find_actual_symmetry_point(roots);
    if(-1 == center_index){
        error_unclassified::out("find symmetry failed. ", " threshold is ",threshold, " aerial slice is ",slice);
        return {vec2<T>(), false};
    }
    if(is_ghost_cd){
        T symmetry_point = (roots[center_index] + roots[center_index + 1]) / 2;
        T half_len = data.measured_cd / 2;
        //== TODO :expand points here
        return {vec2<T>{-half_len, half_len} + symmetry_point, false};
    }
    return {vec2<T>{roots[center_index], roots[center_index + 1]}, true};
}
template<class T> inline void post_calib_analysis(std::vector<cutline_data>& gauges, const std::vector<terms_cutline<T>>& edges, const std::vector<T>& coefs, T threshold, point_dbu step, T dbu)
{
    assert(gauges.size() == edges.size() && coefs.size() == edges.front().size());
    for(size_t i = 0; i < gauges.size(); i++){
        auto resist_cutline = get_resist_slice_image(edges.at(i), coefs);
        auto& data = gauges.at(i);
        int axis = get_cutline_dir(data.cutline);
        auto [cd, cd_found] = get_on_measured_cd_positive_slice<T>(resist_cutline, threshold, data, data.cutline[0][axis], step[axis]);
        T sim_cd = cd[1] - cd[0];
        data.post_calib_results.push_back(sim_cd);
        data.post_calib_results.push_back(dbu_to_um(sim_cd - data.measured_cd, dbu * 1e3));
    }
}
inline void display_cutline(const cutline_data& data, const std::vector<double>& cutline_image,const point_dbu& start_dbu, const point_dbu& step_dbu, float dbu)
{
    auto cutline = data.cutline;
    auto start = dbu_to_um(convert_to<vec2<float>>(start_dbu), dbu);
    auto step = dbu_to_um(convert_to<vec2<float>>(step_dbu), dbu);
    auto center = dbu_to_um((cutline[0] + cutline[1] - 1) / 2, dbu);
    auto [features_in_dbu, dir] = get_feature_pos_from_cutline(cutline, data.measured_cd, start_dbu, step_dbu, dbu);
    auto features = dbu_to_um(convert_to<vec<float, 5>>(features_in_dbu), dbu);
    features += start[dir];
    plot_curves(std::vector<std::vector<double>>{
            cutline_image, std::vector<double>{double(-1 == data.polar ? 0 : 1)}, 
            std::vector<double>{
                linear_interpolate<double>::eval(features_in_dbu[1]/ step_dbu[dir], cutline_image), 
                linear_interpolate<double>::eval(features_in_dbu[2]/ step_dbu[dir], cutline_image)
                // threshold, threshold
            }, 
            std::vector<double>{0.1, 0.1}
        }, 
        {start[dir], features[0], features[1],               features[3]}, 
        {step[dir],  step[dir],   features[2] - features[1], features[4] - features[3]}, 
        {"cutline (um)", "center", "on", "out", "cd found(linear-interpolate)"}, {"b--", "r-x", "g--o", "b--x", "r--x"});  
}
inline void display_cutline_with_cd(const cutline_data& data, const std::vector<double>& cutline_image,const point_dbu& start_dbu,  const point_dbu& step_dbu, float dbu, float threshold = 0.5)
{
    auto cutline = data.cutline;
    auto start = dbu_to_um(convert_to<vec2<float>>(start_dbu), dbu);
    auto step = dbu_to_um(convert_to<vec2<float>>(step_dbu), dbu);
    auto center = dbu_to_um((cutline[0] + cutline[1] - 1) / 2, dbu);
    auto [features_in_dbu, dir] = get_feature_pos_from_cutline(cutline, data.measured_cd, start_dbu, step_dbu, dbu);
    auto features = dbu_to_um(convert_to<vec<float, 5>>(features_in_dbu), dbu);
    features += start[dir];
    auto [cd, found] = get_on_measured_cd_positive_slice<double>(cutline_image, threshold, data, start[dir], step[dir]);
    assert(found);
    plot_curves(std::vector<std::vector<double>>{
            cutline_image, std::vector<double>{double(-1 == data.polar ? 0 : 1)}, 
            std::vector<double>{
                //== TODO : 由于目前只有 mask, 变化太剧烈， 误差较大
                linear_interpolate<double>::eval(features_in_dbu[1]/ step_dbu[dir], cutline_image), 
                linear_interpolate<double>::eval(features_in_dbu[2]/ step_dbu[dir], cutline_image)
                // threshold, threshold
            }, 
            std::vector<double>{0.1, 0.1}, std::vector<double>{threshold, threshold}
        }, 
        {start[dir], features[0], features[1],               features[3],               float(cd[0])        }, 
        {step[dir],  step[dir],   features[2] - features[1], features[4] - features[3], float(cd[1] - cd[0])}, 
        {"cutline (um)", "center", "on", "out", "cd found(linear-interpolate)"}, {"b--", "r-x", "g--o", "b--x", "r--x"});  
}