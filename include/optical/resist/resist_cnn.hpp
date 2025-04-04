#pragma once
#include "resist_common.hpp"
#include <py_helper.hpp>

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
    using Y = std::vector<T>;
    using A = std::vector<T>;
    using X = std::vector<T>;
    static X calib_selected_terms(const std::vector<cutline_data>& gauges, const std::vector<terms_features_intensity<T>>& features, const std::vector<size_t>& term_enable)
    {
        assert(gauges.size() == features.size());
        assert(*std::max_element(term_enable.begin(), term_enable.end()) < features.front().size());
        // [terms coef..., threshold]
        X x(term_enable.size() + 1);
        // [0,0,1,1, 0,0,1,1, 0,0,1,1]
        Y y(features.size() * 4); 
        // [terms intensity..., -1] x2 (on)
        // [terms intensity..., 0] x2  (in - out)
        A a(x.size() * y.size(), -1);
        
        for(size_t i = 0; i < features.size(); i++){
            const terms_features_intensity<T>& feature = features.at(i);
            T weight = T(gauges.at(i).weight);
            weight = 1;
            y.at(4 * i + 2) = gauges.at(i).polar * weight;
            y.at(4 * i + 3) = gauges.at(i).polar * weight;

            T* row = a.data() + i * 4 * x.size();
            for(size_t ix : term_enable){
                const auto [in, on_lhs, on_rhs, out_lhs, out_rhs] = feature.at(ix);
                row[ix]                = weight * (on_lhs);
                row[x.size() + ix]     = weight * (on_rhs);
                row[2 * x.size() + ix] = weight * (in - out_lhs);
                row[3 * x.size() + ix] = weight * (in - out_rhs);
            }
        }
        auto py_a = create_ndarray_from_vector(a, {int(x.size()), int(y.size())});
        auto py_y = create_ndarray_from_vector(y, {int(y.size())});
        auto x_in_py = convert_to<std::vector<T>>(py_plugin::ref()["optimize"]["svd"](py_a, py_y));
        assert(x_in_py.size() == x.size());
        return x_in_py;
    }
    static X calib(const std::vector<cutline_data>& gauges, const std::vector<terms_features_intensity<T>>& features)
    {
        std::vector<size_t> term_enable(features.front().size());
        std::iota(term_enable.begin(), term_enable.end(), 0);
        return calib_selected_terms(gauges, features, term_enable);
    }
    static T calib_optical_threshold(const std::vector<cutline_data>& gauges, const std::vector<terms_features_intensity<T>>& features)
    {
        constexpr size_t index_of_optical_term = 0;
        std::vector<size_t> term_enable{index_of_optical_term};
        X x = calib_selected_terms(gauges, features, term_enable);
        assert(x.size() == 2);
        return x.back() / x.front();
    }
};