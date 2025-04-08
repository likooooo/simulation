#pragma once
#include <kernels/gauss_laguerre.hpp>
#include <fft/conv.hpp>
#include <mekil/mkl_wrapper.hpp>
#include <functional>
#include "resist_common.hpp"

template<class T>struct resist_blackbox
{
    template<class FuncConvWithKernel> static std::vector<std::vector<T>> gauss_laguerre_conv_linear(const vec2<size_t> shape, T sigma, size_t max_associated_order, size_t max_laguerre_order, FuncConvWithKernel&& conv_image_with)
    {
        std::vector<std::vector<T>> linears;
        linears.reserve(max_laguerre_order * max_associated_order);
        for(size_t laguerre_order = 0; laguerre_order < max_laguerre_order; laguerre_order++)
        for(size_t associated_order = 0; associated_order < max_associated_order; associated_order++)
        {
            std::vector<T> kernel(shape[0] * shape[1]);
            kernels::gauss_laguerre<T, 2, false>(kernel.data(), shape, sigma, laguerre_order, associated_order);
            conv_image_with(kernel);
            linears.push_back(std::move(kernel));
        }
        return linears;
    }
    static std::vector<std::vector<T>> gauss_laguerre_conv_quadratic(const std::vector<std::vector<T>>& linear_results)
    {
        std::vector<std::vector<T>> quadratic_result;
        size_t N = linear_results.size();
        quadratic_result.reserve(N * (N + 1) / 2);
        // 对称正定矩阵，只用计算上三角部分
        for(int r = 0; r < N; r++)
        for(int c = r; c < N; c++){
            std::vector<T> prod;
            auto& a = linear_results.at(r);
            auto& b = linear_results.at(c);
            prod.reserve(a.size());
            assert(a.size() == b.size());
            std::transform(a.begin(), a.end(), b.begin(), std::back_insert_iterator(prod), [](T l, T r){return (l * std::conj<T>(r)).real();});
            quadratic_result.push_back(std::move(prod));
        }
        return quadratic_result;
    }
    static std::vector<T> quadratic_coefficients_diagonalize(const std::vector<T>& coef_of_quad, size_t eigen_value_count)
    {
        std::vector<T> eigenval(eigen_value_count);
        LAPACKE_theev(LAPACK_COL_MAJOR, 'V', 'U', eigen_value_count, coef_of_quad.data(), eigen_value_count, eigenval.data());
        return eigenval;
    }
    static std::tuple<std::vector<std::vector<T>>, size_t> gauss_laguerre(const std::vector<T>& input, vec2<size_t> shape, real_t<T> sigma, size_t max_associated_order, size_t max_laguerre_order)
    {
        auto linear = gauss_laguerre_conv_linear(shape, sigma, max_associated_order, max_laguerre_order, [&](std::vector<T>& kernel){
            std::vector<T> k = input;
            k.reserve(k.size() + 2 * shape[1]);
            kernel.reserve(kernel.size() + 2 * shape[1]);
            conv<T, complex_t<T>>(kernel.data(), k.data(), shape[1], shape[0]);
            CenterCornerFlip(kernel.data(), shape[0], shape[1]);
        });
        size_t N = linear.size();
        {
            auto quadratic = gauss_laguerre_conv_quadratic(linear);
            linear.reserve(linear.size() + quadratic.size());
            std::copy(quadratic.begin(), quadratic.end(), std::back_insert_iterator(linear));
        }
        return std::tuple<std::vector<std::vector<T>>, size_t>(std::move(linear), N);
    }

    using Y = std::vector<T>;
    using A = std::vector<T>;
    using X = std::vector<T>;
    static X calib_svd_selected_terms(const std::vector<cutline_data>& gauges, const std::vector<terms_features_intensity<T>>& features, const std::vector<size_t>& term_enable, T threshold)
    {
        assert(gauges.size() == features.size());
        assert(*std::max_element(term_enable.begin(), term_enable.end()) < features.front().size());
        // [terms coef...]
        X x(term_enable.size());
        // [threshold, threshold, 1,1, threshold, threshold,1,1]
        Y y(features.size() * 4); 
        // [terms intensity..., -1] x2 (on)
        // [terms intensity..., 0] x2  (in - out)
        A a(x.size() * y.size(), -1);
        
        for(size_t i = 0; i < features.size(); i++){
            const terms_features_intensity<T>& feature = features.at(i);
            T weight = T(gauges.at(i).weight);
            y.at(4 * i)     = weight * threshold;
            y.at(4 * i + 1) = weight * threshold;
            y.at(4 * i + 2) = gauges.at(i).polar * weight;
            y.at(4 * i + 3) = gauges.at(i).polar * weight;

            T* row = a.data() + i * 4 * x.size();
            for(size_t ix : term_enable){
                const auto [in, on_lhs, on_rhs, out_lhs, out_rhs] = feature.at(ix);
                row[ix]                = weight * (on_lhs);
                row[x.size() + ix]     = weight * (on_rhs);
                //== 如果没有下面两个条件, 解出来的 X 全都是 0, 因为 A*0 = 0
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
    static X calib_svd(const std::vector<cutline_data>& gauges, const std::vector<terms_features_intensity<T>>& features, T threshold)
    {
        std::vector<size_t> term_enable(features.front().size());
        std::iota(term_enable.begin(), term_enable.end(), 0);
        return calib_svd_selected_terms(gauges, features, term_enable, threshold);
    }
    using MatrixP = std::vector<T>;
    static MatrixP vec_outer_product(const std::vector<T>& vec)
    {
        MatrixP p(vec.size() * vec.size());
        size_t index_p = 0;
        for(size_t y = 0; y < vec.size(); y++)
        for(size_t x = 0; x < vec.size(); x++, index_p++)
            p.at(index_p) = vec.at(x) * vec.at(y);
        return p;
    }
    using VectorQ = std::vector<T>;
    static X calib_osqp_selected_terms(const std::vector<cutline_data>& gauges, const std::vector<terms_features_intensity<T>>& features, const std::vector<size_t>& term_enable, T threshold)
    {
        auto error_function = [](const terms_features_intensity<T>& f, size_t feature_index){
            std::vector<T> v(f.size() + 1);
            std::transform(f.begin(), f.end(), [](const auto& n){return n.at(feature_index);});
            v.back() = -1;
            return v;
        };
        MatrixP p = vec_outer_product(error_function(features.front(), 1)) * gauges.front().weight;
        p += (vec_outer_product(error_function(features.front(), 2)) * gauges.front().weight);
        for(size_t i = 1; i < features.size(); i++){
            p += (vec_outer_product(error_function(features.at(i), 1)) * gauges.at(i).weight);
            p += (vec_outer_product(error_function(features.at(i), 2)) * gauges.at(i).weight);
        }
        VectorQ q(size_t(std::sqrt(p.size())));
        for(const auto& feature : features){
            assert(feature.size() + 1 == q.size());
        }        
    }
};