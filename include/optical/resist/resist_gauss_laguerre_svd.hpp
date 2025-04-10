#pragma once
#include <kernels/gauss_laguerre.hpp>
#include <fft/conv.hpp>
#include <mekil/mkl_wrapper.hpp>
#include <functional>
#include "resist_common.hpp"

template<class T>struct resist_gauss_laguerre_svd
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
    static MatrixP cal_matrix_p(const std::vector<terms_features_intensity<T>>& features, const std::vector<T>& weights)
    {
        auto on_point_intensity_should_be_equal_to_threshold = [&](const terms_features_intensity<T>& f, size_t feature_index){
            std::vector<T> vec(f.size() + 1);
            std::transform(f.begin(), f.end(), vec.begin(),[&](const auto& n){return n.at(feature_index);});
            vec.back() = -1;
            return vec;
        };
        auto apply_out_product_to_p = [&](MatrixP& p, const std::vector<T>& vec, T weight){
            size_t index_p = 0;
            for(size_t y = 0; y < vec.size(); y++)
            for(size_t x = 0; x < vec.size(); x++, index_p++)
                p.at(index_p) += (vec.at(x) * vec.at(y) * weight);
        };
        const size_t N = features.front().size() + 1;
        MatrixP p(N * N, 0);
        for(size_t i = 0; i < features.size(); i++){
            apply_out_product_to_p(p, on_point_intensity_should_be_equal_to_threshold(features.at(i), 1), weights.at(i));
            apply_out_product_to_p(p, on_point_intensity_should_be_equal_to_threshold(features.at(i), 2), weights.at(i));
        }
        return p;
    }
    using VectorQ = std::vector<T>;
    static VectorQ cal_vector_q(const std::vector<terms_features_intensity<T>>& features, const std::vector<T>& weights)
    {
        auto in_out_point_intensity_should_be_high_contrast = [&](VectorQ& q, const terms_features_intensity<T>& f, T weight){
            std::vector<T> vec(f.size() + 1);
            for(size_t i = 0; i < f.size(); i++){
                const auto& n = f.at(i);
                q.at(i) += (T(-1) * (n.at(0) - n.at(3) - n.at(4)) * weight);
            }
        };
        const size_t N = features.front().size() + 1;
        VectorQ q(N, 0);
        for(size_t i = 0; i < features.size(); i++)
            in_out_point_intensity_should_be_high_contrast(q, features.at(i), weights.at(i));
        return q;
    }
    static std::tuple<std::vector<T>, std::vector<T>> equality_constraint(const std::vector<terms_features_intensity<T>>& features, T threshold)
    {
        size_t N = features.front().size() + 1;

        std::vector<T> b(3 + 1);
        std::vector<T> A(N * b.size());
        size_t r = 0;
        std::fill(A.begin(), A.begin() + features.front().size(), 0);
        A.at(N - 1) =  1;
        b.at(0) = threshold;
        for(const auto& f : features){
            auto it = A.begin() + N;
            std::transform(f.begin(), f.end(), it, it,[](const vec<T, 5>& v, T n){return n + v.at(0);});
            b.at(1) += 1;
            it = A.begin() + 2 * N;
            std::transform(f.begin(), f.end(), it, it,[](const vec<T, 5>& v, T n){return n + v.at(1) + v.at(2);});
            b.at(2) += (threshold + threshold);
            it = A.begin() + 3 * N;
            std::transform(f.begin(), f.end(), it, it,[](const vec<T, 5>& v, T n){return n + v.at(3) + v.at(4);});
            b.at(3) += 0;
        }

        return {A, b};
    }
    
    static MatrixP cal_matrix_p(const std::vector<terms_dense_features_intensity<T>>& features, const std::vector<T>& weights)
    {
        auto on_point_intensity_should_be_equal_to_threshold = [&](const terms_dense_features_intensity<T>& f, size_t feature_index){
            std::vector<T> vec(f.size() + 1);
            std::transform(f.begin(), f.end(), vec.begin(),[&](const auto& n){return n.at(1).at(feature_index);});
            vec.back() = -1;
            return vec;
        };
        auto apply_out_product_to_p = [&](MatrixP& p, const std::vector<T>& vec, T weight){
            size_t index_p = 0;
            for(size_t y = 0; y < vec.size(); y++)
            for(size_t x = 0; x < vec.size(); x++, index_p++)
                p.at(index_p) += (vec.at(x) * vec.at(y) * weight);
        };
        const size_t N = features.front().size() + 1;
        MatrixP p(N * N, 0);
        for(size_t i = 0; i < features.size(); i++){
            apply_out_product_to_p(p, on_point_intensity_should_be_equal_to_threshold(features.at(i), 0), weights.at(i));
            apply_out_product_to_p(p, on_point_intensity_should_be_equal_to_threshold(features.at(i), 1), weights.at(i));
        }
        return p;
    }
    static VectorQ cal_vector_q(const std::vector<terms_dense_features_intensity<T>>& features, const std::vector<T>& weights)
    {
        //== linear miniminze
        auto in_out_point_intensity_should_be_high_contrast = [&](VectorQ& q, const terms_dense_features_intensity<T>& f, T weight){
            std::vector<T> vec(f.size() + 1);
            for(size_t i = 0; i < f.size(); i++){
                const auto& n = f.at(i);
                T in_data_sum = std::accumulate(n.at(0).begin(), n.at(0).end(), T(0), [](T sum, T a){return sum + a;}) / n.at(0).size();
                T out_data_sum = std::accumulate(n.at(2).begin(), n.at(2).end(), T(0), [](T sum, T a){return sum + a;}) / n.at(2).size();
                q.at(i) += ((out_data_sum - in_data_sum) * weight);
            }
        };
        const size_t N = features.front().size() + 1;
        VectorQ q(N, 0);
        for(size_t i = 0; i < features.size(); i++)
            in_out_point_intensity_should_be_high_contrast(q, features.at(i), weights.at(i));
        return q;
    }
    static std::tuple<std::vector<T>, std::vector<T>> equality_constraint(const std::vector<terms_dense_features_intensity<T>>& features, T threshold)
    {
        size_t N = features.front().size() + 1;

        std::vector<T> b(3 + 1);
        std::vector<T> A(N * b.size());
        size_t r = 0;
        std::fill(A.begin(), A.begin() + features.front().size(), 0);
        A.at(N - 1) =  1;
        b.at(0) = threshold;
        for(const auto& f : features){
            auto it = A.begin() + N;
            std::transform(f.begin(), f.end(), it, it,[](const dense_intensity_feature<T>& v, T n){return n + v.at(0).at(v.at(0).size() / 2);});
            b.at(1) += 1;
            it = A.begin() + 2 * N;
            std::transform(f.begin(), f.end(), it, it,[](const dense_intensity_feature<T>& v, T n){
                    assert(v.at(1).size() == 2);
                    return n + v.at(1).front() + v.at(1).back();
                }
            );
            b.at(2) += (threshold + threshold);
            it = A.begin() + 3 * N;
            std::transform(f.begin(), f.end(), it, it,[](const dense_intensity_feature<T>& v, T n){return n + v.at(2).front() + v.at(2).back();});
            b.at(3) += 0;
        }

        return {A, b};
    }
    static X calib_svd_kkt(const MatrixP& p, const VectorQ& q, const std::vector<T>& A, const std::vector<T>& b)
    {
        auto py_p = create_ndarray_from_vector(p, {int(q.size()), int(q.size())});
        auto py_q = create_ndarray_from_vector(q, {int(q.size())});
        auto py_A = create_ndarray_from_vector(A, {int(A.size()/b.size()), int(b.size())});
        auto py_b = create_ndarray_from_vector(b, {int(b.size())});

        auto x_in_py = convert_to<std::vector<T>>(py_plugin::ref()["optimize"]["optimization_with_kkt"](py_p, py_q, py_A, py_b));
        assert(x_in_py.size() == q.size());
        return x_in_py;
    }
    template<class TFeature>
    static X calib_svd_kkt(const std::vector<cutline_data>& gauges, const TFeature& features, T threshold)
    {
        assert(features.size() == gauges.size());
        std::vector<T> weights(gauges.size());
        std::transform(gauges.begin(), gauges.end(), weights.begin(), [](const cutline_data& data){return T(data.weight);});
        MatrixP p = cal_matrix_p(features, weights);
        VectorQ q = cal_vector_q(features, weights);
        auto [A, b] = equality_constraint(features, threshold);
        return calib_svd_kkt(p, q, A, b);
    }
    
    using VectorL = std::vector<T>;
    using VectorU = VectorL;
    using MatrixA = std::vector<std::vector<T>>;

    // static std::tuple<MatrixA, VectorL, VectorU> cal_constrains(size_t N, vec2<T> threshold_boundary)
    // {
    //     MatrixA ac;
    //     VectorL lb;
    //     VectorU ub;
    //     std::vector<T> constrain(N);
    //     N.back() = 1;

    //     ac.push_back(constrain);
    //     lb.push_back(threshold_boundary[0]);
    //     ub.push_back(threshold_boundary[1]);
    // }
};