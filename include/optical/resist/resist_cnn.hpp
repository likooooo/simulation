#pragma once
#include <kernels/gauss_laguerre.hpp>
#include <fft/conv.hpp>
#include <mekil/mkl_wrapper.hpp>
#include <functional>
#include "resist_common.hpp"

template<class T>
struct resist_cnn
{
    using XN = std::vector<T>;
    using MatrixP = std::vector<XN>;
    using VectorQ = std::vector<T>;
};

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
};