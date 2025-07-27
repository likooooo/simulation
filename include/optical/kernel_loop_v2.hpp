#pragma once
#include <array>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <type_traist_notebook/type_traist.hpp>

namespace kernels_v2
{
    template<class TPos, size_t N, class PixelFunc, int Dim = N - 1>
    constexpr inline void __kernel_loop_impl(const std::array<size_t, N>& shape, const std::array<TPos, N>& center, PixelFunc&& func, std::array<size_t, N>& indices) 
    {
        if constexpr (Dim > -1) 
        {
            for (indices[Dim] = 0; indices[Dim] < shape[Dim]; ++indices[Dim]) 
            {
                __kernel_loop_impl<TPos, N, PixelFunc, Dim - 1>(shape, center, std::forward<PixelFunc>(func), indices);
            }
        } 
        else 
        {
            func(center, indices);
        }
    }

    template<class T, size_t N, class PixelFunc>
    constexpr inline void kernel_loop(const std::array<size_t, N>& shape, PixelFunc&& func) 
    {
        std::array<T, N> center{};
        for(size_t i = 0; i < N; i++) center.at(i) = shape.at(i) / 2;
        
        std::array<size_t, N> indices{};
        __kernel_loop_impl<T>(shape, center, std::forward<PixelFunc>(func), indices);
    }
}