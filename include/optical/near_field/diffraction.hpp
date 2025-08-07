#pragma once
#include <optical/simulation_grid_info.hpp>
#include <optical/kernel_loop_v2.hpp>
#include <optical/source/source.hpp>

template<class rT> std::vector<vec2<int>> get_diffraction_order(const grid_info<rT, 2>& info, vec2<rT> offset_sigmaxy)
{
    auto [start, step] = info.fourier;
    start += offset_sigmaxy;
    vec2<int> lb = convert_to<vec2<int>>((vec2<rT>{-1, -1} + start) / step) - 1;
    vec2<int> ub = convert_to<vec2<int>>((vec2<rT>{1, 1}   + start) / step) + 1;

    std::vector<vec2<int>> orders;
    orders.reserve((ub[1] - ub[0]) * (lb[1] - lb[0]));
    for(int y = lb[1]; y < ub[1]; y++){
        for(int x = lb[0]; x < ub[0]; x++){
            vec2<int> order{x, y};
            if(1 < vector_norm((step * order) - start)) continue;
            orders.push_back(order);
        }
    }
    return orders;
}

template<class rT> void get_diffraction_source_grid(const grid_info<rT, 2>& info, const std::vector<complex_t<rT>>& fourier_spectrum, vec2<rT> offset_sigmaxy)
{
    std::vector<vec2<int>> orders = get_diffraction_order(info, offset_sigmaxy);
    auto [lb, ub] = std::minmax_element(orders.begin(), orders.end(), [](vec2<int> a, vec2<int> b){return full_compare<vec2<int>>::less(a, b);})
    source_grid<rT> sg;
    sg.basis = polarization_basis::TETM;
    sg.step = info.fourier.step;
    
}