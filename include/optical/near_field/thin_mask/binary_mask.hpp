#pragma once
#include <optical/near_field/thin_mask.hpp>

template<class T, class Image = std::vector<complex_t<T>>>
struct binary_mask
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    dbu_grid_start_step<rT> meta;
    constexpr static rT max_sigma = 0;
    binary_mask(cT clear_transmition, cT absorber_transmition)
    {
        grid_info_in_dbu info;
        shapes_dbu shapes;
        auto [mask, grid_info_in_dbu] = thin_mask<cT, Image>::mask_image(info, shapes, 1);
        mask *= (absorber_transmition - clear_transmition);
        mask += clear_transmition;
    }
};