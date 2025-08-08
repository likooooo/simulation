#pragma once
#include <optical/simulation_grid_info.hpp>
#include <optical/kernel_loop_v2.hpp>
#include <optical/source/source.hpp>

template<class rT> std::vector<vec2<int>> get_diffraction_order(const grid_info<rT, 2>& info, vec2<rT> offset_sigmaxy)
{
    auto [start, step] = info.fourier; // Wafer P.O.V
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
template<class rT>int away_from_zero(rT x) {
    if (x > 0)
        return std::ceil(x);   // 向正无穷
    else if (x < 0)
        return std::floor(x);  // 向负无穷
    else
        return 0;
}
template<class rT>
struct diffraction
{
    using cT = real_t<rT>;
    grid_info<rT, 2> gi;
    std::vector<vec2<int>> desired_orders;
    std::vector<source_point<rT>> source_points;
    diffraction(const grid_info<rT, 2>& grid_info) : gi(grid_info), desired_orders(get_diffraction_order(grid_info, {0, 0}))
    {
        source_points.reserve(desired_orders.size());
    }
    const std::vector<source_point<rT>>& update_diffraction_source_points(const std::vector<complex_t<rT>>& fourier_spectrum_center_zero)
    {
        if(fourier_spectrum_center_zero.size() != gi.total_size()){
            std::cerr << "size not match with grid_info" << std::endl; 
            return source_points;
        }
        source_points.clear();
        auto center = gi.tilesize / 2;
        for(auto order : desired_orders){
            auto [ix, iy] = (order + center);
            cT coef = fourier_spectrum_center_zero.at(ix + iy * gi.tilesize.at(1));
            coef /= fourier_spectrum_center_zero.size();
            source_point<rT> sp;
            sp.intensity = coef * std::conj(coef);
            //== TODO : calc DOP here ??

            sp.e_field_direction += std::atan2<rT>(order[1], order[0]); 
            sp.sigmaxy = gi.fourier.step * order;
            source_points.push_back(sp);
        }
        return source_points;
    }
    std::vector<rT> get_imaging_pupil_intensity(vec2<size_t> shape)
    {
        std::vector<rT> image(shape[0] * shape[1], 0);
        vec2<rT> step = 2/convert_to<vec2<rT>>(shape);
        vec2<rT> half_pixel = step / 2;
        for(const auto& sp : source_points){
            vec2<rT> gird_pos = sp.sigmaxy / step + half_pixel;
            vec2<int> from = convert_to<vec2<int>>(gird_pos);
            vec2<int> to{int(away_from_zero(gird_pos[0])), int(away_from_zero(gird_pos[1]))};
            auto [movex, movey] = to - from;
            // int movex = 1, movey = 1;
            assert(1>= std::abs(movex));
            assert(1>= std::abs(movey));
            matrix2x2 coefs = linear_interpolate<rT>::get_coefs<2>(gird_pos - (from * step));
            vec2<size_t> idx = shape / 2 + from;
            image.at(shape[0] * idx[1] + idx[0])                   = coefs[0][0] * sp.intensity;
            image.at(shape[0] * idx[1] + idx[0] + movex)           = coefs[0][1] * sp.intensity;
            image.at(shape[0] * (idx[1] + movey) + idx[0])         = coefs[1][0] * sp.intensity;
            image.at(shape[0] * (idx[1] + movey) + idx[0] + movex) = coefs[1][1] * sp.intensity;
        }
        imshow(image);
        return image;
    }
};