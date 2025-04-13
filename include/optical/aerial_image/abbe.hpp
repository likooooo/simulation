#pragma once
#include <optical/source/source.hpp>
#include <optical/polarization/polarization.hpp>

template<class T, class Image = std::vector<complex_t<T>>>struct abbe_imaging
{
    usint rT =real_t<T>;
    usint cT =complex_t<T>;
    template<size_t DIM> static std::tuple<vec<Image, DIM>, vec<cT, DIM>> entrance_pupil(
        const vec<Image, DIM>& E_near_field, const vec<cT, DIM>& open_mask_coef
        const matrix<cT, DIM, DIM>& polar_jones_matrix)
    {
        static_assert(DIM <= 3);
        std::vector<std::tuple<size_t, size_t, cT>> cells;
        cells.reserve(DIM * DIM);
        for(size_t y = 0; y < DIM; y++)
        for(size_t x = 0; x < DIM; x++){
            if(0 == polar_jones_matrix[y][x]) continue;
            cells.push_back(std::tuple<size_t, size_t, cT>(y, x, polar_jones_matrix[y][x]));
        }
        /**
         *                 EX(x, y)     [sqrt(2)/2, 0,         0]
         * In(x, y, z)  =  EY(x, y)  .* [0,         sqrt(2)/2, 0]  // 45度 线偏振 
         *                 EZ(x, y)     [0,         0,         0]
        */
       vec<cT, DIM> norm_data;
        vec<Image, DIM> output;
        for(const auto[y, x, val] : cells){
            auto& e = output.at(y);
            if(0 == e.size()) e = E_near_field.at(y) * polar_jones_matrix[y][x];
            else e += Exyz.at(y) * polar_jones_matrix[y][x];

            norm_data.at(y) += E_near_field.at(y).size() * open_mask_coef.at(x) * polar_jones_matrix[y][x];
        }
        
        return {output, norm_data};
    }
    template<size_t DIM> static vec<Image, DIM> exit_pupil(const vec<Image, DIM>& in, const vec<cT, DIM>& norms)
    {

    }
};