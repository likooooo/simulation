#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <optical/polynomials.hpp>
#include <kernels/kernel_loop.hpp>
#include <assert.h>

template<class T, size_t M = 4> struct zernike_radial_table
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using radial = std::vector<rT>;
    using debug = debug_print<zernike_radial_table<T, M>>;
    using error = error_print<zernike_radial_table<T, M>>;
    constexpr static size_t size_M = M  + 1;
    constexpr static std::array<size_t, size_M> __init_L() 
    {
        std::array<size_t, size_M> L{0};
        for(size_t i =0; i <= M; i++) L.at(i) = (M -i)/2;
        return L;
    }
    constexpr static vec<size_t, size_M> L = __init_L();
    constexpr static size_t __table_count()
    {
        size_t sum = 0;
        for(size_t n : L) sum += (n + 1);
        return sum;
    } 
    constexpr static size_t table_count = __table_count();
    constexpr static size_t index(size_t m, size_t l)
    {
        size_t index = 0;
        for(size_t i = 0; i < m; i++) index += L.at(i) + 1;
        return index + l;
    }
    constexpr static rT theta(size_t m, rT fx, rT fy)
    {
        return m * std::atan2(fy, fx);
    }    
    constexpr static rT cal_norm(size_t m, size_t l)
    {
        return 2 * (m + 2 * l + 1);
    }

    std::array<radial, table_count> zk_radials;
    std::vector<std::pair<vec2<size_t>, std::vector<cT>>>  zernike_image;
    zernike_radial_table(size_t points) {
        for(auto& v : zk_radials) v.resize(points);
        init_L0();
        init_L1();
        init_rest();
        norm_all();
    }
    std::vector<std::pair<vec2<size_t>, std::vector<cT>>>& init_zernike_image(size_t size) 
    {
        if(0 < zernike_image.size() && size * size == zernike_image.front().second.size()) return zernike_image;

        size_t xsize = size, ysize = size;
        auto& zernike = *this;
        zernike_image.clear();
        zernike_image.reserve(zernike.size_M * zernike.size_M);
        float stepx = rT(2 * zernike.zk_radials.front().size()) / xsize;
        float stepy = rT(2 * zernike.zk_radials.front().size()) / ysize;
        cT center{xsize/2.0f, ysize/2.0f};
        for(size_t m = 0, i = 0; m < zernike.size_M; m++){
            for(size_t l = 0; l <= zernike.L.at(m); l++, i++){
                zernike_image.push_back(std::make_pair(vec2<size_t>{m,l}, std::vector<cT>(xsize * ysize)));
                cT* pData = zernike_image.back().second.data();
                kernels::center_zero_loop_square_r<rT, 2>({ysize, xsize}, {stepy, stepx}, [&](const std::array<rT, 2> pos, rT r){
                    r = std::sqrt(r);
                    if(r < zernike.zk_radials.at(i).size()) {
                        rT rVal = cubic_interpolate<rT>::eval(r, zernike.zk_radials.at(i));
                        rT theta = zernike.theta(m, pos[1], pos[0]);
                        *pData = std::exp(cT(0, theta)) * rVal / std::sqrt(zernike.cal_norm(m, l));
                    }
                    pData++;
                });
            }
        }
        return zernike_image;
    }
    radial apply_aberration_m0_to_pupil(std::vector<matrix2x3<cT>>& pupil_radial, const std::vector<std::tuple<size_t, size_t, cT>>& poly_coefs) const
    {
        radial weights{0};
        for(const auto [m, l, coef] : poly_coefs){
            if(0 != m || 0 == coef) continue;
            if(l >= L.at(m).size()){
                auto type_str = TypeReflection<decltype(*this)>;
                auto L_str = to_string(L);
                debug("%s's L=%s\n", type_str.c_str(), L_str.c_str());
                error(poly_coefs, {"* input args", "m", "l", "coefficients"});
                continue;
            }
            assert(coef.imag() == 0);
            weights += (zk_radials.at(index(m, l)) * coef.real());
        }
        if(pupil_radial.size() == weights.size()){
            for(size_t i = 0; i < pupil_radial.size(); i++){
                pupil_radial.at(i) *= std::exp(cT(0, 1) * weights.at(i));
            }
        }
        else{
            rT step = rT(weights.size()) /  rT(pupil_radial.size());
            for(size_t i = 0; i < pupil_radial.size(); i++){
                pupil_radial.at(i) *= std::exp(cT(0, 1) * cubic_interpolate<T>::eval(step * i, weights));
            }
        }
        return weights;
    }
   
    std::vector<cT> gen_aberration_pupil_image(const vec2<size_t>& shape, const vec2<rT>& step /* NA / (shape -1) */, const std::vector<std::tuple<size_t, size_t, cT>>& poly_coefs) const
    {
        print_table(poly_coefs, {"M", "L", "coef"});
        std::vector<cT> pupil_image(shape[0] * shape[1]);
        cT * p = pupil_image.data();
        kernels::center_zero_loop_square_r<rT, 2>(shape, step, [&](const vec2<rT>& fxy, rT r){
            r = std::sqrt(r);
            if(r > 1.0) {
                p++;
                return;
            }
            cT total_phase = 0;
            for(const auto [m, l, coef] : poly_coefs){
                const auto& zernike_r = zk_radials.at(index(m, l));
                std::cout << zernike_r << std::endl;
                rT r_pixel = cubic_interpolate<rT>::eval(r * zernike_r.size() , zernike_r);
                float theta = this->theta(m, fxy[1], fxy[0]);
                total_phase += (coef.real() * theta * r_pixel + coef.imag() * (0.5_PI - theta) * r_pixel); 

                // cT angle = std::pow(cT(c, s),  int(m));
                // total_phase += (angle.real() * coef.real() * r_pixel + angle.imag() * coef.imag() * r_pixel); 
            }
            *p = std::exp(cT(2_PI_I) * total_phase);
            p++;
        });
        return pupil_image;
    }
    // {
    //     // std::map<size_t, vec2<std::vector<cT>>> radial_intergral_L;
    //     // for(const auto [m, l, coef_cos_sin] : poly_coefs){
    //     //     if(0 == m || 0 == coef_cos_sin) continue;
    //     //     if(l >= L.at(m).size()){
    //     //         auto type_str = TypeReflection<decltype(*this)>;
    //     //         auto L_str = to_string(L);
    //     //         debug("%s's L=%s\n", type_str.c_str(), L_str.c_str());
    //     //         error(poly_coefs, {"* input args", "m", "l", "coefficients"});
    //     //         continue;
    //     //     }
    //     //     const auto& r = zk_radials.at(index(m, l));
    //     //     auto& [sum_cos, sum_sin] = radial_intergral_L[m];
    //     //     auto [cos_wgt, sin_wgt] = coef_cos_sin;
    //     //     if(0 == sum.size()){ 
    //     //         sum_cos = convert_t<std::vector<cT>>(r) * coef_cos_sin;
    //     //         sum_sin = convert_t<std::vector<cT>>(r) * std::conj(coef_cos_sin);
    //     //     }
    //     //     else{
    //     //         sum_cos += convert_t<std::vector<cT>>(r) * coef_cos_sin;
    //     //         sum_sin += convert_t<std::vector<cT>>(r) * std::conj(coef_cos_sin);
    //     //     }
    //     // }

    //     // std::vector<cT> pupil_image(shape[0] * shape[1]);
    //     // kernels::center_zero_loop_square_r<rT, 2>(shape, step, [](const std::array<rT, N>& fxy, const rT r_pow_2){
    //     //     if(r_pow_2 > 1.0) return;
    //     //     const auto [c, s] = fxy / std::sqrt(r_pow_2);
    //     //     cT angle_step(c, s);
    //     //     cT total_phase = 0;
            
    //     //     for(const auto [m, l, coef] : poly_coefs){

    //     //     }
    //     //     for(const auto& [m, intergral_L] : radial_intergral_L){
    //     //         cT a = std::pow(angle_step,  m - 1);
    //     //         const auto& [sum_cos, sum_sin] = intergral_L;
    //     //         total_phase += (sum_cos.at() * std::conj(a) + sum_sin * a) / rT(2);
    //     //     }
    //     // });
    //     // std::vector<cT> weights;
    //     // for(const auto [m, l, coef] : poly_coefs){
    //     //     if(0 == m || 0 == coef) continue;
    //     //     if(l >= L.at(m).size()){
    //     //         auto type_str = TypeReflection<decltype(*this)>;
    //     //         auto L_str = to_string(L);
    //     //         debug("%s's L=%s\n", type_str.c_str(), L_str.c_str());
    //     //         error(poly_coefs, {"* input args", "m", "l", "coefficients"});
    //     //         continue;
    //     //     }
    //     //     auto r = zk_radials.at(index(m, l));
    //     //     std::vector<cT> cos_sin(r.size(), coef);
    //     //     if(0 == weights.size()) weights = cos_sin * r;
    //     //     else weights += (cos_sin * r);
    //     // }
    //     // if(pupil_radial.size() == weights.size()){
    //     //     for(size_t i = 0; i < pupil_radial.size(); i++){
    //     //         pupil_radial.at(i) *= std::exp(cT(2_PI_I) * weights.at(i));
    //     //     }
    //     // }
    //     // else{
    //     //     rT step = rT(weights.size()) /  rT(pupil_radial.size());
    //     //     for(size_t i = 0; i < pupil_radial.size(); i++){
    //     //         pupil_radial.at(i) *= std::exp(cT(2_PI_I) * cubic_interpolate<T>::eval(step * i, weights));
    //     //     }
    //     // }
    //     // return pupil_image;
    //     // if(pupil_radial.size() == weights.size()){
    //     //     for(size_t i = 0; i < pupil_radial.size(); i++){
    //     //         pupil_radial.at(i) *= std::exp(cT(2_PI_I) * weights.at(i));
    //     //     }
    //     // }
    //     // else{
    //     //     rT step = rT(weights.size()) /  rT(pupil_radial.size());
    //     //     for(size_t i = 0; i < pupil_radial.size(); i++){
    //     //         pupil_radial.at(i) *= std::exp(cT(2_PI_I) * cubic_interpolate<T>::eval(step * i, weights));
    //     //     }
    //     // }
    //     return weights;

    // }
private:
    constexpr void init_L0()
    {
        for(size_t m = 0; m <= M; m++){
            radial& v = zk_radials.at(index(m, 0));
            std::iota(v.begin(), v.end(), rT(0));
            v /= rT(v.size() - 1);
            for(auto& r : v) r = std::pow(r, m);
        }
    }    
    constexpr void init_L1()
    {
        for(size_t m = 0; m <= M - 2; m++){
            radial& m0 = zk_radials.at(index(m, 0));
            radial& p0 = zk_radials.at(index(m + 2, 0));
            radial& m1 = zk_radials.at(index(m, 1));
            m1 = p0 * (m + 2) - m0 * (m + 1);
        }
    }
    constexpr void init_rest()
    {
        radial r(zk_radials.at(0).size());
        std::iota(r.begin(), r.end(), rT(0)); r /= rT(r.size() - 1);
        for(size_t m = 0; m <= M; m++)
        for(size_t l = 2; l <= L.at(m); l++){
            int q = m;
            int p = 2 * l + m;
            // std::cout << "1. " << std::make_tuple(q, p)<< std::endl;
            rT rK1 = rT(2.0) / ((p + q) * (p - q) * (p - 2));
            rT k2 = 2 * p * (p - 1) * (p - 2);
            rT k3 = -q * q * (p - 1) - p * (p - 1) * (p - 2);
            rT k4 = rT(-0.5) * (p * (p + q - 2) * (p - q - 2));
            // std::cout << "2. " << std::make_tuple(rK1, k2, k3, k4) << std::endl;

            const radial& m2 = zk_radials.at(index(m, l - 2));
            const radial& m1 = zk_radials.at(index(m, l - 1));
            radial& n = zk_radials.at(index(m, l));
            // std::cout << r << std::endl;
            // std::cout << m1 << std::endl;
            // std::cout << m2 << std::endl;
            n  = ((r * r * k2 + k3) * m1 + m2 * k4) * rK1; 
        }
    }

    constexpr void norm_all()
    {
        for(size_t m = 0; m <= M; m++)
        for(size_t l = 0; l <= L.at(m); l++){
            rT polydeg = cal_norm(m, l);
            zk_radials.at(index(m, l)) *=  std::sqrt(polydeg); 
        }
    }
};