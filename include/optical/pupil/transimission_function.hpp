#pragma once
#include <optical/simulation_grid_info.hpp>

template<class T> struct transimission_function
{
    using rT = real_t<T>;
    using cT = complex_t<T>;

    static rT k_in_corner_zero(size_t idx, size_t N, rT freq_step){return freq_step * rT(int(idx) - int(idx >= N/2) * int(N)) / rT(N);}
    static rT k_in_center_zero(int idx, size_t N, rT freq_step){return freq_step * rT(idx - int(idx >= N/2) * int(N)) / rT(N);}
    static cT dielectric(cT nk){return nk * nk;}
    static cT fz(vec2<rT> fxy, cT dielectric = cT(1))
    {
        return std::sqrt(cT(dielectric - fxy[0] * fxy[0] - fxy[1] * fxy[1]));
    }
    template<class T1 = rT>static T1 delta_phase_angle(rT freq, rT delta, T1 k)
    {
        return 2_PI * freq * delta * k;
    }
    static rT phase_modulate_phase_angle(rT freq, vec2<rT> delta_xy, vec2<rT> fxy)
    {
        return delta_phase_angle<rT>(freq, delta_xy[0], fxy[0]) + 
               delta_phase_angle<rT>(freq, delta_xy[1], fxy[1]);
    }
    static cT propagation_phase_angle(rT freq, rT deltaz, rT fz)
    {
        return delta_phase_angle<rT>(freq, deltaz, fz);
    }
    static cT propagation_phase_angle(rT freq, rT deltaz, cT fz)
    {
        return delta_phase_angle<cT>(freq, deltaz, fz);
    }
    static cT propagation_phase_angle(rT freq, rT deltaz, vec2<rT> fxy, cT dielectric = cT(1))
    {
        return propagation_phase_angle(freq, deltaz, fz(fxy, dielectric));
    }

    //==
    // freq     : NA/lambda / M
    // delta_xy : physical unit. same with lambda
    // fxy      : k(i, N, step)
    // 
    // 如果改用像素单位
    // freq     : = 1
    // delta_xy : 像素单位
    // fxy      : (i, N, 1)
    static cT phase_modulate(rT freq, vec2<rT> delta_xy, vec2<rT> fxy)
    {
        return std::exp(1_I * phase_modulate_phase_angle(freq, delta_xy, fxy))
    }
    static cT propagation(rT freq, rT deltaz, vec2<rT> fxy, cT dielectric = cT(1))
    {
        return std::exp(1_I * propagation_phase_angle(freq, deltaz, fxy, dielectric));
    }
    static cT propagation(rT freq, vec3<rT> delta, vec3<rT> kxyz)
    {
        auto [dx,dy,dz] = delta;
        auto [fx,fy,fz] = kxyz;
        return std::exp(
            1_I * phase_modulate_phase_angle(freq, {dx, dy}, {fx, fy}) + 
            1_I * propagation_phase_angle(freq, dz, fz)
        );
    }
    static cT propagation(rT freq, vec3<rT> delta, vec2<rT> fxy, cT dielectric = cT(1))
    {
        return std::exp(
            1_I * phase_modulate_phase_angle(freq, {delta[0], delta[1]}, fxy) + 
            1_I * propagation_phase_angle(freq, delta[2], fxy, dielectric)
        );
    }
    

    static matrix2x3<cT> sp_polorization_basis(vec2<rT> fxy, cT nk = cT(1))
    {
        auto [n, k] = nk;
        return
        {
            //== s-polor (TE)
            0, 1, 0, 
            //== p-polor (TM)
            fz(fxy, dielectric(nk)) / n , 0, -std::hypot(fxy[0], fxy[1]) / n 
            //== 
            // std::hypot(fxy[0], fxy[1]) / n, 0, fz(fxy, dielectric(nk)) / n
        };
    }
};