#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <kernels/kernel_loop.hpp>
#include <optical/polynomials.hpp>
#include <assert.h>

template<class T>
using dual_vec = vec2<T>;

enum class polar_type
{
    X, Y, XY, TE, TM, LINEAR, UNPOLARIZED
};
template<polar_type p, class rT> inline rT get_polar_angle(rT fx, rT fy, rT linear_pol_angle = NAN)
{
    switch(p)
    {
        case polar_type::X: return 0;
        case polar_type::Y: return 0.5_PI;
        case polar_type::XY:  {
            if((fx - fy < 0 && fx + fy > 0) || (fx - fy > 0 && fx + fy < 0)) return 0.5_PI;
            else if ((fx - fy > 0 && fx + fy > 0) || (fx - fy < 0 && fx + fy < 0)) return 0;
            return 0.25_PI;
        }
        case polar_type::TE: return std::atan2(fy, fx) + 0.5_PI;
        case polar_type::TM: return std::atan2(fy, fx);
        case polar_type::LINEAR: return linear_pol_angle;
        default: return NAN;
    }
}

template<class T>
struct polarization
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using polar_vec = vec3<cT>;

    // 光的偏振态可以用两个正交的偏振方向表示
    static dual_vec<polar_vec> basis(rT polar_angle)
    {
        using std::sin, std::cos;
        return !std::isnan(polar_angle) ? dual_vec<polar_vec>{
            cT( cos(polar_angle)), cT(sin(polar_angle)), cT(0),
            cT(-sin(polar_angle)), cT(cos(polar_angle)), cT(0),
        } : 
        dual_vec<polar_vec>{
            cT(1), cT(0), cT(0),
            cT(0), cT(1), cT(0)
        };
    }
    static dual_vec<rT> componetes(rT intensity, rT degree_of_polarization)
    {
        return linear_interpolate<rT>::get_coef(degree_of_polarization) * intensity;
    }
    static std::pair<dual_vec<polar_vec>, dual_vec<rT>> polarization_state(rT intensity, rT dop, rT polar_angle)
    {
        return {basis(polar_angle, 0), componetes(intensity, dop)};
    }
    // 互强度矩阵描述光场中不同点的相干性，是部分相干光场的统计特性
    union mutual_intensity_matrix
    {
        matrix<cT, 2, 2> m;
        struct {
            cT xx, xy;
            cT yx, yy;
        };
    };

    static mutual_intensity_matrix density_matrix(dual_vec<polar_vec> polar_basis, dual_vec<rT> components)
    {
        mutual_intensity_matrix result{0};
        matrix<cT, 2, 2>& correlation_matrix = result.m;
        for(size_t i = 0; i < polar_basis.size(); i++){
            const auto [x, y, z] = polar_basis.at(i);
            assert(z == cT(0));
            correlation_matrix.at(0) += vec2<cT>{x, y} * (x * components.at(i));
            correlation_matrix.at(1) += vec2<cT>{x, y} * (y * components.at(i));
        }
        return result;
    }
    template<class TPolarFunc> static void source_density_matrix(
        matrix<cT, 2, 2>* output, const rT* pSource, rT src_threshold, 
        vec2<size_t> shape, vec2<rT> step, const dual_vec<rT>& components, // should be {0.5, 0.5} if unpolarized 
        TPolarFunc&& cal_polar_angle
    )
    {
        kernels::center_zero_loop_square_r(shape, step,
            [&](const vec2<rT> k, rT norm_k){
                if(*pSource > src_threshold){
                    const auto [ky, kx] = k;
                    *output = density_matrix(basis(cal_polar_angle(kx, ky)), components).m * (*pSource);
                }
                else{
                    *output = matrix<cT, 2, 2>{0};
                }
                output++;
                pSource++;
            }
        );
    }
    template<class TPolarFunc> static void source_density_matrix(
        matrix<rT, 2, 2>* output, const rT* pSource, rT src_threshold, 
        vec2<size_t> shape, vec2<rT> step, rT degree_of_polarization,
        TPolarFunc&& cal_polar_angle
    )
    {
        source_density_matrix(output, pSource, src_threshold, shape, step, linear_interpolate<rT>::get_coef(degree_of_polarization), std::forward<TPolarFunc>(cal_polar_angle));
    }
};