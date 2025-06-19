#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include "snells_law.hpp"

//== proof of s components :  https://en.wikipedia.org/wiki/Fresnel_equations#s_components
//   proof of p components :  https://en.wikipedia.org/wiki/Fresnel_equations#p_components
template<class T, char polarization = 's'> struct fresnel_equations
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using snell = snells_law<T>;
    // E_input + E_refrecton = E_transmission
    // H_input * cos(theta_input) - H_refrecton * cos(theta_input) = H_transmission * cos(theta_transmission)
    constexpr static cT refrection_coef(cT n1, cT theta1, cT n2, cT theta2)
    {
        using std::cos;
        if constexpr(polarization == 's')
        {
            return (n1 * cos(theta1) - n2 * cos(theta2)) / (n1 * cos(theta1) + n2 * cos(theta2));
        }
        else if constexpr(polarization == 'p')
        {
            return (n2 * cos(theta1) - n1 * cos(theta2)) / (n2 * cos(theta1) + n1 * cos(theta2));
        }
        else
        {
            unreachable_constexpr_if<int>{};
        }
    }
    constexpr static cT transmission_coef(cT n1, cT theta1, cT n2, cT theta2)
    {
        using std::cos;
        if constexpr(polarization == 's')
        {
            return (n1 * cos(theta1) * rT(2)) / (n1 * cos(theta1) + n2 * cos(theta2));
        }
        else if constexpr(polarization == 'p')
        {
            return (n1 * cos(theta1) * rT(2)) / (n2 * cos(theta1) + n1 * cos(theta2));
        }
        else
        {
            unreachable_constexpr_if<int>{};
        }
    }
    constexpr static cT refrection_coef(cT n1, cT theta1, cT n2){return refrection_coef(n1, theta1, n2, snells_law(n1, n2, theta1));}
    constexpr static cT transmission_coef(cT n1, cT theta1, cT n2){return transmission_coef(n1, theta1, n2, snells_law(n1, n2, theta1));}
    constexpr static rT reflected_power(cT r){ return std::norm(r); }
    constexpr static rT transmitted_power(cT t, cT n1, cT theta1, cT n2, cT theta2)
    { 
        using std::cos, std::conj, std::norm;
        if constexpr(polarization == 's')
        {
            return  norm(t) * ((n2 * cos(theta2)).real()) / ((n1 * cos(theta1)).real());
        }
        else if constexpr(polarization == 'p')
        {
            return  norm(t) * ((n2 * conj(cos(theta2))).real()) / ((n1 * conj(cos(theta1))).real());
        }
        else
        {
            unreachable_constexpr_if<int>{};
        }
    }
    //== https://en.wikipedia.org/wiki/Fresnel_equations#Complex_amplitude_reflection_and_transmission_coefficients
    constexpr static bool check_refrection_transmission(cT r, cT t, complex_t<T> n1, complex_t<T> n2, rT eps = 1e-6)
    {
        if constexpr(polarization == 's')
        {
            return std::abs(cT(1) + r - t) < eps;
        }
        else if constexpr(polarization == 'p')
        {
            return std::abs(cT(1) + r - (n2 / n1 * t)) < eps;

        }
    }
    constexpr static bool check_power_refrection_transmission(cT r_power, cT t_power, rT eps = 1e-6)
    {
        // 1 = R + T;
        // std::cout << std::make_tuple(r_power, t_power, std::abs(r_power + t_power - cT(1))) << std::endl;
        return std::abs(r_power + t_power - cT(1)) < eps;
    }
};
// == https://en.wikipedia.org/wiki/Schlick%27s_approximation
template<class T> constexpr static inline complex_t<T> fast_refrection_unpolarized(T theta, complex_t<T> n1, complex_t<T> n2)
{
    complex_t<T> R0 = std::pow((n1-n2)/(n1 + n2), 2);
    return R0 + (1 - R0) * std::pow(1- std::cos(theta), 5);
}
template<class T> constexpr static inline std::tuple<bool, complex_t<T>> is_total_inner_reflection(complex_t<T> theta_in, complex_t<T> n1, complex_t<T> n2)
{
	assert(0 <= std::abs(theta_in) && std::abs(theta_in) <= 0.5_PI);
	if(std::norm(n1) <= std::norm(n2)){
        return {false, complex_t<T>(NAN)};
    }
	complex_t<T> critical_angle = std::asin(std::sin(0.5_PI) * n2/n1);
	return {true, critical_angle};
}

//== https://en.wikipedia.org/wiki/Fresnel_number#Application
template<class T> inline T fresnel_number(T aperture_length, T distance_from_aperture_to_monitor, T lambda){
    return (aperture_length * aperture_length) / (distance_from_aperture_to_monitor * lambda);
} 
template<class T> inline size_t which_diffraction_theory_work (T aperture_length, T distance_from_aperture_to_monitor, T lambda)
{
    T fresnel_n = fresnel_number(aperture_length, distance_from_aperture_to_monitor, lambda);
    size_t i = 0;
    vec3<T> limits{T(0.1), T(10), T(1e6)};
    for(; i < 3; i++){
        if(limits.at(i) > fresnel_n) break;
    }
    // 0 : [0, 0.1)   fraunhofer diffraction 
    // 1 : [0.1, 10)  fresnel diffraction
    // 2 : [10, +inf) rayleigh-sommerfield diffraction
    return i;
}
template<class T> inline T max_dz_of_angular_spectrum(T lambda, T spectrum_step)
{
    //== 借助 layleigh-length, 此处是我的猜想
    return T(1_PI) * spectrum_step * spectrum_step / lambda;
}
//== TODO : fresnel diffraction  
//    https://en.wikipedia.org/wiki/Fresnel_diffraction

//== TODO : fraunhofer diffraction 
//    https://en.wikipedia.org/wiki/Fraunhofer_diffraction