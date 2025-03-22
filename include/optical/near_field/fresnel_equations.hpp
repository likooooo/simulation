#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include "snells_law.hpp"

template<class T, char polarization = 's'>
struct fresnel_equations
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using snell = snells_law<T>;
    constexpr static cT refraction_coef(cT n1, cT theta1, cT n2, cT theta2)
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
    constexpr static cT refraction_coef(cT n1, cT theta1, cT n2){return refraction_coef(n1, theta1, n2, snells_law(n1, n2, theta1));}
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
};