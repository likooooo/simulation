#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <assert.h>
#include <optical/pupil/film_stack_solver.hpp>
#include "fresnel_equations.hpp"

template<class T, char polarization> struct transfer_matrix_method
{
    static_assert(polarization == 's' || polarization == 'p', "light polarization should be 's' or 'p'");
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using essential_matrix = std::vector<matrix2x2<cT>>;
    using meterial = typename film_stack_solver<T>::meterial;
    using fresnel = fresnel_equations<T, polarization>;
    constexpr static std::vector<cT> propagate_direction(const std::vector<meterial>& nk, cT crao)
    {
        std::vector<cT> dir(nk.size());
        for(size_t i = 0; i < nk.size(); i++)
            dir.at(i) = snells_law<T>::refraction_angle(nk.front().nk, nk.at(i).nk, crao); 
        return dir;        
    }
    constexpr static essential_matrix interface_transfer_matrix(const std::vector<meterial>& nk, const std::vector<cT>& propagate_dir)
    {
        const size_t num_layers = nk.size();
        essential_matrix m(num_layers); 
        for(size_t i = 0; i < num_layers - 1; i++)
        {
            cT t = fresnel::transmission_coef(
                nk.at(i).nk,     propagate_dir.at(i),
                nk.at(i + 1).nk, propagate_dir.at(i + 1)
            );
            cT r = fresnel::refraction_coef(
                nk.at(i).nk,     propagate_dir.at(i),
                nk.at(i + 1).nk, propagate_dir.at(i + 1)
            );
            matrix2x2<cT>& continuous_condition = m.at(i);
            continuous_condition = {
                cT(1)/t, r/t,
                r/t,     cT(1)/t
            };
        }
        return m;
    }
    constexpr static std::vector<cT> phase_terms(const std::vector<meterial>& nk, const std::vector<cT>& propagate_dir, rT lambda)
    {
        std::vector<cT> phase(nk.size());
        for(size_t i = 0; i < nk.size(); i++)
            phase.at(i) = cT(2_PI_I) * (nk.at(i).depth / lambda) * std::cos(propagate_dir.at(i));
        return phase;
    }

    constexpr static vec2<rT> get_r_t_power_from_tmm(const matrix2x2<cT>& m, cT n1, cT theta1, cT n2, cT theta2)
    {
        return {fresnel::reflected_power(m[1][0] / m[0][0]), fresnel::transmitted_power(cT(1)/m[0][0], n1, theta1, n2, theta2)}; 
    }

};
using TMM_TEf = transfer_matrix_method< float, 's'>;
using TMM_TMf = transfer_matrix_method< float, 'p'>;
using TMM_TEd = transfer_matrix_method<double, 's'>;
using TMM_TMd = transfer_matrix_method<double, 'p'>;