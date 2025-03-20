#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <assert.h>
#include <snells_law.hpp>
#include <pupil/film_stack_solver.hpp>

template<class T> struct transfer_matrix_method_solver 
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using essential_matrix = std::vector<dual_vec<matrix2x2<cT>>>;
    using inv_essential_matrix = essential_matrix;
    using phase_matrix = std::vector<matrix2x2<cT>> ;
    using meterial = film_stack_solver<T>::meterial;

    constexpr static void eval(const std::vector<meterial>& nk)
    {
        const size_t num_layers = nk.size();
        std::vector<cT> t_list(num_layers * num_layers);
        std::vector<cT> r_list(num_layers * num_layers);
    }
    constexpr static std::vector<cT> init_propagate_direction(const std::vector<cT>& nk, cT crao){
        return snells_law<T>::refraction_angle(nk, crao);
    }
    constexpr static essential_matrix init_essential_matrix(
        const std::vector<cT>& nk, const std::vector<cT>& propagate_direction)
    {
        assert(nk.size() == propagate_direction.size());

        essential_matrix m(nk.size());
        for(size_t i = 0; i < nk.size(); i++)
        {
            m.at(i) = {
                //== essential matrix TE
                1,                                              1,
                nk.at(i) * std::cos(propagate_direction.at(i)), -nk.at(i) * std::cos(propagate_direction.at(i)),
                //== essential matrix TM
                std::cos(propagate_direction.at(i)), std::cos(propagate_direction.at(i))
                nk.at(i)                           , -nk.at(i)
            }
        }
    }
    constexpr static inv_essential_matrix init_inv_essential_matrix(essential_matrix m)
    {
        
        for(auto& [TE, TM] : m)
        {

        }
    }

    constexpr static phase_matrix init_phase_matrix(
        const std::vector<cT>& nk, const std::vector<rT>& thickness, const std::vector<cT>& propagate_direction, rT wave_length
    )
    {
        assert(thickness.size() == propagate_direction.size());
        phase_matrix m(thickness.size());
        for(size_t i = 0; i < m.size(); i++)
        {
            m.at(i) = {
                std::exp(cT(-2_PI_I) * nk.at(i) * std::cos(propagate_direction.at(i)) * (thickness.size() / wave_length)), 0
                0, std::exp(cT(2_PI_I) * nk.at(i) * std::cos(propagate_direction.at(i)) * (thickness.size() / wave_length))
            };
        }
    }
};
using TMM_solverf = transfer_matrix_method_solver<float>;
using TMM_solverd = transfer_matrix_method_solver<double>;