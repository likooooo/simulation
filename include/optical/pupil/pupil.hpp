#pragma once
#include "film_stack_solver.hpp"
#include "jones_pupil.hpp"
#include "zernike.hpp"

template<class T> struct pupil_radial
{
private:
    static std::vector<T> make_kr(T NA, int numPoints)
    {
        std::vector<T> vec(numPoints);
        std::iota(vec.begin(), vec.end(), 0);
        T step = NA / (numPoints - 1); 
        vec *= step;
        return vec; // [0, NA]
    }
public:
    using cT = complex_t<T>;
    using rT = real_t<T>;
    constexpr static cT vacuum_nk{T(1)};
    static std::vector<matrix2x3<cT>> init_anamorphic_pupil_radial(int numPoints, T NA, cT nI)
    {
        std::vector<matrix2x3<cT>> pupil(numPoints);
        std::vector<T> kr = make_kr(NA, numPoints);
        const auto [n, k] = nI;
        std::transform(kr.begin(), kr.end(), pupil.begin(), [&](auto kr){
            return matrix2x3<cT>{
                //== TE
                0, 1, 0, 
                //== TM
                std::sqrt(n * n - kr * kr) / n , 0, -kr/n 
            };
        });
        return pupil;
    }
    static void apply_defocus_to_pupil_radial(std::vector<matrix2x3<cT>>& pupil, cT nI, T deltaZ, T NA, T wavelength, T reduction_ratio = T(1)){
       
        std::vector<T> kr = make_kr(NA, pupil.size());
        for(size_t i = 0; i < kr.size(); i++){
            cT cPhase(0, deltaZ / wavelength);
            T k = kr.at(i) * reduction_ratio;
            cPhase = cPhase * kz<cT>(nI, k) / nI;
            pupil.at(i) *= std::exp(cPhase); 
        }
    }
    static void apply_obliquity_factor_to_pupil_radial(std::vector<matrix2x3<cT>>& pupil, T NA, T reduction_ratio, cT index){
        std::vector<T> kr = make_kr(NA, pupil.size());
        for(size_t i = 0; i < kr.size(); i++)
        {
            if(kr.at(i)  > vacuum_nk.real() || kr.at(i) > index.real())
            {
                pupil.at(i) = {0};
            }
            else
            {
                auto cosO = kz(vacuum_nk.real(), kr.at(i)) / vacuum_nk.real();
                auto cosI = kz(index.real(), kr.at(i)) / index.real();
                pupil.at(i) *= std::sqrt(cosO / cosI);
            }
        }
    }
    static void apply_film_stack_to_pupil_radial(std::vector<matrix2x3<cT>>& pupil, 
        T NA, T lambda, cT index, T depth, 
        const std::vector<typename film_stack_solver<T>::meterial>& meterials)
    {
        film_stack_solver<T> solver(lambda, 1.0);
        const std::vector<T> kr = make_kr(NA, pupil.size());
        for(size_t i = 0; i < kr.size(); i++)
        {
            matrix2x3<cT>& p = pupil.at(i);
            T k = kr.at(i);
            auto E_transfer_matrix = solver.solveTE(meterials, k);
            auto B_transfer_matrix = solver.solveTM(meterials, k);
            vec3<cT> E_at_depth = solver.solve_at_depth(meterials, E_transfer_matrix, k, depth);
            vec3<cT> B_at_depth = solver.solve_at_depth(meterials, B_transfer_matrix, k, depth);
            p.at(0) = solver.solve_field(p.at(0), E_transfer_matrix.at(0).at(0), E_at_depth, B_at_depth, B_at_depth);
            p.at(1) = solver.solve_field(p.at(1), E_transfer_matrix.at(0).at(0), E_at_depth, B_at_depth, B_at_depth);
        }
    }
};

// freq-shift
//