#pragma once
#include "film_stack_solver.hpp"
#include "jones_pupil.hpp"
#include "zernike.hpp"
#include <py_helper.hpp>
#include <optical/optical_numerics.hpp>

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
                //==
                // kr/n, 0, std::sqrt(n * n - kr * kr) / n
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

template<class T> inline T get_pupil_na(T NA)
{
    //== 高 NA 系统, 需要考虑 richards-wolf 矢量衍射理论
    //   NA < 0.4 时，s/p偏振差异 <1% (信息来源 deepseek)
    return T(0.4) <= NA ? T(1) : NA;
}

template<class T> inline void apply_anamorphic_effect(std::vector<matrix2x3<complex_t<T>>>& pupil, vec2<size_t> shape, vec2<T> step,
    T crao, T azimuth,
    T delta_z_mask, T delta_z_imaging,
    T NA, T lambda, complex_t<T> nkIn, complex_t<T> nkOut, T reduction_ratio_x = 1, T reduction_ratio_y = 1)
{
    using cT = complex_t<T>;
    const vec2<T> shift{std::sin(crao) * std::sin(azimuth), std::sin(crao) * std::cos(azimuth)};
    const T pupil_na = get_pupil_na(NA);
    // TODO : pupil shift
    error_unclassified::out("TODO : pupil shift in apply_anamorphic_effect");
    //== apply object-defocus & shift & polarization to pupil
    matrix2x3<cT>* p = pupil.data();
    kernels::center_zero_loop_square_r<T, 2>(shape, step, 
        [&](const vec2<T> fyx, T kr_2){
            const vec2<T> fr = (fyx + shift); 
            T r = std::sqrt(vector_norm(fr));
            if(r <= pupil_na){
                //== 1. apply object-defocus to pupil
                if(0 != delta_z_mask){
                    cT phase_object_defocus = cT(2_PI_I) * delta_z_mask/ lambda * std::sqrt(nkOut * nkOut - r * r);
                    (*p) *= std::exp(phase_object_defocus);
                } 
                T fzr = std::sqrt(T(1) - vector_norm(fr));
                matrix2x3<cT> Hc{
                    fzr, 0, -fr[1],
                    0, fzr, -fr[0],
                };
                assert(is_almost_equal<T>(std::abs((Hc[0] * vec3<cT>{fr[1], fr[0], fzr}) | cT(1)), 0));
                assert(is_almost_equal<T>(std::abs((Hc[1] * vec3<cT>{fr[1], fr[0], fzr}) | cT(1)), 0));
                Hc[0] /= std::sqrt(vector_norm(Hc[0]));
                Hc[1] /= std::sqrt(vector_norm(Hc[1]));
                //== 2.1. rotate xyz to sp
                T theta = r < 1e-6 ? 0 : std::atan2(fr[0], fr[1]);
                using std::sin, std::cos;
                const matrix3x3<cT> rotate_to_sp{
                    sin(theta)      , -cos(theta)     , 0,
                    cos(theta) * fzr, sin(theta) * fzr, -r,
                    fr[1]           , fr[0]           , fzr
                };
                assert(1e-6 > std::abs((((matrix3x3<cT>{
                        fr[0] / r    , -fr[1]/r     ,  0,
                        fr[1] * fzr/r, fr[0] * fzr/r, -r, 
                        fr[1]        , fr[0]        , fzr
                    } - rotate_to_sp) | vec3<cT>{1, 1, 1}) | cT(1)))
                );
                const matrix2x3<cT> Tc{
                    rotate_to_sp | Hc[0],
                    rotate_to_sp | Hc[1]
                };
                //== 2.2. projection polarization to sp cood
                const matrix3x3<cT> Mr{
                    p->at(0)[0], p->at(1)[0], r,
                    p->at(0)[1], p->at(1)[1], 0,
                    p->at(0)[2], p->at(1)[2], fzr
                };
                Hc[0] = Mr | Tc[0];
                Hc[1] = Mr | Tc[1];
                //== 2.3. rotate sp to xyz
                const matrix3x3<cT> rorate_to_xyz{
                    cos(theta) , sin(theta), 0,
                    -sin(theta), cos(theta), 0,
                    0          , 0         , 1
                };
                p->at(0) = rorate_to_xyz | Hc[0];
                p->at(1) = rorate_to_xyz | Hc[1];
            }
            p++;
        }
    );

    //== apply obliquity factor & imaging-defocus to pupil
    std::vector<matrix2x3<cT>> pupil_final(pupil.size());
    p = pupil_final.data();
    kernels::center_zero_loop_square_r<T, 2>(shape, step, 
        [&](const vec2<T> fyx, T kr_2){
            T fr = std::sqrt(kr_2);
            auto [fyR, fxR] = fyx / vec2<T>{reduction_ratio_y, reduction_ratio_x};
            T Or = std::hypot(fxR, fyR);
            cT obliquityFactor = 0;
            if (fr < nkIn.real() && Or <= nkOut.real()) {
                T cosThetaO = std::sqrt((nkOut.real()) * (nkOut.real()) - Or*Or) / (nkOut.real());
                T cosThetaI = std::sqrt((nkIn.real()) * (nkIn.real()) - fr*fr) / (nkIn.real());
                obliquityFactor = std::sqrt(cosThetaO / cosThetaI);
                vec2<T> pos = convert_to<vec2<T>>(shape) * T(0.5) + (fyx / step / vec<T, 2>{reduction_ratio_y, reduction_ratio_x});
                *p = cubic_interpolate<T>:: template eval<2>(pos, pupil, shape);
            }
            if(0 != delta_z_imaging){
                cT phase_imaging_defocus = cT(2_PI_I) * delta_z_imaging / lambda * std::sqrt(nkIn * nkIn - fr * fr);
                obliquityFactor *= std::exp(phase_imaging_defocus);
            } 
           *p *= obliquityFactor;
           p++;
        }
    );
    pupil.swap(pupil_final);
}

template<class rT> inline void pupil_golden_check(const std::vector<matrix2x3<complex_t<rT>>>& pupil_image, vec2<size_t> shape, const vec3<std::string>& golden_path = {})
{
    using cT = complex_t<rT>;
    auto[TE_x, TE_y, TE_z, TM_x, TM_y, TM_z] = decompose_from<matrix2x3<cT>, 
        cT, cT, cT, 
        cT, cT, cT>(pupil_image);
    {
        if(!golden_path[0].empty()) TE_y -= std::get<0>(load_image<std::complex<float>>(golden_path[0]));
        const auto&[real, imag] = decompose_from<cT, rT, rT>(TE_y);
        imshow(real, convert_to<std::vector<size_t>>(shape));
    }
    {
        if(!golden_path[1].empty()) TM_x -= std::get<0>(load_image<std::complex<float>>(golden_path[1]));
        const auto&[real, imag] = decompose_from<cT, rT, rT>(TM_x);
        imshow(real, convert_to<std::vector<size_t>>(shape));
    }
    {
        if(!golden_path[2].empty()) TM_z -= std::get<0>(load_image<std::complex<float>>(golden_path[2]));
        const auto&[real, imag] = decompose_from<cT, rT, rT>(TM_z);
        imshow(real, convert_to<std::vector<size_t>>(shape));
    }
};
template<class T> struct pupil_base
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    //== [-2NA, 2NA]
    vec2<size_t> shape;
    vec2<rT> step;
    rT actual_na;
    rT pupil_na;

    pupil_base(rT actual_na, rT freq_step) 
        : shape({}), step({freq_step, freq_step}), actual_na(actual_na), pupil_na(get_pupil_na(actual_na))
    {
        const size_t N = 2 * (std::ceil(2*actual_na/freq_step) + 1);
        shape = {N, N};
    }
};

template<class T> struct scalar_pupil : pupil_base<T>
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using zk_table = zernike_radial_table<T, 10>;
    std::vector<cT> pupil_image;
    scalar_pupil(rT actual_na, rT freq_step, 
        const std::vector<std::tuple<size_t, size_t, cT>>& poly_coefs = {std::tuple<size_t, size_t, cT>(0, 0, cT(1))}) 
        : pupil_base<T>(actual_na, freq_step)
    {
        pupil_image = zk_table(this->shape[0]).gen_pupil_image_with_zernike(this->shape, this->step, this->pupil_na, poly_coefs);
    }
};
template<class T> struct vector_pupil : public pupil_base<T>
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using zk_table = zernike_radial_table<T, 10>;
    std::vector<matrix2x3<cT>> pupil_image;
    vector_pupil(rT actual_na, rT freq_step, 
        const std::vector<std::tuple<size_t, size_t, cT>>& poly_coefs = {std::tuple<size_t, size_t, cT>(0, 0, cT(1))}) 
        : pupil_base<T>(actual_na, freq_step)
    {
        const std::vector<T> zernike_image = zk_table(this->shape[0]).gen_pupil_image_with_zernike(
            this->shape, this->step, this->pupil_na, poly_coefs);
        const T* pZernike = zernike_image.data();
        // imshow(zernike_image, convert_to<std::vector<size_t>>(shape));
        pupil_image.resize(zernike_image.size());
        matrix2x3<cT>* p = pupil_image.data();
    
        std::vector<matrix2x3<cT>> pupil_radials = pupil_radial<T>::init_anamorphic_pupil_radial(this->shape[0], this->pupil_na, cT(1));
        kernels::center_zero_loop_square_r<rT, 2>(this->shape, this->step, [&](const vec2<rT>& fxy, rT r){
            r = std::sqrt(r);
            matrix2x3<cT> val{0};
            if(r <= this->pupil_na) {
                r /= this->pupil_na;
                val = cubic_interpolate<rT>::eval(r * pupil_radials.size(), pupil_radials);
                val *= std::exp(cT(0, *pZernike));
            } 
            *p = val;
            pZernike++;
            p++;
        });
    }
};







