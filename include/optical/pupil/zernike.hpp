#pragma once
#include <type_traist_notebook/type_traist.hpp>

template<class T, size_t M, size_t N> struct zernike_radial_table
{
    static_assert(M % 2 == 0);
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using radial = vec<rT, N>;
    constexpr static vec<size_t, M + 1> L = init_L();
    std::array<std::vector<radial>,M + 1> zk_matrix;
    zernike_radial_table(){
        for(size_t i = 0; i <= M; i++) zk_matrix.at(i).resize(L.at(i) + 1);
        init_L0();
        init_L1();
        init_rest();
        norm_all();
    }

    radial solveM0(std::vector<cT>& pupil_radial, const std::vector<std::tuple<size_t, size_t, T>>& poly_coefs) const
    {
        radial weights{0};
        for(const auto [m, l, coef] : poly_coefs){
            if(0 != m || 0 == coef) continue;
            if(l >= zk_matrix.at(m).size()){
                auto type_str = TypeReflection<decltype(*this)>;
                auto L_str = to_string(L);
                debug_print<zernike_radial_table<T>>("%s's L=%s\n", type_str.c_str(), L_str.c_str());
                error_print<zernike_radial_table<T>>(poly_coefs, {"* input args", "m", "l", "coefficients"});
                continue;
            }
            weights += (zk_matrix.at(m).at(l) * coef);
        }
        if(pupil_radial.size() == weights.size()){
            for(size_t i = 0; i < weights.size(); i++){
                pupil_radial.at(i) *= std::exp(cT(0, 1) * weights.at(i));
            }
        }
        return weights
    }

private:
    static constexpr vec<size_t, N + 1> init_L()
    {
        vec<size_t, N + 1> L;
        for(size_t i = 0; i <= M; i++){
            L.at(i) = (M - i) / 2;
        }
        return L;
    }
    void init_L0()
    {
        constexpr size_t L = 0;
        for(size_t i = 0; i <= M; i++){
            radial& v = zk_matrix.at(i).at(L);
            std::iota(v.begin(), v.end(), rT(0));
            v /= rT(N - 1);
        }
    }    
    void init_L1()
    {
        for(size_t i = 0; i <= M - 2; i++){
            radial& m0 = zk_matrix.at(i).at(0);
            radial& p0 = zk_matrix.at(i + 2).at(0);
            radial& m1 = zk_matrix.at(i).at(1);
            m1 = p0 * (m + 2) - m0 * (m + 1);
        }
    }
    void init_rest()
    {
        radial r;
        std::iota(r.begin(), r.end(), rT(0)); r /= rT(N - 1);
        for(size_t m = 0; m <= M; m++)
        for(size_t l = 2; l <=L.at(m); l++){
            size_t q = m;
            size_t p = 2 * l + m;

            rT rK1 = rT(2.0) / ((p + q) * (p - q) * (p - 2));
            rT k2 = 2 * p * (p - 1) * (p - 2);
            rT k3 = -q * q * (p - 1) - p * (p - 1) * (p - 2);
            rT k4 = rT(-0.5) * (p * (p + q - 2) * (p - q - 2));

            const radial& m2 = zk_matrix.at(m).at(l - 2);
            const radial& m1 = zk_matrix.at(m).at(l - 1);
            radial& n = zk_matrix.at(m).at(l);
            
            n  = ((r * r * k2 + k3) * m1 + m2 * k4) * rK1; 
        }
    }
    void norm_all()
    {
        for(size_t m = 0; m <= M; m++)
        for(size_t l = 0; l <= L.at(m); l++){
            size_t polydeg = 2 * (m + 2 * l + 1);
            zk_matrix.at(m).at(l) *= std::sqrt<rT>(polydeg);
        }
    }
};