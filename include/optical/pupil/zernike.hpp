#pragma once
#include <type_traist_notebook/type_traist.hpp>

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

    std::array<radial, table_count> zk_radials;
    zernike_radial_table(size_t points) {
        for(auto& v : zk_radials) v.resize(points);
        init_L0();
        init_L1();
        init_rest();
        norm_all();
    }

    // radial solveM0(std::vector<cT>& pupil_radial, const std::vector<std::tuple<size_t, size_t, T>>& poly_coefs) const
    // {
    //     radial weights{0};
    //     for(const auto [m, l, coef] : poly_coefs){
    //         if(0 != m || 0 == coef) continue;
    //         if(l >= zk_matrix.at(m).size()){
    //             auto type_str = TypeReflection<decltype(*this)>;
    //             auto L_str = to_string(L);
    //             debug("%s's L=%s\n", type_str.c_str(), L_str.c_str());
    //             error(poly_coefs, {"* input args", "m", "l", "coefficients"});
    //             continue;
    //         }
    //         weights += (zk_matrix.at(m).at(l) * coef);
    //     }
    //     if(pupil_radial.size() == weights.size()){
    //         for(size_t i = 0; i < weights.size(); i++){
    //             pupil_radial.at(i) *= std::exp(cT(0, 1) * weights.at(i));
    //         }
    //     }
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
    void init_L1()
    {
        for(size_t m = 0; m <= M - 2; m++){
            radial& m0 = zk_radials.at(index(m, 0));
            radial& p0 = zk_radials.at(index(m + 2, 0));
            radial& m1 = zk_radials.at(index(m, 1));
            m1 = p0 * (m + 2) - m0 * (m + 1);
        }
    }
    void init_rest()
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
    void norm_all()
    {
        for(size_t m = 0; m <= M; m++)
        for(size_t l = 0; l <= L.at(m); l++){
            rT polydeg = 2 * (m + 2 * l + 1);
            zk_radials.at(index(m, l)) *=  std::sqrt(polydeg); 
        }
    }
};