#pragma once
#include "parametric_source.hpp"
#include <kernels/kernel_loop.hpp>

template<class T> T k0(T lambda, T NA)
{
    return 2_PI / lambda * NA;
}

template<class T> struct source_point
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    rT intensity{1};
    
    // [0, 1]
    rT sigmax{0}, sigmay{0};
    
    //  0 : unpolarization
    //  1 : full TE-polarization
    // -1 : full TM-polarization
    rT dop{0};
    
    //  0       : linearly polarized
    // -1       : left  circularly polarized
    //  1       : right circularly polarized
    //  (-1, 0) : left  elliptical polarized
    //  ( 0, 1) : right  elliptical polarized
    rT ellipticity{0};

    // degree in rad
    rT e_field_direction{0.5_PI};
    
    bool is_satisfy_bragg_condition() const
    {
        return int(sigmax) == sigmax && int(sigmay) == sigmay;
    }
    vec2<rT> k_vector(rT lambda, rT slow_rate = 1) const
    {
        vec2<rT> r{sigmax, sigmay};
        return k0(lambda, slow_rate) * r;   
    }

    vec2<rT> get_TEM_coef() const
    {
        rT TM = (1 - dop) * 0.5;
        rT TE = 1 - TM;
        return {TE, TM};
    }
    vec2<vec2<cT>> polarization_state() const
    {
        auto cal_polarization_state = [](rT rad, rT ellipticity)
        {
            //== https://en.wikipedia.org/wiki/Polarization_(waves)#Polarization_ellipse
            vec2<cT> polar{std::cos(e_field_direction), std::sin(e_field_direction)};
            polar.at(1) *= std::exp(cT(0, 0.5_PI * ellipticity));
            return polar;
        }
        auto [TE_coef, TM_coef] = get_TEM_coef();
        vec2<cT> TE_ps = 0 == TE_coef ? vec2<cT>() : cal_polarization_state(e_field_direction, ellipticity);
        vec2<cT> TM_ps = 0 == TM_coef ? vec2<cT>() :cal_polarization_state(e_field_direction - 0.5_PI, ellipticity);
        return {TE_ps, TM_ps};
    }

};

template<class T> struct source
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    static constexpr vec3<rT> get_Exyz(cT kxy, rT alpha, rT beta)
    {
        rT gamma = std::sqrt(rT(1) - alpha * alpha - beta * beta);
        const auto [Ex, Ey] = kxy;
        if(gamma < 1e-6) return {Ex, Ey, rT(0)};
        rT Ez = -(Ex * alpha + Ey * beta) / gamma;
        rT norm = std::abs(kxy);
        return {Ex / norm, Ey / norm, Ez / norm};
    }
    static constexpr vec2<rT> get_crao_azimuth(rT alpha, rT beta){
        return {std::asin(std::hypot(alpha, beta)), std::atan2(beta, alpha)};
    }
    static constexpr matrix3x3<rT> rotate_matrix(rT crao = 0, rT azimuth = 0)
    {
        // 1. rotate crao with Y
        // 2. rotate azimuth with Z
        using std::cos, std::sin;
        rT theta = crao, phi = azimuth;
        matrix3x3<rT> m{
            cos(theta) * cos(phi),-sin(phi), sin(theta)*cos(phi),
            cos(theta) * sin(phi), cos(phi), sin(theta)*sin(phi),
                      -sin(theta),    rT(0),          cos(theta)
        };
        return m;
    }
    static constexpr matrix3x3<rT> rotate_inv_matrix(rT crao = 0, rT azimuth = 0)
    {
        using std::cos, std::sin;
        rT theta = crao, phi = azimuth;
        matrix3x3<rT> m{
            cos(theta) * cos(phi), cos(theta)*sin(phi),-sin(theta),
                        -sin(phi),            cos(phi),      rT(0),
            sin(theta) * cos(phi), sin(theta)*sin(phi), cos(theta)
        };
        return m;
    }
    
    struct source_point
    {
        using print_type = std::tuple<vec3<rT>, vec2<rT>, rT>;
        rT intensity;
        vec2<rT> sigmaXY;
        vec3<rT> k;
    };
    static constexpr std::vector<source_point> get_source_points(
        rT* pSourceImage, vec2<size_t> shape, vec2<rT> step, rT threshold, rT crao = 0, rT azimuth = 0)
    {
        std::vector<source_point> points;
        points.reserve(shape[0] * shape[1]);
        kernels::center_zero_loop_square_r(shape, step,
             [&](vec2<rT> k, rT r_square){
                if(*pSourceImage > threshold){
                    source_point p;
                    p.sigmaXY = {k.at(1), k.at(0)};
                    p.k = {k.at(1), k.at(0), rT(1)};
                    p.k = (rotate_matrix(crao, azimuth) | p.k);
                    p.intensity = *pSourceImage; 
                    points.push_back(p);
                }
                pSourceImage++;
            }
        );
        return points;
    }
    static constexpr std::vector<source_point> get_source_points(
        rT* pSourceImage, vec2<size_t> shape, rT NA, rT threshold, rT crao = 0, rT azimuth = 0)
    {
        vec2<rT> step = vec2<rT>{NA * 2, NA * 2} / shape;
        return get_source_points(pSourceImage, shape, step, threshold, crao, azimuth);
    }
};