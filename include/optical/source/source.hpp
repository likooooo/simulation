#pragma once
#include "parametric_source.hpp"
#include <kernels/kernel_loop.hpp>

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