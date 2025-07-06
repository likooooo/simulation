#pragma once
#include <type_traist_notebook/type_traist.hpp>

enum class source_type{
    circular = 0, annular, dipole_fan,quasar, leaf2, leaf4
};
template<class T> struct parametric_source
{
    static_assert(std::is_floating_point_v<T>);

    struct traditional_source_params{T sigma{1}, centerX{0}, centerY{0};};
    constexpr static T source_boundary_epsion =1e-5;

    static bool is_between_angle(const T x, const T y, const T angle1, const T angle2){ 
        T dr = std::hypot(x, y);
        T cos_theta1 = (x * std::cos(angle1) + y * std::sin(angle1)) / dr;
        T cos_theta2 = (x * std::cos(angle2) + y * std::sin(angle2)) / dr;
        T theta1 = std::acos(cos_theta1); // [0, pi]
        T theta2 = std::acos(cos_theta2); // [0, pi]
        
        T delta_angle = std::abs(angle1 - angle2);
        if(is_almost_equal(std::max(theta1, theta2), delta_angle, source_boundary_epsion)) return true;
        else if(std::max(theta1, theta2) > delta_angle) return false;
        return true;
    }

    template<class TFunc> static void foreach_source_pixel(T* pSource, size_t xsize, size_t ysize, TFunc shape_maker)
    {
        const T srcSigmaStepX = 2.0 / (xsize - 1);
        const T srcSigmaStepY = 2.0 / (ysize - 1);
        
        for (int iy = 0; iy < ysize; iy++){
            int indY = iy - ysize/2;
            T sigmaDy = indY * srcSigmaStepY;
            for (int ix = 0; ix < xsize; ix++, pSource++){
                int indX = ix - xsize/2;
                T sigmaDx = indX * srcSigmaStepX;
                shape_maker(*pSource, sigmaDx, sigmaDy);
            }
        }
    }

    static void get_traditional_source(T* pSource, size_t xsize, size_t ysize, const traditional_source_params& params){ 
        const auto [sigma, centerX, centerY] = params;
        auto shape_maker =  [&](T& source, T sigmaDx, T sigmaDy){
            T dr = std::hypot(sigmaDx - centerX, sigmaDy - centerY); 
            if (dr <= sigma || is_almost_equal(dr, sigma, source_boundary_epsion)) source += 1.0;
        };
        foreach_source_pixel(pSource, xsize, ysize, shape_maker);
    }

    struct annular_source_params
    {
        T sigmaOut{1};
        T sigmaIn{0.5};
        T sigmaInnerShiftX{0};
        T sigmaInnerShiftY{0};
        T sigmaShiftX{0};
        T sigmaShiftY{0};
    };
    static void get_annular_source(T* pSource, size_t xsize, size_t ysize, const annular_source_params& params){
        const auto [sigmaOut, sigmaIn, sigmaInnerShiftX,  sigmaInnerShiftY, sigmaShiftX, sigmaShiftY]= params; 
        const auto innerCircCenterX = sigmaShiftX + sigmaInnerShiftX;
        const auto innerCircCenterY = sigmaShiftY + sigmaInnerShiftY;
        auto shape_maker =  [&](T& source, T sigmaDx, T sigmaDy){
            T drOut = std::hypot(sigmaDx - sigmaShiftX, sigmaDy - sigmaShiftY); 
            T drIn = std::hypot(sigmaDx - innerCircCenterX, sigmaDy - innerCircCenterY); 
            if(((drOut < sigmaOut) && (drIn > sigmaIn))|| is_almost_equal(drOut, sigmaOut) || is_almost_equal(drIn, sigmaIn))
                source += 1.0;
        };
        foreach_source_pixel(pSource, xsize, ysize, shape_maker);
    }
    struct dipole_fan_source_params
    {
        annular_source_params args;
        T rotAngle{0}, spanAngle{M_PI_4f};
    };

    static void get_dipole_fan_source(T* pSource, size_t xsize, size_t ysize, const dipole_fan_source_params& params){
        const auto [annular_params, rotAngle, spanAngle] = params;
        get_annular_source(pSource, xsize, ysize, annular_params);
        T fan1_Angle1 = rotAngle + spanAngle * 0.5;
        T fan1_Angle2 = rotAngle - spanAngle * 0.5;
        T fan2_Angle1 = fan1_Angle1 + M_PI;
        T fan2_Angle2 = fan1_Angle2 + M_PI;
        auto shape_maker = [&](T& source, T sigmaDx, T sigmaDy){
            if(is_almost_equal(sigmaDx, 1e-4f) && is_almost_equal(sigmaDy, 1e-4f)) return;
            if(!is_between_angle(sigmaDx, sigmaDy, fan1_Angle1, fan1_Angle2) && 
            !is_between_angle(sigmaDx, sigmaDy, fan2_Angle1, fan2_Angle2)) 
                source = 0.0;
        };

        foreach_source_pixel(pSource, xsize, ysize, shape_maker);
    }
    using quadratic_fan_source_params = dipole_fan_source_params;
    static void get_quadratic_fan_source(T* pSource, size_t xsize, size_t ysize, const quadratic_fan_source_params& params)
    {
        const auto [annular_params, rotAngle, spanAngle] = params;
        get_annular_source(pSource, xsize, ysize, annular_params);
        const T fan1_Angle1 = rotAngle + spanAngle * 0.5;
        const T fan1_Angle2 = rotAngle - spanAngle * 0.5;
        const T fan2_Angle1 = fan1_Angle1 + M_PI * 0.5;
        const T fan2_Angle2 = fan1_Angle2 + M_PI * 0.5;
        const T fan3_Angle1 = fan1_Angle1 + M_PI;
        const T fan3_Angle2 = fan1_Angle2 + M_PI;
        const T fan4_Angle1 = fan1_Angle1 + M_PI * 1.5;
        const T fan4_Angle2 = fan1_Angle2 + M_PI * 1.5;
        auto shape_maker = [&](T& source, T sigmaDx, T sigmaDy){
            sigmaDx -= annular_params.sigmaShiftX;
            sigmaDy -= annular_params.sigmaShiftY;
            if(is_almost_equal(sigmaDx, source_boundary_epsion) && is_almost_equal(sigmaDy, source_boundary_epsion))
                return;
            if(!is_between_angle(sigmaDx, sigmaDy, fan1_Angle1, fan1_Angle2) &&
            !is_between_angle(sigmaDx, sigmaDy, fan2_Angle1, fan2_Angle2) &&
            !is_between_angle(sigmaDx, sigmaDy, fan3_Angle1, fan3_Angle2) &&
            !is_between_angle(sigmaDx, sigmaDy, fan4_Angle1, fan4_Angle2))
                source = 0;
        };
        foreach_source_pixel(pSource, xsize, ysize, shape_maker);
    }

    struct dipole_leaf_source_params{T sigma_D{0.2}, sigma_d{0.3}, rotAngle{0}, sigmaShiftX{0}, sigmaShiftY{0};};
    static void cut_leaf(std::vector<T>& leaf, T target){
        for(auto& n :leaf) {
            if(!is_almost_equal(n, target, source_boundary_epsion)) n = 0;
        }
    };

    static void get_dipole_leaf_source(T* pSource, size_t xsize, size_t ysize, const dipole_leaf_source_params& params)
    {
        auto [sigma_D, sigma_d, rotAngle, sigmaShiftX, sigmaShiftY] = params;
        sigma_D *= 2; // distance between the center of two leaves
        const T leftLeaf_leftCirc_dis = 2.0 - (sigma_d - 0.5 * (sigma_D + sigma_d - 2)); // distance from circle center to (0,0) before rotate and shift, longer
        const T leftLeaf_rightCirc_dis = 0.5 * (sigma_D + sigma_d - 2); // distance from circle center to (0,0) before rotate and shift, shoter

        std::vector<T> leftLeaf(xsize* ysize), rightLeaf(xsize * ysize);
        auto cal_half_leaf = [&](std::vector<T>& leaf, T leftLeaf_leftCirc_dis){
            const T leftLeaf_leftCirc_centerX = -1.0 * leftLeaf_leftCirc_dis * std::cos(rotAngle) + sigmaShiftX;
            const T leftLeaf_leftCirc_centerY = -1.0 * leftLeaf_leftCirc_dis * std::sin(rotAngle) + sigmaShiftY;
            get_traditional_source(leaf.data(), xsize, ysize, traditional_source_params{1.0, leftLeaf_leftCirc_centerX, leftLeaf_leftCirc_centerY});
        };
        cal_half_leaf(leftLeaf, leftLeaf_leftCirc_dis);
        cal_half_leaf(leftLeaf, leftLeaf_rightCirc_dis);
        cal_half_leaf(rightLeaf, -leftLeaf_leftCirc_dis);
        cal_half_leaf(rightLeaf, -leftLeaf_rightCirc_dis);

        cut_leaf(leftLeaf, 2.0);
        cut_leaf(rightLeaf, 2.0);
        std::transform(leftLeaf.begin(), leftLeaf.end(), rightLeaf.begin(), pSource, [](auto a, auto b){
            return 1.0f * static_cast<int>(
                is_almost_equal(a, T(2.0), source_boundary_epsion) ||
                is_almost_equal(b, T(2.0), source_boundary_epsion)
            );
        });
    }
    static void get_leaf_by_four_circles(std::vector<T>& leaf,size_t xsize, size_t ysize, const std::array<std::complex<T>, 4>& CircCenterPt){
        get_traditional_source(leaf.data(), xsize, ysize, {1.0, CircCenterPt.at(0).real(), CircCenterPt.at(0).imag()});
        get_traditional_source(leaf.data(), xsize, ysize, {1.0, CircCenterPt.at(1).real(), CircCenterPt.at(1).imag()});
        get_traditional_source(leaf.data(), xsize, ysize, {1.0, CircCenterPt.at(2).real(), CircCenterPt.at(2).imag()});
        get_traditional_source(leaf.data(), xsize, ysize, {1.0, CircCenterPt.at(3).real(), CircCenterPt.at(3).imag()});
        cut_leaf(leaf, 4.0);
    }

    struct quadratic_leaf_source_params{
        dipole_leaf_source_params params; 
        T aspectRatio{T(0.2)};
    };
    static void get_quadratic_leaf_source(T* pSource, size_t xsize, size_t ysize, const quadratic_leaf_source_params& params)
    {
        const auto [dipole_params, aspectRatio] = params;
        const auto [sigma_D, sigma_d, rotAngle, sigmaShiftX, sigmaShiftY] = dipole_params;
        const T dis1 = 1.0 - 0.5 * sigma_d; // distance from center of right circle to (0,0)
        const T dis2 = 1.0 - 0.5 * sigma_d * aspectRatio; // distance from center of up circle to (0,0)
        const int n = xsize* ysize;
        const T inverse_rotAngle = -1.0 * rotAngle;

        std::vector<T> leftLeaf(n), rightLeaf(n), upLeaf(n), downLeaf(n);

        auto cal_base_leaf_args = [&](T d1, T d2){
            std::complex<T> left(-1.0 * d1 * std::cos(inverse_rotAngle), d1 * std::sin(inverse_rotAngle));
            std::complex<T> top(d2 * std::sin(inverse_rotAngle), d2 * std::cos(inverse_rotAngle));
            return std::array<std::complex<T>, 4>{left, left * T(-1.0), top, top * T(-1.0)};
        };
        auto cal_leaf = [&](const std::array<std::complex<T>, 4>& base, const T a, const T b){
            std::array<std::complex<T>, 4> result = base;
            const std::complex<T> offset(
                a * sigma_D * std::cos(rotAngle) + sigmaShiftX, 
                b * sigma_D * std::sin(rotAngle) + sigmaShiftY
            );
            for(auto& b : result) b += offset;
            return result;
        };
        std::array<std::complex<T>, 4> base, leaf_args;
        base = cal_base_leaf_args(dis1, dis2);
        leaf_args = cal_leaf(base, -1.0, -1.0);
        get_leaf_by_four_circles(leftLeaf, xsize, ysize, leaf_args);
        leaf_args = cal_leaf(base, 1.0, 1.0);
        get_leaf_by_four_circles(rightLeaf, xsize, ysize, leaf_args);
        base = cal_base_leaf_args(dis2, dis1);
        leaf_args = cal_leaf(base, -1.0, 1.0);
        get_leaf_by_four_circles(upLeaf, xsize, ysize, leaf_args);
        leaf_args = cal_leaf(base, 1.0, -1.0);
        get_leaf_by_four_circles(downLeaf, xsize, ysize, leaf_args);
        for(int i =0; i < leftLeaf.size(); i++){
            if((leftLeaf.at(i) + rightLeaf.at(i) + upLeaf.at(i) + downLeaf.at(i)) > 0)pSource[i] = 1;
        }
    }

    struct shift_effect_params{
        T intenTiltX, intenTiltY, intenEllipHV, intenEllipST, shiftX, shiftY;
    };
    static void add_shift_effect_to_source(T* pSource, size_t xsize, size_t ysize, const shift_effect_params& params)
    {
        const auto [intenTiltX, intenTiltY, intenEllipHV, intenEllipST, shiftX, shiftY] = params;
        auto shape_maker =  [&](T& source, T dx, T dy){
            T theta = atan2(dy-shiftY, dx-shiftX); // [-pi pi]
            T distance = std::hypot(dx-shiftX, dy-shiftY);
            T intenShiftCoeff = (intenTiltX*cos(theta) + intenTiltY*sin(theta)) 
                + (intenEllipHV*cos(2*theta) - intenEllipST*sin(2*theta)) * distance * distance;
            source *= (1+intenShiftCoeff);
        };
        foreach_source_pixel(pSource, xsize, ysize, shape_maker);
    }
};
