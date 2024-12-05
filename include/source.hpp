#include <type_traist_notebook.hpp>
#include <cmath>
#include <algorithm>
#include <array>
#include <complex>

namespace std
{
    template<class T> std::complex<T> make_complex(T real, T imag){return std::complex<T>(real, imag);}
}

enum class source_type{
    circular = 0, annular, dipole_fan,quasar, leaf2, leaf4
};
using OLint =int;
using OLreal = float;
using OLbool = bool;
struct traditional_source_params
{
    OLreal sigma, centerX, centerY;
};
constexpr OLreal source_boundary_epsion =1e-5;

template<class T>bool is_almost_equal(T a, T b, T epsion = std::numeric_limits<T>::epsilon())
{
    static_assert(std::is_floating_point_v<T>);
    return epsion > std::fabs(a-b);
}
OLbool is_between_angle(const OLreal x, const OLreal y, const OLreal angle1, const OLreal angle2){ 
    OLreal dr = std::hypot(x, y);
    OLreal cos_theta1 = (x * std::cos(angle1) + y * std::sin(angle1)) / dr;
    OLreal cos_theta2 = (x * std::cos(angle2) + y * std::sin(angle2)) / dr;
    OLreal theta1 = std::acos(cos_theta1); // [0, pi]
    OLreal theta2 = std::acos(cos_theta2); // [0, pi]
    
    OLreal delta_angle = std::abs(angle1 - angle2);
    if(is_almost_equal(std::max(theta1, theta2), delta_angle, source_boundary_epsion)) return true;
    else if(std::max(theta1, theta2) > delta_angle) return false;
    return true;
}

template<class TFunc> void foreach_source_pixel(OLreal* pSource, int xsize, int ysize, TFunc shape_maker)
{
    const OLreal srcSigmaStepX = 2.0 / xsize;
    const OLreal srcSigmaStepY = 2.0 / ysize;
    for (OLint iy = 0; iy < ysize; iy++){
        OLint indY = iy - ysize/2;
        OLreal sigmaDy = indY * srcSigmaStepY;
        for (OLint ix = 0; ix < xsize; ix++, pSource++){
            OLint indX = ix - xsize/2;
            OLreal sigmaDx = indX * srcSigmaStepX;
            shape_maker(*pSource, sigmaDx, sigmaDy);
        }
    }
}

void get_traditional_source(OLreal* pSource, int xsize, int ysize, const traditional_source_params& params){ 
    const auto [sigma, centerX, centerY] = params;
    auto shape_maker =  [&](OLreal& source, OLreal sigmaDx, OLreal sigmaDy){
        OLreal dr = std::hypot(sigmaDx - centerX, sigmaDy - centerY); 
        if (dr < sigma || is_almost_equal(dr, sigma, source_boundary_epsion)) source += 1.0;
    };
    foreach_source_pixel(pSource, xsize, ysize, shape_maker);
}

struct annular_source_params
{
    OLreal sigmaOut;
    OLreal sigmaIn;
    OLreal sigmaInnerShiftX;
    OLreal sigmaInnerShiftY;
    OLreal sigmaShiftX;
    OLreal sigmaShiftY;
};
void get_annular_source(OLreal* pSource, int xsize, int ysize, const annular_source_params& params){
    const auto [sigmaOut, sigmaIn, sigmaInnerShiftX,  sigmaInnerShiftY, sigmaShiftX, sigmaShiftY]= params; 
    const auto innerCircCenterX = sigmaShiftX + sigmaInnerShiftX;
    const auto innerCircCenterY = sigmaShiftY + sigmaInnerShiftY;
    auto shape_maker =  [&](OLreal& source, OLreal sigmaDx, OLreal sigmaDy){
        OLreal drOut = std::hypot(sigmaDx - sigmaShiftX, sigmaDy - sigmaShiftY); 
        OLreal drIn = std::hypot(sigmaDx - innerCircCenterX, sigmaDy - innerCircCenterY); 
        if(((drOut < sigmaOut) && (drIn > sigmaIn))|| is_almost_equal(drOut, sigmaOut) || is_almost_equal(drIn, sigmaIn))
            source += 1.0;
    };
    foreach_source_pixel(pSource, xsize, ysize, shape_maker);
}
struct dipole_fan_source_params
{
    annular_source_params args;
    OLreal rotAngle, spanAngle;
};

void get_dipole_fan_source(OLreal* pSource, int xsize, int ysize, const dipole_fan_source_params& params){
    const auto [annular_params, rotAngle, spanAngle] = params;
    get_annular_source(pSource, xsize, ysize, annular_params);
    OLreal fan1_Angle1 = rotAngle + spanAngle * 0.5;
    OLreal fan1_Angle2 = rotAngle - spanAngle * 0.5;
    OLreal fan2_Angle1 = fan1_Angle1 + M_PI;
    OLreal fan2_Angle2 = fan1_Angle2 + M_PI;
    auto shape_maker = [&](OLreal& source, OLreal sigmaDx, OLreal sigmaDy){
        if(is_almost_equal(sigmaDx, 1e-4f) && is_almost_equal(sigmaDy, 1e-4f)) return;
        if(!is_between_angle(sigmaDx, sigmaDy, fan1_Angle1, fan1_Angle2) && 
           !is_between_angle(sigmaDx, sigmaDy, fan2_Angle1, fan2_Angle2)) 
            source = 0.0;
    };

    foreach_source_pixel(pSource, xsize, ysize, shape_maker);
}
using quadratic_fan_source_params = dipole_fan_source_params;
void get_quadratic_fan_source(OLreal* pSource, int xsize, int ysize, const quadratic_fan_source_params& params)
{
    const auto [annular_params, rotAngle, spanAngle] = params;
    get_annular_source(pSource, xsize, ysize, annular_params);
    const OLreal fan1_Angle1 = rotAngle + spanAngle * 0.5;
    const OLreal fan1_Angle2 = rotAngle - spanAngle * 0.5;
    const OLreal fan2_Angle1 = fan1_Angle1 + M_PI * 0.5;
    const OLreal fan2_Angle2 = fan1_Angle2 + M_PI * 0.5;
    const OLreal fan3_Angle1 = fan1_Angle1 + M_PI;
    const OLreal fan3_Angle2 = fan1_Angle2 + M_PI;
    const OLreal fan4_Angle1 = fan1_Angle1 + M_PI * 1.5;
    const OLreal fan4_Angle2 = fan1_Angle2 + M_PI * 1.5;
    auto shape_maker = [&](OLreal& source, OLreal sigmaDx, OLreal sigmaDy){
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

struct dipole_leaf_source_params{OLreal sigma_D, sigma_d, rotAngle, sigmaShiftX, sigmaShiftY;};
void cut_leaf(std::vector<OLreal>& leaf, OLreal target){
    for(auto& n :leaf) {
        if(!is_almost_equal(n, target, source_boundary_epsion)) n = 0;
    }
};

void get_dipole_leaf_source(OLreal* pSource, int xsize, int ysize, const dipole_leaf_source_params& params)
{
    auto [sigma_D, sigma_d, rotAngle, sigmaShiftX, sigmaShiftY] = params;
    sigma_D *= 2; // distance between the center of two leaves
    const OLreal leftLeaf_leftCirc_dis = 2.0 - (sigma_d - 0.5 * (sigma_D + sigma_d - 2)); // distance from circle center to (0,0) before rotate and shift, longer
    const OLreal leftLeaf_rightCirc_dis = 0.5 * (sigma_D + sigma_d - 2); // distance from circle center to (0,0) before rotate and shift, shoter

    std::vector<OLreal> leftLeaf(xsize* ysize), rightLeaf(xsize * ysize);
    auto cal_half_leaf = [&](std::vector<OLreal>& leaf, OLreal leftLeaf_leftCirc_dis){
        const OLreal leftLeaf_leftCirc_centerX = -1.0 * leftLeaf_leftCirc_dis * std::cos(rotAngle) + sigmaShiftX;
        const OLreal leftLeaf_leftCirc_centerY = -1.0 * leftLeaf_leftCirc_dis * std::sin(rotAngle) + sigmaShiftY;
        get_traditional_source(leaf.data(), xsize, ysize, traditional_source_params{1.0, leftLeaf_leftCirc_centerX, leftLeaf_leftCirc_centerY});
    };
    cal_half_leaf(leftLeaf, leftLeaf_leftCirc_dis);
    cal_half_leaf(leftLeaf, leftLeaf_rightCirc_dis);
    cal_half_leaf(rightLeaf, -leftLeaf_leftCirc_dis);
    cal_half_leaf(rightLeaf, -leftLeaf_rightCirc_dis);

    cut_leaf(leftLeaf, 2.0);
    cut_leaf(rightLeaf, 2.0);
    std::transform(leftLeaf.begin(), leftLeaf.end(), rightLeaf.begin(), pSource, [](auto a, auto b){
        return 1.0f * static_cast<OLint>(
            is_almost_equal(a, OLreal(2.0), source_boundary_epsion) ||
            is_almost_equal(b, OLreal(2.0), source_boundary_epsion)
        );
    });
}
void get_leaf_by_four_circles(std::vector<OLreal>& leaf,int xsize, int ysize, const std::array<std::complex<OLreal>, 4>& CircCenterPt){
    get_traditional_source(leaf.data(), xsize, ysize, {1.0, CircCenterPt.at(0).real(), CircCenterPt.at(0).imag()});
    get_traditional_source(leaf.data(), xsize, ysize, {1.0, CircCenterPt.at(1).real(), CircCenterPt.at(1).imag()});
    get_traditional_source(leaf.data(), xsize, ysize, {1.0, CircCenterPt.at(2).real(), CircCenterPt.at(2).imag()});
    get_traditional_source(leaf.data(), xsize, ysize, {1.0, CircCenterPt.at(3).real(), CircCenterPt.at(3).imag()});
    cut_leaf(leaf, 4.0);
}

struct quadratic_leaf_source_params{
    dipole_leaf_source_params params; 
    OLreal aspectRatio;
};
void get_quadratic_leaf_source(OLreal* pSource, int xsize, int ysize, const quadratic_leaf_source_params& params)
{
    const auto [dipole_params, aspectRatio] = params;
    const auto [sigma_D, sigma_d, rotAngle, sigmaShiftX, sigmaShiftY] = dipole_params;
    const OLreal dis1 = 1.0 - 0.5 * sigma_d; // distance from center of right circle to (0,0)
    const OLreal dis2 = 1.0 - 0.5 * sigma_d * aspectRatio; // distance from center of up circle to (0,0)
    const int n = xsize* ysize;
    const OLreal inverse_rotAngle = -1.0 * rotAngle;

    std::vector<OLreal> leftLeaf(n), rightLeaf(n), upLeaf(n), downLeaf(n);

    auto cal_base_leaf_args = [&](OLreal d1, OLreal d2){
        std::complex<OLreal> left(-1.0 * d1 * std::cos(inverse_rotAngle), d1 * std::sin(inverse_rotAngle));
        std::complex<OLreal> top(d2 * std::sin(inverse_rotAngle), d2 * std::cos(inverse_rotAngle));
        return std::array<std::complex<OLreal>, 4>{left, left * OLreal(-1.0), top, top * OLreal(-1.0)};
    };
    auto cal_leaf = [&](const std::array<std::complex<OLreal>, 4>& base, const OLreal a, const OLreal b){
        std::array<std::complex<OLreal>, 4> result = base;
        const std::complex<OLreal> offset(
            a * sigma_D * std::cos(rotAngle) + sigmaShiftX, 
            b * sigma_D * std::sin(rotAngle) + sigmaShiftY
        );
        for(auto& b : result) b += offset;
        return result;
    };
    std::array<std::complex<OLreal>, 4> base, leaf_args;
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

//   1973            if (srcShape == "CIRCULAR"){                                                                     │
// │   1974                GetNonSmoothedTraditionalSrc(illumSource, srcPara);                                          │
// │   1975            }else if(srcShape == "ANNULAR"){                                                                 │
// │  >1976                c(illumSource, srcPara);                                              │
// │   1977            }else if(srcShape == "DIPOLE-FAN"){                                                              │
// │   1978                GetNonSmoothedDipoleFanSrc(illumSource, srcPara);                                            │
// │   1979            }else if(srcShape == "QUASAR"){                                                                  │
// │   1980                GetNonSmoothedQuadraFanSrc(illumSource, srcPara);                                            │
// │   1981            }else if(srcShape == "2-LEAF"){                                                                  │
// │   1982                GetNonSmoothedDipoleLeafSrc(illumSource, srcPara);                                           │
// │   1983            }else if(srcShape == "4-LEAF"){                                                                  │
// │   1984                GetNonSmoothedQuadraLeafSrc(illumSource, srcPara);                                           │
// │   1985            }  

//     sourcePara_.param.clear();
//     std::string srcShape = GetPyDictValue<std::string>(srcPara, "type");
//     if (srcShape == "CIRCULAR"){ //圆形
//         OLreal sigma = GetPyDictValue<OLreal>(srcPara, "sigma");
//         OLreal ShiftX = GetPyDictValue<OLreal>(srcPara, "shift_x");
//         OLreal ShiftY = GetPyDictValue<OLreal>(srcPara, "shift_y");
//         sourcePara_.type = "CIRCULAR";
//         sourcePara_.param.push_back(sigma);
//         sourcePara_.param.push_back(ShiftX);
//         sourcePara_.param.push_back(ShiftY);
//     }else if(srcShape == "ANNULAR"){ //环形
//         OLreal sigmaOut = GetPyDictValue<OLreal>(srcPara, "sigma_out");
//         OLreal sigmaIn = GetPyDictValue<OLreal>(srcPara, "sigma_in");
//         OLreal innerShiftX = GetPyDictValue<OLreal>(srcPara, "inner_shift_x");
//         OLreal innerShiftY = GetPyDictValue<OLreal>(srcPara, "inner_shift_y");
//         OLreal ShiftX = GetPyDictValue<OLreal>(srcPara, "shift_x");
//         OLreal ShiftY = GetPyDictValue<OLreal>(srcPara, "shift_y");
//         sourcePara_.type = "ANNULAR";
//         sourcePara_.param.push_back(sigmaOut);
//         sourcePara_.param.push_back(sigmaIn);
//         sourcePara_.param.push_back(innerShiftX);
//         sourcePara_.param.push_back(innerShiftY);
//         sourcePara_.param.push_back(ShiftX);
//         sourcePara_.param.push_back(ShiftY);
//     }else if(srcShape == "DIPOLE-FAN"){ //二级扇形
//         OLreal sigmaOut = GetPyDictValue<OLreal>(srcPara, "sigma_out");
//         OLreal sigmaIn = GetPyDictValue<OLreal>(srcPara, "sigma_in");
//         OLreal rotAngle = GetPyDictValue<OLreal>(srcPara, "rot_angle");
//         OLreal spanAngle = GetPyDictValue<OLreal>(srcPara, "span_angle");
//         OLreal innerShiftX = GetPyDictValue<OLreal>(srcPara, "inner_shift_x");
//         OLreal innerShiftY = GetPyDictValue<OLreal>(srcPara, "inner_shift_y");
//         OLreal ShiftX = GetPyDictValue<OLreal>(srcPara, "shift_x");
//         OLreal ShiftY = GetPyDictValue<OLreal>(srcPara, "shift_y");
//         sourcePara_.type = "DIPOLE-FAN";
//         sourcePara_.param.push_back(sigmaOut);
//         sourcePara_.param.push_back(sigmaIn);
//         sourcePara_.param.push_back(rotAngle);
//         sourcePara_.param.push_back(spanAngle);
//         sourcePara_.param.push_back(innerShiftX);
//         sourcePara_.param.push_back(innerShiftY);
//         sourcePara_.param.push_back(ShiftX);
//         sourcePara_.param.push_back(ShiftY);
//     }else if(srcShape == "QUASAR"){ //四级扇形
//         OLreal sigmaOut = GetPyDictValue<OLreal>(srcPara,"sigma_out");
//         OLreal sigmaIn = GetPyDictValue<OLreal>(srcPara,"sigma_in");
//         OLreal rotAngle = GetPyDictValue<OLreal>(srcPara,"rot_angle");
//         OLreal spanAngle = GetPyDictValue<OLreal>(srcPara,"span_angle");
//         OLreal innerShiftX = GetPyDictValue<OLreal>(srcPara,"inner_shift_x");
//         OLreal innerShiftY = GetPyDictValue<OLreal>(srcPara,"inner_shift_y");
//         OLreal ShiftX = GetPyDictValue<OLreal>(srcPara, "shift_x");
//         OLreal ShiftY = GetPyDictValue<OLreal>(srcPara, "shift_y");
//         sourcePara_.type = "QUASAR";
//         sourcePara_.param.push_back(sigmaOut);
//         sourcePara_.param.push_back(sigmaIn);
//         sourcePara_.param.push_back(rotAngle);
//         sourcePara_.param.push_back(spanAngle);
//         sourcePara_.param.push_back(innerShiftX);
//         sourcePara_.param.push_back(innerShiftY);
//         sourcePara_.param.push_back(ShiftX);
//         sourcePara_.param.push_back(ShiftY);
//     }else if(srcShape == "2-LEAF"){ //二级叶形
//         OLreal D = GetPyDictValue<OLreal>(srcPara,"sigma_out");
//         OLreal d = GetPyDictValue<OLreal>(srcPara,"sigma_in");
//         OLreal rotAngle = GetPyDictValue<OLreal>(srcPara,"rot_angle"); 
//         OLreal ShiftX = GetPyDictValue<OLreal>(srcPara, "shift_x");
//         OLreal ShiftY = GetPyDictValue<OLreal>(srcPara, "shift_y");  
//         sourcePara_.type = "2-LEAF";
//         sourcePara_.param.push_back(D);
//         sourcePara_.param.push_back(d);
//         sourcePara_.param.push_back(rotAngle);
//         sourcePara_.param.push_back(ShiftX);
//         sourcePara_.param.push_back(ShiftY);
//     }else if(srcShape == "4-LEAF"){ //四级叶形
//         OLreal D = GetPyDictValue<OLreal>(srcPara,"sigma_out");
//         OLreal d = GetPyDictValue<OLreal>(srcPara,"sigma_in");
//         OLreal aspectRatio = GetPyDictValue<OLreal>(srcPara,"aspect_ratio");
//         OLreal rotAngle = GetPyDictValue<OLreal>(srcPara,"rot_angle");
//         OLreal ShiftX = GetPyDictValue<OLreal>(srcPara, "shift_x");
//         OLreal ShiftY = GetPyDictValue<OLreal>(srcPara, "shift_y");
//         sourcePara_.type = "4-LEAF";
//         sourcePara_.param.push_back(D);
//         sourcePara_.param.push_back(d);
//         sourcePara_.param.push_back(aspectRatio);
//         sourcePara_.param.push_back(rotAngle);
//         sourcePara_.param.push_back(ShiftX);
//         sourcePara_.param.push_back(ShiftY);
//     }
// }