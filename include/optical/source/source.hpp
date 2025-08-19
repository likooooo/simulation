#pragma once
#include "parametric_source.hpp"
#include <optical/simulation_grid_info.hpp>
#include <kernels/kernel_loop.hpp>
#include <assert.h>
#include <py_helper.hpp>

template<class T> T k0(T lambda, T NA_with_M)
{
    return 2_PI / lambda * NA_with_M;
}
template<class T> bool check_na_valid(T NA, T immersion_index)
{
    return immersion_index >= NA;
}
    
template<class T> struct source_point
{
    //== reference : panoramictech/UserGuide/Aerial/Illuminations/PupilMapIllumination.html
    using rT = real_t<T>;
    using cT = complex_t<T>;

    //== 强度正比于幅值的平方
    // \( I \propto A^2 \)
    rT intensity{1};

    // degree in rad
    rT e_field_direction{0.5_PI};

    // degree of polarization
    rT DOP{0};

    //  0       : linearly polarized
    // -1       : left  circularly polarized
    //  1       : right circularly polarized
    //  (-1, 0) : left  elliptical polarized
    //  ( 0, 1) : right  elliptical polarized
    rT ellipticity{0};

    // [0, 1]
    vec2<rT> sigmaxy{0, 0};

    bool is_valid() const
    {
        return 0 != intensity;
    }
    bool is_satisfy_bragg_condition() const
    {
        return int(sigmaxy.at(0)) == sigmaxy.at(0) && int(sigmaxy.at(0)) == sigmaxy.at(0);
    }
    vec2<rT> k_vector(rT lambda, rT NA_with_M = 1) const
    {
        return k0(lambda, NA_with_M) * sigmaxy;   
    }
    vec2<cT> polarization_state() const
    {
        //== https://en.wikipedia.org/wiki/Polarization_(waves)#Polarization_ellipse
        vec2<cT> polar{std::cos(e_field_direction), std::sin(e_field_direction)};
        //==
        // [1          0] [x]
        // [0 e^(i*pi/2)] [y]
        polar.at(1) *= std::exp(cT(0, 0.5_PI * ellipticity));
        return polar;
    }
    vec2<vec3<cT>> sp_polarization_state() const
    {
        vec3<cT> TE{0, 1, 0};
        vec3<cT> TM{std::exp(cT(0, 0.5_PI * ellipticity)), 0, 0};

        vec2<vec3<cT>> local_polarized_basis; 
        local_polarized_basis.at(0) = get_current_polarized_dir(TE);
        local_polarized_basis.at(1) = get_current_polarized_dir(TM);
        return local_polarized_basis;
    }
    //  0 : unpolarization
    //  1 : full TE-polarization
    // -1 : full TM-polarization
    static vec2<rT> get_TEM_coef(rT dop = 0)
    {
        rT TM = (1 - dop) * 0.5;
        rT TE = 1 - TM;
        return {TE, TM};
    }
    source_point decompose_unpolarized_component()
    {
        auto [te_component, tm_component] = (get_TEM_coef(DOP) * intensity);
        this->DOP = 1;
        this->intensity = te_component;

        source_point sp = *this;
        sp.intensity = tm_component;
        sp.e_field_direction = this->e_field_direction - 0.5_PI;
        return sp;
    }
    vec2<source_point> decompose_componets() const
    {
        source_point te = *this;
        source_point tm = te.decompose_unpolarized_component();
        return {te, tm};
    }
    vec2<rT> get_crao_azimuth() const
    {
        auto [alpha, beta] = sigmaxy;
        return {std::asin(std::hypot(alpha, beta)), std::atan2(beta, alpha)};
    }
    vec3<cT> get_current_polarized_dir(vec3<cT> normal_incidence_polor_dir = {0, 1, 0}) const
    {
        auto [crao, azimuth] = get_crao_azimuth();
        return convert_to<matrix3x3<cT>>(rotate_matrix(crao, azimuth)) | normal_incidence_polor_dir;
    }
    constexpr static matrix3x3<rT> rotate_matrix(rT crao = 0, rT azimuth = 0)
    {
        using std::cos, std::sin;
        rT theta = crao, phi = azimuth;
        matrix3x3<rT> m{
            cos(theta) * cos(phi),-sin(phi), sin(theta)*cos(phi),
            cos(theta) * sin(phi), cos(phi), sin(theta)*sin(phi),
                      -sin(theta),    rT(0),          cos(theta)
        };
        return m;
    }
    bool operator == (const source_point& o) const
    {
        return 0 == std::memcmp(&o, this, sizeof(source_point));
    }
};
template<class T> std::ostream& operator<<(std::ostream& s, const source_point<T> & sp) 
{
    return s << sp.sigmaxy[0] << sp.sigmaxy[1] << sp.intensity << sp.e_field_direction << sp.DOP << sp.ellipticity;
}
template<class T> std::ostream& operator<<(std::ostream& s, const std::vector<source_point<T>> & a) 
{
    using print_type = std::tuple<size_t, T, T, T, T, T, T>;
    std::vector<print_type> lines; 
    lines.reserve(a.size());
    for(size_t i = 0; i < a.size(); i++){
        const source_point<T>& sp = a.at(i);
        auto dir = sp.e_field_direction * real_t<T>(180)/1_PI;
        lines.push_back(print_type(
            i, sp.sigmaxy[0], sp.sigmaxy[1], sp.intensity, 
            dir, sp.DOP, sp.ellipticity
        ));
    }
    print_table(s, lines, {"#", "Sigma-X", "Sigma-Y", "Intensity", "E-field-direction", "DOP", "Ellipticity"}, -1);
    return s;
}

template<class rT> std::vector<vec2<int>> get_diffraction_order(const grid_info<rT, 2>& info, vec2<rT> offset_sigmaxy)
{
    auto [start, step] = info.fourier; // Wafer P.O.V
    start += offset_sigmaxy;
    vec2<int> lb = convert_to<vec2<int>>((vec2<rT>{-1, -1} + start) / step) - 1;
    vec2<int> ub = convert_to<vec2<int>>((vec2<rT>{1, 1}   + start) / step) + 1;

    std::vector<vec2<int>> orders;
    orders.reserve((ub[1] - ub[0]) * (lb[1] - lb[0]));
    for(int y = lb[1]; y < ub[1]; y++){
        for(int x = lb[0]; x < ub[0]; x++){
            vec2<int> order{x, y};
            if(1 < vector_norm((step * order) - start)) continue;
            orders.push_back(order);
        }
    }
    return orders;
}
enum polarization_basis
{
    Descartes = 0, TETM, X_Y_Zone, count
};
template<class T> struct source_grid
{
    using sp_t = source_point<T>;
    using rT = typename sp_t::rT;
    using cT = typename sp_t::cT;
    using parametric_source_t = parametric_source<T>;
    std::vector<sp_t> source_points;
    vec2<size_t> shape; 
    vec2<rT> step; 
    polarization_basis basis{polarization_basis::Descartes};
    grid_info<rT, 2> get_grid_info(rT lambda) const
    {
        grid_info<rT, 2> gi;
        gi.tilesize = shape;
        gi.fourier.step  = step;
        gi.fourier.start = step * (shape / 2 - 1);
        vec2<rT> pitch   = lambda / gi.fourier.step;
        gi.spatial.step  = pitch /(shape - 1);
        gi.spatial.start = gi.spatial.step * (shape / 2 - 1);
        gi.is_coord_at_center = true;
        return gi;
    }
    source_grid(size_t sample_size_odd, polarization_basis pb = polarization_basis::Descartes) : 
        source_points(std::vector<sp_t>(sample_size_odd * sample_size_odd)), 
        shape({sample_size_odd, sample_size_odd}), 
        step({rT(2) / (sample_size_odd - 1), rT(2) / (sample_size_odd - 1)}),
        basis(pb)
    {
        assert(sample_size_odd % 2 == 1);
    }
    source_grid(const std::vector<vec2<int>>& diffrac_orders, const grid_info<rT, 2>& info) :
        shape(info.tilesize), step(info.fourier.step), basis(polarization_basis::TETM)
    {
        source_points.reserve(diffrac_orders.size());
        for(vec2<int> ixy : diffrac_orders){
            sp_t sp;
            sp.sigmaxy   = step * ixy;
            sp.e_field_direction += std::atan2(ixy[1], ixy[0]);
        }
    }
    source_grid() = default;
    template<class TP> void init(const TP& params, rT e_field_direction = 0.5_PI, rT ellipticity = 0, size_t polarization = polarization_basis::Descartes)
    {
        size_t size = shape[0];
        source_grid& sg = *this;
        std::vector<rT> intensity_image(sg.source_points.size(), 0);
        if constexpr(std::is_same_v<TP, typename parametric_source_t::traditional_source_params>)
            parametric_source_t::get_traditional_source(intensity_image.data(), size, size, params); 
        else if constexpr(std::is_same_v<TP, typename parametric_source_t::annular_source_params>)
            parametric_source_t::get_annular_source(intensity_image.data(), size, size, params); 
        else if constexpr(std::is_same_v<TP, typename parametric_source_t::dipole_fan_source_params>)
            parametric_source_t::get_dipole_fan_source(intensity_image.data(), size, size, params); 
        else if constexpr(std::is_same_v<TP, typename parametric_source_t::quadratic_fan_source_params>)
            parametric_source_t::get_quadratic_fan_source(intensity_image.data(), size, size, params); 
        else if constexpr(std::is_same_v<TP, typename parametric_source_t::dipole_leaf_source_params>)
            parametric_source_t::get_dipole_leaf_source(intensity_image.data(), size, size, params); 
        else if constexpr(std::is_same_v<TP, typename parametric_source_t::quadratic_leaf_source_params>)
            parametric_source_t::get_quadratic_leaf_source(intensity_image.data(), size, size, params); 
        else 
            unreachable_constexpr_if<TP>();

        for(int y = 0; y < size; y++){
            for(int x = 0; x < size; x++){
                auto& sp = sg.source_points.at(y * size + x);
                sp.intensity = intensity_image.at(y * size + x); 
                sp.ellipticity = ellipticity;
                sp.e_field_direction = e_field_direction;
                sp.sigmaxy = sg.step * vec2<int>{x - int(size /2), y - int(size /2)};
                if(sg.basis == polarization_basis::TETM) sp.e_field_direction += std::atan2(sp.sigmaxy[1], sp.sigmaxy[0]);
            }
        }
    }
    template<class TP> static source_grid create(size_t size, const TP& params, rT e_field_direction = 0.5_PI, rT ellipticity = 0, size_t polarization = polarization_basis::Descartes)
    {
        source_grid sg(size, (polarization_basis)polarization);
        sg.init<TP>(params, e_field_direction, ellipticity, polarization);
        return sg;
    }
    
    static vec2<rT> get_dc_from_chief_ray(rT theta, rT phi)
    {
        using std::sin, std::cos;
        return {sin(theta) * cos(phi), sin(theta) * sin(phi)};
    }
    void shift_dc(const vec2<rT> DC_location = {0, 0})
    {
        if(DC_location == vec2<rT>{0, 0}) return;
        for(auto& sp : source_points) sp.sigmaxy -= DC_location;
    }
    void shift_dc(rT theta, rT phi, rT maskNA)
    {
        shift_dc(get_dc_from_chief_ray(theta, phi) / maskNA);
    }

    void clear_invalid_source_points(rT intensity_threshold = 1e-2)
    {
        source_points.erase(std::remove_if(source_points.begin(), source_points.end(), 
            [intensity_threshold](const sp_t& sp) { return sp.intensity < intensity_threshold;}), source_points.end()
        );
    }
    void decompose_polarized_components()
    {
        source_grid& sg = *this;
        std::vector<sp_t>& polarized_components = sg.source_points;
        sg.clear_invalid_source_points();
        polarized_components.reserve(polarized_components.size() * 2);
        const size_t n = polarized_components.size();
        for(size_t i = 0; i < n; i++)
        {
            sp_t sp = polarized_components.at(i).decompose_unpolarized_component();
            if(sp.is_valid()) polarized_components.push_back(sp);
        }
    }
    rT calc_max_sigma(const vec2<rT> DC_location = {0, 0}) const
    {
        std::vector<rT> sigma_r;
        sigma_r.reserve(source_points.size());
        std::transform(source_points.begin(), source_points.end(), std::back_insert_iterator(sigma_r), 
            [&DC_location](const sp_t& sp){
                if(!sp.is_valid()) return rT(0);
                vec2<rT> sigmaxy = sp.sigmaxy - DC_location;
                return std::hypot(sigmaxy.at(0), sigmaxy.at(1));
            }
        );
        return *std::max_element(sigma_r.begin(), sigma_r.end());
    }
    std::vector<rT> get_intensity_bitmap() const
    {
        std::vector<rT> intensity_image;
        intensity_image.reserve(source_points.size());
        for(size_t i = 0; i < source_points.size(); i++) intensity_image.push_back(source_points.at(i).intensity);
        return intensity_image;
    }
    std::vector<rT> get_e_field_direction() const
    {
        std::vector<rT> dir;
        dir.reserve(source_points.size());
        for(size_t i = 0; i < source_points.size(); i++) dir.push_back(source_points.at(i).e_field_direction);
        return dir;
    }
    std::vector<rT> get_ellipticity() const
    {
        std::vector<rT> ellipticity;
        ellipticity.reserve(source_points.size());
        for(size_t i = 0; i < source_points.size(); i++) ellipticity.push_back(source_points.at(i).ellipticity);
        return ellipticity;
    }
    std::vector<vec2<rT>> get_sigmaxy_bitmap() const
    {
        std::vector<vec2<rT>> sigmaxy;
        sigmaxy.reserve(source_points.size());
        for(size_t i = 0; i < source_points.size(); i++) sigmaxy.push_back(source_points.at(i).sigmaxy);
        return sigmaxy;
    }
    void plot(const std::string& title, const grid_info<rT>& mask_grid_info, 
        rT theta = 0, rT phi = 0, rT wavelength = 0, rT NA = 0, rT M = 0
    ) const
    {
        bool is_wafer_pov = (wavelength == 0 || NA == 0 || M == 0);

        source_grid sg = *this;
        //== shift mask-POV
        vec2<rT> shift_by_angle{0};
        if(!is_wafer_pov){
            sg.clear_invalid_source_points();
            shift_by_angle = get_dc_from_chief_ray(theta, phi) / (NA / M);
            vec2<rT> shift_sigma_xy = mask_grid_info.fourier.start + shift_by_angle;
            sg.shift_dc(shift_sigma_xy);
        }
        sg.decompose_polarized_components();

        std::vector<vec2<rT>> sigmaxy = sg.get_sigmaxy_bitmap();
        std::vector<rT> ellipticity = sg.get_ellipticity();
        std::vector<std::string> color;
        {
            color.reserve(sg.source_points.size());
            std::transform(ellipticity.begin(), ellipticity.end(), std::back_insert_iterator(color), 
                [](rT ellip){return ellip >= 0 ? "red" : "blue";}
            );
        }
        std::vector<rT> intensity = sg.get_intensity_bitmap();
        std::vector<rT> dir = sg.get_e_field_direction();

        for(vec2<int> ixy : get_diffraction_order(mask_grid_info, shift_by_angle)){
            vec2<rT> xy = mask_grid_info.fourier.step * ixy;
            sigmaxy.push_back(xy);
            dir.push_back(0);
            intensity.push_back(rT(1));
            ellipticity.push_back(0);
            color.push_back("green");

            sigmaxy.push_back(xy);
            dir.push_back(0.5_PI);
            intensity.push_back(rT(1));
            ellipticity.push_back(0);
            color.push_back("green");
        }

        if(!is_wafer_pov){
            for(auto& xy : sigmaxy) xy *= k0(wavelength, NA / M);
        }
        rT line_length_scalar = (0.5 / (shape[0] - 1));
        intensity *= line_length_scalar;
        plot_field<rT>(sigmaxy, dir, intensity, ellipticity, color, title);
    }
    void plot_wafer_pov(const grid_info<rT>& info = grid_info<rT>()) const
    {
        plot("wafer P.O.V", info);
    }
    void plot_mask_pov(const grid_info<rT>& info, rT theta, rT phi, rT wavelength, rT NA, rT M) const
    {
        plot("mask P.O.V", info, theta, phi, wavelength, NA, M);
    }
};
template<class T> std::ostream& operator<<(std::ostream& s, const source_grid<T> & a) 
{
    constexpr vec3<const char*>basis_str{"X-Y", "S-P", "TODO"};
    return s << " shape=" << a.shape << "\n polarization=" << 
       basis_str.at(size_t(a.basis)) << std::endl << a.source_points;
}

// 近场插值 : file:///G:/Document/YuWei/hyper%20litho/panoramictech/UserGuide/Articles/NonConstantScatteringCoefficients/default.html
// 坐标约定 : file:///G:/Document/YuWei/hyper%20litho/panoramictech/UserGuide/Articles/CoordinateSystemConventions/default.html
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