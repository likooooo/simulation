#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <type_traist_notebook/uca/backend.hpp>
#include <py_helper.hpp>
#include <assert.h>

template<class T> uca::backend<T>& get_math_backend();

template<class T, class TTo> struct rebind
{
    using type = TTo;
};
template<class T, class TTo> using rebind_t = typename rebind<T, TTo>::type;
template<class T, size_t N, class TTo> struct rebind<std::array<T, N>, TTo>
{
    using type = std::array<rebind_t<T, TTo>, N>;
};
template<class T, class TTo> struct rebind<std::vector<T>, TTo>
{
    using type = std::vector<rebind_t<T, TTo>>;
};

using dbu_t = int64_t;
template<size_t N> using point_nd_dbu = vec<dbu_t, N>;

template<class T, class T1> rebind_t<T, T1> dbu_to_physical(const T& t, T1 dbu)
{
    return convert_to<rebind_t<T, T1>>(t) * dbu;
}
template<class T, class T1> rebind_t<T, dbu_t> physical_to_dbu(const T& t, T1 dbu)
{
    return convert_to<rebind_t<T, dbu_t>>(t / dbu);
}
// template<class T, class T1> T dbu_to_physical(const T& t, T1 dbu){return t * dbu;}
// template<class T, class T1> T physical_to_dbu(const T& t, T1 dbu){return t / dbu;}
template<class rT> inline rT get_lens_cutoff_frequency(rT lambda, rT NA)
{
    return 2 * NA / lambda;
}
template<class rT> inline rT get_simulation_system_cutoff_frequency(rT lambda, rT NA)
{
    //== 从 near field 到 aerial image 的 intensity 会有一个平方, 带宽需要 *2
    constexpr size_t doubling_bandwith_for_squaring_intensity = 1;
    return get_lens_cutoff_frequency(lambda, NA) * doubling_bandwith_for_squaring_intensity;
}
template<class rT, size_t _dim = 2> struct grid_info
{
    constexpr static size_t dim = _dim;
    using point_dbu_t = point_nd_dbu<dim>;
    using point_physical_t = vec<rT, dim>;

    rT dbu;
    struct {point_dbu_t start, step;} spatial;
    struct {point_physical_t start, step;} fourier;
    vec<size_t, dim> tilesize;
    //== is_coord_at_center = true 离散点的值代表区间的中点 
    bool is_coord_at_center {false};
    point_dbu_t shape_in_dbu() const
    {
        return spatial.step * tilesize;
    }
    point_dbu_t operator [](const point_dbu_t& pos_in_dbu) const
    {
        point_dbu_t idx = (pos_in_dbu - spatial.start) / spatial.step;
        assert(!(tilesize > idx));
        return idx;
    }
    vec<size_t, dim> stride() const
    {
        vec<size_t, dim> s;
        s.fill(1);
        for(size_t i = 1; i < dim; i++){
            for(size_t j = 0; j < i; j++){
                s.at(i) *= tilesize.at(j);
            }
        }
        // std::reverse(s.begin(), s.end());
        return s;
    }
    size_t operator()(const vec<size_t, dim>& indices_which_last_dim_change_fastest) const
    {
        vec<size_t, dim> s = stride();
        s *= indices_which_last_dim_change_fastest;
        return std::accumulate(s.begin(), s.end(), size_t(0));
    }
    size_t total_size() const
    {
        size_t n = 1;
        for(auto i : tilesize) n *= i;
        return n;
    }
    vec2<point_dbu_t> get_spatial_roi() const
    {
        return {spatial.start, spatial.start + spatial.step * tilesize};
    }
    vec2<point_physical_t> get_fourier_roi() const
    {
        return {fourier.start, fourier.start + fourier.step * tilesize};
    }
    template<class T>void display(const std::vector<T>& rowdata) const
    {
        imshow(rowdata, convert_to<std::vector<size_t>>(tilesize));
    }

    // //== user definition
    // template<class UnaryOperation> point_dbu_t physical_to_dbu(const point_physical_t& pos_in_physical, 
    //     UnaryOperation op = [](rT n){return typename point_dbu_t::value_type(std::floor(n));}
    // )
    // {
    //     point_dbu_t n;
    //     std::transform(pos_in_physical.begin(), pos_in_physical.end(), n.begin(), op);
    //     return n;
    // }
    static grid_info create_grid_info(vec<size_t, dim> shape, rT lambda, rT sigma , rT NA, vec2<point_physical_t> roi, rT dbu)
    {
        point_physical_t pitch  = roi.at(1) - roi.at(0);
        grid_info info;
        info.dbu = dbu;
        info.tilesize = shape;
        info.fourier.start = {0};
        info.fourier.step  = lambda / (NA + sigma) / pitch; 
        info.spatial.start = physical_to_dbu(roi.at(0), dbu);
        info.spatial.step  = physical_to_dbu(pitch / shape, dbu);
        return info;
    }
    static vec<size_t, dim> nyquist_sampling_rate(point_physical_t pitch, rT lambda, rT sigma , rT NA, size_t USF = 8)
    {
        rT cutoff_freq = get_lens_cutoff_frequency(lambda, NA + sigma);
        vec<size_t, dim> shape;
        for(size_t i = 0; i < dim; i++) shape.at(i) = USF * std::ceil(pitch.at(i) * 2 * cutoff_freq);
        return shape;
    }
    //== shape 可以改变
    static grid_info create_grid_info_bloch_mode(vec<size_t, dim> min_shape, rT lambda, rT sigma , rT NA, vec2<point_physical_t> roi, rT dbu)
    {
        grid_info info;
        info.dbu = dbu;
        point_physical_t pitch  = roi.at(1) - roi.at(0);
        //== [0, 2pi)
        info.fourier.start = {0}; 
        info.fourier.step  = lambda / (NA + sigma) / pitch;
        rT cutoff_freq = get_lens_cutoff_frequency(lambda, NA + sigma);
        for(size_t i = 0; i < dim; i++){
            dbu_t size = min_shape.at(i);
            info.tilesize.at(i) = std::max<dbu_t>(std::ceil(pitch.at(i) * 2 * cutoff_freq /* -cutoff_freq ~ +cutoff_freq */), size);
            if(size != info.tilesize.at(i)) {
                error_unclassified::out("create grid info size too small in dim-", i, 
                    ". reset from ", size, " to ", info.tilesize.at(i)
                );
            }
        }
        info.spatial.start = physical_to_dbu(roi.at(0), dbu);
        info.spatial.step  = physical_to_dbu(pitch / info.tilesize, dbu);
        return info;
    }
    //== roi 可以改变, shape 固定
    static grid_info create_grid_info_opc_mode(vec<size_t, dim> shape, rT lambda, rT sigma, rT NA, vec2<point_physical_t> roi, rT dbu)
    {
        grid_info info;
        info.dbu = dbu;
        info.tilesize = shape;
        rT cutoff_freq = get_lens_cutoff_frequency(lambda, NA + sigma);
        info.spatial.step.fill(std::ceil(rT(1) / cutoff_freq / dbu));
        point_dbu_t pitch_in_dbu = info.spatial.step * shape;
        info.spatial.start = (physical_to_dbu(roi[0] + roi[1], dbu) - pitch_in_dbu) / 2;
        info.fourier.start = {0};
        info.fourier.step = lambda / (NA + sigma) / dbu_to_physical(pitch_in_dbu, dbu);
        return info;
    }
    static grid_info create_grid_info_opc_mode(point_dbu_t shape, rT lambda, rT sigma, rT NA, point_physical_t center, rT dbu)
    {
        return create_grid_info_opc_mode(shape, lambda, sigma, NA, {center, center}, dbu);
    }
    
    static rT get_demangnification(rT sigma, rT NA)
    {
        return (sigma + NA) / NA;
    }
    //== without multiply k0. (2pi *NA/lambda)
    static vec2<vec3<rT>> mask_pov_k_space_boundary(rT NA, rT sigma, vec2<rT> shift_effict_by_chief_ray = {0, 0})
    { 
        //== k-space boundary without lens effect 
        vec2<rT> center{0, 0};
        rT radius = rT(1) / NA;
        vec3<rT> gray_boundary{center[0], center[1], radius};
        
        //== k-space boundary with lens effect 
        center += shift_effict_by_chief_ray;
        radius = rT(1) / get_demangnification(sigma, NA);
        vec3<rT> blue_boundary{center[0], center[1], radius};
        return {gray_boundary, blue_boundary};
    }
    static vec2<vec3<rT>> wafer_pov_k_space_boundary(rT NA)
    {
        vec2<rT> center{0, 0};
        rT radius = rT(1) / NA;
        vec3<rT> gray_boundary{center[0], center[1], radius};
        vec3<rT> blue_boundary{center[0], center[1], 1};
        return {gray_boundary, blue_boundary};
    }
};
template<class rT, size_t dim = 2> inline std::ostream& operator << (std::ostream& stream, const grid_info<rT, dim>& in) 
{
    std::vector<std::tuple<std::string, std::string>> msg{
        std::make_tuple(std::string("tilesize      :"), to_string(in.tilesize)),
        std::make_tuple(std::string("spatial start :"), to_string(dbu_to_physical(in.spatial.start, in.dbu))),
        std::make_tuple(std::string("spatial step  :"), to_string(dbu_to_physical(in.spatial.step, in.dbu))),
        std::make_tuple(std::string("fourier start :"), to_string(in.fourier.start)),
        std::make_tuple(std::string("fourier step  :"), to_string(in.fourier.step)), 
        std::make_tuple(std::string("coords place  :"), std::string(vec2<const char*>{"|*--|corner", "|-*-|center"}[size_t(in.is_coord_at_center)])), 
    };
    print_table(stream, msg, {"* grid info", ""}, -1);
    return stream;
}
// template<class T, size_t DIM=2>
// struct grid_start_step
// {
//     using value_type = T;
//     constexpr static size_t dim = DIM;
//     template<size_t DIM1> using rebind_t = grid_start_step<T, DIM1>; 

//     struct StartStep{vec<T, DIM> start, step;};
//     StartStep spatial, fourier;
//     vec<size_t, DIM> tilesize;
// };
// template<class TMeta, size_t N> inline auto change_dim(const TMeta& in)
// {
//     typename TMeta::rebind<N> meta{0};
//     for(size_t i = 0; i < TMeta::dim; i++){
//         meta.tilesize[i]      = in.tilesize[i];
//         meta.spatial.start[i] = in.spatial.start[i];
//         meta.spatial.step[i]  = in.spatial.step[i];
//         meta.fourier.start[i] = in.fourier.start[i];
//         meta.fourier.step[i]  = in.fourier.step[i];
//     }
//     return meta;
// }
// template<class TMeta, class print_to = debug_unclassified> inline void print_grid_start_step(const TMeta& in, const std::string& tag = "") 
// {
//     std::vector<std::tuple<std::string, std::string>> msg{
//         std::make_tuple(std::string("tilesize      :"), to_string(in.tilesize)),
//         std::make_tuple(std::string("spatial start :"), to_string(in.spatial.start)),
//         std::make_tuple(std::string("spatial step  :"), to_string(in.spatial.step)),
//         std::make_tuple(std::string("fourier start :"), to_string(in.fourier.start)),
//         std::make_tuple(std::string("fourier step  :"), to_string(in.fourier.step)), 
//     };
//     print_to(msg, {"* " + tag, TypeReflection<TMeta>()});
// }

// template<class T>using point = vec2<T>;
// template<class T>using rectangle = vec2<point<T>>;


// template<class T> inline grid_start_step<T> optical_numerics(const rectangle<T>& roi, const vec2<T>& ambit, const vec2<size_t>& tilesize, T maxNA, T lambda)
// {
//     grid_start_step<T> grid;
//     grid.tilesize = tilesize;

//     //== 1. spatial step
//     vec2<T> nyquist_freq;
//     nyquist_freq.fill(((T(2) * maxNA) / lambda) * doubling_bandwith_for_squaring_signal_system);
//     grid.spatial.step = vec2<T>{T(1), T(1)} / nyquist_freq;

//     //== 2. spatial start
//     vec2<T> spatial_domain;
//     for(int i = 0; i < spatial_domain.size(); i++)
//     spatial_domain.at(i) = tilesize.at(i) / nyquist_freq.at(i); 
//     const auto& [from, to] = roi;
//     grid.spatial.start = (from + to - spatial_domain) * T(0.5); 
//     if((spatial_domain - (ambit * 2)) < (to - from)) // check
//     {
//         print_table(std::cerr, std::vector<std::tuple<std::string, vec2<T>>>{
//             std::make_tuple(std::string("roi from       :"), from), 
//             std::make_tuple(std::string("roi to         :"), to),
//             std::make_tuple(std::string("ambit          :"), ambit),
//             std::make_tuple(std::string("spatial domain :"), spatial_domain),
//             std::make_tuple(std::string("need more tile :"), ((to - from) - (spatial_domain - (ambit * 2))) * nyquist_freq),
//         }, {"* roi or ambit is too large", ""});
//     }

//     //== 3. fourier start & step
//     grid.fourier.step = vec2<T>{lambda, lambda} / spatial_domain;
//     grid.fourier.start = vec2<T>{0, 0}; // use DC-corner
//     return grid;
// }
// template<class T> inline grid_start_step<T> bloch_optical_numerics(rectangle<T> spatial_domain, vec2<size_t> tilesize, T maxNA, T lambda)
// {
//     grid_start_step<T> grid;

//     //== 1. fourier start & step
//     grid.fourier.step = vec2<T>{lambda, lambda} / (spatial_domain[1] - spatial_domain[0]);
//     grid.fourier.start = vec2<T>{0, 0}; 
    
//     //== 2. tilesize
//     grid.tilesize[0] = std::max(tilesize[0], doubling_bandwith_for_squaring_signal_system * size_t(std::ceil(2 * maxNA / grid.fourier.step[0])));
//     grid.tilesize[1] = std::max(tilesize[1], doubling_bandwith_for_squaring_signal_system * size_t(std::ceil(2 * maxNA / grid.fourier.step[1])));
    
//     //== 3. spatial start & step
//     grid.spatial.start = spatial_domain.at(0);
//     grid.spatial.step = (spatial_domain[1] - spatial_domain[0]) / grid.tilesize;
//     return grid;
// }

// #include "geometry.hpp"
// using grid_info_in_dbu = dbu_grid_start_step<double>;
// template<class TUserConfig> dbu_grid_start_step<typename TUserConfig::value_type> optical_numerics_in_dbu(line_dbu cutline, const TUserConfig& config)
// {
//     using rT = typename TUserConfig::value_type;
//     auto roi = convert<line_dbu, rectangle<rT>>{}(cutline);
//     dbu_to_physical(roi, config.dbu);
//     auto grid_in_um = optical_numerics<rT>(roi, config.ambit, config.tilesize, config.maxNA, config.wavelength); 
//     // print_grid_start_step(grid_in_um, "   origin grid-um");

//     //== 1. spatial step
//     physical_to_dbu(grid_in_um.spatial.step, config.dbu);
//     vec2<int64_t> spatial_step_in_dbu = {int64_t(std::floor(grid_in_um.spatial.step[0])), int64_t(std::floor(grid_in_um.spatial.step[1]))};
//     //== 2. spatial domain & start
//     auto tile_size = convert<vec2<size_t>, vec2<int64_t>>{}(config.tilesize);
//     vec2<int64_t> spatial_domain_in_dbu = tile_size * spatial_step_in_dbu;
//     const auto& [from_in_dbu, to_in_dbu] = cutline;
//     vec2<int64_t> spatial_start_in_dbu = ((from_in_dbu + to_in_dbu - spatial_domain_in_dbu) + 1) / 2; 
    
//     //== check tilesize
//     vec2<rT> ambit = config.ambit; physical_to_dbu(ambit, config.dbu);
//     vec2<int64_t> ambit_in_dbu = {int64_t(std::ceil(ambit[0])), int64_t(std::ceil(ambit[1]))};
//     vec2<int64_t> span_in_dbu = to_in_dbu - from_in_dbu;
//     int64_t length_in_dbu = std::max(span_in_dbu[0], span_in_dbu[1]);
//     span_in_dbu = {length_in_dbu, length_in_dbu};
//     if((spatial_domain_in_dbu - (ambit_in_dbu * 2)) < span_in_dbu)
//     {
//         print_table(std::cerr, std::vector<std::tuple<std::string, vec2<int64_t>>>{
//             std::make_tuple(std::string("roi from       :"), from_in_dbu), 
//             std::make_tuple(std::string("roi to         :"), to_in_dbu),
//             std::make_tuple(std::string("ambit          :"), ambit_in_dbu),
//             std::make_tuple(std::string("spatial domain :"), spatial_domain_in_dbu),
//             std::make_tuple(std::string("need more tile :"), ((to_in_dbu - from_in_dbu) - (spatial_domain_in_dbu - (ambit_in_dbu * 2))) / spatial_step_in_dbu),
//         }, {"*    roi or ambit is too large", ""});
//     }
//     grid_info_in_dbu grid_in_dbu;
//     // grid_start_step<rT> grid_in_dbu;
//     grid_in_dbu.tilesize = config.tilesize;
//     grid_in_dbu.spatial.start = spatial_start_in_dbu;
//     grid_in_dbu.spatial.step  = spatial_step_in_dbu; 
//     grid_in_dbu.dbu = config.dbu;
//     rT lambda_in_dbu      = config.wavelength; 
//     physical_to_dbu(lambda_in_dbu, config.dbu);
//     grid_in_dbu.fourier.step  = vec2<rT>{lambda_in_dbu, lambda_in_dbu} / spatial_domain_in_dbu;
//     grid_in_dbu.fourier.start = vec2<rT>{0, 0}; 
    
//     if(debug_unclassified::verbose())
//     {
//         grid_start_step<rT> grid_to_um;
//         grid_to_um.tilesize      = grid_in_dbu.tilesize; 
//         grid_to_um.fourier       = grid_in_dbu.fourier;
//         grid_to_um.spatial.start = convert_to<vec2<rT>>(grid_in_dbu.spatial.start); 
//         grid_to_um.spatial.step  = convert_to<vec2<rT>>(grid_in_dbu.spatial.step); 
//         dbu_to_physical(grid_to_um.spatial.start, config.dbu);
//         dbu_to_physical(grid_to_um.spatial.step, config.dbu);
//         static bool alread_printted = false;
//         if(!alread_printted){
//             alread_printted = true;
//             print_grid_start_step(grid_to_um, "    simulation domain(um)");
//         }
//         else{
//             debug_unclassified::out("*    simulation domain(um)\nspatial start :              ", grid_to_um.spatial.start);
//         }
//     }
//     return grid_in_dbu;
// }
