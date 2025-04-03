#pragma once
#include <type_traist_notebook/type_traist.hpp>

template<class T, size_t DIM=2>
struct grid_start_step
{
    using value_type = T;
    constexpr static size_t dim = DIM;
    template<size_t DIM1> using rebind_t = grid_start_step<T, DIM1>; 

    struct StartStep{vec<T, DIM> start, step;};
    StartStep spatial, fourier;
    vec<size_t, DIM> tilesize;
};
template<class TMeta, size_t N> inline auto change_dim(const TMeta& in)
{
    typename TMeta::rebind<N> meta{0};
    for(size_t i = 0; i < TMeta::dim; i++){
        meta.tilesize[i]      = in.tilesize[i];
        meta.spatial.start[i] = in.spatial.start[i];
        meta.spatial.step[i]  = in.spatial.step[i];
        meta.fourier.start[i] = in.fourier.start[i];
        meta.fourier.step[i]  = in.fourier.step[i];
    }
    return meta;
}
template<class TMeta, class print_to = debug_unclassified> inline void print_grid_start_step(const TMeta& in, const std::string& tag = "") 
{
    std::vector<std::tuple<std::string, std::string>> msg{
        std::make_tuple(std::string("tilesize      :"), to_string(in.tilesize)),
        std::make_tuple(std::string("spatial start :"), to_string(in.spatial.start)),
        std::make_tuple(std::string("spatial step  :"), to_string(in.spatial.step)),
        std::make_tuple(std::string("fourier start :"), to_string(in.fourier.start)),
        std::make_tuple(std::string("fourier step  :"), to_string(in.fourier.step)), 
    };
    print_to(msg, {"* " + tag, TypeReflection<TMeta>()});
}

template<class T>using point = vec2<T>;
template<class T>using rectangle = vec2<point<T>>;


constexpr size_t doubling_bandwith_for_squaring_signal_system = 2;
template<class T> inline grid_start_step<T> optical_numerics(const rectangle<T>& roi, const vec2<T>& ambit, const vec2<size_t>& tilesize, T maxNA, T lambda)
{
    grid_start_step<T> grid;
    grid.tilesize = tilesize;

    //== 1. spatial step
    vec2<T> nyquist_freq;
    nyquist_freq.fill(((T(2) * maxNA) / lambda) * doubling_bandwith_for_squaring_signal_system);
    grid.spatial.step = vec2<T>{T(1), T(1)} / nyquist_freq;

    //== 2. spatial start
    vec2<T> spatial_domain;
    for(int i = 0; i < spatial_domain.size(); i++)
    spatial_domain.at(i) = tilesize.at(i) / nyquist_freq.at(i); 
    const auto& [from, to] = roi;
    grid.spatial.start = (from + to - spatial_domain) * T(0.5); 
    if((spatial_domain - (ambit * 2)) < (to - from)) // check
    {
        print_table(std::cerr, std::vector<std::tuple<std::string, vec2<T>>>{
            std::make_tuple(std::string("roi from       :"), from), 
            std::make_tuple(std::string("roi to         :"), to),
            std::make_tuple(std::string("ambit          :"), ambit),
            std::make_tuple(std::string("spatial domain :"), spatial_domain),
            std::make_tuple(std::string("need more tile :"), ((to - from) - (spatial_domain - (ambit * 2))) * nyquist_freq),
        }, {"* roi or ambit is too large", ""});
    }

    //== 3. fourier start & step
    grid.fourier.step = vec2<T>{lambda, lambda} / spatial_domain;
    grid.fourier.start = vec2<T>{0, 0}; // use DC-corner
    return grid;
}
template<class T> inline grid_start_step<T> bloch_optical_numerics(rectangle<T> spatial_domain, vec2<size_t> tilesize, T maxNA, T lambda)
{
    grid_start_step<T> grid;

    //== 1. fourier start & step
    grid.fourier.step = vec2<T>{lambda, lambda} / (spatial_domain[1] - spatial_domain[0]);
    grid.fourier.start = vec2<T>{0, 0}; 
    
    //== 2. tilesize
    grid.tilesize[0] = std::max(tilesize[0], doubling_bandwith_for_squaring_signal_system * size_t(std::ceil(2 * maxNA / grid.fourier.step[0])));
    grid.tilesize[1] = std::max(tilesize[1], doubling_bandwith_for_squaring_signal_system * size_t(std::ceil(2 * maxNA / grid.fourier.step[1])));
    
    //== 3. spatial start & step
    grid.spatial.start = spatial_domain.at(0);
    grid.spatial.step = (spatial_domain[1] - spatial_domain[0]) / grid.tilesize;
    return grid;
}

#include "geometry.hpp"
template<class T> struct dbu_grid_start_step
{
    using value_type = int64_t;
    constexpr static size_t dim = 2;
    template<size_t DIM1> using rebind_t = grid_start_step<int64_t, DIM1>; 
    struct StartStep{point_dbu start, step;}spatial;
    typename grid_start_step<T, dim>::StartStep fourier;
    vec<size_t, dim> tilesize;
};
using grid_info_in_dbu = dbu_grid_start_step<double>;
template<class TUserConfig> grid_info_in_dbu optical_numerics_in_dbu(cutline_dbu cutline, const TUserConfig& config)
{
    auto roi = convert<cutline_dbu, rectangle<double>>{}(cutline);
    dbu_to_um(roi, config.dbu);
    auto grid_in_um = optical_numerics<double>(roi, config.ambit, config.tilesize, config.maxNA, config.wavelength); 
    print_grid_start_step(grid_in_um, "origin grid-um");

    //== 1. spatial step
    um_to_dbu(grid_in_um.spatial.step, config.dbu);
    vec2<int64_t> spatial_step_in_dbu = {int64_t(std::floor(grid_in_um.spatial.step[0])), int64_t(std::floor(grid_in_um.spatial.step[1]))};
    //== 2. spatial domain & start
    auto tile_size = convert<vec2<size_t>, vec2<int64_t>>{}(config.tilesize);
    vec2<int64_t> spatial_domain_in_dbu = tile_size * spatial_step_in_dbu;
    const auto& [from_in_dbu, to_in_dbu] = cutline;
    vec2<int64_t> spatial_start_in_dbu = ((from_in_dbu + to_in_dbu - spatial_domain_in_dbu) + 1) / 2; 
    
    //== check tilesize
    vec2<double> ambit = config.ambit; um_to_dbu(ambit, config.dbu);
    vec2<int64_t> ambit_in_dbu = {int64_t(std::ceil(ambit[0])), int64_t(std::ceil(ambit[1]))};
    vec2<int64_t> span_in_dbu = to_in_dbu - from_in_dbu;
    int64_t length_in_dbu = std::max(span_in_dbu[0], span_in_dbu[1]);
    span_in_dbu = {length_in_dbu, length_in_dbu};
    if((spatial_domain_in_dbu - (ambit_in_dbu * 2)) < span_in_dbu)
    {
        print_table(std::cerr, std::vector<std::tuple<std::string, vec2<int64_t>>>{
            std::make_tuple(std::string("roi from       :"), from_in_dbu), 
            std::make_tuple(std::string("roi to         :"), to_in_dbu),
            std::make_tuple(std::string("ambit          :"), ambit_in_dbu),
            std::make_tuple(std::string("spatial domain :"), spatial_domain_in_dbu),
            std::make_tuple(std::string("need more tile :"), ((to_in_dbu - from_in_dbu) - (spatial_domain_in_dbu - (ambit_in_dbu * 2))) / spatial_step_in_dbu),
        }, {"* roi or ambit is too large", ""});
    }
    grid_info_in_dbu grid_in_dbu;
    // grid_start_step<double> grid_in_dbu;
    grid_in_dbu.tilesize = config.tilesize;
    grid_in_dbu.spatial.start = spatial_start_in_dbu;
    grid_in_dbu.spatial.step  = spatial_step_in_dbu; 
    double lambda_in_dbu      = config.wavelength; 
    um_to_dbu(lambda_in_dbu, config.dbu);
    grid_in_dbu.fourier.step  = vec2<double>{lambda_in_dbu, lambda_in_dbu} / spatial_domain_in_dbu;
    grid_in_dbu.fourier.start = vec2<double>{0, 0}; 
    
    if(debug_unclassified::verbose())
    {
        grid_start_step<double> grid_to_um;
        grid_to_um.tilesize      = grid_in_dbu.tilesize; 
        grid_to_um.fourier       = grid_in_dbu.fourier;
        grid_to_um.spatial.start = convert_to<vec2<double>>(grid_in_dbu.spatial.start); 
        grid_to_um.spatial.step  = convert_to<vec2<double>>(grid_in_dbu.spatial.step); 
        dbu_to_um(grid_to_um.spatial.start, config.dbu);
        dbu_to_um(grid_to_um.spatial.step, config.dbu);
        print_grid_start_step(grid_to_um, "grid-dbu-aligined to grid-um");
    }
    return grid_in_dbu;
}
