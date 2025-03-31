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
template<class TMeta> inline void print_grid_start_step(const TMeta& in, const std::string& tag = "") 
{
    std::vector<std::tuple<std::string, std::string>> msg{
        std::make_tuple(std::string("tilesize      :"), to_string(in.tilesize)),
        std::make_tuple(std::string("spatial start :"), to_string(in.spatial.start)),
        std::make_tuple(std::string("spatial step  :"), to_string(in.spatial.step)),
        std::make_tuple(std::string("fourier start :"), to_string(in.fourier.start)),
        std::make_tuple(std::string("fourier step  :"), to_string(in.fourier.step)), 
    };
    print_table(msg, {"* " + tag, TypeReflection<TMeta>()});
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