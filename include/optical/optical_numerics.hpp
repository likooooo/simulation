#pragma once
#include <type_traist_notebook/type_traist.hpp>

template<class T>
struct grid_start_step
{
    struct StartStep{vec2<T> start, step;};
    StartStep spatial, fourier;
    vec2<size_t> tilesize;
    void print() const
    {
        std::vector<std::tuple<std::string, vec2<T>>> msg{
            std::make_tuple(std::string("spatial start :"), spatial.start),
            std::make_tuple(std::string("spatial step  :"), spatial.step),
            std::make_tuple(std::string("fourier start :"), fourier.start),
            std::make_tuple(std::string("fourier step  :"), fourier.step), 
        };
        print_table(msg, {"* optical numerics with tilesize=" + to_string(tilesize), ""});
    }
};

constexpr size_t doubling_bandwith_for_squaring_signal_system = 2;
template<class T> inline grid_start_step<T> optical_numerics(const vec2<vec2<T>>& roi, const vec2<T>& ambit, const vec2<size_t>& tilesize, T maxNA, T lambda)
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
template<class T> inline grid_start_step<T> bloch_optical_numerics(vec2<vec2<T>> spatial_domain, vec2<size_t> tilesize, T maxNA, T lambda)
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