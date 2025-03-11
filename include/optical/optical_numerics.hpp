#pragma once
#include <type_traist_notebook/type_traist.hpp>

template<class T>
struct GridStartStep
{
    struct StartStep{vec2<T> start, step;};
    StartStep spatial, fourier;
    void print() const
    {
        std::vector<std::tuple<std::string, vec2<T>>> msg{
            std::make_tuple(std::string("spatial start :"), spatial.start),
            std::make_tuple(std::string("spatial step  :"), spatial.step),
            std::make_tuple(std::string("fourier start :"), fourier.start),
            std::make_tuple(std::string("fourier step  :"), fourier.step), 
        };
        print_table(msg, {"* optical numerics", ""});
    }
};

template<class T> inline GridStartStep<T> optical_numerics(vec2<vec2<T>> roi, vec2<T> ambit, vec2<size_t> tilesize, T maxNA, T lambda)
{
    GridStartStep<T> grid;

    //== 1. spatial step
    constexpr T doubling_bandwith_for_squaring_signal_system = T(2);
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