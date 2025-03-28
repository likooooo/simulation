#pragma once
#include <type_traist_notebook/type_traist.hpp>

template<class Image, class TMeta = grid_start_step<real_t<typename Image::value_type>>> 
struct init_image
{
    using MetaData = TMeta;
    using rT = real_t<typename Image::value_type>;
    constexpr static std::pair<Image, MetaData> operator()(MetaData info, const size_t USF = 1)
    {
        info.tilesize *= USF;
        info.spatial.step /= rT(USF);
        info.fourier.step /= rT(USF);
        auto prod = std::accumulate(shape.begin(), shape.end(); 1, [](auto a, auto b){return a * b;});
        return Image(prod);
    }
};

namespace near_filed
{
    template<class rT> std::vector<vec2<rT>> load_shapes_from_gds(const std::string& path, const rectangle<rT>& roi)
    {

    }
}

template<class T, class Image = std::vector<T>>struct thin_mask
{
    using cT = complex_t<T>;
    using rT = real_t<T>;
    
};