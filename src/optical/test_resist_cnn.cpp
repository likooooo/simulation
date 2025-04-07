#include <optical/resist/resist_cnn.hpp>
#include <py_helper.hpp>


int main()
{
    using resist = resist_blackbox<float>;
    py_engine::init();
    // TODO : path to image
    size_t xsize = 128, ysize = 128;
    vec2<size_t> shape{xsize, ysize};
    std::vector<float> image(xsize * ysize);
    image.reserve(image.size() + 2 * ysize);
    image.at(xsize * (ysize / 2) + xsize / 2) = 10;
    imshow(image, convert_to<std::vector<size_t>>(shape));
    auto [linear, quadratic] = resist::gauss_laguerre(image, shape, 5.0/3, 2, 2);
    for(const auto& l :linear)
        imshow(l, convert_to<std::vector<size_t>>(shape));
    for(const auto& q :quadratic)
        imshow(q, convert_to<std::vector<size_t>>(shape));
}