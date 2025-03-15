#include <optical/pupil/zernike.hpp>
#include <py_helper.hpp>

int main()
{
    using zk_table = zernike_radial_table<float, 4>;
    zk_table zernike(5);
    std::cout << TypeReflection<zk_table>() << " (max_M, max_L) = "<< std::make_tuple(zernike.size_M - 1, zernike.L) << std::endl
              << "  |-  Zernike(0, l)  = layler series expansion" << std::endl
              << "  |-  Zernike others = legendre polynomial expansion" << std::endl << std::endl;
    for(size_t m = 0, i = 0; m < zernike.size_M; m++)
    {
        for(size_t l = 0; l <= zernike.L.at(m); l++, i++)
        {
            printf("*Z(%zu, %zu)\n", m, l);
            std::cout << zernike.zk_radials.at(i) << std::endl;
        }
    }
    size_t xsize = 10, ysize = 10;
    std::vector<std::vector<float>> images;
    images.reserve(zernike.size_M * zernike.size_M);
    float step = float(2 * zernike.zk_radials.front().size()) / xsize;
    complex_t<float> center{xsize/2.0f, ysize/2.0f};
    for(size_t m = 0, i = 0; m < zernike.size_M; m++)
    {
        for(size_t l = 0; l <= zernike.L.at(m); l++, i++)
        {
            if(m){
                images.emplace_back(xsize * ysize);
                images.emplace_back(xsize * ysize);
            }
            else{
                images.emplace_back(xsize * ysize);
            }
            for(size_t y = 0; y < ysize; y++)
            for(size_t x = 0; x < xsize; x++)
            {
                float r = std::abs(complex_t<float>(x, y) - center) * step;
                if(m){
                    float rVal = cubic_interpolate<float>::eval(r, zernike.zk_radials.at(i));
                    float theta = zk_table::theta(m, x - center.real(), y - center.imag());
                    images.back().at(y * xsize + x) = rVal * std::cos(theta);
                    (images.rbegin() + 1)->at(y * xsize + x) = rVal * std::cos(theta);
                }else{
                    images.back().at(y * xsize + x) = cubic_interpolate<float>::eval(r, zernike.zk_radials.at(i));
                }
            }
        }
    }
    std::cout << images.size() << std::endl;
    py_loader::init();
    py_plot::get_default_visualizer_dir() = "/usr/local/bin";
    for(const auto& im : images)
    {
        imshow(im, {ysize, xsize});
    }
}