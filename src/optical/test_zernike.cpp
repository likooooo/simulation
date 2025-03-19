#include <optical/pupil/zernike.hpp>
#include <py_helper.hpp>
#include <kernels/kernel_loop.hpp>
int main(int argc, char *argv[])
{
    using zk_table = zernike_radial_table<float, 10>;
    zk_table zernike(50);
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
    size_t xsize = 100, ysize = 100;
    std::vector<std::pair<vec2<size_t>, std::vector<complex_t<float>>>> images;
    images.reserve(zernike.size_M * zernike.size_M);
    float step = float(2 * zernike.zk_radials.front().size()) / xsize;
    complex_t<float> center{xsize/2.0f, ysize/2.0f};
    for(size_t m = 0, i = 0; m < zernike.size_M; m++)
    {
        for(size_t l = 0; l <= zernike.L.at(m); l++, i++)
        {
            images.push_back(std::make_pair(vec2<size_t>{m,l}, std::vector<complex_t<float>>(xsize * ysize)));
            complex_t<float>* pData = images.back().second.data();
            kernels::center_zero_loop_square_r<float, 2>({ysize, xsize}, {step, step}, [&](const std::array<float, 2> pos, float r){
                r = std::sqrt(r);
                if(r < zernike.zk_radials.at(i).size()){
                    float rVal = cubic_interpolate<float>::eval(r, zernike.zk_radials.at(i));
                    float theta = zk_table::theta(m, pos[0], pos[1]);
                    *pData = std::exp(complex_t<float>(0, theta)) * rVal;
                }
                pData++;
            });
        }
    }
    std::cout << images.size() << std::endl;
    py_loader::init();
    py_plot::get_default_visualizer_dir() = "/usr/local/bin";
    
    int select_m = -1, select_l = -1;
    for(int i = 1; i < argc; i++){
        if('m' == argv[i][0]){
            i++;
            select_m =  argv[i][0] - '0';
        }
        if('l' == argv[i][0]){
            i++;
            select_l =  argv[i][0] - '0';
        }
    }
    for(const auto& [info, im] : images)
    {
        auto [m, l] = info;
        if(!((l == select_l || select_l == -1) && (m == select_m || select_m == -1))) continue;
        std::cout << "Z(m,l)=" << info << std::endl;
        std::vector<float> temp(im.size());
        std::transform(im.begin(), im.end(), temp.begin(),[](complex_t<float> c){return c.real();});
        imshow(temp, {ysize, xsize});
    }
}