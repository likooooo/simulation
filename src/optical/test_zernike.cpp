#include <optical/pupil/zernike.hpp>
#include <py_helper.hpp>
#include <kernels/kernel_loop.hpp>
#include <kernels/bloch.hpp>

void display_zernike_image(int select_m, int select_l, size_t size, int flag);
void display_aberration_pupil_image(int m, int l);
int main(int argc, char *argv[])
{
    py_engine::init();
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
    // display_aberration_pupil_image(std::max(select_m, 0), std::max(select_l, 0));
    display_zernike_image(select_m, select_l, 100, 1);
}
void display_zernike_image(int select_m, int select_l, size_t size, int flag)
{
    using zk_table = zernike_radial_table<float, 10>;
    static zk_table zernike(50);
    static bool print_enable = true;
    if(print_enable){
        print_enable = false;
        std::cout << TypeReflection<zk_table>() << " (max_M, max_L) = "<< std::make_tuple(zernike.size_M - 1, zernike.L) << std::endl
                  << "  |-  Zernike(0, l)  = layler series expansion" << std::endl
                  << "  |-  Zernike others = legendre polynomial expansion" << std::endl << std::endl;
    }

    for(const auto& [info, im] : zernike.init_zernike_image(size)){
        auto [m, l] = info;
        if(!((l == select_l || select_l == -1) && (m == select_m || select_m == -1))) continue;
        std::cout << "Z(m,l)=" << info << std::endl;
        std::vector<float> temp(im.size());
        std::transform(im.begin(), im.end(), temp.begin(),[flag](complex_t<float> c){return flag > 0 ? c.real() : c.imag();});
        std::cout << "    " << zernike.zk_radials.at(zernike.index(m, l)) << std::endl;
        imshow(temp, {size, size});
    }
}
void display_aberration_pupil_image(int m, int l)
{
    using cT =complex_t<float>;
    using zk_table = zernike_radial_table<float, 10>;
    zk_table zernike(50);
    size_t size = 100;
    float step = 2.0/(size - 1);
    float freq = 1;
    float dz = 0.1;
    cT norm = cT(2_PI * dz * freq);
    auto im = zernike.gen_aberration_pupil_image({size, size}, {step, step}, 1, {
        // std::tuple<size_t, size_t, cT>(m, l, cT(1, 1))
        // std::tuple<size_t, size_t, cT>(0, 0, norm),
        std::tuple<size_t, size_t, cT>(0, 2, norm / (-2.0f)),
        // std::tuple<size_t, size_t, cT>(0, 4, norm / (-8.0f)),
        // std::tuple<size_t, size_t, cT>(0, 6, norm / (-16.0f)),
        // std::tuple<size_t, size_t, cT>(0, 8, norm * 5.0f / (-128.0f))
    });
    // imshow(im, {size, size});
    // std::transform(im.begin(), im.end(), temp.begin(),[](float c){return std::exp(cT(0, c)).real();});    
    {
        // auto [it_min, it_max] = std::minmax_element(im.begin(), im.end());
        // std::cout << (*it_max / * it_min) << std::endl;
    }
    imshow(im, {size, size});

    std::vector<float> temp(im.size());
    std::vector<complex_t<float>> vec(size * size);
    kernels::free_propagation<float, 2>(vec.data(), {size, size}, {step, step}, freq, dz);
    std::transform(vec.begin(), vec.end(), temp.begin(),[](complex_t<float> c){return std::arg(c);});
    {
        // auto [it_min, it_max] = std::minmax_element(temp.begin(), temp.end());
        // std::cout << (*it_max / * it_min) << std::endl;
    }
    imshow(temp, {size, size});
}

BOOST_PYTHON_MODULE(lib_test_zernike) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("display_zernike_image", &display_zernike_image);
    py::def("display_aberration_pupil_image", &display_aberration_pupil_image);
}