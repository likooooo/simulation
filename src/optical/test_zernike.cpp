#include <optical/pupil/zernike.hpp>
#include <py_helper.hpp>
#include <kernels/kernel_loop.hpp>

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
    auto im = zernike.gen_aberration_pupil_image({size, size}, {step, step}, {
        std::tuple<size_t, size_t, cT>(m, l, cT(1.0 / 2_PI / std::sqrt(zk_table::cal_norm(m, l))))
    });
    std::vector<float> temp(im.size());
    std::transform(im.begin(), im.end(), temp.begin(),[](complex_t<float> c){return c.real();});
    imshow(temp, {size, size});
}

BOOST_PYTHON_MODULE(lib_test_zernike) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("display_zernike_image", &display_zernike_image);
    py::def("display_aberration_pupil_image", &display_aberration_pupil_image);
}