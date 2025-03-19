#include <optical/polarization/polarization.hpp>

int main()
{
    using polar = polarization<float>;
    std::vector<std::tuple<float, dual_vec<typename polar::polar_vec>>> basis_result;
    for(float theta : {0_PI, 0.25_PI, 0.5_PI, 0.75_PI, 1_PI, 1.25_PI, 1.5_PI, 2_PI})
        basis_result.push_back(std::make_tuple(theta / M_PI, polar::basis(theta)));
    print_table(basis_result, {"polar angle(PI)", "polar basis"}, 1024);

    using density_result_type = std::tuple<float, dual_vec<float>, matrix<std::complex<float>, 2, 2>>;
    std::vector<density_result_type> density_matrix_rows;
    for(float theta : {0_PI, 0.25_PI, 0.75_PI, 1.25_PI})
    for(dual_vec<float> c : {dual_vec<float>{1, 0}, dual_vec<float>{0, 1}, dual_vec<float>{0.5, 0.5}})
    {
        density_matrix_rows.push_back(density_result_type(theta / M_PI, c,  polar::density_matrix(polar::basis(theta), c).m));
    }
    print_table(density_matrix_rows, {"polar angle(PI)", "components","density matrix"}, 1024);

    size_t xsize = 3, ysize = 3;
    std::vector<float> src(xsize * ysize, 1);src.at(xsize * (ysize / 2)+ xsize/2) = 0;
    std::vector<matrix<std::complex<float>, 2, 2>> SDM(xsize * ysize);
    polar::source_density_matrix(
        SDM.data(), src.data(), 1e-6, {ysize, xsize}, {0.25f, 0.25f}, {1.0f, 0.0f},
        [](float fx, float fy){return get_polar_angle<polar_type::TM>(fx, fy);}
    );
    for(size_t y =0; y < ysize; y++)
    {
        for(size_t x =0; x < xsize; x++)
            std::cout << SDM.at(y * xsize + x) <<  " ";
        std::cout << std::endl;
    }
}