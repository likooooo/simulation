#include <optical/pupil/zernike.hpp>

int main()
{
    zernike_radial_table<float, 28> zernike(5);
    // std::cout << "(size_M, L) = "<< std::make_tuple(zernike.size_M, zernike.L) << std::endl;
    for(size_t m = 0, i = 0; m < zernike.size_M; m++)
    {
        for(size_t l = 0; l <= zernike.L.at(m); l++, i++)
        {
            printf("*Z(%zu, %zu)\n", m, l);
            std::cout << zernike.zk_radials.at(i) << std::endl;
        }
    }
}