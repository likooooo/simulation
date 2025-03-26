
#include <optical/optical_numerics.hpp>

template<class T>void sim_config_test()
{
    T sigma = 0;
    T maxNA = 0.5 + sigma;
    T lambda = 13.5;
    vec2<T> from{0, -200};
    vec2<T> to{400, 200};
    vec2<T> ambit{0, 0};
    vec2<size_t> tilesize{160, 160};
    optical_numerics({from, to}, ambit, tilesize, maxNA, lambda).print();
    bloch_optical_numerics({from, to}, tilesize, maxNA, lambda).print();
}
int main()
{
    sim_config_test<float>();
    // sim_config_test<double>();
}
