#include <optical/polynomials.hpp>

template<class TIntgrate, class T>
void test_gauss_function(T sigma, size_t N, T step)
{
    T norm = T(1) /  (sigma * std::sqrt(2_PI));
    auto gauss = [&](T x)
    {
        return norm * std::exp(-x* x / (2 * sigma *sigma));
    };
    T sum = 0;
    T x = 0;
    for(size_t i = 0; i < N; i++)
    {
        sum += gauss(x);
        sum += gauss(-x);
        x += step;
    }
    sum *= step;
    T intgrate_result = TIntgrate::eval(gauss, -step * N, step * N);
    std::cout << "(sum, intgrate) result = " << std::make_tuple(sum, intgrate_result) << std::endl;
}
int main()
{
    vec<double, 5> golden_result = {NAN, 0.682689492, 0.95449974, 0.9973002, 0.99993666};
    using intgrate = gauss_integrate<double, 5>;

    // TODO : 分析 sigma, N 和误差的关系
    // sigma 越小， 越需要更大的 N

    for(auto sigma : {5.0, 50.0})
    for(auto N : {1024})
    for(auto c : {1, 2, 3, 4})
    {
        double step = c * sigma/ N;
        std::cout << "test-" << c << " with (sigma, resolution, n-sigma)="<< std::make_tuple(sigma, N, c) << std::endl;
        std::cout << "golden result = " << golden_result.at(c) << std::endl;
        test_gauss_function<intgrate, double>(sigma, N, step);
        std::cout << std::endl;
    }
    // std::cout << newton_raphson<legendre_poly<double, 5>>::max_tolerance << std::endl;
    // std::cout << std::setprecision(10) << std::setiosflags(std::ios::fixed);
    // std::cout << intgrate::nodes << std::endl;
    // std::cout << intgrate::weights << std::endl;
}