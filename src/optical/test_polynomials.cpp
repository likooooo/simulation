#include <optical/polynomials.hpp>
#include <py_helper.hpp>
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
void test()
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
template<class TInterpolate, size_t DIM> void interoplate_test(const matrix<float, DIM, DIM>& on_grid_value)
{
    constexpr matrix<float, DIM, DIM> compile_time_coef = TInterpolate:: template get_coefs<2>({0.5, 0.5}); 
    std::cout << compile_time_coef << std::endl;
    std::cout << TInterpolate:: template get_coefs<2>({0, 0}) << std::endl;
    std::cout << TInterpolate:: template get_coefs<2>({0, 1}) << std::endl;
    std::cout << TInterpolate:: template get_coefs<2>({1, 0}) << std::endl;
    std::cout << TInterpolate:: template get_coefs<2>({1, 1}) << std::endl;

    size_t N = 100;
    float step = float(DIM - 1)/(N -1);
    std::vector<float> values(N * N);
    for(size_t y = 0; y < N; y++)
    for(size_t x = 0; x < N; x++){
        values.at(y * N + x) = TInterpolate:: template eval<float, 2>({x *step, y * step}, on_grid_value);
    }
    imshow(values, {N, N});

}
int main()
{
    py_engine::init();
    // test();

    // //== max value should be sqrt(2)
    // interoplate_test<linear_interpolate<float>, 2>({
    //     0.0f, 1.0f,
    //     1.0f, std::sqrt(2.0f)
    // });
    auto f_pow3 = [](float x){
        return std::pow<float>(x, 3);
    };
    // //== max value should be pow(3 * sqrt(2), 3)
    // interoplate_test<cubic_interpolate<float>, 4>({
    //     f_pow3(0), f_pow3(1),                f_pow3(2),                f_pow3(3),
    //     f_pow3(1), f_pow3(std::sqrt(1 + 1)), f_pow3(std::sqrt(4 + 1)), f_pow3(std::sqrt(9 + 1)), 
    //     f_pow3(2), f_pow3(std::sqrt(1 + 4)), f_pow3(std::sqrt(4 + 4)), f_pow3(std::sqrt(9 + 4)),
    //     f_pow3(3), f_pow3(std::sqrt(1 + 9)), f_pow3(std::sqrt(4 + 9)), f_pow3(std::sqrt(9 + 9))
    // });

    size_t N = 100;
    auto cal_grid = [](size_t N, float step){
        std::vector<float> on_grid(N * N);
        for(size_t y = 0; y < N; y++)
        for(size_t x = 0; x < N; x++)
            on_grid.at(y * N + x) = std::pow(std::hypot<float, float>(float(x) - N/2, float(y) - N/2) * step, 3);
        return on_grid;
    };
    auto on_grid = cal_grid(N, 1.0f/(N - 1));
    const size_t scalar = 3;
    // imshow(cal_grid(N, 1.0f/(N - 1)), {N, N});
    N *= scalar;
    auto golden = cal_grid(N, 1.0f/(N - 1));
    // imshow(golden, {N, N});
    
    std::vector<float> test(N*N);
    float step = float(N / scalar - 1) / (N - 1);
    for(size_t y = 0; y < N; y++)
    for(size_t x = 0; x < N; x++)
        test.at(y * N + x) =  cubic_interpolate<float>:: template eval<2>(vec2<float>{float(y), float(x)} * step, on_grid, {N / scalar, N/ scalar});
    auto error = golden - test;
    float rms = 0;
    for(auto& r : error) {
        r = std::abs(r);
        rms += r * r;
    }
    std::cout << "RMS=" << std::sqrt(rms /= error.size()) << std::endl;
    imshow(error, {N, N});
}