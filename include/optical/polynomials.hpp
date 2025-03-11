#pragma once
#include <type_traist_notebook/type_traist.hpp>
// template<class T>
// struct  


// template <class TYPE> OLreal CSvpExpanderT<TYPE>::Denormalize(const OLreal x)
// {
//     return varMin_ + (varMax_ - varMin_) * (x + 1) / 2;
// }

// template <class TYPE> OLreal CSvpExpanderT<TYPE>::NormalizeIt(const OLreal X)
// {
//     if (varMax_ > varMin_)
//         return 2 * (X - varMin_) / (varMax_ - varMin_) - 1;
//     else
//         return -1.0;
// }


template<class T, size_t n>
struct legendre_poly
{
    using value_type = T;
    constexpr static T eval(T x)
    {
        //== {\displaystyle (n+1)P_{n+1}=(2n+1)xP_{n}-nP_{n-1}\}
        if constexpr(N == 0) return T(1);
        if constexpr(N == 1) return x;
        return ((2*N - 1) * x * legendre_poly<N-1>::eval>(x) - (N - 1) * legendre_poly<N-2>::eval(x)) / N;
    }
    constexpr static T eval_derived(T x)
    {
        return N * (x * legendre_poly<N>::eval(x) - legendre_poly<N-1>::eval(x)) / (x * x - 1.0);
    }
    constexpr static T eval_initial_guess(size_t root_index){
        return std::cos(M_PI * (root_index + 0.75) / (N + 0.5));
    }
};
template<class TPoly>
struct newton_raphson
{
    // https://en.wikipedia.org/wiki/Newton's_method
    using T = typename TPoly::value_type;
    T optimize(T initial_guess, T tolerance = 1e-12)
    {
        T x = initial_guess;
        T delta;
        do{
            //== (f(x) - f(root)) / delta = f'(x) . f(root) = 0
            delta = (TPoly::eval(x) - 0) / TPoly::eval_derived(x);
            x -= delta;
        } while (std::abs(delta) > tolerance);
        return x;    
    }
};

template<class T, size_t N>
struct gauss_integrate
{
    using poly = legendre_poly<T, 10>;
    constexpr vec<T, N> nodes = generate_gauss_weights().first;
    constexpr vec<T, N> weights = generate_gauss_weights().second;
private:
    constexpr static std::pair<vec<T, N>, vec<T, N>> generate_gauss_weights()
    {
        vec<T, N> nodes;
        vec<T, N> weights;
        // 对称性：仅计算正根，负根对称
        for(size_t i = 0; i < (N + 1)/2; i++)
        {
            T initial_guess = poly::eval_initial_guess(i);
            T root = newton_raphson<poly>::optimize(initial_guess);
            nodes.at(i) = -root;
            nodes.at(N - 1 - i) = root;
            T deriv = poly::eval_derived(root);
            weights.at(i) = weights.at(N - 1 -i) = T(2) / ((T(1) - root  root) * deriv * deriv);
        }
        return std::make_pair(nodes, weights);
    }
};


// // 动态生成的高斯积分
// double DynamicGaussIntegrate(
//     const std::function<double(double)>& func,
//     double a,
//     double b,
//     int n = 5
// ) {
//     std::vector<double> nodes, weights;
//     GenerateGaussWeights(n, nodes, weights);
    
//     double sum = 0.0;
//     double scale = (b - a) / 2.0;
//     double shift = (a + b) / 2.0;
    
//     for (int i = 0; i < n; ++i) {
//         double x_mapped = scale * nodes[i] + shift;
//         sum += weights[i] * func(x_mapped);
//     }
//     return sum * scale;
// }