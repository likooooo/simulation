#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <kernels/kernel_loop.hpp>

template <class T, size_t N, size_t Order> struct nested_array {using type = vec<typename nested_array<T, N, Order - 1>::type, N>;};
template <class T, size_t N> struct nested_array<T, N, 1> {using type = vec<T, N>;};
template <class T, size_t N, size_t Order> using tensor = typename nested_array<T, N, Order>::type;
template<size_t DIM, size_t I, class TTensor> constexpr auto& get(TTensor& t, const vec<size_t, DIM>& indexs) 
{
    if constexpr(I != DIM - 1){
        return get<DIM, I + 1>(t.at(indexs.at(I)), indexs);
    } else{
        return t.at(indexs.at(I));
    }
}
template<class T_coef, size_t Order>
struct lagrange_interpolate
{
    constexpr static size_t N = Order + 1;
    //== dx = pos - std::floor(pos)
    constexpr static std::array<T_coef, N> get_coef(const T_coef dx)
    {
        std::array<T_coef, N> coefs{0};
        for(size_t n = 0; n < N; n++){
            coefs.at(n) = get_Ln_numerator(dx, n) / get_Ln_denominator(n);
        }
        return coefs;
    }
    template<size_t tensor_order, size_t ...Is> constexpr static tensor<T_coef, N, tensor_order> __get_coefs_impl(const vec<T_coef, tensor_order>& dxyzn, 
        std::index_sequence<Is...>)
    {
        vec<std::array<T_coef, N>, tensor_order> coefs{0};
        for(size_t i = 0; i < tensor_order; i++) coefs.at(i) = get_coef(dxyzn.at(i));
        tensor<T_coef, N, tensor_order> t{0};
        vec<size_t, tensor_order> shape{0};
        for(size_t i = 0; i <tensor_order; i++) shape.at(i) = N;
        kernels::kernel_loop<size_t, tensor_order>(shape, [&](const vec<size_t, tensor_order>& unused_center, const vec<size_t, tensor_order>& indices){
            T_coef& n = get<tensor_order, 0>(t, indices);
            T_coef prod = 1;
            for(size_t i = 0; i < tensor_order; i++) prod *= coefs.at(i).at(indices.at(i));
            n = prod;
        }); 
        return t;
    }
    template<size_t tensor_order> constexpr static tensor<T_coef, N, tensor_order> get_coefs(const vec<T_coef, tensor_order>& dxyzn)
    {
        return __get_coefs_impl(dxyzn, std::make_index_sequence<tensor_order>{});
    }
    constexpr static std::pair<size_t, std::array<T_coef, N>> interpolate_info(const T_coef x, const size_t unit_count){
        size_t index = std::floor(x);
        index = std::min(index, unit_count -1 - Order);
        return std::make_pair(index, get_coef(x - index));
    } 
    template<class TContainer>constexpr static typename TContainer::value_type eval(const T_coef x, const TContainer& vec){
        auto [index, info] = interpolate_info(x, vec.size());
        typename TContainer::value_type result{0};
        for(size_t i = 0; i < N; i++){
            result += vec.at(index + i) * info.at(i);
        }
        return result;
    } 
    template<class TPixel, size_t tensor_order> constexpr static TPixel eval(const vec<TPixel, tensor_order>& dxyzn, const  tensor<TPixel, N, tensor_order>& on_grid_value)
    {
        auto prod = (on_grid_value * get_coefs<tensor_order>(dxyzn));
        vec<size_t, tensor_order> shape{0};
        for(size_t i = 0; i <tensor_order; i++) shape.at(i) = N;
        TPixel sum{0};
        kernels::kernel_loop<size_t, tensor_order>(shape, 
            [&](const vec<size_t, tensor_order>& unused_center, const vec<size_t, tensor_order>& indices){
                sum += get<tensor_order, 0>(prod, indices);
            }
        ); 
        return sum;
    }
private:
    constexpr static std::array<T_coef, N> get_xn(T_coef dx)
    {
        std::array<T_coef, N> xn{0};
        for(size_t i = 0; i < xn.size(); i++) xn.at(i) = dx - i;
        return xn;
    }
    constexpr static T_coef get_Ln_denominator(const size_t n)
    {
        const std::array<T_coef,N> xn = get_xn(n);
        T_coef Ln = 1;
        for(size_t i = 0; i < N; i++){
	        if(n == i) continue;
            Ln *= xn.at(i);
        }
	    return Ln;
    }
    constexpr static T_coef get_Ln_numerator(const T_coef dx, const size_t n)
    {
        const std::array<T_coef, N> xn = get_xn(dx);
        T_coef Ln = 1;
        for(size_t i = 0; i < N; i++){
	        if(n == i) continue;
            Ln *= xn.at(i);
        }
	    return Ln;
    }
};
template<class T> using linear_interpolate = lagrange_interpolate<T, 1>;
template<class T> using cubic_interpolate = lagrange_interpolate<T, 3>;

template<class T, size_t N>
struct legendre_poly
{
    static_assert(N >=0);
    using value_type = T;
    constexpr static T eval(T x)
    {
        //== {\displaystyle (n+1)P_{n+1}=(2n+1)xP_{n}-nP_{n-1}\}
        if constexpr(N == 0) 
            return T(1);
        else if constexpr(N == 1) 
            return x;
        else 
            return ((2*N - 1) * x * legendre_poly<T, N-1>::eval(x) - (N - 1) * legendre_poly<T, N-2>::eval(x)) / N;
    }
    constexpr static T eval_derived(T x)
    {
        return N * (x * legendre_poly<T, N>::eval(x) - legendre_poly<T, N-1>::eval(x)) / (x * x - 1.0);
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
    constexpr static T max_tolerance = std::numeric_limits<real_t<T>>::epsilon();
    constexpr static T optimize(T initial_guess, T tolerance = max_tolerance, int max_iter = 0)
    {
        T x = initial_guess;
        T delta = 0;
        int i = -1;
        do{
            //== (f(x) - f(root)) / delta = f'(x) . f(root) = 0
            delta = (TPoly::eval(x) - 0) / TPoly::eval_derived(x);
            x -= delta;
            if(0 != max_iter) i++;
        } while (std::abs(delta) > tolerance && i < max_iter);
        return x;    
    }
};

template<class T, size_t N>
struct gauss_integrate
{
    using poly = legendre_poly<T, N>;
    
    constexpr static std::pair<vec<T, N>, vec<T, N>> generate_gauss_weights()
    {
        vec<T, N> nodes{0};
        vec<T, N> weights{0};
        // 对称性：仅计算正根，负根对称
        for(size_t i = 0; i < (N + 1)/2; i++)
        {
            T initial_guess = poly::eval_initial_guess(i);
            T root = newton_raphson<poly>::optimize(initial_guess);
            nodes.at(i) = -root;
            nodes.at(N - 1 - i) = root;
            T deriv = poly::eval_derived(root);
            weights.at(i) = weights.at(N - 1 -i) = T(2) / ((T(1) - root * root) * deriv * deriv);
        }
        return std::make_pair(nodes, weights);
    }
    constexpr static vec<T, N> nodes = generate_gauss_weights().first;
    constexpr static vec<T, N> weights = generate_gauss_weights().second;

    template<class TCallback> static T eval(TCallback&& func, const T a, const T b, 
        const vec<T, N>& nodes = gauss_integrate::nodes, const vec<T, N>& weights = gauss_integrate::weights)
    {
        const T scale = (b - a) / 2.0;
        const T shift = (a + b) / 2.0;

        T sum = 0.0;
        for (size_t i = 0; i < N; ++i)
            sum += weights[i] * func(scale * nodes[i] + shift);
        return sum * scale;
    }
};
