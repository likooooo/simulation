
#pragma once
#include <mekil/cpu_backend.hpp>
#ifdef GPU_BACKEND_ENABLE
#   include <cpp_cuda/gpu_backend.hpp>
#endif

template<class T> uca::backend<T>& get_math_backend_impl()
{
    static uca::backend<T>& ref = uca::cpu<T>::ref();
    return ref;
}
template<> uca::backend<float>& get_math_backend<float>(){return get_math_backend_impl<float>();}
template<> uca::backend<double>& get_math_backend<double>(){return get_math_backend_impl<double>();}
template<> uca::backend<complex_t<float>>& get_math_backend<complex_t<float>>(){return get_math_backend_impl<complex_t<float>>();}
template<> uca::backend<complex_t<double>>& get_math_backend<complex_t<double>>(){return get_math_backend_impl<complex_t<double>>();}