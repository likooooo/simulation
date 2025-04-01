#include <type_traist_notebook/type_traist.hpp>

namespace uca
{    
    template<class T, class vec = std::vector<T>>struct backend_todo
    {
        void (*VtAdd)(vec& y, const vec& x);
        int (*integral_x)(vec2<size_t>, vec& im);
    };
    
    template<class T>struct backend
    {
        // y += x
        void (*VtAdd)(const int n, const T *x, T *y);
        int (*integral_x)(vec2<size_t>, T*);
    };
    
    template<class T> struct gpu_backend : backend<T>
    {
        using value_type = T;
        using alloc_type = std::allocator<T>;
        gpu_backend();
        static gpu_backend& ref()
        {
            static gpu_backend gpu;
            return gpu;
        }
    };
    template<class T>struct cpu_backend : backend<T>
    {
        using value_type = T;
        using alloc_type = std::allocator<T>;
        cpu_backend();
        static cpu_backend& ref()
        {
            static cpu_backend cpu;
            return cpu;
        }
    };
    template<class T> using gpu = gpu_backend<T>;
    template<class T> using cpu = cpu_backend<T>;
}