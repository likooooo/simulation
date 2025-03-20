#pragma once
#include <type_traist_notebook/type_traist.hpp>

template<class T>
struct snells_law
{
    using cT = complex_t<T>;
    constexpr static cT refraction_angle(cT n1, cT n2, cT input = 0_PI){ std::asin(n1 * std::sin(input) / n2); }
    constexpr static std::vector<cT> refraction_angle(const std::vector<cT> nk, cT input = 0_PI)
    { 
        std::vector<cT> results(nk.size());results.at(0) = input;
        for(size_t i = 0; i < nk.size() - 1; i++)
            results.at(i + 1) = refraction_angle(nk.at(i), nk.at(1 + 1), results.at(i));
        return results;
    }
};
