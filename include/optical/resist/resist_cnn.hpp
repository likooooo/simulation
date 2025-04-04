#pragma once
#include "resist_common.hpp"
template<class T>
struct resist_cnn
{
    using XN = std::vector<T>;
    using MatrixP = std::vector<XN>;
    using VectorQ = std::vector<T>;
};

template<class T>
struct resist_least_squares
{

};