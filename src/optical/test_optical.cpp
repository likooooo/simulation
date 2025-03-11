
#include <optical/optical_numerics.hpp>
int main()
{
    optical_numerics<float>(
        vec2<vec2<float>>{vec2<float>{0.f, 0.f}, vec2<float>{100.f, 100.f}}, 
        vec2<float>{1.f, 1.f}, vec2<size_t>{10, 10}, 0.8, 13.5
    ).print();
}
