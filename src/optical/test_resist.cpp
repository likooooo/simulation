#include <optical/resist/resist_common.hpp>

void get_on_measured_cd_positive_slice_test(std::vector<double> cutline, double threshold, vec2<double> result, bool found)
{
    auto [cd, is_found] = get_on_measured_cd_positive_slice(cutline, threshold, cutline_data());
    assert(found == is_found);
    assert(abs(result[0] - cd[0]) < 1e-4); 
    assert(abs(result[1] - cd[1]) < 1e-4); 
}

void actual_cutline_in_nm_from_thin_mask_test(std::vector<double> cutline, double threshold, double result)
{
    cutline_data data;
    auto [cd, found] = get_on_measured_cd_positive_slice<double>(cutline, threshold, data, 0, 5);
    assert(found);
    std::cout << "actually cd is " << result << "(nm) cd found is " << (cd[1] - cd[0]) << "(nm) " << cd << std::endl;

    //== mask image is 5nm/pixel
    assert(std::abs(result - (cd[1] - cd[0])) < 1); 
}

void edge_pixelization_test()
{
    actual_cutline_in_nm_from_thin_mask_test(
        //== L48P110 integral_x only  
        {1.22125e-15,1.22125e-15,1.22125e-15,1.22125e-15,1.22125e-15,1.22125e-15,0.2,0.45,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.45,0.2,1.33227e-15,1.33227e-15,1.33227e-15,1.33227e-15,1.33227e-15,1.33227e-15,0},
        0.25, 48.0
    );
    actual_cutline_in_nm_from_thin_mask_test(
        //== L48P110 integral_x + integral_y
        {1.22125e-15,1.22125e-15,1.22125e-15,1.22125e-15,1.22125e-15,1.22125e-15,0.285,0.78125,0.99625,1,1,1,1,1,1,0.94375,0.51875,0.075,1.33227e-15,1.33227e-15,1.33227e-15,1.33227e-15,1.33227e-15,0},
        0.5, 48.0
    );
    actual_cutline_in_nm_from_thin_mask_test(
        //== sIL62P116 integral_x + integral_y 
        {1,1,1,1,1,0.9325,0.445,0.06,2.22045e-16,2.22045e-16,2.22045e-16,2.22045e-16,2.22045e-16,2.22045e-16,2.22045e-16,2.22045e-16,2.22045e-16,2.22045e-16,0.342,0.8235,0.997,1,1,1,1},
        0.5, 62.0
    );
}

int main()
{
    get_on_measured_cd_positive_slice_test({0, 1, 0}, 0.5, {0.5, 1.5}, true);
    get_on_measured_cd_positive_slice_test({0.5, 1, 1, 0.5}, 0.5, {0, 3}, true);
    edge_pixelization_test();
}