#include <optical/pupil/pupil.hpp>
#include <optical/pupil/zernike.hpp>
#include <kernels/kernel_loop.hpp>
#include <assert.h>
#include <random>
#include <py_helper.hpp>
#include <fft/spectrum_analysis.hpp>
template<class T> struct test
{
    using pupil_print = debug_print<pupil_radial<T>>;
    using cT = complex_t<T>;
    using rT = real_t<T>;
    using meterial = typename film_stack_solver<T>::meterial;

    constexpr static auto vacuum_nk = pupil_radial<T>::vacuum_nk;
    T NA = 0.8, lambda = 13.5;
    pycallback_update_frame display = py_plot::create_callback_simulation_fram_done();
    test()
    {
        pupil_print::print_to() << std::fixed << std::setprecision(5) << std::endl;
        pupil_print::verbose() = false;
    }
    void get_anamorphic_pupil()
    {
        auto print_radials = [&](const std::string& str, const std::vector<matrix2x3<complex_t<T>>>& pupil_radials)
        {
            const auto&[TE_x, TE_y, TE_z, TM_x, TM_y, TM_z] = decompose_from<matrix2x3<complex_t<T>>, complex_t<T>, complex_t<T>, complex_t<T>, complex_t<T>, complex_t<T>, complex_t<T>>(pupil_radials);
            // auto abs = [](decltype(TE_x) n){return n;};
            pupil_print(std::vector<std::tuple<const char*, decltype(abs(TE_x))>>
                {
                    std::make_tuple("TE_x", abs(TE_x)), std::make_tuple("TE_y", abs(TE_y)), std::make_tuple("TE_z", abs(TE_z)), 
                    std::make_tuple("TM_x", abs(TM_x)), std::make_tuple("TM_y", abs(TM_y)), std::make_tuple("TM_z", abs(TM_z)),
                },
                {
                    str, "NA = " + to_string(NA)
                }, 100
            );
            pupil_print("\n\n");
        };
        
        std::vector<matrix2x3<complex_t<T>>> pupil_radials = pupil_radial<T>::init_anamorphic_pupil_radial(5, NA, vacuum_nk);
        print_radials("* init pupil_radial", pupil_radials);
        pupil_radial<T>::apply_defocus_to_pupil_radial(pupil_radials, vacuum_nk, lambda / 4, NA, lambda);
        print_radials("* defocus", pupil_radials);

        pupil_radial<T>::apply_obliquity_factor_to_pupil_radial(pupil_radials, NA, 
            T(1.0 / NA * (1.0 + 1.0 / pupil_radials.size())),   // 增大放大倍率，使得部分高频信息丢失
            vacuum_nk.real()                                    // 增大背景的折射率，增加衰减项
        );
        print_radials("* obliquity", pupil_radials);

        std::vector<meterial> meterials = {
            meterial{1.0, 0},                                  // background
            meterial{1.0, lambda*10},                          // top
            meterial{complex_t<T>(0.926, 0.044), lambda * 10}, // absorber
            meterial{1.0, lambda*10},                          // substrate
            meterial{1.0, 0},                                  // substrate
        };
        T depth = std::accumulate(meterials.begin(), meterials.end() - 1, T(0), [](T sum, meterial t){return sum + t.depth;});
        size_t N = 1024;
        T step = depth / N;
        std::vector<rT> wave(N);
        for(size_t i = 0; i < N; i++)
        {
            pupil_radials = pupil_radial<T>::init_anamorphic_pupil_radial(5, NA, vacuum_nk);
            pupil_radial<T>::apply_film_stack_to_pupil_radial(pupil_radials, NA, lambda, vacuum_nk, 
                step * i, meterials);
            wave.at(i) = pupil_radials.at(0).at(0)[1].real();
            // wave.at(i) = std::abs(pupil_radials.at(0).at(0)[1]);
            // print_radials(to_string(i), pupil_radials);
        }
        meterial::print(meterials) << std::endl;
        printf("After passing through the absorption layer, the energy attenuation is observed...\n");
        printf("\n");
        plot_curves(std::vector<std::vector<rT>>{wave}, {0}, {float(step)}, {"wave"}, {"b--"});  
    }

    void test_single_layer() {
        film_stack_solver<T> solver(lambda, 1);

        uniform_random<T> rSubstrate(0.1, 100);
        uniform_random<T> rDepth(0, 1e3);
        T top = rSubstrate();
        std::vector<meterial> meterials = {
            {cT(top), 0},
            {cT(top), rDepth()},    // 空气层
            {cT(rSubstrate()), 1e9} // 基底
        };
        
        meterial::print(meterials) << std::endl;
        
        //== 垂直入射下， 多层膜和菲涅尔公式结果一致
        T kr = 0; 
        auto te_matrix = solver.solveTE(meterials, kr);
        auto tm_matrix = solver.solveTM(meterials, kr);

        //== 菲涅尔反射
        T expected_r_te = std::abs((meterials[2].nk - meterials[1].nk) / (meterials[2].nk + meterials[1].nk));
        auto eps = std::sqrt(vector_norm(te_matrix[1][1])) - expected_r_te;
        assert(std::abs(eps) < 1e-6);
        printf("    test success. reflection is same behavior with Fresnel-Formula\n");

        //== 菲涅尔透射
        T expected_t_te = std::abs(T(2) * meterials[1].nk / (meterials[2].nk + meterials[1].nk));
        eps = expected_t_te - std::sqrt(vector_norm(te_matrix[2][0]));
        assert(std::abs(eps) < 1e-6);
        printf("    test success. transmition is same behavior with Fresnel-Formula\n");

        //== 我无意中发现的等式
        expected_t_te = std::abs(T(2) * meterials[2].nk / (meterials[2].nk + meterials[1].nk));
        eps = 1 - (std::sqrt(vector_norm(te_matrix[2][0])) + expected_t_te) / 2;
        assert(std::abs(eps) < 1e-6);
        printf("    test success. transmition is same behavior with UNKNOW-Formula\n");

        //== 能量守恒
        auto t = [](const std::string& s, const std::vector<matrix2x3<cT>>& m)
        {
            auto input_energe = std::sqrt(vector_norm(m[1][0]));
            auto transmition_of_substrate_energe = std::sqrt(vector_norm(m[2][0]));
            auto reflection_of_substrate_energe = std::sqrt(vector_norm(m[1][1]));
            // std::cout << std::make_tuple(input_energe - reflection_of_substrate_energe - transmition_of_substrate_energe, 
            //     input_energe, transmition_of_substrate_energe, reflection_of_substrate_energe) << std::endl;
            if(input_energe > transmition_of_substrate_energe){
                assert(is_almost_equal<T>(input_energe, transmition_of_substrate_energe + reflection_of_substrate_energe, 1e-5));
            }
            else{
                assert(is_almost_equal<T>(input_energe, transmition_of_substrate_energe - reflection_of_substrate_energe, 1e-5));
            }
            assert(is_almost_equal<T>(input_energe, 1, 1e-6));
            printf("    test success. %s conservation of energy.\n", s.c_str());
        };
        t("TE", te_matrix);
        t("TM", tm_matrix);
    }
};
void test_single_layer()
{
    test<float> t;
    t.get_anamorphic_pupil();
    for(size_t i = 0; i < 100; i++){
        t.test_single_layer();
        printf("test-%zu end\n\n", i);
    }
}

template<class T> std::vector<matrix2x3<complex_t<T>>> gen_pupil_array(T lambda, T defocus, T NA, complex_t<T> nk = complex_t<T>(1), T freq_step = 0.001, T reduction_ratio = 5,
    const std::vector<std::tuple<size_t, size_t, complex_t<T>>>& poly_coefs= {
            // std::tuple<size_t, size_t, complex_t<T>>(0, 0, complex_t<T>(2_PI)), 
            // std::tuple<size_t, size_t, complex_t<T>>(2, 0, complex_t<T>(2_PI, 0))
            std::tuple<size_t, size_t, complex_t<T>>(0, 0, complex_t<T>(1)) 
    })
{
    using rT = real_t<T>;
    using cT = complex_t<T>;
    using print_type = std::vector<std::tuple<vec3<cT>, vec3<cT>>>;
    using zk_table = zernike_radial_table<T, 10>;
    assert(0 <= NA && NA <= 1);
    //== [-2NA, 2NA]
    const size_t N = 2 * (std::ceil(2*NA/freq_step) + 1);
    const vec2<size_t> shape{N, N};
    const vec2<T> step{freq_step, freq_step};
    const T pupil_na = get_pupil_na(NA);
    std::vector<matrix2x3<cT>> pupil_radials = pupil_radial<T>::init_anamorphic_pupil_radial(N, pupil_na, nk);
    if(false){
        //== test code
        std::vector<matrix2x3<complex_t<float>>> test = pupil_radial<float>::init_anamorphic_pupil_radial(1000, pupil_na, nk);
        std::vector<float>& temp = reinterpret_cast<std::vector<float>&>(test);
        temp -= std::get<0>(load_image<float>("pupil_init.bin", {12, 1000}));
        imshow(temp, {12, 1000});
        auto it = std::max_element(temp.begin(), temp.end());
        debug_unclassified::verbose() = true;
        debug_unclassified::out("pupil init max-error=", *it, " at ", std::distance(temp.begin(), it));
        exit(0);
    }
    const std::vector<T> zernike_image = zk_table(pupil_radials.size()).gen_pupil_image_with_zernike(shape, step, pupil_na, poly_coefs);
    const T* pZernike = zernike_image.data();
    // imshow(zernike_image, convert_to<std::vector<size_t>>(shape));
    std::vector<matrix2x3<cT>> pupil_image(zernike_image.size());
    matrix2x3<cT>* p = pupil_image.data();
    kernels::center_zero_loop_square_r<rT, 2>(shape, step, [&](const vec2<rT>& fxy, rT r){
        r = std::sqrt(r);
        matrix2x3<cT> val{0};
        if(r <= pupil_na) {
            r /= pupil_na;
            val = cubic_interpolate<rT>::eval(r * pupil_radials.size(), pupil_radials);
            val *= std::exp(cT(0, *pZernike));
        } 
        *p = val;
        pZernike++;
        p++;
    });
    // pupil<T> p(NA, freq_step, poly_coefs);
    // auto pupil_image = p.vector_pupil;
    //== pupil with zernike ON
    // pupil_golden_check<rT>(pupil_image, shape, {"pupil_zernike_TE_y.bin", "pupil_zernike_TM_x.bin", "pupil_zernike_TM_z.bin"});
    
    //== pupil with zernike OFF
    // pupil_golden_check({"pupil_TE_y.bin", "pupil_TM_x.bin", "pupil_TM_z.bin"});
    
    // apply_anamorphic_effect<rT>(pupil_image, shape, step, 9_PI/180, 45_PI/180, 0, 0, NA, lambda, cT(1), nk, reduction_ratio, reduction_ratio);
    
    //== pupil with zernike OFF, projection ON
    // pupil_golden_check({"pupil_sp_projection_TE_y.bin", "pupil_sp_projection_TM_x.bin", "pupil_sp_projection_TM_z.bin"});
    
    //== pupil with obliquity factor
    // pupil_golden_check<rT>(pupil_image, shape, {"pupil_final_TE_y.bin", "pupil_final_TM_x.bin", "pupil_final_TM_z.bin"});
    pupil_golden_check<rT>(pupil_image, shape);
    return pupil_image;
}

int main()
{
    py_engine::init();
    gen_pupil_array<float>(13.5, 0, 0.8);
    return 0;
}