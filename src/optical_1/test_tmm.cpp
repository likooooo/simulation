#include <optical/near_field/transfer_matrix_methord.hpp>
#include <py_helper.hpp>

void reflectance_test(const std::string& name, complex_t<float> target)
{
    std::vector<meterialf> nk{
        meterialf{complex_t<float>(1.0, 0), 0},
        meterialf{target, 0},
    };
    std::cout << "* reflectance test " << name << std::endl;
    meterialf::print(nk);
    size_t N = 100;
    std::vector<float> Rp(N);
    std::vector<float> Rs(N);
    float step = 0.5_PI /(N - 1);
    for(size_t i = 0; i < N; i++)
    {
        float crao = i * step;
        auto dir = TMM_TEf::propagate_direction(nk, crao);
        auto tmm = TMM_TEf::interface_transfer_matrix(nk, dir);
        auto [rs, ts] = TMM_TEf::get_r_t_power_from_tmm(tmm.at(0), nk.at(0).nk, dir.at(0), nk.at(1).nk, dir.at(1));
        assert(is_almost_equal<float>(std::norm(rs + ts), 1.0, 1e-6));

        tmm = TMM_TMf::interface_transfer_matrix(nk, dir);
        auto [rp, tp] = TMM_TMf::get_r_t_power_from_tmm(tmm.at(0), nk.at(0).nk, dir.at(0), nk.at(1).nk, dir.at(1));
        assert(is_almost_equal<float>(std::norm(rp + tp), 1.0, 1e-6));
        Rp.at(i) = rp;
        Rs.at(i) = rs;
    }
    auto display =  py_plot::create_callback_simulation_fram_done();
    size_t max_frame = 50;
    std::cout << "    plot of reflectans (p polarized)\n";
    for(size_t i = 0; i < max_frame; i++) display(create_ndarray_from_vector(Rp, {1, int(N)}));
    std::cout << "    plot of reflectans (s polarized)\n";
    for(size_t i = 0; i < max_frame; i++) display(create_ndarray_from_vector(Rs, {1, int(N)}));
    std::vector<float> non_polarized = (Rp + Rs) * float(0.5);
    std::cout << "    plot of reflectans (non-polarized)\n";
    for(size_t i = 0; i < max_frame; i++) display(create_ndarray_from_vector(non_polarized, {1, int(N)}));
    std::cout << name <<" reflectance test success. " << std::endl << std::endl;
} 

int main()
{
    py_engine::init();
    reflectance_test("SiO2", {1.4649, 0.0017953});   // lambda = 587.6 nm
    reflectance_test("Ta2O5", {2.1306, 0.00091143}); // lambda = 587.6 nm
}