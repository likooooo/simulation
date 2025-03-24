#include <optical/near_field/transfer_matrix_methord.hpp>
#include <py_helper.hpp>
#include <algorithm>
#include <map>
#include <functional>

void fresnel_test(const std::string& name, const std::vector<meterialf>& nk, float assert_eps = 1e-6)
{
    std::cout << "* " << name << std::endl;
    size_t N = 100;
    std::vector<float> rs_power(N);
    std::vector<float> ts_power(N);
    std::vector<float> rp_power(N);
    std::vector<float> tp_power(N);
    std::vector<float> rs_real_part(N);
    std::vector<float> ts_real_part(N);
    std::vector<float> rp_real_part(N);
    std::vector<float> tp_real_part(N);
    
    float step = 0.5_PI /(N - 1);
    for(size_t i = 0; i < N; i++)
    {
        float crao = i * step;
        auto dir = TMM_TEf::propagate_direction(nk, crao);

        auto tmm = TMM_TEf::interface_transfer_matrix(nk, dir);
        auto [Rs, Ts] = TMM_TEf::get_r_t_power_from_tmm(tmm.at(0), nk.at(0).nk, dir.at(0), nk.at(1).nk, dir.at(1));
        assert(TMM_TEf::fresnel::check_power_refrection_transmission(Rs, Ts, assert_eps));
        rs_power.at(i) = Rs;
        ts_power.at(i) = Ts;
        auto [rs, ts] = TMM_TMf::get_r_t_from_tmm(tmm.at(0));
        assert(TMM_TEf::fresnel::check_refrection_transmission(rs, ts, nk.at(0).nk, nk.at(1).nk));
        rs_real_part.at(i) = rs.real();
        ts_real_part.at(i) = ts.real();
        // rs_real_part.at(i) = rs.real() / std::abs(rs.real()) * std::abs(rs);
        // ts_real_part.at(i) = ts.real() / std::abs(ts.real()) * std::abs(ts);

        tmm = TMM_TMf::interface_transfer_matrix(nk, dir);
        auto [Rp, Tp] = TMM_TMf::get_r_t_power_from_tmm(tmm.at(0), nk.at(0).nk, dir.at(0), nk.at(1).nk, dir.at(1));
        assert(TMM_TMf::fresnel::check_power_refrection_transmission(Rp, Tp, assert_eps));
        rp_power.at(i) = Rp;
        tp_power.at(i) = Tp;
        auto [rp, tp] = TMM_TMf::get_r_t_from_tmm(tmm.at(0));
        assert(TMM_TMf::fresnel::check_refrection_transmission(rp, tp, nk.at(0).nk, nk.at(1).nk));
        rp_real_part.at(i) = rp.real();
        tp_real_part.at(i) = tp.real();
        // rp_real_part.at(i) = rp.real() / std::abs(rp.real()) * std::abs(rp);
        // tp_real_part.at(i) = tp.real() / std::abs(tp.real()) * std::abs(tp);

    }
    std::vector<float> non_polarized = (rp_power + rs_power) * float(0.5);
    std::vector<std::vector<float>> curves{
 
    };
    auto it_min = std::min_element(rp_power.begin(), rp_power.end());
    float brewster_angle = 90.f/(N - 1) * std::distance(rp_power.begin(), it_min);
    std::cout <<"    brewster's angle = " << brewster_angle << "°" << " at Rp_power = " << *it_min << std::endl;
    float brewster_angle_analytic_solution = TMM_TMf::fresnel::snell::brewster_angle(nk.at(0).nk, nk.at(1).nk) / 1_PI * 180;
    std::cout <<"    brewster's angle(analytic solution) = " << brewster_angle_analytic_solution << "°" << std::endl;
    assert(is_almost_equal(brewster_angle,  brewster_angle_analytic_solution, 0.09f));

    auto it_critical_angle = std::find_if(rp_power.begin(), rp_power.end(), [assert_eps](float rp){return 1.0f - rp < assert_eps;});
    if(rp_power.end() != it_critical_angle)
    {
        float critical_angle = 90.f/(N - 1) * std::distance(rp_power.begin(), it_critical_angle);
        std::cout <<"    total internal reflection = " << critical_angle << "°" << " at Rp_power = " << *it_critical_angle << std::endl;
    }
    else
    {
        std::cout <<"    NO total internal reflection was found" << std::endl;
    }
    auto unpolarized_reflection   = (rs_power + rp_power) * 0.5f;
    auto unpolarized_transimition = (ts_power + tp_power) * 0.5f;
    plot_curves(
        std::vector<std::vector<float>>{
            rs_power, ts_power, rp_power, tp_power, unpolarized_reflection, unpolarized_transimition
        },
        std::vector<float>(6, 0),
        std::vector<float>(6, 90.f/(N - 1)),
        std::vector<std::string>{
            "R_s", "T_s", "R_p", "T_p", "R_unpolarized", "T_unpolarized"
            //  "Rs", "Ts", "Rp", "Tp"
        },
        std::vector<std::string>{
            "b--", "b-", "r--", "r-", "c.", "y."
        }
    );
    plot_curves(
        std::vector<std::vector<float>>{
            rs_real_part, ts_real_part, rp_real_part, tp_real_part
        },
        std::vector<float>(4, 0),
        std::vector<float>(4, 90.f/(N - 1)),
        std::vector<std::string>{
             "rs", "ts", "rp", "tp"
        },
        std::vector<std::string>{
            "b--", "b-", "r--", "r-", 
        }
    );
    std::cout << "---> test success. " << name << std::endl << std::endl;
} 


int main(int argc, char** argv)
{
    py_loader::init();
    py_plot::get_default_visualizer_dir() = "/usr/local/bin";
    std::vector<meterialf> glass_to_air{
        meterialf{complex_t<float>(1.5, 0), 0},
        meterialf{complex_t<float>(1.0, 0), 0},
    };std::vector<meterialf> air_to_glass{
        meterialf{complex_t<float>(1.0, 0), 0},
        meterialf{complex_t<float>(1.5, 0), 0},
    };
    using callback =  std::function<void()>;
    std::map<std::string, callback> callbacks{
        {"air_to_glass", callback([&](){fresnel_test("fresnel equations(air to glass)", air_to_glass);})},
        {"glass_to_air", callback([&](){fresnel_test("fresnel equations(glass to air)", glass_to_air);})}
    };
    std::string key = argc == 2 ? std::string(argv[1]) : "glass_to_air";
    auto it_call = callbacks.find(key);
    assert(callbacks.end() != it_call);
    (it_call->second)();   
}