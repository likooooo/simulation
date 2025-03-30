#pragma once
#include "geometry.hpp"

struct cutline_jobs
{
    std::string clip_dir;
    py::object workspace;
    np::array2df cutlines_in_um;
    np::array2df mid_points_in_um;
    size_t size() const{
        return cutlines_in_um.shape(0);
    }
    cutline_jobs() = default;
    std::filesystem::path clip_path(size_t n) const{
        return std::filesystem::path(clip_dir) / (std::to_string(n) + ".oas");
    }
    double file_dbu() const{
        static double n = 0;
        if(0 == n) n = py_plugin::call<double>("klayout_op", "get_dbu", clip_path(0).c_str());
        return n;
    }
    struct user_config{
        double wavelength;
        double maxNA;
        double maxSigma;
        vec2<size_t> tilesize;
        vec2<double> ambit;

        std::string gauge_file;
        std::string oas_file;
        std::string cell_name;
        int layer_id;
        double dbu;
        using print_type = std::tuple<double, int, std::string, std::string, std::string, vec2<double>, vec2<size_t>, double, double, double>;
    }user;

    static user_config get_user_config(const std::string& path)
    {
        using T = inverse_tuple_t<user_config::print_type>;
        auto workspace = py_plugin::exec({path});

        std::array<std::string, 10> params{
            "wavelength", "maxNA", "maxSigma", 
            "tilesize", "ambit", 
            "gauge_file", "oas_file", 
            "cell_name", "layer_id", "dbu"
        };
        std::reverse(params.begin(), params.end());
        auto config_in_py = get_dict_values<user_config::print_type>(std::vector<std::string>(params.begin(), params.end()), workspace);
        user_config cfg = reinterpret_cast<user_config&>(config_in_py);

        print_table(std::vector<user_config::print_type>{config_in_py}, params, 70);
        return cfg;
    }
    static cutline_jobs cutline_clip_flow(const cutline_jobs::user_config& config, vec2<double> shape_in_um)
    {
        std::vector<std::string> cmd {"./core_plugins/gauge_io.py",
            config.gauge_file, config.oas_file, config.cell_name, std::to_string(config.layer_id),
            "--shape", std::to_string(shape_in_um[0]) + ", " + std::to_string(shape_in_um[1]), 
            //"--verbose", "3"
        };
        auto workspace = py_plugin::exec(cmd);
        std::string clip_dir = py::extract<std::string>(workspace["workdir"]);
        np::array2df cutlines_in_um = np::array(workspace["cutlines_in_um"]);
        np::array2df mid_points_in_um = np::array(workspace["mid_points_in_um"]);
        std::cout << "clip_dir, job-count = " << std::make_tuple(clip_dir, cutlines_in_um.shape(0)) << std::endl;
        return cutline_jobs{clip_dir, workspace, cutlines_in_um, mid_points_in_um, config};
    }
    
};
