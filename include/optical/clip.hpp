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
        int verbose;
        using print_type = std::tuple<int, double, int, std::string, std::string, std::string, vec2<double>, vec2<size_t>, double, double, double>;
    }user;

    static std::pair<user_config, py::dict> get_user_config(const std::string& path)
    {
        using T = inverse_tuple_t<user_config::print_type>;
        auto workspace = py_plugin::exec({path});

        std::array<std::string, 11> params{
            "wavelength", "maxNA", "maxSigma", 
            "tilesize", "ambit", 
            "gauge_file", "oas_file", 
            "cell_name", "layer_id", "dbu", "verbose"
        };
        std::reverse(params.begin(), params.end());
        auto config_in_py = get_dict_values<user_config::print_type>(std::vector<std::string>(params.begin(), params.end()), workspace);
        user_config cfg = reinterpret_cast<user_config&>(config_in_py);
        debug_unclassified::verbose() =  (-1 < cfg.verbose); 
        debug_unclassified(std::vector<user_config::print_type>{config_in_py}, params, 70);
        return {cfg, convert_to<py::dict>(workspace)};
    }
    static cutline_jobs cutline_clip_flow(const cutline_jobs::user_config& config, vec2<double> shape_in_um)
    {
        std::vector<std::string> cmd {"./core_plugins/gauge_io.py",
            config.gauge_file, config.oas_file, config.cell_name, std::to_string(config.layer_id),
            "--shape", std::to_string(shape_in_um[0]) + ", " + std::to_string(shape_in_um[1]), 
            "--verbose", std::to_string(config.verbose)
        };
        auto workspace = py_plugin::exec(cmd);
        std::string clip_dir = py::extract<std::string>(workspace["workdir"]);
        np::array2df cutlines_in_um = np::array(workspace["cutlines_in_um"]);
        np::array2df mid_points_in_um = np::array(workspace["mid_points_in_um"]);
        std::cout << "clip_dir, job-count = " << std::make_tuple(clip_dir, cutlines_in_um.shape(0)) << std::endl;
        return cutline_jobs{clip_dir, workspace, cutlines_in_um, mid_points_in_um, config};
    }
    
    static void cutline_clip_flow_v1(const cutline_jobs& config, vec2<double> shape_in_um)
    {
        const std::array<std::string, 2> suffix{";", ""};
        std::string start_points;
        auto [pLines, line_size] = ndarray_ref_no_padding<rectangle<double>>(config.cutlines_in_um);
        for(size_t i = 0; i < line_size; i++, pLines++){
            const auto [from, to] = *pLines;
            vec2<double> start = (to + from - shape_in_um)/2;
            start_points += std::to_string(start[0]) + ", " + std::to_string(start[1]) + suffix.at(i == (line_size - 1));
        }
        
        std::vector<std::string> cmd {"./core_plugins/klayout_op.py",
            config.user.oas_file, config.user.cell_name, std::to_string(config.user.layer_id),
            "--start-points", start_points,
            "--shape", std::to_string(shape_in_um[0]) + ", " + std::to_string(shape_in_um[1]), 
        };
        
        auto workspace = py_plugin::exec(cmd);
    }
};
