#pragma once
#include "geometry.hpp"

struct cutline_data
{
    std::string pattern_name;
    cutline_dbu cutline;
    int polar;
    double measured_cd;
    double weight;
    static inline std::vector<std::string> ref_names;
    std::vector<double> post_calib_results;
    cutline_data() = default;
    using print_type = std::tuple<std::string, cutline_dbu, int , double, double, std::vector<double>>;
    print_type to_tuple() const
    {
        return std::make_tuple(pattern_name, cutline, polar, measured_cd, weight, post_calib_results);
    }
    static void print(const std::vector<cutline_data>& lines)
    {
        std::vector<print_type> rows; rows.reserve(lines.size());
        std::transform(lines.begin(), lines.end(), std::back_insert_iterator(rows), [](const auto& l){return l.to_tuple();});
        debug_unclassified(rows, {"pattern-name", "cutline(dbu)", "polar", "measured-cd(dbu)", "weight", "post-calib-results"});
        debug_unclassified::out("    post-calib-results are ", ref_names);
    }
    static void regist_geometry_pyclass()
    {
        py::class_<cutline_data>("cutline_data").def(py::init<>())       
            .def_readwrite("pattern_name",& cutline_data::pattern_name)
            .def_readwrite("cutline",& cutline_data::cutline)
            .def_readwrite("polar",& cutline_data::polar)
            .def_readwrite("measured_cd",& cutline_data::measured_cd)
            .def_readwrite("weight",& cutline_data::weight)
            .def_readwrite("ref_names",& cutline_data::ref_names)
            .def_readwrite("post_calib_results",& cutline_data::post_calib_results);
    }
    
    static std::vector<cutline_data> load_gauge_file(const std::string& path, double dbu)
    {
        py::object gauge_table = py_plugin::ref()["gauge_io"]["read_gauge_file"](path);
        auto lines = convert_to<py::list>(py_plugin::ref()["gauge_io"]["get_cutline_datas"](gauge_table, dbu));
        size_t N = len(lines);
        std::vector<cutline_data> gg_basic_table;
        gg_basic_table.reserve(N);
        for(size_t i = 0; i < N; i++)
        {
            pyobject_wrapper line(lines[i]);
            cutline_data temp;
            temp.pattern_name = convert_to<std::string>(line["pattern_name"]);
            temp.cutline = convert_to<cutline_dbu>(line["cutline"]);
            temp.polar = convert_to<int>(line["polar"]);
            temp.measured_cd = convert_to<double>(line["measured_cd"]);
            temp.weight = convert_to<double>(line["weight"]);
            gg_basic_table.push_back(temp);
        }    
        cutline_data::ref_names = convert_to<std::vector<std::string>>(py_plugin::ref()["gauge_io"]["post_calib_options"]);
        auto max_weight = std::max_element(gg_basic_table.begin(), gg_basic_table.end(), [](const auto& a, const auto& b){return a.weight < b.weight;})->weight;
        for(size_t i = 0; i < N; i++) gg_basic_table.at(i).weight /= max_weight;
        return gg_basic_table;
    }
};
struct user_config
{
    using print_type = std::tuple<int, double, int, std::string, std::string, std::string, vec2<double>, vec2<size_t>, double, double, double>;
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
    int vb;
    bool verbose() const {return -1 < vb;}

    static std::pair<user_config, py::dict> load_form_file(const std::string& path)
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
        debug_unclassified::verbose() =  cfg.verbose(); 
        debug_unclassified(std::vector<user_config::print_type>{config_in_py}, params, 70);
        return {cfg, convert_to<py::dict>(workspace)};
    }
};
struct clip_data
{
    std::string workdir;
    clip_data() = default;
    std::filesystem::path clip_path(size_t n) const{
        return std::filesystem::path(workdir) / (std::to_string(n) + ".oas");
    }
    static clip_data cutline_clip_flow(const user_config& config, const std::vector<cutline_data>& table, vec2<double> shape_in_um)
    {
        std::vector<double> cutlines_in_um;
        cutlines_in_um.reserve(4 * table.size());
        for(const auto& t : table)
        {
            auto [from, to] = t. cutline;
            cutlines_in_um.push_back(config.dbu * from[0]);
            cutlines_in_um.push_back(config.dbu * from[1]);
            cutlines_in_um.push_back(config.dbu * to[0]);
            cutlines_in_um.push_back(config.dbu * to[1]);
        }

        py::object str = py_plugin::ref()["gauge_io"]["clip_flow"](
            create_ndarray_from_vector(cutlines_in_um, {4, int(table.size())}), 
            config.oas_file, config.cell_name, config.layer_id, shape_in_um, config.vb
        );
        return clip_data{py::extract<std::string>(str)};
    }
    std::filesystem::path clip_workdir(size_t n) const{
        std::filesystem::path dir(workdir);
        dir /= std::to_string(n);
        if (!std::filesystem::exists(dir)) assert(std::filesystem::create_directories(dir));
        return dir;
    }
};
