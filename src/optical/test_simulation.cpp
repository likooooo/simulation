#include <py_helper.hpp>
#include <optical/optical_numerics.hpp>
#include <type_traist_notebook/type_traist.hpp>
#include <filesystem>
#include <optical/clip.hpp>
#include <optical/near_field/thin_mask.hpp>
#include <optical/resist/resist_cnn.hpp>
// #include <cpp_cuda/cuda_vector.hpp>
#include <type_traist_notebook/uca/backend.hpp>
#ifdef CPU_BACKEND_ENABLE
#   include <mekil/cpu_backend.hpp>
#endif
#ifdef GPU_BACKEND_ENABLE
#   include <cpp_cuda/gpu_backend.hpp>
#endif

void regist_simulation_pyclass();
void simulation_flow(const std::string&);

bool regist_py = py_engine::regist_py_custom(cutline_data::regist_geometry_pyclass) && py_engine::regist_py_custom(regist_simulation_pyclass);
int main()
{
    py_engine::init();
    catch_py_error(simulation_flow("/home/like/repos/simulation/config/calib_193.py"));
    py_engine::dispose();
    std::cout << "simulation end" << std::endl;
}
BOOST_PYTHON_MODULE(lib_test_simulation) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("simulation_flow", &simulation_flow);
}


// const uca::backend<double>& backend = uca::cpu<double>::ref();
const uca::backend<double>& backend = uca::cpu<double>::ref();
std::tuple<std::vector<double>, grid_info_in_dbu> do_cutline_job(const std::string& oas_path, const cutline_data& data, const user_config& user_config, const py::object params_optional)
{
    const cutline_dbu& cutline = data.cutline;
    auto startstep_in_dbu = optical_numerics_in_dbu(cutline, user_config);
    //== load subclip
    using thin_mask_solver = thin_mask<double>;
    shapes_dbu shapes = near_filed::load_shapes_from_file(oas_path.c_str(), user_config.cell_name, user_config.layer_id);
    auto [x, y, mask_info] = thin_mask_solver::edge_pixelization(startstep_in_dbu, shapes, convert_to<size_t>(params_optional["mask_USF"]),  convert_to<double>(params_optional["mask_edge_dissect_coef"]));
    print_grid_start_step<grid_info_in_dbu, debug_print<thin_mask_solver>>(mask_info, "    intergral image");
    backend.integral_y(mask_info.tilesize, x.data());
    backend.integral_x(mask_info.tilesize, y.data());
    backend.VtAdd(x.size(), x.data(), y.data());
    auto [edge_image, cutline_meta] = thin_mask_solver::get_edge_from_rasterization(mask_info, y, cutline);
    using print_type = std::tuple<std::string, std::string, vec2<double>, double, std::vector<double>>;
    std::vector<print_type> rows{
        print_type(
            std::filesystem::path(oas_path).filename(), 
            data.pattern_name, 
            dbu_to_um(convert_to<vec2<double>>(cutline_meta.spatial.start), user_config.dbu), 
            data.weight,
            edge_image 
        )
    };
    debug_unclassified(rows, {"path", "pattern-name", "start(um)", "weight", "intensity"}, 100);
    auto cutline_debug_list = convert_to<std::vector<int>>(params_optional["cutline_debug_list"]);
    if(std::find(cutline_debug_list.begin(), cutline_debug_list.end(), std::stoi(std::filesystem::path(oas_path).stem().string())) != cutline_debug_list.end())
    {
        display_cutline(data, edge_image, cutline_meta.spatial.start, cutline_meta.spatial.step, user_config.dbu, convert_to<double>(params_optional["threshold_guess"]));
        imshow(y, convert_to<std::vector<size_t>>(mask_info.tilesize));
        //== TODO : extract contour
        // py_plugin::ref()["extract_contours"]["show_contuor"](create_ndarray_from_vector(y, convert_to<std::vector<int>>(mask_info.tilesize)),  convert_to<double>(params_optional["threshold_guess"]));
    }

    return {edge_image, cutline_meta};
}

std::tuple<std::vector<terms_cutline<double>>, std::vector<terms_features_intensity<double>>, grid_info_in_dbu>
gather_cutline_and_features(const user_config& config,  py::dict user_config_table, const std::vector<cutline_data>& gg_table, const clip_data& clip)
{
    std::vector<terms_cutline<double>> edges;
    std::vector<terms_features_intensity<double>> cutline_features;
    edges.reserve(gg_table.size());
    cutline_features.reserve(gg_table.size());
    grid_info_in_dbu edge_meta;
    for(size_t i = 0; i < gg_table.size(); i++)
    {
        if(gg_table.at(i).weight == 0) continue;
        const auto& [edge_image, meta] = do_cutline_job(clip.clip_path(i), gg_table.at(i), config, user_config_table);
        edge_meta = meta;
        edges.push_back(terms_cutline<double>{edge_image});
        auto [start, step] = meta.spatial;
        assert(gg_table.at(i).cutline[0] == start);

        auto features = get_feature_intensity_from_cutline<double>(
            gg_table.at(i).cutline, 
            gg_table.at(i).measured_cd,
            edges.back(), start, step, config.dbu
        );
        cutline_features.push_back(features);
        debug_unclassified::out("    cutline features is ", features);
    }
    return {edges, cutline_features, edge_meta};
}

struct simulation{
    user_config config;
    py::dict user_config_table;
    std::vector<cutline_data> gg_table;
    clip_data clip;
    std::vector<terms_cutline<double>> edges;
    std::vector<terms_features_intensity<double>> cutline_features;
    grid_info_in_dbu edge_meta;
    double threshold;

    py::dict load_user_config(const std::string& config_path, bool verbose = false)
    {
        verbose_guard<debug_unclassified> vb(debug_unclassified::verbose() || verbose);

        auto results = user_config::load_form_file(config_path);
        config = std::get<0>(results);
        user_config_table = std::get<1>(results);
        vb.backup = config.verbose(); 
        return user_config_table;
    }
    void load_gauge_file(bool verbose = false)
    {
        verbose_guard<debug_unclassified> vb(debug_unclassified::verbose() || verbose);
        gg_table = cutline_data::load_gauge_file(config.gauge_file, config.dbu);
        cutline_data::print(gg_table);
    }
    void clip_cutline(bool verbose = false)
    {
        verbose_guard<debug_unclassified> vb(debug_unclassified::verbose() || verbose);
        clip = clip_data::cutline_clip_flow(config, gg_table);
    }
    void gather_cutline_and_features(bool verbose = false)
    {
        verbose_guard<debug_unclassified> vb(debug_unclassified::verbose() || verbose);
        auto result = ::gather_cutline_and_features(config, user_config_table, gg_table, clip);
        edges = std::get<0>(result);
        cutline_features = std::get<1>(result);
        edge_meta = std::get<2>(result);
    }
    void calib_optical_threshold(bool verbose = false)
    {
        verbose_guard<debug_unclassified> vb(debug_unclassified::verbose() || verbose);
        threshold = resist_least_squares<double>::calib_optical_threshold(gg_table, cutline_features);
        post_calib_analysis(gg_table, edges, {1}, threshold, edge_meta.spatial.step, config.dbu);
        cutline_data::print(gg_table);
        std::cout << "    optical threshold is " << threshold << std::endl;
    }
};


void regist_simulation_pyclass()
{
    py::class_<simulation>("simulation").def(py::init<>())       
        .def("load_user_config", &simulation::load_user_config)
        .def("load_gauge_file", &simulation::load_gauge_file)
        .def("clip_cutline", &simulation::clip_cutline)
        .def("gather_cutline_and_features", &simulation::gather_cutline_and_features)
        .def("calib_optical_threshold", &simulation::calib_optical_threshold)
    ;
}

void simulation_flow(const std::string& config_path)
{
    simulation s;
    s.load_user_config(config_path);
    s.load_gauge_file();
    s.clip_cutline();
    s.gather_cutline_and_features(false);
    s.calib_optical_threshold(true);
}