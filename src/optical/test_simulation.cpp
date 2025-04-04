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

bool verbose = false;

void regist_geometry_pyclass();
void simulation_flow(const std::string&);
int main()
{
    py_engine::init();
    cutline_data::regist_geometry_pyclass();
    add_path_to_sys_path("core_plugins");
    catch_py_error(simulation_flow("/home/like/repos/simulation/config/calib_193.py"));
    py_engine::dispose();
    printf("simulation end\n");
}
BOOST_PYTHON_MODULE(lib_test_simulation) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("simulation", &simulation_flow);
}


// const uca::backend<double>& backend = uca::cpu<double>::ref();
const uca::backend<double>& backend = uca::cpu<double>::ref();
std::tuple<std::vector<double>, grid_info_in_dbu> do_cutline_job(const std::string& oas_path, const cutline_data& data, const user_config& user_config, const py::object params_optional)
{
    const cutline_dbu& cutline = data.cutline;
    auto startstep_in_dbu = optical_numerics_in_dbu(cutline, user_config);
   
    //== cutline subclip
    auto shape = convert_to<vec2<double>>((startstep_in_dbu.spatial.step * startstep_in_dbu.tilesize)) * user_config.dbu;
    
    //== load subclip
    using thin_mask_solver = thin_mask<double>;
    shapes_dbu shapes = near_filed::load_shapes_from_file(oas_path.c_str(), user_config.cell_name, user_config.layer_id);
    // debug_print<thin_mask_solver>::verbose() = user_config.verbose();
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
    debug_unclassified(rows, {"path", "pattern-name", "start(um)", "weight", "intensity"}, -1);
    // if(std::filesystem::path(oas_path).stem().string() == std::to_string(user_config.vb))
    auto cutline_debug_list = convert_to<std::vector<int>>(params_optional["cutline_debug_list"]);
    if(std::find(cutline_debug_list.begin(), cutline_debug_list.end(), std::stoi(std::filesystem::path(oas_path).stem().string())) != cutline_debug_list.end())
    {
        display_cutline(data, edge_image, cutline_meta.spatial.start, cutline_meta.spatial.step, user_config.dbu,  convert_to<double>(params_optional["threshold_guess"]));
        imshow(y, convert_to<std::vector<size_t>>(mask_info.tilesize));
    }
    return {edge_image, cutline_meta};
}

void simulation_flow(const std::string& config_path)
{

    //== user settings
    auto [user_config, user_config_table] = user_config::load_form_file(config_path);
    
    //== load gauge file
    auto gg_table = cutline_data::load_gauge_file(user_config.gauge_file, user_config.dbu);
    cutline_data::print(gg_table);

    //== calc startstep
    auto startstep_in_dbu = optical_numerics_in_dbu(gg_table.at(0).cutline, user_config);
   
    //== cutline subclip
    auto clip = clip_data::cutline_clip_flow(user_config, gg_table,
        convert_to<vec2<double>>((startstep_in_dbu.spatial.step * startstep_in_dbu.tilesize)) * user_config.dbu
    );
   
    debug_unclassified::verbose() = true;
    //== calc cutline
    std::vector<terms_cutline<double>> edges;
    edges.reserve(gg_table.size());
    std::vector<std::vector<vec<double, 5>>> cutline_features;
    for(size_t i = 0; i < gg_table.size(); i++)
    {
        if(gg_table.at(i).weight == 0) continue;
        auto [edge_image, meta] = do_cutline_job(clip.clip_path(i), gg_table.at(i), user_config, user_config_table);
        edges.push_back(terms_cutline<double>{edge_image});
        auto [start, step] = meta.spatial;
        auto features = get_feature_intensity_from_cutline<double>(
            gg_table.at(i).cutline, 
            um_to_dbu((const double)gg_table.at(i).measured_cd, user_config.dbu),
            edges.back(), start, step
        );
        cutline_features.push_back(features);
        debug_unclassified::out("    features is ", features);
    }
    std::cout << "    optical threshold is " << resist_least_squares<double>::calib_optical_threshold(gg_table, cutline_features) << std::endl;
}