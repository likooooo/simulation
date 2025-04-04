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
    debug_print<thin_mask_solver>::verbose() = user_config.verbose();
    auto [x, y, mask_info] = thin_mask_solver::edge_pixelization(startstep_in_dbu, shapes, convert_to<size_t>(params_optional["mask_USF"]),  convert_to<double>(params_optional["mask_edge_dissect_coef"]));
    backend.integral_y(mask_info.tilesize, x.data());
    backend.integral_x(mask_info.tilesize, y.data());
    backend.VtAdd(x.size(), x.data(), y.data());

    auto [edge_image, cutline_meta] = thin_mask_solver::get_edge_from_rasterization(mask_info, y, cutline);
    if(std::filesystem::path(oas_path).stem().string() == std::to_string(user_config.vb))
    {
        imshow(y, convert_to<std::vector<size_t>>(mask_info.tilesize));
        auto start = dbu_to_um(convert_to<vec2<float>>(cutline_meta.spatial.start), user_config.dbu);
        auto step = dbu_to_um(convert_to<vec2<float>>(cutline_meta.spatial.step), user_config.dbu);
        // plot_curves(std::vector<std::vector<double>>{cutline}, {start[0]}, {float(step[0])}, {"cutline (um)"}, {"b--"});  
        auto center = dbu_to_um((cutline[0] + cutline[1]) / 2, user_config.dbu);
        auto [features_in_dbu, dir] = get_feature_pos_from_cutline(cutline, data.measured_cd, cutline_meta.spatial.start, cutline_meta.spatial.step);
        auto features = dbu_to_um(convert_to<vec<float, 5>>(features_in_dbu), user_config.dbu);
        features += start[dir];
        plot_curves(std::vector<std::vector<double>>{
            edge_image, std::vector<double>{1}, std::vector<double>{0.5, 0.5}, std::vector<double>{0.25, 0.25}}, 
            {start[dir], features[0], features[1], features[3]}, 
            {step[dir], step[dir] , features[2] - features[1], features[4] - features[3]}, 
            {"cutline (um)", "center", "on", "out"}, {"b--", "r-x", "g--x", "r--x"});  
    }
    std::cout << "cutline of " << oas_path << " is " << edge_image << std::endl;
    return {edge_image, cutline_meta};
}

void simulation_flow(const std::string& config_path)
{
    debug_unclassified::verbose() = true;

    //== user settings
    auto [user_config, user_config_table] = user_config::load_form_file(config_path);
    
    //== load gauge file
    auto gg_table = cutline_data::load_gauge_file(user_config.gauge_file, user_config.dbu);
   
    //== calc startstep
    auto startstep_in_dbu = optical_numerics_in_dbu(gg_table.at(0).cutline, user_config);
   
    //== cutline subclip
    auto clip = clip_data::cutline_clip_flow(user_config, 
        convert_to<vec2<double>>((startstep_in_dbu.spatial.step * startstep_in_dbu.tilesize)) * user_config.dbu
    );
   
    //== calc cutline
    std::vector<std::vector<std::vector<double>>> edges;
    edges.reserve(gg_table.size());
    for(size_t i = 0; i < gg_table.size(); i++)
    {
        auto [edge_image, meta] = do_cutline_job(clip.clip_path(i), gg_table.at(i), user_config, user_config_table);
        edges.push_back(std::vector<std::vector<double>>{edge_image});
        auto [start, step] = meta.spatial;
        std::cout << "features=" << get_feature_intensity_from_cutline<double>(
            gg_table.at(i).cutline, 
            gg_table.at(i).measured_cd,
            edges.back(), start, step
        ) << std::endl;
    }
}