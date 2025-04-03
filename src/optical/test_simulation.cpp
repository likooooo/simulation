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

void simulation_flow(const std::string&);
int main()
{
    py_engine::init();
    add_path_to_sys_path("core_plugins");
    catch_py_error
    (simulation_flow("/home/like/repos/simulation/config/calib_193.py"));
    py_engine::dispose();
}
BOOST_PYTHON_MODULE(lib_test_simulation) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("simulation", &simulation_flow);
}

struct cutline_data
{
    cutline_dbu cl;
    int64_t design_cd;
    double weight;
    std::shared_ptr<std::vector<std::string>> ref_names;
    std::vector<double> others;
};

std::vector<cutline_data> load_gauge_file(const std::string& path)
{
    // read_gauge_file
    py::list_tuple gauge_table = convert_to<py::list_tuple>(py_plugin::ref()["gauge_io"]["read_gauge_file"](path));
    np::array2df cutlines_in_dbu = convert_to<np::array2df>(py_plugin::ref()["gauge_io"]["get_culine_from_gauge_table"](gauge_table, 1));
    np::array2di design_cd = convert_to<np::array2df>(py_plugin::ref()["gauge_io"]["get_design_cd"](gauge_table));
    auto [pRec, n] = ndarray_ref_no_padding<vec<double, 4>>(cutlines_in_dbu);
    auto [pCD, n1] = ndarray_ref_no_padding<vec<int64_t, 1>>(design_cd);
    assert(n == n1);
    std::vector<cutline_data> gg_table; 
    gg_table.reserve(n);
    for(size_t i = 0; i < n; i++, pRec++, pCD++){
        const auto& rec = *pRec;
        cutline_data data;
        data.cl = cutline_dbu{
            static_cast<int64_t>(rec[0]), static_cast<int64_t>(rec[1]), 
            static_cast<int64_t>(rec[2]), static_cast<int64_t>(rec[3])
        };
        data.design_cd = (*pCD)[0];
        gg_table.push_back(data);
    }
    return gg_table;
}

// const uca::backend<double>& backend = uca::cpu<double>::ref();
const uca::backend<double>& backend = uca::cpu<double>::ref();
std::tuple<std::vector<double>, grid_info_in_dbu> do_cutline_job(const std::string& oas_path, const cutline_dbu& cutline, const cutline_jobs::user_config& user_config ,const py::object params_optional)
{
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
    if(user_config.verbose())
    {
        imshow(y, convert_to<std::vector<size_t>>(mask_info.tilesize));
        auto start = dbu_to_um(convert_to<vec2<float>>(cutline_meta.spatial.start), user_config.dbu);
        auto step = dbu_to_um(convert_to<vec2<float>>(cutline_meta.spatial.step), user_config.dbu);
        // plot_curves(std::vector<std::vector<double>>{cutline}, {start[0]}, {float(step[0])}, {"cutline (um)"}, {"b--"});  
        plot_curves(std::vector<std::vector<double>>{edge_image}, {0}, {float(1)}, {"cutline (pixel)"}, {"b--"});  
    }
    std::cout << "cutline of " << oas_path << " is " << edge_image << std::endl;
    return {edge_image, cutline_meta};
}

void simulation_flow(const std::string& config_path)
{
    auto [user_config, params] = cutline_jobs::get_user_config(config_path);
    //== load gauge file & calc startstep
    auto gg_table = load_gauge_file(user_config.gauge_file);

    auto startstep_in_dbu = optical_numerics_in_dbu(gg_table.at(0).cl, user_config);
   
    //== cutline subclip
    auto shape = convert_to<vec2<double>>((startstep_in_dbu.spatial.step * startstep_in_dbu.tilesize)) * user_config.dbu;
    cutline_jobs jobs = cutline_jobs::cutline_clip_flow(user_config, shape);
    std::vector<std::vector<std::vector<double>>> edges;
    edges.reserve(gg_table.size());
    for(size_t i = 0; i < gg_table.size(); i++)
    {
        auto [edge_image, meta] = do_cutline_job(jobs.clip_path(i), gg_table.at(i).cl, user_config, params);
        edges.push_back(std::vector<std::vector<double>>{edge_image});
        auto [start, step] = meta.spatial;
        std::cout << "features=" << get_feature_intensity_from_cutline<double>(gg_table.at(i).cl, gg_table.at(i).design_cd, edges.back(), start, step) << std::endl;
    }
}

// void a(cutline_dbu cutline, int64_t design_cd, const std::vector<std::vector<double>>& yArray, point_dbu start, point_dbu step)
// {
//     get_feature_intensity_from_cutline<double>(cutline, design_cd, yArray, start, step);
// }
