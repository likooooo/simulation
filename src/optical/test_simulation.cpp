#include <filesystem>
#include <optical/resist/resist_least_squares.hpp>
#include <optical/gather_cutline.hpp>

void regist_simulation_pyclass();
void simulation_flow(const std::string&);

bool regist_py = py_engine::regist_py_custom(cutline_data::regist_geometry_pyclass) && py_engine::regist_py_custom(regist_simulation_pyclass);
const uca::backend<double>& backend = uca::cpu<double>::ref();
// const uca::backend<double>& backend = uca::cpu<double>::ref();

int main()
{
    py_engine::init();
    catch_py_error(simulation_flow("/home/like/repos/simulation/config/calib_193.py"));
    py_engine::dispose();
    std::cout << "thin-mask simulation end" << std::endl;
}
BOOST_PYTHON_MODULE(lib_test_simulation) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("simulation_flow", &simulation_flow);
}


inline std::tuple<terms_cutline<double>, grid_info_in_dbu> get_thin_mask_cutline(const std::string& oas_path, const cutline_data& data, const user_config& user_config, const py::object params_optional)
{
    const cutline_dbu& cutline = data.cutline;
    auto startstep_in_dbu = optical_numerics_in_dbu(cutline, user_config);
    //== load subclip
    using thin_mask_solver = thin_mask<double>;
    shapes_dbu shapes = near_filed::load_shapes_from_file(oas_path.c_str(), user_config.cell_name, user_config.layer_id);
    auto [mask, mask_info] = thin_mask_solver::mask_image(startstep_in_dbu, shapes, convert_to<size_t>(params_optional["mask_USF"]),  convert_to<double>(params_optional["mask_edge_dissect_coef"]));
    auto [edge_image, cutline_meta] = thin_mask_solver::get_edge_from_rasterization(mask_info, mask, cutline);
    //== debug cutline
    check_cutline<double>(edge_image, cutline_meta, oas_path, data, user_config, params_optional, 
        [&](const cutline_data& data, const std::vector<double>& cutline_image,const point_dbu& start_dbu,  const point_dbu& step_dbu, float dbu){
            display_cutline_with_cd(data, cutline_image, start_dbu, step_dbu, dbu,  convert_to<double>(params_optional["threshold_guess"]));
            imshow(mask, convert_to<std::vector<size_t>>(mask_info.tilesize));
            // py_plugin::ref()["extract_contours"]["show_contuor"](create_ndarray_from_vector(mask, convert_to<std::vector<int>>(mask_info.tilesize)),  convert_to<double>(params_optional["threshold_guess"]));
        }
    );
    return {terms_cutline<double>{edge_image}, cutline_meta};
}

struct thin_mask_simulation : simulation_common{
    double threshold;
    void gather_cutline_and_features(bool verbose = false)
    {
        verbose_guard<debug_unclassified> vb(debug_unclassified::verbose() || verbose);
        auto result = ::gather_cutline_and_features(config, user_config_table, gg_table, clip, get_thin_mask_cutline);
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
    py::class_<thin_mask_simulation>("thin_mask_simulation").def(py::init<>())      
        .def_readwrite("gauge_table",& thin_mask_simulation::gg_table)
        .def("load_user_config", &thin_mask_simulation::load_user_config)
        .def("load_gauge_file", &thin_mask_simulation::load_gauge_file)
        .def("clip_cutline", &thin_mask_simulation::clip_cutline)
        .def("gather_cutline_and_features", &thin_mask_simulation::gather_cutline_and_features)
        .def("calib_optical_threshold", &thin_mask_simulation::calib_optical_threshold)
    ;
}

void simulation_flow(const std::string& config_path)
{
    thin_mask_simulation s;
    s.load_user_config(config_path);
    s.load_gauge_file();
    s.clip_cutline();
    s.gather_cutline_and_features(false);
    s.calib_optical_threshold(true);
}