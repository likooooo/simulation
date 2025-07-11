#include <optical/resist/resist_gauss_laguerre_svd.hpp>
#include <optical/gather_cutline.hpp>
#include <optical/resist/resist_least_squares.hpp>

using resist = resist_gauss_laguerre_svd<double>;
void simulation_flow(const std::string& config_path);
void regist_simulation_pyclass();
bool regist_py = py_engine::regist_py_custom(cutline_data::regist_geometry_pyclass) && py_engine::regist_py_custom(regist_simulation_pyclass);
int main()
{
    py_engine::init();
    catch_py_error(simulation_flow("/home/like/repos/simulation/config/calib_193.py"));
    py_engine::dispose();
    std::cout << "gauss-laguerre simulation end" << std::endl;
}
BOOST_PYTHON_MODULE(lib_test_resist_cnn) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
}

const uca::backend<double>& backend = uca::cpu<double>::ref();
std::tuple<terms_cutline<double>, grid_info_in_dbu> get_resist_cutline(const std::string& oas_path, const cutline_data& data, const user_config& user_config, const py::object params_optional)
{
    const cutline_dbu& cutline = data.cutline;
    auto startstep_in_dbu = optical_numerics_in_dbu(cutline, user_config);
    //== load subclip
    using thin_mask_solver = thin_mask<double>;
    shapes_dbu shapes = near_filed::load_shapes_from_file(oas_path.c_str(), user_config.cell_name, user_config.layer_id);
    auto [mask, mask_info] = thin_mask_solver::mask_image(startstep_in_dbu, shapes, convert_to<size_t>(params_optional["mask_USF"]),  convert_to<double>(params_optional["mask_edge_dissect_coef"]));
    assert(mask_info.spatial.step[0] ==  mask_info.spatial.step[1]);
    auto [images, N] = resist::gauss_laguerre(mask, mask_info.tilesize, 
        convert_to<double>(params_optional["gauss_laguerre_sigma_in_dbu"]) / mask_info.spatial.step[0], 
        convert_to<size_t>(params_optional["associate_order"]), 
        convert_to<size_t>(params_optional["laguerre_order"])
    );
    terms_cutline<double> terms;
    terms.reserve(images.size() + 1);
    {
        auto [term, cutline_meta_all_the_same] = thin_mask_solver::get_edge_from_rasterization(mask_info, mask, cutline);
        terms.push_back(term);
    }
    grid_info_in_dbu cutline_meta;
    for(const auto& l : images){
        auto [term, cutline_meta_all_the_same] = thin_mask_solver::get_edge_from_rasterization(mask_info, l, cutline);
        terms.push_back(term);
        cutline_meta = cutline_meta_all_the_same;
        //== debug cutline
        check_cutline<double>(term, cutline_meta, oas_path, data, user_config, params_optional, 
            [&](const cutline_data& data, const std::vector<double>& cutline_image, const point_dbu& start_dbu,  const point_dbu& step_dbu, float dbu){
                display_cutline(data, cutline_image, start_dbu, step_dbu, dbu);
                imshow(l, convert_to<std::vector<size_t>>(mask_info.tilesize));
            }
        );
    }
    auto resist_coefficients = convert_to<std::vector<double>>(params_optional["resist_coefficients"]);
    if(resist_coefficients.size())
    {
        auto resist = mask;
        VecScala<double>(resist.size(), resist_coefficients.front(), resist.data());

        for(size_t i = 0; i < images.size(); i++)
        {
            VecScala<double>(images.at(i).size(), resist_coefficients.at(i + 1), images.at(i).data());
            backend.VtAdd(resist.size(), images.at(i).data(), resist.data());
        }
        // py_plugin::ref()["extract_contours"]["find_and_plot_contours"](create_ndarray_from_vector(resist, convert_to<std::vector<int>>(mask_info.tilesize)),  convert_to<double>(params_optional["threshold_guess"]));
        auto [optical_cutline, cutline_meta_all_the_same] = thin_mask_solver::get_edge_from_rasterization(mask_info, resist, cutline);
        display_cutline_with_cd(data, optical_cutline,  cutline_meta_all_the_same.spatial.start, cutline_meta_all_the_same.spatial.step, user_config.dbu, convert_to<double>(params_optional["threshold_guess"]));
        imshow(resist, convert_to<std::vector<size_t>>(mask_info.tilesize));
    }
    return {terms, cutline_meta};
}

struct resis_simulation : simulation_common{
    double threshold;
    std::vector<terms_dense_features_intensity<double>> cutline_features;
    void gather_cutline_and_features(bool verbose = false)
    {
        verbose_guard<debug_unclassified> vb(debug_unclassified::verbose() || verbose);
        auto result = gather_dense_feature_from_cutline(config, user_config_table, gg_table, clip, get_resist_cutline);
        edges = std::get<0>(result);
        cutline_features = std::get<1>(result);
        edge_meta = std::get<2>(result);
    }
    void calib_optical_threshold(bool verbose = false)
    {
        threshold = convert_to<double>(user_config_table["threshold_guess"]);
        verbose_guard<debug_unclassified> vb(debug_unclassified::verbose() || verbose);
        auto x = resist::calib_osqp_with_equality_constraints(gg_table, cutline_features, threshold);
        threshold = x.back();
        x.pop_back();
        post_calib_analysis(gg_table, edges, x, threshold, edge_meta.spatial.step, config.dbu);
        cutline_data::print(gg_table);
        using print_type = std::tuple<std::string, std::string>;
        std::vector<print_type> rows{
            print_type("threshold", to_string(threshold)),
            print_type("coefficients", to_string(x)),
            print_type("resist-format", to_string(vec3<double>{
                convert_to<double>(user_config_table["gauss_laguerre_sigma_in_dbu"]), 
                convert_to<double>(user_config_table["associate_order"]), 
                convert_to<double>(user_config_table["laguerre_order"])
            })),
        };
        print_table(rows, {"resist model", ""}, -1);
        auto [it_min, it_max] = std::minmax_element(gg_table.begin(), gg_table.end(), [](const cutline_data& a, const cutline_data& b){
            return a.post_calib_results.at(1) < b.post_calib_results.at(1); 
        });
        double start = it_min->post_calib_results.at(1);
        double end = it_max->post_calib_results.at(1);
        const size_t N = 100;
        double step = (end - start) / (N - 1);

        std::vector<double> errors(N, 0);
        for(const auto& data : gg_table){
            errors.at(0 != step ? std::floor((data.post_calib_results.at(1) - start)/step) : 0) += 1;
        }
        plot_curves(std::vector<std::vector<double>>{errors}, {float(start)},{float(step)}, {"error distribution"}, { "g-x"});  
    }
};

void regist_simulation_pyclass()
{
    py::class_<resis_simulation>("resis_simulation").def(py::init<>())       
        .def_readwrite("gauge_table",& resis_simulation::gg_table)
        .def("load_user_config", &resis_simulation::load_user_config)
        .def("load_gauge_file", &resis_simulation::load_gauge_file)
        .def("clip_cutline", &resis_simulation::clip_cutline)
        .def("gather_cutline_and_features", &resis_simulation::gather_cutline_and_features)
        .def("calib_optical_threshold", &resis_simulation::calib_optical_threshold)
    ;
}

void simulation_flow(const std::string& config_path)
{
    resis_simulation s;
    s.load_user_config(config_path);
    s.load_gauge_file();
    s.clip_cutline();
    s.gather_cutline_and_features();
    s.calib_optical_threshold(true);
}