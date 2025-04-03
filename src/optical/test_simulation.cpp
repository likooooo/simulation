#include <py_helper.hpp>
#include <optical/optical_numerics.hpp>
#include <type_traist_notebook/type_traist.hpp>
#include <filesystem>
#include <optical/geometry.hpp>
#include <optical/clip.hpp>
#include <optical/near_field/thin_mask.hpp>
// #include <cpp_cuda/cuda_vector.hpp>
#include <type_traist_notebook/uca/backend.hpp>
#include <mekil/cpu_backend.hpp>
#include <cpp_cuda/gpu_backend.hpp>

bool verbose = false;

void simulation_flow(const std::string&);
int main()
{
    py_engine::init();
    add_path_to_sys_path("core_plugins");
    // catch_py_error
    (simulation_flow("/home/like/repos/simulation/config/calib_193.py"));
    py_engine::dispose();
}
BOOST_PYTHON_MODULE(lib_test_simulation) {
    py_engine::init();
    py_engine::init_exception_for_pycall();
    py::def("simulation", &simulation_flow);
}

std::vector<cutline_dbu> load_gauge_file(const std::string& path)
{
    // read_gauge_file
    py::list_tuple gauge_table = convert_to<py::list_tuple>(py_plugin::ref()["gauge_io"]["read_gauge_file"](path));
    np::array2df cutlines_in_um = convert_to<np::array2df>(py_plugin::ref()["gauge_io"]["get_culine_from_gauge_table"](gauge_table, 1));
    auto [pRec, n] = ndarray_ref_no_padding<vec<double, 4>>(cutlines_in_um);
    std::vector<cutline_dbu> cutlines; 
    cutlines.reserve(n);
    for(size_t i = 0; i < n; i++, pRec++){
        const auto& rec = *pRec;
        cutlines.push_back(cutline_dbu{
            static_cast<int64_t>(rec[0]), static_cast<int64_t>(rec[1]), 
            static_cast<int64_t>(rec[2]), static_cast<int64_t>(rec[3])}
        );
    }
    return cutlines;
}

template<class T> struct dbu_grid_start_step
{
    using value_type = int64_t;
    constexpr static size_t dim = 2;
    template<size_t DIM1> using rebind_t = grid_start_step<int64_t, DIM1>; 
    struct StartStep{point_dbu start, step;}spatial;
    typename grid_start_step<T, dim>::StartStep fourier;
    vec<size_t, dim> tilesize;
};

template<class TUserConfig> dbu_grid_start_step<double> optical_numerics_in_dbu(cutline_dbu cutline, const TUserConfig& config)
{
    auto roi = convert<cutline_dbu, rectangle<double>>{}(cutline);
    dbu_to_um(roi, config.dbu);
    auto grid_in_um = optical_numerics<double>(roi, config.ambit, config.tilesize, config.maxNA, config.wavelength); 
    print_grid_start_step(grid_in_um, "origin grid-um");

    //== 1. spatial step
    um_to_dbu(grid_in_um.spatial.step, config.dbu);
    vec2<int64_t> spatial_step_in_dbu = {int64_t(std::floor(grid_in_um.spatial.step[0])), int64_t(std::floor(grid_in_um.spatial.step[1]))};
    //== 2. spatial domain & start
    auto tile_size = convert<vec2<size_t>, vec2<int64_t>>{}(config.tilesize);
    vec2<int64_t> spatial_domain_in_dbu = tile_size * spatial_step_in_dbu;
    const auto& [from_in_dbu, to_in_dbu] = cutline;
    vec2<int64_t> spatial_start_in_dbu = ((from_in_dbu + to_in_dbu - spatial_domain_in_dbu) + 1) / 2; 
    
    //== check tilesize
    vec2<double> ambit = config.ambit; um_to_dbu(ambit, config.dbu);
    vec2<int64_t> ambit_in_dbu = {int64_t(std::ceil(ambit[0])), int64_t(std::ceil(ambit[1]))};
    vec2<int64_t> span_in_dbu = to_in_dbu - from_in_dbu;
    int64_t length_in_dbu = std::max(span_in_dbu[0], span_in_dbu[1]);
    span_in_dbu = {length_in_dbu, length_in_dbu};
    if((spatial_domain_in_dbu - (ambit_in_dbu * 2)) < span_in_dbu)
    {
        print_table(std::cerr, std::vector<std::tuple<std::string, vec2<int64_t>>>{
            std::make_tuple(std::string("roi from       :"), from_in_dbu), 
            std::make_tuple(std::string("roi to         :"), to_in_dbu),
            std::make_tuple(std::string("ambit          :"), ambit_in_dbu),
            std::make_tuple(std::string("spatial domain :"), spatial_domain_in_dbu),
            std::make_tuple(std::string("need more tile :"), ((to_in_dbu - from_in_dbu) - (spatial_domain_in_dbu - (ambit_in_dbu * 2))) / spatial_step_in_dbu),
        }, {"* roi or ambit is too large", ""});
    }
    dbu_grid_start_step<double> grid_in_dbu;
    // grid_start_step<double> grid_in_dbu;
    grid_in_dbu.tilesize = config.tilesize;
    grid_in_dbu.spatial.start = spatial_start_in_dbu;
    grid_in_dbu.spatial.step  = spatial_step_in_dbu; 
    double lambda_in_dbu      = config.wavelength; 
    um_to_dbu(lambda_in_dbu, config.dbu);
    grid_in_dbu.fourier.step  = vec2<double>{lambda_in_dbu, lambda_in_dbu} / spatial_domain_in_dbu;
    grid_in_dbu.fourier.start = vec2<double>{0, 0}; 
    
    if(debug_unclassified::verbose())
    {
        grid_start_step<double> grid_to_um;
        grid_to_um.tilesize      = grid_in_dbu.tilesize; 
        grid_to_um.fourier       = grid_in_dbu.fourier;
        grid_to_um.spatial.start = convert_to<vec2<double>>(grid_in_dbu.spatial.start); 
        grid_to_um.spatial.step  = convert_to<vec2<double>>(grid_in_dbu.spatial.step); 
        dbu_to_um(grid_to_um.spatial.start, config.dbu);
        dbu_to_um(grid_to_um.spatial.step, config.dbu);
        print_grid_start_step(grid_to_um, "grid-dbu-aligined to grid-um");
    }
    return grid_in_dbu;
}


void simulation_flow(const std::string& config_path)
{
    // const uca::backend<double>& backend = uca::cpu<double>::ref();
    const uca::backend<double>& backend = uca::cpu<double>::ref();

    auto [user_config, params] = cutline_jobs::get_user_config(config_path);
    //== load gauge file & calc startstep
    auto cutlines = load_gauge_file(user_config.gauge_file);
    auto startstep_in_dbu = optical_numerics_in_dbu(cutlines.at(0), user_config);
   
    //== cutline subclip
    auto shape = convert_to<vec2<double>>((startstep_in_dbu.spatial.step * startstep_in_dbu.tilesize)) * user_config.dbu;
    cutline_jobs jobs = cutline_jobs::cutline_clip_flow(user_config, shape);
    
    //== load subclip
    using thin_mask_solver = thin_mask<double, dbu_grid_start_step<double>>;
    shapes_dbu shapes = near_filed::load_shapes_from_file(jobs.clip_path(0).c_str(), user_config.cell_name, user_config.layer_id);
    debug_print<thin_mask_solver>::verbose() = -1 < user_config.verbose;
    auto [x, y, mask_info] = thin_mask_solver::edge_pixelization(startstep_in_dbu, shapes, convert_to<size_t>(params["mask_USF"]),  convert_to<double>(params["mask_edge_dissect_coef"]));
    // imshow(x, convert_to<std::vector<size_t>>(mask_info.tilesize));
    backend.integral_y(mask_info.tilesize, x.data());
    // imshow(x, convert_to<std::vector<size_t>>(mask_info.tilesize));

    // imshow(y, convert_to<std::vector<size_t>>(mask_info.tilesize));
    backend.integral_x(mask_info.tilesize, y.data());
    // imshow(y, convert_to<std::vector<size_t>>(mask_info.tilesize));

    backend.VtAdd(x.size(), x.data(), y.data());

    //== gpu backend
    // cuda::device_vector<double> cx, cy; cx << x; cy << y;
    // uca::gpu<double>::ref().VtAdd(x.size(), cx.data(), cy.data());
    // y <<cy;


    //== compare dissect coef
    // auto [x1, y1, mask_info1] = thin_mask_solver::edge_pixelization(startstep_in_dbu, shapes, 0.5);
    // uca::cpu<double>::ref().VtAdd(x1.size(), x1.data(), y1.data());
    // y -= y1;
    // std::cout << *std::max_element(y.begin(), y.end()) << std::endl;
    imshow(y, convert_to<std::vector<size_t>>(mask_info.tilesize));

    auto [cutline, cutline_meta] = thin_mask_solver::get_edge_from_rasterization(mask_info, y, cutlines.at(0));
    {
        auto start = dbu_to_um(convert_to<vec2<float>>(cutline_meta.spatial.start), user_config.dbu);
        auto step = dbu_to_um(convert_to<vec2<float>>(cutline_meta.spatial.step), user_config.dbu);
        // plot_curves(std::vector<std::vector<double>>{cutline}, {start[0]}, {float(step[0])}, {"cutline (um)"}, {"b--"});  
        plot_curves(std::vector<std::vector<double>>{cutline}, {0}, {float(1)}, {"cutline (pixel)"}, {"b--"});  

    }
}