#include <py_helper.hpp>
#include <optical/optical_numerics.hpp>
#include <type_traist_notebook/type_traist.hpp>
#include <filesystem>
#include <optical/geometry.hpp>
#include <optical/clip.hpp>

void simulation_flow();
int main()
{
    py_engine::init();
    add_path_to_sys_path("core_plugins");
    simulation_flow();
    py_engine::dispose();
}

std::vector<cutlinei> load_gauge_file(const std::string& path)
{
    // read_gauge_file
    py::list_tuple gauge_table = py::extract<py::list_tuple>(py_plugin::ref()["gauge_io"]["read_gauge_file"](path));
    np::array2df cutlines_in_um = py::extract<np::array2df>(py_plugin::ref()["gauge_io"]["get_culine_from_gauge_table"](gauge_table, 1));
    auto [pRec, n] = ndarray_ref_no_padding<vec<double, 4>>(cutlines_in_um);
    std::vector<cutlinei> cutlines; 
    cutlines.reserve(n);
    for(size_t i = 0; i < n; i++, pRec++){
        const auto& rec = *pRec;
        cutlines.push_back(cutlinei{
            static_cast<int64_t>(rec[0]), static_cast<int64_t>(rec[1]), 
            static_cast<int64_t>(rec[2]), static_cast<int64_t>(rec[3])}
        );
    }
    return cutlines;
}

template<class TUserConfig> grid_start_step<double> optical_numerics_in_dbu(cutlinei cutline, const TUserConfig& config)
{
    auto roi = convert<cutlinei, rectangle<double>>{}(cutline);
    dbu_to_um(roi, config.dbu);
    auto grid_in_um = optical_numerics<double>(roi, config.ambit, config.tilesize, config.maxNA, config.wavelength); 
    // grid_in_um.print();

    //== 1. spatial step
    um_to_dbu(grid_in_um.spatial.step, config.dbu);
    vec2<int64_t> spatial_step_in_dbu = {int64_t(std::floor(grid_in_um.spatial.step[0])), int64_t(std::floor(grid_in_um.spatial.step[1]))};
    //== 2. spatial domain & start
    auto tile_size = convert<vec2<size_t>, vec2<int64_t>>{}(config.tilesize);
    vec2<int64_t> spatial_domain_in_dbu = tile_size * spatial_step_in_dbu;
    const auto& [from_in_dbu, to_in_dbu] = cutline;
    vec2<int64_t> spatial_start_in_dbu = ((from_in_dbu + to_in_dbu - spatial_step_in_dbu) + 1) / 2; 
    
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
    grid_start_step<double> grid_in_dbu;
    grid_in_dbu.tilesize = config.tilesize;
    grid_in_dbu.spatial.start = convert<vec2<int64_t>, vec2<double>>{}(spatial_start_in_dbu);
    grid_in_dbu.spatial.step  = convert<vec2<int64_t>, vec2<double>>{}(spatial_step_in_dbu); 
    double lambda_in_dbu = config.wavelength; 
    um_to_dbu(lambda_in_dbu, config.dbu);
    grid_in_dbu.fourier.step  = vec2<double>{lambda_in_dbu, lambda_in_dbu} / spatial_domain_in_dbu;
    grid_in_dbu.fourier.start = vec2<double>{0, 0}; 
    // grid_in_dbu.print();
    return grid_in_dbu;
}

void simulation_flow()
{
    auto user_config = cutline_jobs::get_user_config("/home/like/repos/simulation/config/calib_193.py");
    //== load gauge file & calc startstep
    auto cutlines = load_gauge_file(user_config.gauge_file);
    auto startstep_in_dbu = optical_numerics_in_dbu(cutlines.at(0), user_config);
    //== cutline subclip
    auto shape =  (startstep_in_dbu.spatial.step * startstep_in_dbu.tilesize) * user_config.dbu;
    std::cout << shape << std::endl;
    cutline_jobs jobs = cutline_jobs::cutline_clip_flow(user_config, shape);
    auto path = jobs.clip_path(0);

    //== load subclip
    py::tuple poly_and_holes = py_plugin::call<py::tuple>("klayout_op", "load_oas_vertexs", 
        jobs.clip_path(0).c_str(), user_config.cell_name, user_config.layer_id
    );
    np::list_array2d polys = py::extract<np::list_array2d>(poly_and_holes[0]);
    np::list_array2d holes = py::extract<np::list_array2d>(poly_and_holes[1]);
    foreach_shape_lines(polys, [](point_dbu from, point_dbu to){
        auto unit = unit_vector(from, to);
        auto norm = norm_vector(from, to);
        std::cout << from << to << point_in_domain_clamp({0, 0}, from, to)  << std::endl;
    });
    

    // std::cout << len(polys) << std::endl;
    // std::cout << len(holes) << std::endl;

    // np::array2df poly = py::extract<np::array2df>(polys[0]);
    // auto [pPoly, poly_point_size] = ndarray_ref_no_padding<point_dbu>(poly);
    // std::cout << "[";
    // for(size_t i = 0; i < poly_point_size; i++)
    //     std::cout << pPoly[i] << vec2<const char*>({"\n", "]\n"}).at((i + 1) == poly_point_size);

        
    // np::array2df hole = py::extract<np::array2df>(holes[0]);
    // auto [phole, hole_point_size] = ndarray_ref_no_padding<point_dbu>(hole);
    // std::cout << "[";
    // for(size_t i = 0; i < hole_point_size; i++)
    //     std::cout << phole[i] << vec2<const char*>({"\n", "]\n"}).at((i + 1) == hole_point_size);
}