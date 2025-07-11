#pragma once
#include <type_traist_notebook/uca/backend.hpp>
#include <optical/near_field/thin_mask.hpp>
#include <optical/resist/resist_common.hpp>
#ifdef CPU_BACKEND_ENABLE
#   include <mekil/cpu_backend.hpp>
#endif
#ifdef GPU_BACKEND_ENABLE
#   include <cpp_cuda/gpu_backend.hpp>
#endif

template<class TGetherCutline> inline std::tuple<std::vector<terms_cutline<double>>, std::vector<terms_features_intensity<double>>, grid_info_in_dbu> gather_feature_from_cutline(
    const user_config& config, py::dict user_config_table, const std::vector<cutline_data>& gg_table, const clip_data& clip, TGetherCutline&& callback)
{
    std::vector<terms_cutline<double>> edges;
    std::vector<terms_features_intensity<double>> cutline_features;
    edges.reserve(gg_table.size());
    cutline_features.reserve(gg_table.size());
    grid_info_in_dbu edge_meta;
    for(size_t i = 0; i < gg_table.size(); i++)
    {
        if(gg_table.at(i).weight == 0) continue;
        const auto& [terms, meta] = callback(clip.clip_path(i), gg_table.at(i), config, user_config_table);
        edge_meta = meta;
        edges.push_back(terms);
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
template<class TGetherCutline> inline std::tuple<std::vector<terms_cutline<double>>, std::vector<terms_dense_features_intensity<double>>, grid_info_in_dbu> gather_dense_feature_from_cutline(
    const user_config& config, py::dict user_config_table, const std::vector<cutline_data>& gg_table, const clip_data& clip, TGetherCutline&& callback)
{
    std::vector<terms_cutline<double>> edges;
    std::vector<terms_dense_features_intensity<double>> cutline_features;
    edges.reserve(gg_table.size());
    cutline_features.reserve(gg_table.size());
    grid_info_in_dbu edge_meta;
    for(size_t i = 0; i < gg_table.size(); i++)
    {
        if(gg_table.at(i).weight == 0) continue;
        const auto& [terms, meta] = callback(clip.clip_path(i), gg_table.at(i), config, user_config_table);
        edge_meta = meta;
        edges.push_back(terms);
        auto [start, step] = meta.spatial;
        assert(gg_table.at(i).cutline[0] == start);

        auto features = get_dense_feature_intensity_from_cutline<double>(
            gg_table.at(i).cutline, 
            gg_table.at(i).measured_cd,
            edges.back(), start, step, config.dbu
        );
        cutline_features.push_back(features);
        debug_unclassified::out("    cutline dense-features is ", features);
    }
    return {edges, cutline_features, edge_meta};
}
template<class T, class TDisplay> void check_cutline(const std::vector<T>& edge_image, const grid_info_in_dbu& cutline_meta, 
    const std::string& oas_path, const cutline_data& data, const user_config& user_config, const py::object params_optional,
    TDisplay&& display)
{
    using print_type = std::tuple<std::string, std::string, vec2<double>, double, std::vector<double>>;
    std::vector<print_type> rows{
        print_type(
            std::filesystem::path(oas_path).filename(), 
            data.pattern_name, 
            dbu_to_physical(convert_to<vec2<double>>(cutline_meta.spatial.start), user_config.dbu), 
            data.weight,
            edge_image 
        )
    };
    debug_unclassified(rows, {"path", "pattern-name", "start(um)", "weight", "intensity"}, 100);
    auto cutline_debug_list = convert_to<std::vector<int>>(params_optional["cutline_debug_list"]);
    if(std::find(cutline_debug_list.begin(), cutline_debug_list.end(), std::stoi(std::filesystem::path(oas_path).stem().string())) != cutline_debug_list.end())
    {
        display(data, edge_image, cutline_meta.spatial.start, cutline_meta.spatial.step, user_config.dbu);
    }
}

struct simulation_common{
    user_config config;
    py::dict user_config_table;
    std::vector<cutline_data> gg_table;
    clip_data clip;
    std::vector<terms_cutline<double>> edges;
    std::vector<terms_features_intensity<double>> cutline_features;
    grid_info_in_dbu edge_meta;

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
};
