#pragma once
#include <type_traist_notebook/type_traist.hpp>
#include <variant>
#include <assert.h>
#include <optical/simulation_grid_info.hpp>
namespace litho_lab_f
{
    using rT = float;
    using cT = complex_t<rT>;
    using scalar_field_r = std::vector<rT>;
    using scalar_field_c = std::vector<cT>;
    using vec_field2_r = std::vector<vec2<rT>>;
    using vec_field2_c = std::vector<vec2<cT>>;
    using vec_field3_r = std::vector<vec3<rT>>;

    struct polarization
    {
        virtual vec_field2_c get_field_xy(const scalar_field_r&) const;
        virtual vec_field2_c get_field_xy(const scalar_field_c&) const;
        virtual vec_field2_c get_field_xy(const vec_field2_r&) const;
        virtual vec_field2_c get_field_xy(const vec_field2_c&) const;
    };
    struct global_polarization : polarization
    {
        vec2<cT> jones_vector;
        template<class T> vec_field2_c get_field_xy_impl(const std::vector<T>& f) const
        {
            vec_field2_c n(f.size());
            std::transform(f.begin(), f.end(), n.begin(),[&](T r){return jones_vector * r;});
            return n;
        }
        vec_field2_c get_field_xy(const scalar_field_r& f) const override{return get_field_xy_impl(f);}
        vec_field2_c get_field_xy(const scalar_field_c& f) const override{return get_field_xy_impl(f);}
        vec_field2_c get_field_xy(const vec_field2_r& f) const override{return get_field_xy_impl(f);}
        vec_field2_c get_field_xy(const vec_field2_c& f) const override{return get_field_xy_impl(f);}
    };
    struct local_polarization : polarization
    {
        vec_field2_c jones_vector_field;
        template<class T> vec_field2_c get_field_xy_impl(const std::vector<T>& f) const
        {
            assert(jones_vector.size() == f.size());
            return jones_vector_field * f;
        }
        vec_field2_c get_field_xy(const scalar_field_r& f) const override{return get_field_xy_impl(f);}
        vec_field2_c get_field_xy(const scalar_field_c& f) const override{return get_field_xy_impl(f);}
        vec_field2_c get_field_xy(const vec_field2_r& f) const override{return get_field_xy_impl(f);}
        vec_field2_c get_field_xy(const vec_field2_c& f) const override{return get_field_xy_impl(f);}
    };
    
    struct field_data
    {
        std::variant<scalar_field_r, scalar_field_c, vec_field2_r, vec_field2_c, vec_field3_r> d;
        grid_start_step<rT, 3> meta{0};
    };
    struct phase_radius
    {
        scalar_field_r phase;
    };
    struct gauss_wave
    {
        scalar_field_r depth;
        scalar_field_r value;
        grid_start_step<rT, 3> meta{0};
        
    };
};