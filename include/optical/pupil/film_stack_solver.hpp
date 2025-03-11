#pragma once
#include <type_traist_notebook/type_traist.hpp>

template<class T> inline T kz(T refraction_index, T kr){
    return T(2_PI) * std::sqrt(refraction_index * refraction_index - kr * kr);
}
template<class rT>struct film_stack_solver{
    using cT = complex_t<rT>;
    /**
     *              nTop (E,Er)
     * ---------------------- 
     *           ^ 
     *  distance | 
     *           |  nBottom (E0, Er0)
     * ---------------------- 
    */
    static std::array<cT, 2> thin_mask_transfer_E(cT E0, cT Er0, cT nTop, cT nBottom, rT distance, cT kr, rT lambda,cT omega = cT{1.0})
    {
        // std::cout << "input :" << std::make_tuple(E0, Er0) << std::make_tuple(nTop, nBottom) << std::endl;
        cT kz_bottom = kz<cT>(nBottom * omega, kr);
        E0 *= std::exp(cT(0, -1) * kz_bottom * distance / lambda);
        Er0 *= std::exp(cT(0, 1) * kz_bottom * distance / lambda);
        
        cT kz_top = kz<cT>(nTop * omega, kr);
        cT E = ((kz_top + kz_bottom) * E0 + (kz_top - kz_bottom) * Er0) / (kz_top * rT(2.0)); 
        cT Er = ((kz_top - kz_bottom) * E0 + (kz_top + kz_bottom) * Er0) / (kz_top * rT(2.0)); 
        // std::cout << "output : " << std::make_tuple(E, Er) << std::endl;
        return {E, Er};
    }
    static matrix2x3<cT> TE(std::array<cT, 2> d){auto [E, Er] = d; return TE(E, Er);}
    static matrix2x3<cT> TE(cT E, cT Er){
        matrix2x3<cT> m;
        m.at(0) = vec3<cT>{0, E, 0};
        m.at(1) = vec3<cT>{0, Er, 0};
        return m;
    }
    static std::array<cT, 2> thin_mask_transfer_B(cT B0, cT Br0, cT nTop, cT nBottom, rT distance, cT kr, rT lambda, cT omega = rT(1.0))
    {
        cT kz_bottom = kz(nBottom * omega, kr);
        B0 *= std::exp(cT(0, -1) * kz_bottom * distance / lambda);
        Br0 *= std::exp(cT(0, 1)*kz_bottom * distance / lambda);
        
        cT kz_top = kz(nTop * omega, kr);
        nBottom *= nBottom;
        nTop *= nTop;
        cT B = ((nBottom *kz_top + nTop * kz_bottom) * B0 +
            (nBottom *kz_top - nTop * kz_bottom) * Br0) / (nBottom * kz_top * rT(2));
        cT Br = ((nBottom *kz_top - nTop * kz_bottom) * B0 +
            (nBottom *kz_top + nTop * kz_bottom) * Br0) / (nBottom * kz_top * rT(2));  
        return {B, Br};
    }
    static matrix2x3<cT> TM(std::array<cT, 2> d, cT nTop, rT kr, rT omega = rT(1.0)){auto [B, Br] = d; return TM(B, Br, nTop, kr, omega);}
    static matrix2x3<cT> TM(cT B, cT Br, cT nTop, rT kr, rT omega = rT(1.0)){
        matrix2x3<cT> m{0};
        cT z = kz<cT>(nTop * omega, kr);
        omega *= rT(2_PI);
        m.at(0) = vec3<cT>{cT(z), cT(0), cT(-kr)};
        m.at(0) *= (B/(omega * nTop * nTop));
        m.at(1) = vec3<cT>{cT(-z), cT(0), cT(-kr)};
        m.at(1) *= (Br/(omega * nTop * nTop));
        return m;
    }

    static void normalize_transfer_matrix(std::vector<matrix2x3<cT>>& film_stack_matrix){
        auto& top = film_stack_matrix.front().at(0);
        rT N = std::sqrt(vector_norm(top));
        for(auto& m : film_stack_matrix) m /= N;
    }
    struct meterial{
        cT nk; rT depth;
        static std::ostream& print(const std::vector<meterial>& meterials)
        {
            debug_print<>(
                reinterpret_cast<const std::vector<std::tuple<rT, rT, rT>>&>(meterials),
                std::array<std::string, 3>{"thickness", "absorption coefficient", "refractive index"}
            );
            return debug_print<>::print_to();
        }
    };
    rT lambda, omega;
    film_stack_solver(rT l, rT w) : lambda(l), omega(w){};
    std::vector<matrix2x3<cT>> solveTE(const std::vector<meterial>& meterials, rT kr, cT E0 = cT(1), cT Er0 = cT(0))
    {
        std::vector<matrix2x3<cT>> m(meterials.size());
        m.back() = TE(E0, Er0);
        if(1 < meterials.size()){
            std::transform(meterials.rbegin() + 1, meterials.rend(), meterials.rbegin(), m.rbegin() +1, 
            [&](const meterial& top, const meterial& bottom){
                auto[E, Er] = thin_mask_transfer_E(E0, Er0, top.nk, bottom.nk, bottom.depth, kr, lambda, omega);
                E0 =E; Er0 = Er;
                return TE(E, Er);
            }
        ); 
        } 
        normalize_transfer_matrix(m);
        return m;
    }
    std::vector<matrix2x3<cT>> solveTM(const std::vector<meterial>& meterials, rT kr, cT B0 = cT(1), cT Br0 = cT(0))
    {
        std::vector<matrix2x3<cT>> m(meterials.size());
        m.back() = TM(B0, Br0, meterials.back().nk, kr, omega);
        if(1 < meterials.size()){
                std::transform(meterials.rbegin() + 1, meterials.rend(), meterials.rbegin(), m.rbegin() +1, 
                [&](const meterial& top, const meterial& bottom){
                    auto[B, Br] = thin_mask_transfer_B(B0, Br0, top.nk, bottom.nk, bottom.depth, kr, lambda, omega);
                    B0 =B; Br0 = Br;
                    return TM(B, Br, top.nk, kr, omega);
                }
            );
        }
        normalize_transfer_matrix(m);
        return m;
    }
    vec3<cT> solve_at_depth(const std::vector<meterial>& meterials, const std::vector<matrix2x3<cT>>& m, cT kr, rT depth)
    {
        int i = 0;
        for(; i < meterials.size(); i++){
            if(depth < meterials.at(i).depth){
                depth = meterials.at(i).depth - depth;      
                break;
            }
            depth -= meterials.at(i).depth;
        }
        const auto& transmition = m.at(i).at(0);
        const auto& reflection  = m.at(i).at(1);
        cT phase_term = kz(meterials.at(i).nk * omega, kr) * depth / lambda;
        if(i == meterials.size() - 1){
            debug_print("depth is inside of substrate. %f\n", depth);
            return transmition * std::exp(cT(0, 1) * phase_term);

        }else{
            return transmition * std::exp(cT(0, -1) * phase_term) + reflection  * std::exp(cT(0, 1) * phase_term);
        }
    }
    vec3<cT> solve_field(const vec3<cT>& in, const vec3<cT>& e_top, const vec3<cT>& e_at_depth, const vec3<cT>& b_top,const vec3<cT>& b_at_depth)
    {
        return in * conj(e_top) * e_at_depth + in * conj(b_top) * b_at_depth;
    }
};
