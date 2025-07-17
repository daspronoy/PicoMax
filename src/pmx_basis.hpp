/*
    pmx_basis.cpp | includes basis-related functinalities of picomax

    * kpoint generation for BZ integration
    * gpoint generation for planewave basis set
    * high-symmetry point
    * transverse vector

    Author: Jungho Mun
    Date: June 6, 2024
*/

#pragma once
#ifndef _PMX_BASIS_HPP_
#define _PMX_BASIS_HPP_

/*
    header files
*/
#include <Eigen/Core>
#include "pmx.hpp"

/*
    function declaration
*/
namespace pmx {

    // generate a gpoint grid
    void generate_gpt (lattice &lat);
    void generate_gpt_rect (lattice &lat);
    void generate_loci (lattice &lat);

    // import kpoint file
    void import_kpt (env &dat);

    // generate a kpoint grid inside the 1st BZ
    void generate_kpt_bz (lattice &lat);
    void generate_kpt_fcc (lattice &lat);
    void generate_kpt_hex (lattice &lat);
    void generate_kpt_hex2 (lattice &lat);
    void generate_kpt_fcc2 (lattice &lat);

    void generate_kpt_cohen_chadi (lattice &lat);
    inline bool kroneckerDelta(int m, int n);
    void generate_kpt_monkhorst_pack (env &dat);

    // check if kpoint is inside the 1st BZ
    bool checkinbz_fcc (Eigen::Vector3d K);
    bool checkinbz_hex(Eigen::Vector3d K, double f);
    bool checkinbz_hex2(Eigen::Vector3d K);
    bool checkonbz_fcc (Eigen::Vector3d K);

    // high-symmetry points
    Eigen::Vector3d highsympoint_fcc (std::string s);
    Eigen::Vector3d highsympoint_hex (std::string s, double f);
    Eigen::Vector3d highsympoint_hex2d (std::string s);

    // polarization unit vectors
    inline Eigen::Vector3d uvec_L (Eigen::Vector3d v);
    inline Eigen::Vector3d uvec_T (Eigen::Vector3d v);
    std::vector<Eigen::Vector3cd> uvec_LT (Eigen::Vector3d v);
    std::vector<Eigen::Vector3cd> ovec_LT (Eigen::Vector3cd v1, Eigen::Vector3cd v2);

    // Eigen::Vector3d normvec(Eigen::Vector3d v);
    // Eigen::Vector3d perpvec(Eigen::Vector3d v);
    // Eigen::Vector3d tanvec(Eigen::Vector3d v);
    // Eigen::Vector3d tanvecq(Eigen::Vector3d v, Eigen::Vector3d q);
    // std::vector<Eigen::Vector3cd> unitvec (Eigen::Vector3d nvec);
    // std::vector<Eigen::Vector3cd> unitvecq (Eigen::Vector3d nvec, Eigen::Vector3d qvec);
    // Eigen::Vector3cd **uvecqm (Eigen::Vector3d Q, std::vector<Eigen::Vector3d> G);

} /* namespace pmx */


#endif /* _PMX_BASIS_HPP_ */