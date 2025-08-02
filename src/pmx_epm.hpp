#pragma once
#ifndef _PMX_EPM_HPP_
#define _PMX_EPM_HPP_


/*
    header files
*/
#include "pmx.hpp"


/*
    function declaration
*/
namespace pmx {

    double get_soc_form_factor(const std::map<double, double>& so_map, double q_sq_norm_dimless);
    Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<Eigen::Vector3d> T, pmx::mater mat);
    std::complex<double> pseudopotential(Eigen::Vector3d G, std::vector<Eigen::Vector3d> T, 
                                        pmx::mater mat);    void empiricalpseudopotentialmethod(pmx::env &dat);
    // Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<double> v);
    // std::complex<double> pseudopotential(Eigen::Vector3d Kvec, std::vector<double> v);
    std::complex<double> nonlocalpseudopotential(Eigen::Vector3d Kkvecm, Eigen::Vector3d Kkvecn, std::vector<double>v);
    double sphbesj(int n, double r);
    void setRefEnergy(pmx::env &dat);
    void calculateSpinProjection(pmx::env &dat);
    double calculate_spin_z(const Eigen::VectorXcd& eigenvector);

} /* namespace pmx */

#endif /* _PMX_EPM_HPP_ */
