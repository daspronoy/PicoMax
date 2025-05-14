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

    Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<Eigen::Vector3d> T, mater mat);
    std::complex<double> pseudopotential(Eigen::Vector3d G, std::vector<Eigen::Vector3d> T, 
                                        mater mat);
    void empiricalpseudopotentialmethod(env &dat);
    // Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<double> v);
    // std::complex<double> pseudopotential(Eigen::Vector3d Kvec, std::vector<double> v);
    std::complex<double> nonlocalpseudopotential(Eigen::Vector3d Kkvecm, Eigen::Vector3d Kkvecn, std::vector<double>v);
    double sphbesj(int n, double r);
    void setRefEnergy(env &dat);

} /* namespace pmx */

#endif /* _PMX_EPM_HPP_ */