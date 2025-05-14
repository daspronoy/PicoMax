/*
    pmx_tbg.cpp | 

    Author: Jungho Mun
    Date: November 5, 2024
*/

/*
    header files
*/
#include "pmx_tbg.hpp"


/*
    function definitions
*/
namespace pmx {

void tbg_eband(env &dat){

    tbg t;
    t.theta = dat.tbgv[0];
    // double n = 31;
    // t.theta = acos((3*pow(n,2)+3*n+0.5)/(3*pow(n,2)+3*n+1));
    t.u = dat.tbgv[1];
    t.up = dat.tbgv[2];
    t.s = dat.tbgv[3];
    t.ngmax = dat.tbgv[4];
    t.kcut = dat.tbgv[5];

    t.initialize();

    dat.eband = new double*[NQ];
    // NBAND = 2*NPW;

    // #pragma omp parallel for
    for (int q=0; q<NQ; q++){
        // Eigen::Vector3d qq = dat.lat.Q[q];
        double qx = dat.lat.Q[q](0);
        double qy = dat.lat.Q[q](1);


        // Hamiltonian
        Eigen::MatrixXcd H = t.gen_hamiltonian(qx,qy);

        // diagonalize the Hamiltonian
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver(H);
        if (eigsolver.info() != Eigen::Success){ std::cout << "ERROR: Eigensolver failed to solve the matrix." << std::endl; abort(); }

        // truncate to NBAND and apply the energy offset
        dat.eband[q] = new double[NBAND];
        for (int n=0; n<NBAND; n++){
            dat.eband[q][n] = double(eigsolver.eigenvalues()(NPW-dat.nvalence+n));
        }
    }
    // #pragma omp barrier

    return;
}





}/* namespace pmx */