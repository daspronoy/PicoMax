#pragma once
#ifndef PMX_H_
#define PMX_H_


/*
    libraries
*/
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <iomanip> //For controlling precision in display
#include <omp.h>
// #include "stdafx.h"
#include <stdlib.h>
// #include "interpolation.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/LU>
#include <Eigen/Eigenvalues> 
// #include <unsupported/Eigen/Splines>
#include <Eigen/Dense>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
// #define EIGEN_USE_MKL_ALL



/*
	Physical constants:
	m_o*c^2 (eV), hbar_c (ev*m), xscale:1 Angstrom, vac_epsilon, rydberg, echarge, cscale, xl0
*/
#define pi 3.14159265358979323846
#define hbar    1.054571817e-27 // [erg*s] (CGS)
#define e       4.803204712570263e-10 // [esu] (CGS)
#define e_m     9.1093837015e-28 // [g] (CGS)
#define eV      1.602176634e-12 // [erg] (CGS)
#define Ry      2.179872361103589e-11 // [erg] (CGS)
#define Ry2eV   13.605699169854621 // Ry/eV
#define angstrom  1e-8 // [cm] (CGS)
#define KCON    3.809982111485962 // hbar^2/(2*m_e*eV*angstrom^2)

/*
    global variables (defined in io.cpp)
*/
extern int NPW;
extern int MAXPW;
extern double a;
extern int NEPS;
extern int NBAND;
extern int NKPT;
extern int NFREQ;
extern int NQ;
extern int DELTA;
extern double EPSILON;
extern int ORIGINSHIFT;
extern bool EPM_NONLOCAL;
extern int NVALENCE;



namespace pmx {

// environment variables
struct env {
    
    // general options
    int debug = 0; // debug level
    char *outputfile; // filename of outputfile
    char *wdir; // working directory
    char *kpointfile; // filename of kpointfile
    int logfile = 0; // whether to save the logfile

    // planewave grids
    double encut = 200; // plane-wave energy cutoff [eV]
    double gcut = 0; // plane-wave cutoff [2*pi/a]; determined from encut, ENCUT = KCON*(2*pi/a)^2*GCUT^2
    std::vector<Eigen::Vector3d> G; // G-vector grids [2*pi/a]
    int **loci; // local-field plane-wave index
    int kpoint = 0; // kpoint basic option; 0: full BZ, 1: Monkhorst-Pack grids, 2: Cohen-Chadi grids
    int kpointorder = 3; //
    double **K; // k-vector grids [2*pi/a]
    double *KW; // kpoint weights

    // frequency
    double *freq;// frequency [eV]
    double dfreq = 0;// frequency interval [eV]
    int kk = 0;// use kramers-kronig transformation to obtain the real parts of dielectric constant
    double epsilon = 0.05; // delta function epsilon
    // int delta = 0; // delta function model; 0: lorentzian, 1: gaussian

    // q-vector
    std::vector<Eigen::Vector3d> Q;// q-vector grids [2*pi/a]
    std::vector<int> qnum;// number of q-vectors in each path
    std::vector<std::string> qpath;// path of q-vector
    Eigen::Vector3d gammazero = {1e-6,0,0}; // shift for the gamma-point to avoid singularity

    // material parameter
    std::vector<double> epm;// empirical pseudopotential parameters [rydberg] {vs3,vs8,vs11,va3,va4,va11,a0,b0,A2,rl0,rl2}















    // double GCUT = 3.5;
    // int dim = 3;//
    // int bandscheme = 0; // 0: empirical pseudopotential method, 1: epm with nonlocal
    // char *OUTPUT_FILE;
    // int KPOINT_SCHEME = 0; // 0: Full-BZ; 1: Monkhorst-Pack; 2: Cohen-Chadi; 
    // int KPOINT_ORDER = 3; // Cohen-Chadi grid order
    // // kgrid and kweight [k]
    // int gridscheme = 0;
    // int gridorder = 3;
    // char *kgridfile;
    // // reciprocal lattice vectors {b1, b2, b3} in units of 2*pi/a
    // std::vector<Eigen::Vector3d> b {
    //         Eigen::Vector3d(3), 
    //         Eigen::Vector3d(3), 
    //         Eigen::Vector3d(3)
    //     };

    // q-vector
    
    // int nq;// total number of q-vectors

    // material parameters
    // double a;// lattice constant [angstroms]
    

    double energyoffset;
    double eigerror;

    // lattice permittivity matrix [q][m][n][f]
    double ****realepsL;
    double ****imagepsL;
    double ****realepsT;
    double ****imagepsT;

    // electronic bandstructure [q][n]
    double **band;
    
    // // polarization?
    // double **realpolL;
    // double **imagpolL;
    // double **realpolT;
    // double **imagpolT;


    

    // //recycle bin
    // double **points;
    // double *weights;
    // Eigen::Vector3d qvector = Eigen::Vector3d(3); // q-vector for dielectric function calculations
    // int numomega = 1; // number of freq points
    // double omega;// frequency [eV]
    // double vs3, vs4, vs8, vs11, va3, va4, va8, va11;// pseudopotentials in rydberg units
    // double a0, b0, A2, rl0, rl2;// nonlocal psuedopotential parameters in rydberg units

};

// parse the input option strings
class inputParser {
    public:

        // constructor
        inputParser (int &argc, char **argv){
            for (int i=1; i<argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }

        // check if the option is saved
        bool existOption(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }

        // return the value of the option
        const std::string& valueOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr = std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }

        // return vector of a comma-separated option value
        const std::vector<std::string> vectorOption(const std::string &option) const{
            std::string str = valueOption(option);
            std::stringstream ss(str);
            std::vector<std::string> v;
            while (ss.good()){
                std::string substr;
                getline(ss,substr,',');
                v.push_back(substr);
            }
            return v;
        }

    private:
        std::vector<std::string> tokens;
};

/*
    header files
*/
// io.cpp
void doc();
void ver();
void progmsg();
void inputProcess(inputParser inp, env &dat);
void printInput(env &dat);
void printBand(env &dat);
void printEpsilon(env &dat);
void writeEpsilon(env &dat);

// basis.cpp
void generate_planewaves (env &dat);
std::vector<Eigen::Vector3d> find_G (int g2, int n);
void find_locGi (env &dat);
void importkpointfile (env &dat);
void genkpoint_cohen_chadi (env &dat);
inline bool kroneckerDelta(int m, int n);
void genkpoint_monkhorst_pack (env &dat);
void importkpointfile (env &dat);

// epm.cpp
void empiricalpseudopotentialmethod(env &dat);
Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<double> v);
std::complex<double> pseudopotential(Eigen::Vector3d Kvec, std::vector<double> v);
std::complex<double> nonlocalpseudopotential(Eigen::Vector3d Kkvecm, Eigen::Vector3d Kkvecn, std::vector<double>v);
double sphbesj(int n, double r);
void setRefEnergy(env &dat);

// chi.cpp
void permittivity(env &dat);
Eigen::Vector3d perpvec(Eigen::Vector3d v);
Eigen::Vector3d tanvec(Eigen::Vector3d v, Eigen::Vector3d q);
Eigen::Vector3d normvec(Eigen::Vector3d v);
extern std::complex<double> (*diracdelta) (std::complex<double> value);
inline std::complex<double> diracdelta_lorentzian(std::complex<double> x);
inline std::complex<double> diracdelta_gaussian(std::complex<double> x);
void kramerskronigtransform(double *ReX, double *ImX, double *w, double dw);
} /* namespace pmx */
#endif /* PMX_H_ */