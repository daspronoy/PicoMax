#pragma once
#ifndef _PMX_H_
#define _PMX_H_

/*
    pmx.hpp | main header file for picomax

    Author: Jungho Mun
    Date: June 6, 2024
*/


/*
    libraries
*/
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <complex>
#include <chrono>
#include <ctime>
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
#include <thread>

/*
    header files
*/
// #include "pmx_io.hpp"
// #include "pmx_basis.hpp"
// #include "pmx_epm.hpp"
// #include "pmx_chi.hpp"


// #include <hdf5.h>
// #include "H5Cpp.h"
// 
// #include <highfive/highfive.hpp>

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
    global variables (defined in pmx.cpp)
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
extern int NTSR;
extern int NATOM;
extern std::complex<double> im;



namespace pmx {

// crystal lattice parameters
struct lattice {
    int dim = 3;// dimension (2 or 3)
    std::string name;// {zb,wz}
    std::vector<Eigen::Vector3d> lattice;// lattice vectors in units of [a]
    std::vector<Eigen::Vector3d> reciprocal;// reciprocal lattice vectors in units of [2*pi/a]
    std::vector<Eigen::Vector3d> atomic;// atomic basis vectors in units of [a]
    Eigen::Vector3d offset;// 
    
    std::vector<Eigen::Vector3d> G;// g-vector
    std::vector<Eigen::Vector3d> K;// k-vector
    std::vector<double> KW;// k-vector weights
    int **loci; // local-field plane-wave index

    double a; // lattice constant [angstrom]
    double f = sqrt(8.0/3.0); // (c/a) fraction
    double u = 3.0/8.0;

    // planewave grids
    double encut = 200; // plane-wave energy cutoff [eV]
    double gcut = 0; // plane-wave cutoff [2*pi/a]; determined from encut, ENCUT = KCON*(2*pi/a)^2*GCUT^2
    int gsym = 1;// point group symmetry of g-vector grids? (enforced by truncation)

    int kpointorder = 3;//
    double bz_volume = 4;// bz volume [(2*pi/a)^3]

    // q-vector
    std::vector<Eigen::Vector3d> Q;// q-vector grids [2*pi/a]
    std::vector<int> qnum;// number of q-vectors in each path
    std::vector<std::string> qpath;// path of q-vector
    Eigen::Vector3d gammazero = {0,0,0}; // shift for the gamma-point to avoid singularity
};

// material parameters
struct mater {
    int NV;// number of vg values
    double g2max;// maximum |g|^2
    bool nonlocal = false;

    std::vector<double> g2;// |g_i|^2 [NV]
    double **vg;// v_{|g_i|^2}^{a} [NV x NATOM]
};

// environment variables
struct env {

    mater mat;
    lattice lat;

    // general options
    float version = 0.1;
    int debug = 0;// debug level
    char *outputfile;// filename of outputfile
    char *wdir;// working directory
    char *kpointfile;// filename of kpointfile
    int logfile = 0;// whether to save the logfile
    int kpoint = -1; // kpoint basic option; 0: full BZ, 1: Monkhorst-Pack grids, 2: Cohen-Chadi grids

    // frequency
    double *freq;// frequency [eV]
    double dfreq = 0;// frequency interval [eV]
    int kk = 0;// use kramers-kronig transformation to obtain the real parts of dielectric constant
    // double epsilon = 0.05; // delta function epsilon
    int delta = 0; // delta function model; 0: lorentzian, 1: gaussian

    int nvalence = 4;// number of valence electrons
    double energyoffset;// energy reference [eV]
    double eigerror;//
    Eigen::Vector3d refpoint = {0,0,0};// reference kpt where the energy reference is determined; gamma-point by default

    // lattice susceptibility tensor matrix [q][i][j][m][n][f]
    double ******ReXij;
    double ******ImXij;

    // electronic bandstructure [q][n]
    double **eband;
};

/*
    parser input class for parsing the input option strings
*/
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

        // fill in integer-valued option
        void fill_int(inputParser inp, const std::string option, int &value){
            if (inp.existOption(option)){
                value = std::stoi(inp.valueOption(option));
            }
        }

        // fill in double-valued option
        void fill_double(inputParser inp, const std::string option, double &value){
            if (inp.existOption(option)){
                value = std::stod(inp.valueOption(option));
            }
        }

    private:
        std::vector<std::string> tokens;
};

} /* namespace pmx */
#endif /* PMX_H_ */