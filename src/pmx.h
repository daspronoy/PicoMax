

//======================================================//
//                      libriaries                      //
//======================================================//
// #define EIGEN_USE_MKL_ALL
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
#include <iomanip> //For controlling precision in display
#include <omp.h>
// #include "stdafx.h"
#include <stdlib.h>
// #include "interpolation.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/LU>
#include <Eigen/Eigenvalues> 
#include <unsupported/Eigen/Splines>
#include <Eigen/Dense>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

// #include "kPointLatticeGenerator.h"


int NPW; // number of plane waves
int MAXPW;
double a; // lattice constant [angstrom]
int NEPS = 1; // epsilon cutoff
int NBAND = 8; // number of bands
int NKPT; // number of kpoint grids
int NFREQ = 1; // number of frequency









//======================================================//
//                      constants                       //
//======================================================//

#ifndef pi
#define pi 3.14159265358979323846
#endif

/*
	Physical constants:
	m_o*c^2 (eV), hbar_c (ev*m), xscale:1 Angstrom, vac_epsilon, rydberg, echarge, cscale, xl0
*/

#define hbar    1.054571817e-27 // [erg*s] (CGS)
#define e       4.803204712570263e-10 // [esu] (CGS)
#define e_m     9.1093837015e-28 // [g] (CGS)
#define eV      1.602176634e-12 // [erg] (CGS)
#define Ry      2.179872361103589e-11 // [erg] (CGS)
#define Ry2eV   13.605699169854621 // Ry/eV
#define angstrom  1e-8 // [cm] (CGS)
#define KCON    3.809982111485962 // hbar^2/(2*m_e*eV*angstrom^2)


// #define xm0     0.5109990615e+6
// #define hbarc   1.9732705359//10^-7 eVm
// #define emass 	9.1093837015//10^-31 kg
// #define hbar	1.054571817//10^-34 Js 
// #define vac_epsilon 8.854187817 //10^-12 in Farad/meter
// #define rydberg	13.60569312299426 // in eV	
// #define echarge	1.602176634 //in 10^-19 eV


// //	cscale = hbar^2/(2*m*xl0^2) 
// #define cscale  3.809984038754390
// #define xl0     1.0e-8
// #define xscale  1.0e-8

// #define PI      3.14159265358979323846
// #define HBAR    1.054571817e-34 // [J*s]
// #define HBAR_eV 6.582119569e-16 // [eV*s]
// #define EPS0    8.8541878128e-12 // [F/m]
// #define Rydberg 13.60569312299426 // [eV]
// #define E_MASS  9.1093837015e-31 // [kg]
// #define E_CHARGE 1.602176634e-19 // [C]   



namespace pmx {

//======================================================//
//                      parameters                      //
//======================================================//

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

// print help
void printHelp(){
    std::cout << "============================================================\n";
    std::cout << "  picomax-test0.1.0\n";
    std::cout << "Usage:\n";
    std::cout << "  picomax [-option] [value]\n";
    std::cout << "============================================================\n";
    std::cout << "-h, -help : print help\n";
    std::cout << "-v, -version : print PicoMax version\n";
    std::cout << "-debug : debug level {0}\n";
    std::cout << "\n";
    std::cout << "-numband : number of electronic bands\n";
    std::cout << "\n";
    std::cout << "-pathq : path of q-vectors\n";
    std::cout << "         {L,G,X,W,K,G}\n";
    std::cout << "-numq : number of q-vectors in each path\n";
    std::cout << "         {20,20,10,10,20}\n";
    std::cout << "\n";
    std::cout << "-numfreq : number of frequency points\n";
    std::cout << "-dfreq : frequency interval [eV]\n";
    std::cout << "\n";
    std::cout << "-gridscheme : {0}\n";
    std::cout << "-gridorder : {2}\n";
    std::cout << "\n";
    std::cout << "-a : lattice constant\n";
    std::cout << "-b1/b2/b3 : reciprocal lattice constants\n";
    std::cout << "       ex) -b1 1,1,-1\n";
    std::cout << "\n";
    std::cout << "============================================================\n";
    // std::cout << "\n";
    // std::cout << "\n";
    // std::cout << "\n";
    return;
}

// print version
void printVersion(){
    std::cout << "picomax-0.10\n";
    std::cout << "  Date: 2023.10.23\n";
    std::cout << "  Author: Sathwik Bharadwaj\n";
    std::cout << "  Email: sdaralag@purdue.edu\n";
    std::cout << "  https://electrodynamics.org\n";
    return;
}

//
struct env {

    double GCUT = 3.5;
    int **loci; // local-field plane-wave index

    // plane waves?
    std::vector<Eigen::Vector3d> G;

     
    int debug = 0; // debug level
    int dim = 3;//
    int kk = 0;// use kramers-kronig relation to obtain real
    int bandscheme = 0; // 0: empirical pseudopotential method, 1: epm with nonlocal


    // kgrid and kweight [k]
    int gridscheme = 0;
    int gridorder = 3;
    char *kgridfile;
    double **kgrid;
    double *kweight;
    // int numk; // number of k-points
    // int numk0; // total number of k-points

    double *freq;// frequency [eV]
    double dfreq;// frequency interval [eV]
    // int numfreq = 1;// number of frequency points

    // reciprocal lattice vectors {b1, b2, b3} in units of 2*pi/a
    std::vector<Eigen::Vector3d> b {
            Eigen::Vector3d(3), 
            Eigen::Vector3d(3), 
            Eigen::Vector3d(3)
        };

    // q-vector
    std::vector<Eigen::Vector3d> q;// q-vector
    std::vector<int> numq;// number of q-vectors in each path
    std::vector<std::string> pathq;// path of q-vector
    int nq;// total number of q-vectors

    // material parameters
    // double a;// lattice constant [angstroms]
    std::vector<double> epm;// nonlocal pseudopotential parameters [rydberg] {vs3,vs8,vs11,va3,va4,va11,a0,b0,A2,rl0,rl2}

    double energyoffset;
    double eigerror;

    // lattice permittivity tensors [q][m][n][f]
    double ****realepsL;
    double ****imagepsL;
    double ****realepsT;
    double ****imagepsT;

    // electronic bandstructure [q][n]
    double **band;
    
    // polarization?
    double **realpolL;
    double **imagpolL;
    double **realpolT;
    double **imagpolT;


    double epsilon = 0.05;

    // //recycle bin
    // double **points;
    // double *weights;
    // Eigen::Vector3d qvector = Eigen::Vector3d(3); // q-vector for dielectric function calculations
    // int numomega = 1; // number of freq points
    // double omega;// frequency [eV]
    // double vs3, vs4, vs8, vs11, va3, va4, va8, va11;// pseudopotentials in rydberg units
    // double a0, b0, A2, rl0, rl2;// nonlocal psuedopotential parameters in rydberg units

};

// preprocess input options into dat
void inputProcess(inputParser inp, env &dat){

    // debug level
    if (inp.existOption("-debug")){
        dat.debug = std::stoi(inp.valueOption("-debug"));
    }

    // GCUT
    if (inp.existOption("-GCUT")){
        dat.GCUT = std::stod(inp.valueOption("-GCUT"));
    }


    // reciprocal lattice vectors b1,b2,b3
    if (inp.existOption("-b1")){
        std::vector<std::string> v;
        v = inp.vectorOption("-b1");
        for (int i=0; i<v.size(); i++){
            dat.b[0](i) = std::stod(v[i]);
        }
    }
    if (inp.existOption("-b2")){
        std::vector<std::string> v;
        v = inp.vectorOption("-b2");
        for (int i=0; i<v.size(); i++){
            dat.b[1](i) = std::stod(v[i]);
        }
    }
    if (inp.existOption("-b3")){
        std::vector<std::string> v;
        v = inp.vectorOption("-b3");
        for (int i=0; i<v.size(); i++){
            dat.b[2](i) = std::stod(v[i]);
        }
    }

    // // number of plane waves
    // if (inp.existOption("-npw")){
    //     NPW = std::stoi(inp.valueOption("-npw"));
    // }

    // k-grids generation
    if (inp.existOption("-gridscheme")){
        dat.gridscheme = std::stoi(inp.valueOption("-gridscheme"));
    }
    if (inp.existOption("-gridorder")){
        dat.gridorder = std::stoi(inp.valueOption("-gridorder"));
    }
    if (inp.existOption("-kgridfile")){
        const char *filename = inp.valueOption("-kgridfile").c_str();
        dat.kgridfile = new char[strlen(filename)+1];
        strcpy(dat.kgridfile,filename);
    }

    // number of bands to consider
    if (inp.existOption("-nband")){
        NBAND = std::stoi(inp.valueOption("-nband"));
    }

    // dimension of dielectric matrix
    if (inp.existOption("-neps")){
        NEPS = std::stoi(inp.valueOption("-neps"));
    }

    // frequency
    if (inp.existOption("-nfreq") && inp.existOption("-dfreq")){
        NFREQ = std::stoi(inp.valueOption("-nfreq"));
        dat.dfreq = std::stod(inp.valueOption("-dfreq"));
        dat.freq = new double [NFREQ];
        for (int f=0; f<NFREQ; f++){
            dat.freq[f] = dat.dfreq*f;
        }
    }
    // switch for kk
    if (inp.existOption("-kk")){
        dat.kk = std::stoi(inp.valueOption("-kk"));
    }

    // q-vector
    if (inp.existOption("-q")){
        std::vector<std::string> v;
        v = inp.vectorOption("-q");
        Eigen::Vector3d vecq;
        for (int i=0; i<v.size(); i++){
            vecq(i) = std::stod(v[i]);
        }
        dat.q.push_back(vecq);
        dat.nq = 1;
    }else if ((inp.existOption("-pathq") && !inp.existOption("-numq"))
            || (!inp.existOption("-pathq") && inp.existOption("-numq"))){
        std::cout << "specify both -pathq and -numq"; abort();
    }else if (inp.existOption("-pathq") && inp.existOption("-numq")){
        std::vector<std::string> v;
        v = inp.vectorOption("-pathq");
        for (int i=0; i<v.size(); i++){
            dat.pathq.push_back(v[i]);
        }
        
        v = inp.vectorOption("-numq");
        for (int i=0; i<v.size(); i++){
            dat.numq.push_back(std::stoi(v[i]));
        }

        // construct q-grids
        dat.nq = std::accumulate(dat.numq.begin(),dat.numq.end(),0)+1;
        const static std::unordered_map<std::string,int> bandind{
            {"G",0},{"X",1},{"L",2},{"W",3},{"U",4},{"K",5}
        };
        Eigen::Vector3d v0,v1,dv;
        for (int i=0; i<dat.numq.size(); i++){
            switch (bandind.at(dat.pathq[i])){
                case 0:
                    v0 << 0.0,0.0,0.0; break;
                case 1:
                    v0 << 1.0,0.0,0.0; break;
                case 2:
                    v0 << 0.5,0.5,0.5; break;
                case 3:
                    v0 << 1.0,0.5,0.0; break;
                case 4:
                    v0 << 1.0,0.25,0.25; break;
                case 5:
                    v0 << 0.75,0.75,0.0; break;
            };
            switch (bandind.at(dat.pathq[i+1])){
                case 0:
                    v1 << 0.0,0.0,0.0; break;
                case 1:
                    v1 << 1.0,0.0,0.0; break;
                case 2:
                    v1 << 0.5,0.5,0.5; break;
                case 3:
                    v1 << 1.0,0.5,0.0; break;
                case 4:
                    v1 << 1.0,0.25,0.25; break;
                case 5:
                    v1 << 0.75,0.75,0.0; break;
            };
            dv = v1-v0;
            int nqi = dat.numq[i];
            for (int j=0; j<nqi; j++){
                dat.q.push_back(v0 + double(j)/double(nqi)*dv);
                if(dat.q.back().norm()==0) dat.q.back()(0) += 1e-6;
            }
        }
        dat.q.push_back(v1);
        if(dat.q.back().norm()==0) dat.q.back()(0) += 1e-6;
    }// q-vector

    // lattice constant
    if (inp.existOption("-a")){
        a = std::stod(inp.valueOption("-a"));
    }

    // parameters for empirical pseudopotential method
    if (inp.existOption("-epm")){
        std::vector<std::string> v;
        v = inp.vectorOption("-epm");
        for (int i=0; i<v.size(); i++){
            dat.epm.push_back(std::stod(v[i]));
        }
        if (v.size()==6) dat.bandscheme = 0;
        else if (v.size()==11) dat.bandscheme = 1;
        else {
            std::cout << "check epm parameters\n";
            abort();
        }
    }
    
    // delta function epsilon
    if (inp.existOption("-epsilon")){
        dat.epsilon = std::stod(inp.valueOption("-epsilon"));
    }

    return;
}

// progress messaging
void progmsg(){
    
}

// print input (for debugging purposes)
void printInput(env &dat){
    std::cout << "Number of plane waves: " << NPW << std::endl;
    std::cout << "Number of bands: " << NBAND << std::endl;
    std::cout << "Dimension of dielectric matrix: " << NEPS << std::endl;
    std::cout << "Material properties" << std::endl;
    std::cout << "Lattice constant: " << a << std::endl;
    std::cout << "Empirical pseudopotential method parameters:" << std::endl;
    std::cout << "vs3\tvs8\tvs11\tva3\tva4\tva11\ta0\tb0\tA2\trl0\trl2" << std::endl;
    for(int i=0; i<dat.epm.size(); i++){
        std::cout << dat.epm[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "Lattice vectors" << std::endl;
    std::cout << "b1\tb2\tb3" << std::endl;
    for(int i=0; i<3; i++){
        std::cout << dat.b[0](i) << "\t" << dat.b[1](i) << "\t" << dat.b[2](i) << std::endl;
    }
    std::cout << "Number of k-grids: " << NKPT << std::endl;
    std::cout << "Order of k-grids: " << dat.gridorder << std::endl;
    std::cout << "Number of q-grids: " << dat.nq << std::endl;
    for(int i=0; i<3; i++){
        for(int q=0; q<dat.nq; q++){
            std::cout << dat.q[q](i) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Number of G-grids: " << NPW << std::endl;
    for(int i=0; i<3; i++){
        for(int g=0; g<NPW; g++){
            std::cout << dat.G[g](i) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Number of K-grids: " << NKPT << std::endl;
    for(int i=0; i<3; i++){
        for(int k=0; k<NKPT; k++){
            std::cout << dat.kgrid[k][i] << "\t";
        }
        std::cout << std::endl;
    }



    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    return;
}

//======================================================//
//                      header files                    //
//======================================================//


// 


void genkgrid_cohen_chadi (env &dat);

void genkgrid_MP (env &dat);

bool kroneckerDelta(int m, int n);

void empiricalpseudopotentialmethod(env &dat);

// inline std::vector<double> millerindex(int n, int N);

// inline Eigen::Vector3d recvec(int n, int N, std::vector<Eigen::Vector3d> b);

// std::complex<double> matrixelement(int m, int n, int N, Eigen::Vector3d kvec, std::vector<Eigen::Vector3d> b, std::vector<double> v, double a);

std::complex<double> pseudopotential(Eigen::Vector3d Kvec, std::vector<double> v);

std::complex<double> nonlocalpseudopotential(Eigen::Vector3d Kkvecm, Eigen::Vector3d Kkvecn, 
                    std::vector<double>v);

double sphbesj(int n, double r);

void printBand(env &dat);



void permittivity(env &dat);

void setRefEnergy(env &dat);

Eigen::Vector3d perpvec(Eigen::Vector3d v);

std::complex<double> diracdelta(std::complex<double> value, double epsilon);

void kramerskronig(env &dat);
void kramerskronigtransform(double *ReX, double *ImX, double *w, double dw);


void printEpsilon(env &dat);
void writeEpsilon(env &dat);




void generate_planewaves (env &dat);
Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G, Eigen::Vector3d K, std::vector<double> v);

std::vector<Eigen::Vector3d> find_G (int g2, int n);
void find_locGi (env &dat);
}



