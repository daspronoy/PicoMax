#pragma once
#ifndef _PMX_TBG_H_
#define _PMX_TBG_H_

/*
    pmx_tbg.cpp | includes twisted bilayer graphene into PicoMax

    Authors: Jungho Mun, Pronoy Das
    Date: November 5, 2024
*/


/*
    header files
*/
#include "pmx.hpp"


/*
    function declaration
*/
namespace pmx {

void tbg_eband(env &dat);

class tbg;

} /* namespace pmx */

class pmx::tbg {
    public:
        std::vector<Eigen::Vector2d> a_vec;
        std::vector<Eigen::Vector2d> b_vec;
        std::vector<Eigen::Vector2d> aM_vec;
        std::vector<Eigen::Vector2d> bM_vec;
        std::vector<Eigen::Vector2d> q_vec;
        std::vector<Eigen::Vector2d> gpt;
        std::vector<int> gi;
        std::vector<std::vector<std::array<int,2>>> g_nn;

        double kcut = 4;
        double kt;
        int ngmax = 20;
        double a = 2.46;// lattice constant [angstorm]
        double u = 0.11;
        double up = 0.11;

        double theta = 0;
        double phi = 0;
        double s = 1;

        /*
            constructor
        */
        tbg () {

            // graphene lattice vectors
            this->a_vec.reserve(2);
            this->a_vec.push_back({1,0});
            this->a_vec.push_back({0.5,0.5/sqrt(3)});
            this->b_vec.reserve(2);
            this->b_vec.push_back({1,-1/sqrt(3)});
            this->b_vec.push_back({0,2/sqrt(3)});
        }

        /*
        
        */
        void initialize() {
            
            // rotation matrix
            Eigen::Matrix2d Rp = this->R(-0.5*(this->theta));
            Eigen::Matrix2d Rn = this->R(0.5*(this->theta));

            // rotated lattice vectors
            Eigen::Vector2d b1p, b1n, b2p, b2n; // a1p, a2n, a2p, a2n, 
            b1p = Rp * (this->b_vec[0]);
            b2p = Rp * (this->b_vec[1]);
            b1n = Rn * (this->b_vec[0]);
            b2n = Rn * (this->b_vec[1]);

            // moire lattice vectors
            this->bM_vec.reserve(2);
            this->bM_vec.push_back(b1p-b1n);
            this->bM_vec.push_back(b2p-b2n);

            // translation vectors
            this->q_vec.reserve(3);
            this->q_vec.push_back(-2.0/3.0*this->bM_vec[1]-1.0/3.0*this->bM_vec[0]);
            this->q_vec.push_back(1.0/3.0*this->bM_vec[1]+2.0/3.0*this->bM_vec[0]);
            this->q_vec.push_back(1.0/3.0*this->bM_vec[1]-1.0/3.0*this->bM_vec[0]);
            this->kt = this->q_vec[0].norm();

            // construct g-vector grid
            int N = this->ngmax;
            double gcut = sqrt(3)*this->kt*this->kcut;
            for (int i=-N; i<=N; i++){for (int j=-N; j<=N; j++){for (int l=0; l<=1; l++){
                Eigen::Vector2d gg = i*this->bM_vec[0]+j*this->bM_vec[1]+l*this->q_vec[0];
                if (gg.norm()<=gcut){
                    this->gpt.push_back(gg);
                    this->gi.push_back(l);
                }
            }}}
            NPW = this->gpt.size();

            // find the nearest-neighbors on the g-vector lattice
            for (int i=0; i<NPW; i++){
                std::vector<std::array<int,2>> g_nni;
                
                if (this->gi[i]==1){// we only consider layer 1 (lower)
                    for (int q=0; q<3; q++){
                        // find g-point in upper layer (K+), such that $K^+ == K^- - Q$
                        for (int j=0; j<NPW; j++){
                            Eigen::Vector2d gq = this->gpt[j] - this->gpt[i] + this->q_vec[q];
                            if (gq.norm()<1e-6){
                                g_nni.push_back({j,q});
                            }
                        }
                    }
                }
                this->g_nn.push_back(g_nni);
            }
        }

        /*
            s : valley index {+1,-1}
        */
        Eigen::MatrixXcd gen_hamiltonian(double qx, double qy){
            Eigen::Vector2d q = {qx,qy};

            s = -(this->s);

            Eigen::Matrix2cd Sx, Sy, U1, U2, U3;

            Sx = Eigen::Matrix2cd {{0,1.0},{1.0,0}};
            Sy = Eigen::Matrix2cd {{0,-im},{im,0}};

            double hv = 2.0*pi*2.1354; // 2.1354 eV

            Eigen::MatrixXcd H_diag(2*NPW,2*NPW);
            for (int i=0; i<NPW; i++){
                int t = this->gi[i];
                int l = sigi(2*t-1);
                Eigen::Matrix2d Ri = R(-0.5*l*s*this->theta);
                Eigen::Vector2d h = -hv*Ri*(q - this->gpt[i]);
                Eigen::Matrix2cd hx, hy;
                hx = h[0]*s*Sx;
                hy = h[1]*Sy;
                H_diag({2*i,2*i+1},{2*i,2*i+1}) = hx+hy;
            }

            std::complex<double> w = exp(im*2.0/3.0*pi);
            U1 = Eigen::Matrix2cd {{this->u,this->up},{this->up,this->u}};
            U2 = Eigen::Matrix2cd {{this->u,this->up * pow(w,-s)},{this->up * pow(w,s),this->u}};
            U3 = Eigen::Matrix2cd {{this->u,this->up * pow(w,s)},{this->up * pow(w,-s),this->u}};

            Eigen::MatrixXcd H_offdiag(2*NPW,2*NPW);
            for (int i=0; i<NPW; i++){
                std::vector<std::array<int,2>> nn = this->g_nn[i];
                for (int j=0; j<nn.size(); j++){
                    int i1 = nn[j][0];// j
                    int iq = nn[j][1];// iq = {0,1,2}

                    H_offdiag({2*i,2*i+1},{2*i1,2*i1+1}) = (iq==1)*U1+(iq==2)*U2+(iq==0)*U3;
                }
            }

            //
            return H_diag + H_offdiag + H_offdiag.adjoint();
        }

        /*
        
        */
        void gen_gpt(){

        }

        /*

        */
        Eigen::Matrix2d R (double x){
            return Eigen::Matrix2d {{cos(x),-sin(x)},{sin(x),cos(x)}};
        }

        /*
        
        */
       int sigi (int x){
            if (x>=0){
                return 1;
            }else{
                return -1;
            }
       }
};



#endif /* _PMX_TBG_H_ */