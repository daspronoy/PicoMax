/*
    basis.cpp




    Author:
*/
#include "pmx.hpp"

namespace pmx{

// generate plane waves G for fcc in rectangular fractional coordinate
void gen_G_fcc_f (env &dat){
    int N = 2;
    NPW = pow(2*N+1,3);
    dat.G.reserve(NPW);

    Eigen::Vector3d b1, b2, b3;
    b1 << -1.0,1.0,1.0;
    b2 << 1.0,-1.0,1.0;
    b3 << 1.0,1.0,-1.0;

    for (int h=-N; h<N; h++){ for (int k=-N; k<N; k++){ for (int l=-N; l<N; l++){
        dat.G.push_back(h*b1+k*b2+l*b3);
    }}}
    return;
}

// generate plane waves G for fcc in spherical...
void generate_planewaves (env &dat){
    std::vector<int> g2 {0,3,4,8,11,12,16,19,20,24,27,32,35,36,40,43,44,48,51,52,56,59,64,67,68,72,75,76,80,83,84,88,91,96,99,100};
    std::vector<int> nnn {1,8,6,12,24,8,6,24,24,24,32,12,48,30,24,24,24,8,48,24,48,72,6,24,48,36,56,24,24,72,48,24,48,24,72,30};
    std::vector<int>::iterator itr = std::lower_bound(g2.begin(), g2.end(), pow(dat.gcut,2));
    int N = int(itr - g2.begin());
    std::vector<int> cnn(N);

    NPW = 0;
    for (int i=0; i<N; i++){
        NPW += nnn[i];
        if (i==0) cnn[i] = 1;
        else cnn[i] = cnn[i-1]+nnn[i];
    }

    dat.G.reserve(NPW);
    for (int i=0; i<NPW; i++){
        dat.G.push_back({0,0,0});
    }

    #pragma omp parallel for
    for (int i=0; i<N; i++){
        std::vector<Eigen::Vector3d> iG = find_G(g2[i],nnn[i]);
        int in;
        if (i==0) in = 0;
        else in = cnn[i-1];
        for (int ic=0; ic<nnn[i]; ic++){
            dat.G[in+ic] = iG[ic];
        }
    }
    #pragma omp barrier

    // find local-field index
    find_locGi(dat);

    return;
}

// find n reciprocal lattice vectors with squared magnitude g2
std::vector<Eigen::Vector3d> find_G (int g2, int n){
    std::vector<Eigen::Vector3d> G(n);

    int ia = int(sqrt(double(g2)));
    int gx2, gy2, gz2, gxy2, gxyz2;
    int ic = 0;

    for (int gx=-ia; gx<=ia; gx++){
        gx2 = gx*gx;
        for (int gy=-ia; gy<=ia; gy++){
            gy2 = gy*gy;
            gxy2 = gx2+gy2;
            for (int gz=-ia; gz<=ia; gz++){
                gz2 = gz*gz;
                gxyz2 = gxy2+gz2;
                if (gxyz2==g2){
                    G[ic] << double(gx), double(gy), double(gz);
                    ic +=1;
                    if (ic==n) return G;
                }
            }
        }
    }
    return G;
}

// find local-field index
void find_locGi (env &dat){

    // initialize
    dat.loci = new int *[NEPS];
    for (int m=0; m<NEPS; m++){
        dat.loci[m] = new int [NPW];
    }

    // g(0)
    for (int i=0; i<NPW; i++){
        dat.loci[0][i] = i;
    }
    #pragma omp parallel for
    for (int m=1; m<NEPS; m++){

        // local-field vector Gm
        Eigen::Vector3d Gm = dat.G[m];
        for (int i=0; i<NPW; i++){
            // bra vector Gi
            Eigen::Vector3d Gi = dat.G[i];

            // if no such vector found, -1
            dat.loci[m][i] = -1;
            for (int j=0; j<NPW; j++){
                // ket vector Gj
                Eigen::Vector3d Gj = dat.G[j];

                // find j, such that Gi = Gj-Gm
                if ((Gi-Gj+Gm).norm()<1e-10){
                    dat.loci[m][i] = j;
                    break;
                }
            }
        }
    }
    #pragma omp barrier
    return;
}




void genkpoint_cohen_chadi (env &dat) {
    //Here we will use Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    int i, j, k, l; //dummy indices
    int N = dat.kpointorder;
    NKPT = int(pow(2.0,N-1)*(pow(2.0,N)+1.0)*(pow(2.0,N-1.0)+1.0)/3.0);
    dat.K = new double*[NKPT];// xi [N x 3]
    dat.KW = new double[NKPT];// wi [N]
    for (i=0; i<NKPT; i++){
        dat.K[i] = new double[3];
    }
    
    //2-point Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    if (N==1){
        // 0 point
        dat.KW[0] = 0.7500; dat.K[0][0] = 0.75; dat.K[0][1] = 0.25; dat.K[0][2] = 0.25;
        // 1 point
        dat.KW[1] = 0.2500; dat.K[1][0] = 0.25; dat.K[1][1] = 0.25; dat.K[1][2] = 0.25;
    }

    //10-point Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    else if (N==2){
        // 0 point
        dat.KW[0] = 0.18750; dat.K[0][0] = 0.875; dat.K[0][1] = 0.375; dat.K[0][2] = 0.125;
        // 1 point
        dat.KW[1] = 0.09375; dat.K[1][0] = 0.875; dat.K[1][1] = 0.125; dat.K[1][2] = 0.125;
        // 2 point
        dat.KW[2] = 0.09375; dat.K[2][0] = 0.625; dat.K[2][1] = 0.625; dat.K[2][2] = 0.125;
        // 3 point
        dat.KW[3] = 0.09375; dat.K[3][0] = 0.625; dat.K[3][1] = 0.375; dat.K[3][2] = 0.375;
        // 4 point
        dat.KW[4] = 0.18750; dat.K[4][0] = 0.625; dat.K[4][1] = 0.375; dat.K[4][2] = 0.125;
        // 5 point
        dat.KW[5] = 0.09375; dat.K[5][0] = 0.625; dat.K[5][1] = 0.125; dat.K[5][2] = 0.125;
        // 6 point
        dat.KW[6] = 0.03125; dat.K[6][0] = 0.375; dat.K[6][1] = 0.375; dat.K[6][2] = 0.375;
        // 7 point
        dat.KW[7] = 0.09375; dat.K[7][0] = 0.375; dat.K[7][1] = 0.375; dat.K[7][2] = 0.125;
        // 8 point
        dat.KW[8] = 0.09375; dat.K[8][0] = 0.375; dat.K[8][1] = 0.125; dat.K[8][2] = 0.125;
        // 9 point
        dat.KW[9] = 0.03125; dat.K[9][0] = 0.125; dat.K[9][1] = 0.125; dat.K[9][2] = 0.125;
    }

    //Higher-order Cohen and Chadi scheme: PHILOSOPHICAL MAGAZINE B 81, NO. 6, 551--559 (2001).
    else if (N>2){
        k = 0;
        for (i=1; i<=pow(2,N-1); i++)	for (j =1; j<=i; j++)	for (l=1; l<=j; l++){
            //std::cout << "Loop 1\n" << std::endl;
            dat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.K[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.K[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.K[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }	
        for (i=pow(2,N-1)+2; i<=pow(2,N); i+=2)	for (j=1; j<=3*pow(2,N-2)-int(i*0.5); j++)	for (l=1; l<=j; l++){
            //std::cout << "Loop 2\n" << std::endl;
            dat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.K[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.K[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.K[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }	
        for (i=pow(2,N-1)+2; i<=3*pow(2,N-2); i+=2)	for (j=3*pow(2,N-2)-int(i*0.5)+1; j<=i; j++)	for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 3\n" << std::endl;
            dat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.K[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.K[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.K[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
        for (i=3*pow(2,N-2)+2; i<=pow(2,N); i+=2)	for (j=3*pow(2,N-2)-int(i*0.5)+1; j<=3*pow(2,N-1)-i; j++)	for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 4\n" << std::endl;
            dat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.K[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.K[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.K[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
        for (i=pow(2,N-1)+1; i<=pow(2,N)-1; i+=2)	for (j=1; j<=3*pow(2,N-2)-int((i-1)*0.5); j++)	for (l=1; l<=j; l++){
            //std::cout << "Loop 5\n" << std::endl;
            dat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.K[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.K[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.K[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
        for (i=pow(2,N-1)+1; i<=3*pow(2,N-2)-1; i+=2)	for (j=3*pow(2,N-2)-int((i-1)*0.5)+1; j<=i; j++)	for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 6\n" << std::endl;
            dat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.K[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.K[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.K[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
        for (i=3*pow(2,N-2)+1; i<=pow(2,N)-1; i+=2)	for (j=3*pow(2,N-2)-int((i-1)*0.5)+1; j<=3*pow(2,N-1)-i; j++)	for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 7\n" << std::endl;
            dat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            dat.K[k][0] = (2.0*double(i)-1)/pow(2,N+1);	dat.K[k][1] = (2.0*double(j)-1)/pow(2,N+1);	dat.K[k][2] = (2.0*double(l)-1)/pow(2,N+1);
            k++;
        }
    }
    return;
}

// Krnoecker Delta Function
inline bool kroneckerDelta(int m, int n) {
	return (m == n);
}

void genkpoint_monkhorst_pack (env &dat) {

}

void importkpointfile (env &dat) {
    std::ifstream kdata(dat.kpointfile);
    int i = 0; //dummy index
    int NKPT0;

    std::string line;
    if (kdata.is_open())
    {
        while(getline (kdata,line)){
        if (!(line[0]=='#' || line=="")){
            //Line 0: Number of integration k-points in entire BZ
            if (i == 0){
                std::stringstream(line) >> NKPT0;
                i++;
                continue;	
            }
            //Line 1: Number of irreducible k-points in BZ
            else if (i==1){
                std::stringstream(line) >> NKPT;  
                dat.K = new double*[NKPT];
                dat.KW = new double[NKPT];
                for (int j=0; j<NKPT; j++){
                    dat.K[j] = new double[3];
                }//end of for j
                i++;
                continue;
            }
            else if (i>1 && i-2<NKPT){
                //Line 1: Number of plane waves
                std::stringstream(line) >> dat.K[i-2][0] >>  dat.K[i-2][1] >>  dat.K[i-2][2] >>  dat.KW[i-2];  
                dat.KW[i-2] = dat.KW[i-2]/NKPT0;
                i++;
                continue; 
            
            }//else if i
            }//end of if line
        }//end of while 
    }//end of if kdata open
    kdata.close();
    return;
}

}