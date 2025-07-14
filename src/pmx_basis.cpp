/*
    pmx_basis.cpp | includes basis-related functinalities of picomax

    * kpoint generation for BZ integration
    * gpoint generation for planewave basis set
    * high-symmetry point
    * polarization unit vector

    Authors: Jungho Mun, Pronoy Das
*/

/*
    header file
*/
#include "pmx_basis.hpp"

/*
    function declaration
*/
namespace pmx{

//=======================================================//
//                  GPOINT generation                    //
//=======================================================//

/*
    generate G-vectors, truncated by `gcut` determined from `encut`
*/
void generate_gpt (lattice &lat){
    double tol = 1e-6;
    int imax = ceil(lat.gcut);
    int N = pow(2*imax+1,lat.dim);

    // generate G-vectors
    std::vector<Eigen::Vector3d> G;
    G.reserve(N);
    if (lat.dim==3){
        Eigen::Vector3d b1 = lat.reciprocal[0];
        Eigen::Vector3d b2 = lat.reciprocal[1];
        Eigen::Vector3d b3 = lat.reciprocal[2];
        for (int i1=-imax; i1<=imax; i1++){for (int i2=-imax; i2<=imax; i2++){for (int i3=-imax; i3<=imax; i3++){
            G.push_back(i1*b1+i2*b2+i3*b3);
        }}}
    }else if (lat.dim==2){
        Eigen::Vector3d b1 = lat.reciprocal[0];
        Eigen::Vector3d b2 = lat.reciprocal[1];
        for (int i1=-imax; i1<=imax; i1++){for (int i2=-imax; i2<=imax; i2++){
            G.push_back(i1*b1+i2*b2);
        }}
    }
    
    std::vector<double> Gnorm;
    Gnorm.reserve(N);
    for (int i=0; i<N; i++){
        Gnorm.push_back(G[i].norm());
    }

    std::vector<int> index(N);
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), 
        [&](int A, int B) -> bool {
            return Gnorm[A] < Gnorm[B];
        });
    std::sort(Gnorm.begin(), Gnorm.end(), std::less());

    NPW = N;
    for (int i=0; i<N; i++){
        if (Gnorm[i]>lat.gcut+tol){
            NPW = i;
            break;
        }
    }
    
    lat.G.reserve(NPW);
    for (int i=0; i<NPW; i++){
        lat.G.push_back(G[index[i]]);
    }

    // generate index shift by lattice vectors
    generate_loci(lat);
}

/*
    generate G-vectors rectangular in the fractional coordinate [hkl]
*/
void generate_gpt_rect (lattice &lat){
    double tol = 1e-6;
    int imax = ceil(lat.gcut);
    int N = pow(2*imax+1,lat.dim);
    NPW = N;

    // generate G-vectors
    std::vector<Eigen::Vector3d> G;
    G.reserve(N);
    if (lat.dim==3){
        Eigen::Vector3d b1 = lat.reciprocal[0];
        Eigen::Vector3d b2 = lat.reciprocal[1];
        Eigen::Vector3d b3 = lat.reciprocal[2];
        for (int i1=-imax; i1<=imax; i1++){for (int i2=-imax; i2<=imax; i2++){for (int i3=-imax; i3<=imax; i3++){
            G.push_back(i1*b1+i2*b2+i3*b3);
        }}}
    }else if (lat.dim==2){
        Eigen::Vector3d b1 = lat.reciprocal[0];
        Eigen::Vector3d b2 = lat.reciprocal[1];
        for (int i1=-imax; i1<=imax; i1++){for (int i2=-imax; i2<=imax; i2++){
            G.push_back(i1*b1+i2*b2);
        }}
    }

    // sort G-vectors
    std::vector<double> Gnorm;
    Gnorm.reserve(N);
    for (int i=0; i<N; i++){
        Gnorm.push_back(G[i].norm());
    }

    std::vector<int> index(N);
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), 
        [&](int A, int B) -> bool {
            return Gnorm[A] < Gnorm[B];
        });
    std::sort(Gnorm.begin(), Gnorm.end(), std::less());

    lat.G.reserve(NPW);
    for (int i=0; i<NPW; i++){
        lat.G.push_back(G[index[i]]);
    }

    // generate index shift by lattice vectors
    generate_loci(lat);
}

/*
    find local-field index i, such that Gi=Gp+Gm

    idx_offset [NEPS x NPW]

    index p: planewave summation index
    Gm: offset wavevector (local-field effect)
*/ 
void generate_loci (lattice &lat){

    // initialize
    lat.loci = new int *[NEPS];
    for (int m=0; m<NEPS; m++){
        lat.loci[m] = new int [NPW];
    }

    // g(0)
    for (int p=0; p<NPW; p++){
        lat.loci[0][p] = p;
    }
    #pragma omp parallel for
    for (int m=1; m<NEPS; m++){

        // local-field vector Gm
        Eigen::Vector3d Gm = lat.G[m];
        for (int p=0; p<NPW; p++){
            // bra vector Gp
            Eigen::Vector3d Gp = lat.G[p];

            // if no such vector found, -1
            lat.loci[m][p] = -1;
            for (int i=0; i<NPW; i++){
                // ket vector Gi
                Eigen::Vector3d Gi = lat.G[i];

                // find i, such that Gi = Gp+Gm
                if ((Gp-Gi+Gm).norm()<1e-10){
                    lat.loci[m][p] = i;
                    break;
                }
            }
        }
    }
    #pragma omp barrier
    return;
}




//=======================================================//
//                  KPOINT generation                    //
//=======================================================//

/*
    generate kpoint using monkhorst-pack scheme

    not implemented yet...
*/
void generate_kpt_monkhorst_pack (env &dat) {

}

/*
    import k-vectors from an external file
*/
void import_kpt (env &dat) {
    std::ifstream kdata(dat.kpointfile);
    int i = 0; //dummy index
    int NKPT0;

    std::string line;
    if (kdata.is_open())
    {
        while(getline(kdata,line)){
            if (!(line[0]=='#' || line=="")){
                //Line 0: Number of integration k-points in entire BZ
                if (i==0){
                    std::stringstream(line) >> NKPT0;
                    i++;
                    continue;	
                }
                //Line 1: Number of irreducible k-points in BZ
                else if (i==1){
                    std::stringstream(line) >> NKPT;

                    dat.lat.K.reserve(NKPT);
                    dat.lat.KW.reserve(NKPT);
                    for (int j=0; j<NKPT; j++){
                        dat.lat.K.push_back({0,0,0});
                        dat.lat.KW.push_back(0);
                    }
                    i++;
                    continue;
                }
                else if (i>1 && i-2<NKPT){
                    //Line 1: Number of plane waves
                    std::stringstream(line) >> dat.lat.K[i-2](0) >> dat.lat.K[i-2](1) >> dat.lat.K[i-2](2) >> dat.lat.KW[i-2];
                    dat.lat.KW[i-2] /= NKPT0;// normalize by the total kpoint number

                    // std::stringstream(line) >> dat.K[i-2][0] >>  dat.K[i-2][1] >>  dat.K[i-2][2] >>  dat.KW[i-2];  
                    // dat.KW[i-2] = dat.KW[i-2]/NKPT0;
                    i++;
                    continue; 
                
                }//else if i
            }//end of if line
        }//end of while 
    }//end of if kdata open
    kdata.close();
    return;
}


/*
    generate kpoint using Cohen-Chadi scheme

    applicable for only fcc lattice?

    [1] Phys. Rev. B 8, 5757 (1973)
    [2] PHILOSOPHICAL MAGAZINE B 81, NO. 6, 551--559 (2001)
*/
void generate_kpt_cohen_chadi (lattice &lat) {
    //Here we will use Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    int i, j, k, l; //dummy indices
    int N = lat.kpointorder;
    NKPT = int(pow(2.0,N-1)*(pow(2.0,N)+1.0)*(pow(2.0,N-1.0)+1.0)/3.0);
    lat.K.reserve(NKPT);
    lat.KW.reserve(NKPT);
    for (int k=0; k<NKPT; k++){
        lat.K.push_back({0,0,0});
        lat.KW.push_back(0);
    }
    
    //2-point Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    if (N==1){
        lat.KW[0] = 0.7500; lat.K[0][0] = 0.75; lat.K[0][1] = 0.25; lat.K[0][2] = 0.25;
        lat.KW[1] = 0.2500; lat.K[1][0] = 0.25; lat.K[1][1] = 0.25; lat.K[1][2] = 0.25;
    }

    //10-point Cohen and Chadi scheme: Phys. Rev. B 8, 5757 (1973).
    else if (N==2){
        lat.KW[0] = 0.18750; lat.K[0][0] = 0.875; lat.K[0][1] = 0.375; lat.K[0][2] = 0.125;
        lat.KW[1] = 0.09375; lat.K[1][0] = 0.875; lat.K[1][1] = 0.125; lat.K[1][2] = 0.125;
        lat.KW[2] = 0.09375; lat.K[2][0] = 0.625; lat.K[2][1] = 0.625; lat.K[2][2] = 0.125;
        lat.KW[3] = 0.09375; lat.K[3][0] = 0.625; lat.K[3][1] = 0.375; lat.K[3][2] = 0.375;
        lat.KW[4] = 0.18750; lat.K[4][0] = 0.625; lat.K[4][1] = 0.375; lat.K[4][2] = 0.125;
        lat.KW[5] = 0.09375; lat.K[5][0] = 0.625; lat.K[5][1] = 0.125; lat.K[5][2] = 0.125;
        lat.KW[6] = 0.03125; lat.K[6][0] = 0.375; lat.K[6][1] = 0.375; lat.K[6][2] = 0.375;
        lat.KW[7] = 0.09375; lat.K[7][0] = 0.375; lat.K[7][1] = 0.375; lat.K[7][2] = 0.125;
        lat.KW[8] = 0.09375; lat.K[8][0] = 0.375; lat.K[8][1] = 0.125; lat.K[8][2] = 0.125;
        lat.KW[9] = 0.03125; lat.K[9][0] = 0.125; lat.K[9][1] = 0.125; lat.K[9][2] = 0.125;
    }

    //Higher-order Cohen and Chadi scheme: PHILOSOPHICAL MAGAZINE B 81, NO. 6, 551--559 (2001).
    else if (N>2){
        k = 0;
        for (i=1; i<=pow(2,N-1); i++) for (j=1; j<=i; j++) for (l=1; l<=j; l++){
            //std::cout << "Loop 1\n" << std::endl;
            lat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            lat.K[k] = {(2.0*double(i)-1)/pow(2,N+1), (2.0*double(j)-1)/pow(2,N+1), (2.0*double(l)-1)/pow(2,N+1)};
            k++;
        }	
        for (i=pow(2,N-1)+2; i<=pow(2,N); i+=2) for (j=1; j<=3*pow(2,N-2)-int(i*0.5); j++) for (l=1; l<=j; l++){
            //std::cout << "Loop 2\n" << std::endl;
            lat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            lat.K[k] = {(2.0*double(i)-1)/pow(2,N+1), (2.0*double(j)-1)/pow(2,N+1), (2.0*double(l)-1)/pow(2,N+1)};
            k++;
        }	
        for (i=pow(2,N-1)+2; i<=3*pow(2,N-2); i+=2) for (j=3*pow(2,N-2)-int(i*0.5)+1; j<=i; j++) for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 3\n" << std::endl;
            lat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            lat.K[k] = {(2.0*double(i)-1)/pow(2,N+1), (2.0*double(j)-1)/pow(2,N+1), (2.0*double(l)-1)/pow(2,N+1)};
            k++;
        }
        for (i=3*pow(2,N-2)+2; i<=pow(2,N); i+=2) for (j=3*pow(2,N-2)-int(i*0.5)+1; j<=3*pow(2,N-1)-i; j++) for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 4\n" << std::endl;
            lat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            lat.K[k] = {(2.0*double(i)-1)/pow(2,N+1), (2.0*double(j)-1)/pow(2,N+1), (2.0*double(l)-1)/pow(2,N+1)};
            k++;
        }
        for (i=pow(2,N-1)+1; i<=pow(2,N)-1; i+=2) for (j=1; j<=3*pow(2,N-2)-int((i-1)*0.5); j++) for (l=1; l<=j; l++){
            //std::cout << "Loop 5\n" << std::endl;
            lat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            lat.K[k] = {(2.0*double(i)-1)/pow(2,N+1), (2.0*double(j)-1)/pow(2,N+1), (2.0*double(l)-1)/pow(2,N+1)};
            k++;
        }
        for (i=pow(2,N-1)+1; i<=3*pow(2,N-2)-1; i+=2) for (j=3*pow(2,N-2)-int((i-1)*0.5)+1; j<=i; j++) for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 6\n" << std::endl;
            lat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            lat.K[k] = {(2.0*double(i)-1)/pow(2,N+1), (2.0*double(j)-1)/pow(2,N+1), (2.0*double(l)-1)/pow(2,N+1)};
            k++;
        }
        for (i=3*pow(2,N-2)+1; i<=pow(2,N)-1; i+=2) for (j=3*pow(2,N-2)-int((i-1)*0.5)+1; j<=3*pow(2,N-1)-i; j++) for (l=1; l<=3*pow(2,N-1)+1-i-j; l++){
            //std::cout << "Loop 7\n" << std::endl;
            lat.KW[k] = (3.0*(1.0-kroneckerDelta(i,j))+3.0*(1.0-kroneckerDelta(j,l))+kroneckerDelta(i,l))/(4.0*pow(8.0,N-1));
            lat.K[k] = {(2.0*double(i)-1)/pow(2,N+1), (2.0*double(j)-1)/pow(2,N+1), (2.0*double(l)-1)/pow(2,N+1)};
            k++;
        }
    }
    return;
}

/*
    Kronecker delta function
*/
inline bool kroneckerDelta(int m, int n) {
	return (m==n);
}

/*
    generate kpoint in the full 1st BZ
*/
void generate_kpt_bz (lattice &lat){
    if (lat.name.compare("zb")==0 || lat.name.compare("rs")==0){// fcc
        generate_kpt_fcc2(lat);
    }else if (lat.name.compare("wz")==0 || lat.name.compare("te")==0){// hex
        generate_kpt_hex(lat);
    }else if (lat.name.compare("hex")==0){// 2d hexagonal
        generate_kpt_hex2(lat);
    }
}

/*
    generate k-vectors for fcc crystal
*/
void generate_kpt_fcc (lattice &lat){
    double tol = 1e-6;
    int N = lat.kpointorder;

    // reciprocal lattice vectors
    Eigen::Vector3d b1 = lat.reciprocal[0];
    Eigen::Vector3d b2 = lat.reciprocal[1];
    Eigen::Vector3d b3 = lat.reciprocal[2];

    // generate fractional coordinate grids
    double df = 1.0/N;
    NKPT = pow(N,3);

    // shift index
    int shiftindex [27][3];
    int cnt = 0;
    for (int i1=-1; i1<=1; i1++){for (int i2=-1; i2<=1;i2++){for (int i3=-1; i3<=1; i3++){
        shiftindex[cnt][0] = i1;
        shiftindex[cnt][1] = i2;
        shiftindex[cnt][2] = i3;
        cnt += 1;
    }}}

    lat.K.reserve(NKPT);
    for (int i1=0; i1<N; i1++){for (int i2=0; i2<N; i2++){for (int i3=0; i3<N; i3++){
        lat.K.push_back(df*i1*b1+df*i2*b2+df*i3*b3);
    }}}
    for (int i=0; i<NKPT; i++){
        if (!checkinbz_fcc(lat.K[i])){
            for (int c=0; c<27; c++){
                if (checkinbz_fcc(lat.K[i]+shiftindex[c][0]*b1+shiftindex[c][1]*b2+shiftindex[c][2]*b3)){
                    lat.K[i] += shiftindex[c][0]*b1+shiftindex[c][1]*b2+shiftindex[c][2]*b3;
                    break;
                }
            }
        }
    }

    // uniform weight of k-points
    lat.KW.reserve(NKPT);
    for (int k=0; k<NKPT; k++){
        lat.KW[k] = 1.0/NKPT;
    }
}

/*
    generate k-vectors in the full 1st BZ of hexagonal crystal structure
*/

void symmetrize_kpoints_xy_hex(lattice &lat) {
    std::vector<Eigen::Vector3d> K_sym;
    std::vector<double> KW_sym;
    
    for (size_t i = 0; i < lat.K.size(); i++) {
        Eigen::Vector3d k = lat.K[i];
        double w = lat.KW[i];
        
        // For tellurium, the main symmetry breaking x↔y
        std::vector<Eigen::Vector3d> symmetric_points;
        std::vector<double> weights;
        
        // Original point
        symmetric_points.push_back(k);
        weights.push_back(w);
        
        // x↔y swapped point
        Eigen::Vector3d k_xy_swap = {k(1), k(0), k(2)};
        if (checkinbz_hex(k_xy_swap, lat.f) && (k - k_xy_swap).norm() > 1e-10) {
            symmetric_points.push_back(k_xy_swap);
            weights.push_back(w);
        }
        
        // Average the weights
        double total_weight = 0;
        for (double wt : weights) total_weight += wt;
        double avg_weight = total_weight / symmetric_points.size();
        
        // Add all symmetric points with equal weights
        for (const auto& sym_k : symmetric_points) {
            K_sym.push_back(sym_k);
            KW_sym.push_back(avg_weight);
        }
    }
    
    lat.K = K_sym;
    lat.KW = KW_sym;
    NKPT = lat.K.size();
}
}


void generate_kpt_hex (lattice &lat){
    double tol = 1e-6;
    int N = lat.kpointorder;

    // reciprocal lattice vectors
    Eigen::Vector3d b1 = lat.reciprocal[0];
    Eigen::Vector3d b2 = lat.reciprocal[1];
    Eigen::Vector3d b3 = lat.reciprocal[2];

    // generate fractional coordinate grids
    double df = 1.0/N;
    NKPT = pow(N,3);

    // shift index
    int shiftindex [27][3];
    int cnt = 0;
    for (int i1=-1; i1<=1; i1++){for (int i2=-1; i2<=1;i2++){for (int i3=-1; i3<=1; i3++){
        shiftindex[cnt][0] = i1;
        shiftindex[cnt][1] = i2;
        shiftindex[cnt][2] = i3;
        cnt += 1;
    }}}

    lat.K.reserve(NKPT);
    for (int i1=0; i1<N; i1++){for (int i2=0; i2<N; i2++){for (int i3=0; i3<N; i3++){
        lat.K.push_back(df*i1*b1+df*i2*b2+df*i3*b3);
    }}}
    for (int i=0; i<NKPT; i++){
        if (!checkinbz_hex(lat.K[i],lat.f)){
            for (int c=0; c<27; c++){
                if (checkinbz_hex(lat.K[i]+shiftindex[c][0]*b1+shiftindex[c][1]*b2+shiftindex[c][2]*b3,lat.f)){
                    lat.K[i] += shiftindex[c][0]*b1+shiftindex[c][1]*b2+shiftindex[c][2]*b3;
                    break;
                }
            }
        }
    }

    // uniform weight of k-points
    lat.KW.reserve(NKPT);
    for (int k=0; k<NKPT; k++){
        lat.KW[k] = 1.0/NKPT;
    }

    // ADD SYMMETRIZATION HERE
    symmetrize_kpoints_xy_hex(lat);
}





/*
    generate k-vectors in the full 1st BZ of hexagonal crystal structure
*/
void generate_kpt_hex2 (lattice &lat){
    double tol = 1e-6;
    int N = lat.kpointorder;

    // reciprocal lattice vectors
    Eigen::Vector3d b1 = lat.reciprocal[0];
    Eigen::Vector3d b2 = lat.reciprocal[1];

    // generate fractional coordinate grids
    double df = 1.0/N;
    NKPT = pow(N,2);

    // shift index
    int shiftindex [9][2];
    int cnt = 0;
    for (int i1=-1; i1<=1; i1++){for (int i2=-1; i2<=1;i2++){
        shiftindex[cnt][0] = i1;
        shiftindex[cnt][1] = i2;
        cnt += 1;
    }}

    lat.K.reserve(NKPT);
    for (int i1=0; i1<N; i1++){for (int i2=0; i2<N; i2++){
        lat.K.push_back(df*i1*b1+df*i2*b2);
    }}
    for (int i=0; i<NKPT; i++){
        if (!checkinbz_hex2(lat.K[i])){
            for (int c=0; c<27; c++){
                if (checkinbz_hex2(lat.K[i]+shiftindex[c][0]*b1+shiftindex[c][1]*b2)){
                    lat.K[i] += shiftindex[c][0]*b1+shiftindex[c][1]*b2;
                    break;
                }
            }
        }
    }

    // uniform weight of k-points
    lat.KW.reserve(NKPT);
    for (int k=0; k<NKPT; k++){
        lat.KW[k] = 1.0/NKPT;
    }
}

/*
    generate k-vectors for fcc crystal
*/
void generate_kpt_fcc2 (lattice &lat){
    double tol = 1e-6;
    int N = lat.kpointorder;

    // reciprocal lattice vectors
    Eigen::Vector3d b1 = lat.reciprocal[0];
    Eigen::Vector3d b2 = lat.reciprocal[1];
    Eigen::Vector3d b3 = lat.reciprocal[2];

    // generate fractional coordinate grids
    double df = 1.0/N;
    // NKPT = pow(N,3);
    int NNN = pow(2*N+1,3);
    std::vector<Eigen::Vector3d> K;
    K.reserve(NNN);
    for (int i1=-N; i1<=N; i1++){for (int i2=-N; i2<=N; i2++){for (int i3=-N; i3<=N; i3++){
        K.push_back(df*i1*b1+df*i2*b2+df*i3*b3);
    }}}

    std::vector<double> KW;
    KW.reserve(NNN);

    double NKPT0 = 0;
    NKPT = 0;
    for (int i=0; i<NNN; i++){
        if (!checkinbz_fcc(K[i])){
            KW[i] = 0;
        }else if (checkonbz_fcc(K[i])){
            KW[i] = 0.5;
        }else{
            KW[i] = 1;
        }
        if (KW[i]!=0){
            NKPT += 1;
            NKPT0 += KW[i];
        }
    }

    lat.K.reserve(NKPT);
    lat.KW.reserve(NKPT);
    for (int i=0; i<NNN; i++){
        if (KW[i]!=0){
            lat.K.push_back(K[i]);
            lat.KW.push_back(KW[i]/NKPT0);
        }
    }
}

/*
    check if K is inside the 1st BZ of fcc
*/
bool checkinbz_fcc (Eigen::Vector3d K){
    double tol = 1e-6;

    if (abs(K(0))>1.0+tol || abs(K(1))>1.0+tol || abs(K(2))>1.0+tol || abs(K(0))+abs(K(1))+abs(K(2))>1.5+tol){
        return false;
    }else{
        return true;
    }

}

/*
    check if K is inside the 1st BZ of hex
*/
bool checkinbz_hex(Eigen::Vector3d K, double f){
    double tol = 1e-6;

    if (abs(K(2))>0.5/f+tol || abs(K(1))>1.0/sqrt(3)+tol || abs(K(0))+abs(K(1))/sqrt(3)>1.0/sqrt(3)+tol){
        return false;
    }else{
        return true;
    }
}

/*
    check if K is inside the 1st BZ of 2D hex
*/
bool checkinbz_hex2(Eigen::Vector3d K){
    double tol = 1e-6;

    if (abs(K(1))>1.0/sqrt(3)+tol || abs(K(0))+abs(K(1))/sqrt(3)>1.0/sqrt(3)+tol){
        return false;
    }else{
        return true;
    }
}

/*
    check if K is on the edge of the 1st BZ (fcc lattice)
*/
bool checkonbz_fcc (Eigen::Vector3d K){
    double tol = 1e-6;

    if (abs(abs(K(0))-1.0)<tol || abs(abs(K(1))-1.0)<tol || abs(abs(K(2))-1.0)<tol 
                                    || abs(abs(K(0))+abs(K(1))+abs(K(2))-1.5)<tol){
        return true;
    }else{
        return false;
    }
}

/*
    check if K is on the edge of the 1st BZ (hexagonal lattice)
*/
bool checkonbz_hex(Eigen::Vector3d K, double f){
    double tol = 1e-6;

    if (abs(abs(K(2))-0.5/f)<tol || abs(abs(K(1))-1.0/sqrt(3))<tol || abs(abs(K(0))+abs(K(1))/sqrt(3)-1.0/sqrt(3))<tol){
        return true;
    }else{
        return false;
    }
}

/*
    returns high-symmetry point vector from symbol (fcc lattice)
*/
Eigen::Vector3d highsympoint_fcc (std::string s){
    const static std::unordered_map<std::string,int> bandind{
        {"G",0},{"X",1},{"L",2},{"W",3},{"U",4},{"K",5}
    };
    Eigen::Vector3d k;
    switch (bandind.at(s)){
        case 0:// Gamma-point
            k <<    0,  0,  0; break;
        case 1:// X
            k <<    1,  0,  0; break;
        case 2:// L
            k <<    0.5,0.5,0.5; break;
        case 3:// W
            k <<    1,  0.5,0; break;
        case 4:// U
            k <<    1,  0.25,0.25; break;
        case 5:// K
            k <<    0.75,0.75,0; break;
    };
    return k;
}

/*
    returns high-symmetry point vector from symbol (hexagonal lattice)
*/
Eigen::Vector3d highsympoint_hex (std::string s, double f){
    const static std::unordered_map<std::string,int> bandind{
        {"G",0},{"A",1},{"H",2},{"K",3},{"L",4},{"M",5}
    };
    Eigen::Vector3d k;
    switch (bandind.at(s)){
        case 0:// Gamma-point
            k <<    0,  0,  0; break;
        case 1:// A
            k <<    0,  0,  0.5/f; break;
        case 2:// H
            k <<    2.0/3.0,  0,  0.5/f; break;
        case 3:// K
            k <<    2.0/3.0,  0,  0; break;
        case 4:// L
            k <<    0.5,-0.5/sqrt(3),0.5/f; break;
        case 5:// M
            k <<    0.5,-0.5/sqrt(3),0; break;
    };
    return k;
}

/*
    returns high-symmetry point vector from symbol (hex2d)
*/
Eigen::Vector3d highsympoint_hex2d (std::string s){
    const static std::unordered_map<std::string,int> bandind{
        {"G",0},{"M",1},{"K",2},{"Kp",3}
    };
    Eigen::Vector3d k;
    switch (bandind.at(s)){
        case 0:// Gamma-point
            k <<    0,  0,  0; break;
        case 1:// M
            k <<    1.0/3.0,  0,  0; break;
        case 2:// K
            k <<    1.0/3.0,  1.0/3.0/sqrt(3),  0; break;
        case 3:// Kp
            k <<    1.0/3.0,  -1.0/3.0/sqrt(3),  0; break;
    };
    return k;
}


//==============================================================//
//                  polarization unit vector                    //
//==============================================================//

/*
    returns the longitudinal polarization unit vector parallel to the momentum
*/
inline Eigen::Vector3d uvec_L (Eigen::Vector3d v){
    if (v.norm()==0){
        v << 1.0, 0.0, 0.0;
    }else{
        v = v/v.norm();
    }
    return v;
}

/*
    returns the transverse polarization unit vector tangential to the momentum
    analogous to the typical polarization unit vector in free space (of tangential fields)
*/
inline Eigen::Vector3d uvec_T (Eigen::Vector3d v){
    Eigen::Vector3d p;
    if (v(1)!=0){
        p << -v(1), v(0), 0;
    }else if (v(1)==0 && v(2)!=0){
        p << v(2), 0, -v(0);
    }else if (v(1)==0 && v(2)==0){
        p << 0, 1, 0;
    }
    p = p/p.norm();

    return p;
}

/*
    returns the longitudinal and transverse polarization unit vectors
*/
std::vector<Eigen::Vector3cd> uvec_LT (Eigen::Vector3d v){
    std::vector<Eigen::Vector3cd> uvec;
    uvec.reserve(3);

    // longitudinal polarization unit vector
    uvec[0] = uvec_L(v);
    
    // transverse polarization unit vectors
    Eigen::Vector3d tvec, pvec;
    tvec = uvec_T(v);
    pvec = (uvec[0].real()).cross(tvec);

    bool COMPLEX = true;

    if (COMPLEX){
        uvec[1] = (tvec+im*pvec)/sqrt(2);
        uvec[2] = (tvec-im*pvec)/sqrt(2);
    }else{
        uvec[1] = tvec;
        uvec[2] = pvec;
    }

    return uvec;
}

/*
    returns the Cartesian polarization unit vectors (x,y,z)
*/
// std::vector<Eigen::Vector3cd> uvec_xyz(){
//     std::vector<Eigen::Vector3cd> uvec;
//     uvec.reserve(3);
    
// }













// inline Eigen::Vector3d normvec(Eigen::Vector3d v){
//     if (v.norm()==0){
//         v << 0.0, 0.0, 1.0;
//     }else{
//         v = v/v.norm();
//     }
//     return v;
// }









// // Original?
// inline Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(1)!=0){
//         p << -v(1), v(0), 0;
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)!=0){
//         p << v(2), 0, -v(0);
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 1, 0;
//         // p = p/p.norm();
//     }
//     return p;
// }

// // G-X (100) direction!!!
// inline Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(2)!=0){
//         p << 0, -v(2), v(1);
//     }else if (v(2)==0 && v(1)!=0){
//         p << -v(1), v(0), 0;
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 0, 1;
//     }
//     p = p/p.norm();
//     return p;
// }




// // Original?
// inline Eigen::Vector3d tanvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(1)!=0){
//         p << -v(1), v(0), 0;
//     }else if (v(1)==0 && v(2)!=0){
//         p << v(2), 0, -v(0);
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 1, 0;
//     }
//     p = p/p.norm();

//     return p;
// }






// inline Eigen::Vector3d tanvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(2)!=0){
//         p << 0, -v(2), v(1);
//     }else if (v(2)==0 && v(1)!=0){
//         p << -v(1), v(0), 0;
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 1, 1;
//     }
//     p = p/p.norm();

//     return p;
// }



// // y-axis
// inline Eigen::Vector3d tanvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(2)!=0){
//         p << v(2), 0, -v(0);
//     }else if (v(2)==0 && v(1)!=0){
//         p << -v(1), v(0), 0;
//     }else if (v(2)==0 && v(1)==0){
//         p << 0, 1, 0;
//     }
//     p = p/p.norm();
//     return p;
// }


// // test
// inline Eigen::Vector3d tanvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     double tx,ty,tz;
//     if (v(2)==0 && v(2)==0){
//         p << 0,1,0;
//     }else{
//         tx = -(pow(v(1),2)+pow(v(2),2));
//         ty = v(1)*v(0);
//         tz = v(2)*v(0);
//         p << tx,ty,tz;
//         if (v(2)<0) p = -p;
//     }
//     p = p/p.norm();

//     return p;
// }






// // tangential vector
// Eigen::Vector3d tanvec(Eigen::Vector3d v, Eigen::Vector3d q){
//     Eigen::Vector3d p;
//     double tol = 1e-7;

//     if (abs(q(1))<tol && abs(q(2))<tol){ // [100]
//         if (v(2)!=0){
//             p << 0, -v(2), v(1); 
//         }else if (v(2)==0 && v(1)!=0){
//             p << -v(1), v(0), 0;
//         }else if (v(1)==0 && v(2)==0){
//             p << 0, 0, 1;
//         }
//     }else if (abs(q(0)-q(1))<tol && abs(q(1)-q(2))<tol && abs(q(2)-q(0))<tol){ // [111]
//         if (v(1)!=0){
//             p << -v(1), v(0), 0;
//         }else if (v(1)==0 && v(2)!=0){
//             p << v(2), 0, -v(0);
//         }else if (v(1)==0 && v(2)==0){
//             p << 0, 1, 0;
//         }
//     }else if (abs(q(0)-q(1))<tol && abs(q(2))<tol){ // [110]
//         if (v(1)!=0){
//             p << -v(1), v(0), 0;
//         }else if (v(1)==0 && v(2)!=0){
//             p << v(2), 0, -v(0);
//         }else if (v(1)==0 && v(2)==0){
//             p << 0, 1, 0;
//         }
//     }else{
//         if (v(1)!=0){
//             p << -v(1), v(0), 0;
//         }else if (v(1)==0 && v(2)!=0){
//             p << v(2), 0, -v(0);
//         }else if (v(1)==0 && v(2)==0){
//             p << 0, 1, 0;
//         }
//     }
//     p = p/p.norm();
//     return p;
// }

// // tangential vector
// Eigen::Vector3d tanvecq(Eigen::Vector3d nvec, Eigen::Vector3d qvec){
//     Eigen::Vector3d tvec;
//     Eigen::Vector3d tmp;
//     double tol = 1e-7;
//     double tx, ty, tz, K;

//     if (abs(qvec(1))<tol && abs(qvec(2))<tol){// (100)
//         if (abs(nvec(1))<tol && abs(nvec(2))<tol){
//             tvec << 0,0,1;
//         }else{
//             tmp << 1,0,0;
//             tvec = nvec.cross(tmp);
//         }
//     }else if (abs(qvec(0)-qvec(1))<tol && abs(qvec(0)-qvec(2))<tol){// (111)
//         if (abs(nvec(1)-nvec(2))<tol && abs(nvec(0)-nvec(2))<tol){
//             tmp << 1,0,0;
//             tvec << nvec.cross(tmp);
//         }else{
//             tmp << 1,1,1;
//             tvec = nvec.cross(tmp);
//         }
//     }else if (abs(qvec(0)-qvec(1))<tol && abs(qvec(2))<tol){// (110)
//         if (abs(nvec(0)-nvec(1))<tol && abs(nvec(2))<tol){
//             tvec << 0,0,1;
//         }else{
//             tmp << 1,1,0;
//             tvec = nvec.cross(tmp);
//         }
//     }else{
//         // if (v(1)!=0){
//         //     p << -v(1), v(0), 0;
//         // }else if (v(1)==0 && v(2)!=0){
//         //     p << v(2), 0, -v(0);
//         // }else if (v(1)==0 && v(2)==0){
//         //     p << 0, 1, 0;
//         // }
//         tvec = tanvec(nvec);
//     }

//     // if (abs(v.dot(p))>tol){
//     //     std::cout << "check tangential vector"; abort();
//     // }


//     tvec /= tvec.norm();

//     return tvec;
// }




// // return a perpendicular vector
// inline Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (abs(v(1))<1e-10 && abs(v(2))<1e-10){
//         p << 0, 1, 1;
//         p = p/p.norm();
//     }else if (abs(v(1)-v(2))<1e-10){
//         p << 1, -v(0)/(v(1)+v(2)), -v(0)/(v(1)+v(2));
//         p = p/p.norm();
//     }else if (abs(v(1)+v(2))<1e-10){
//         p << 1, -v(0)/(v(1)-v(2)), v(0)/(v(1)-v(2));
//         p = p/p.norm();
//     }else if (v(2)!=0){
//         p << 0, -v(2), v(1);
//         p = p/p.norm();
//     }else if (v(2)==0 && v(1)!=0){
//         p << -v(1), v(0), 0;
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 0, 1;
//         p = p/p.norm();
//     }
//     return p;
// }



// Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     Eigen::Vector3d p;
//     if (v(2)!=0){
//         p << -v(2), 0, v(0);
//         p = p/p.norm();
//     }else if (v(2)==0 && v(1)!=0){
//         p << v(1), -v(0), 0;
//         p = p/p.norm();
//     }else if (v(1)==0 && v(2)==0){
//         p << 0, 1, 0;
//         // p = p/p.norm();
//     }
//     return p;
// }




// // 2024.03.09 (100) direction
// inline Eigen::Vector3d perpvec(Eigen::Vector3d v){
//     v = v/v.norm(); // normalize
//     Eigen::Vector3d p;
//     double tx, ty, tz, K;
//     if (v(1)==0 && v(2)==0){
//         p << 0, 0, 1;
//     }else {
//         tx = 1;
//         ty = -0.5/v(2)-0.5*v(0)/v(1);
//         tz = 0.5/v(1)-0.5*v(0)/v(2);
//         // if (v(0)>0){
//         //     tx = 1;
//         //     ty = -0.5/v(2)-0.5*v(0)/v(1);
//         //     tz = 0.5/v(1)-0.5*v(0)/v(2);
//         // }else{
//         //     tx = -1;
//         //     ty = 0.5/v(2)+0.5*v(0)/v(1);
//         //     tz = -0.5/v(1)+0.5*v(0)/v(2);
//         // }
//         p << tx, ty, tz;
//     }
//     p = p/p.norm();
//     // p = p.cross(v);
//     return p;
// }


// std::vector<Eigen::Vector3cd> unitvec (Eigen::Vector3d nvec){
//     std::vector<Eigen::Vector3cd> uvec;
//     uvec.reserve(3);

//     // normal vector
//     if (nvec.norm()==0){
//         nvec << 1.0,0.0,0.0;
//         uvec[0] = nvec;
//     }else{
//         nvec = nvec/nvec.norm();
//         uvec[0] = nvec;
//     }
    
//     // tangential vectors
//     Eigen::Vector3d tvec, pvec;
//     tvec = tanvec(nvec);
//     pvec = nvec.cross(tvec);

//     // uvec[1] = tvec;
//     // uvec[2] = pvec;

//     // std::complex<double> im = std::complex<double>(0.0,1.0);
//     uvec[1] = (tvec+im*pvec)/sqrt(2);
//     uvec[2] = (tvec-im*pvec)/sqrt(2);

//     // double tol = 1e-10;
//     // if (abs(nvec.dot(tvec))>tol || abs(nvec.dot(pvec))>tol){
//     //     std::cout << "unnormal vector";
//     //     abort();
//     // }
    
//     return uvec;
// }


// std::vector<Eigen::Vector3cd> unitvecq (Eigen::Vector3d nvec, Eigen::Vector3d qvec){
//     std::vector<Eigen::Vector3cd> uvec;
//     uvec.reserve(3);

//     // normal vector
//     if (nvec.norm()==0){
//         nvec << 1.0,0.0,0.0;
//         uvec[0] = nvec;
//     }else{
//         nvec = nvec/nvec.norm();
//         uvec[0] = nvec;
//     }
    
//     // tangential vectors
//     Eigen::Vector3d tvec, pvec;
//     pvec = tanvecq(nvec,qvec); // ~ y
//     tvec = pvec.cross(nvec); // x = cross(y,z)

//     bool REAL = false;
//     if (REAL){
//         // uvec[1] = tvec;
//         // uvec[2] = pvec;

//         uvec[1] = (tvec-pvec)/sqrt(2);
//         uvec[2] = (tvec+pvec)/sqrt(2);
//     }else{
//         // std::complex<double> im = std::complex<double>(0.0,1.0);
//         uvec[1] = (tvec+im*pvec)/sqrt(2);
//         uvec[2] = (tvec-im*pvec)/sqrt(2);

//         // uvec[1] = (tvec-pvec)/2+im*(tvec+pvec)/2;
//         // uvec[2] = (tvec-pvec)/2-im*(tvec+pvec)/2;
//     }
//     return uvec;
// }

// Eigen::Vector3cd **uvecqm (Eigen::Vector3d Q, std::vector<Eigen::Vector3d> G){
//     Eigen::Vector3cd **uvec;
//     uvec = new Eigen::Vector3cd *[NEPS];
//     for (int m=0; m<NEPS; m++) uvec[m] = new Eigen::Vector3cd [3];
    
//     // normal vector
//     for (int m=0; m<NEPS; m++){
//         uvec[m][0] << Q+G[m];
//         if (uvec[m][0].norm()==0){
//             uvec[m][0] << 1.0,0.0,0.0;
//         }else{
//             uvec[m][0] /= uvec[m][0].norm(); // normalize
//         }
//     }

//     // tangential vectors
//     double tol = 1e-8;
//     double tx, ty, tz, K;
//     double s;
//     if (abs(Q(1))<tol && abs(Q(2))<tol){ // [100]
//         for (int m=0; m<NEPS; m++){
//             if (abs(uvec[m][0](1))<tol && abs(uvec[m][0](2))<tol){
//                 uvec[m][1] << 0,1,0;
//             }else{
//                 uvec[m][1](0) = -(pow(uvec[m][0](1),2)+pow(uvec[m][0](2),2));
//                 uvec[m][1](1) = uvec[m][0](1)*uvec[m][0](0);
//                 uvec[m][1](2) = uvec[m][0](2)*uvec[m][0](0);
//                 uvec[m][1] /= uvec[m][1].norm(); // normalize
//             }
//         }
//     }else if (abs(Q(0)-Q(1))<tol && abs(Q(1)-Q(2))<tol){ // [111]
//         for (int m=0; m<NEPS; m++){
//             if (m==0){
//                 uvec[m][1] << 0.0,-1.0,1.0;
//                 uvec[m][1] /= uvec[m][1].norm();
//             }else if (m==1){
//                 uvec[m][1] << 1.0,0.0,-1.0;
//                 uvec[m][1] /= uvec[m][1].norm();
//             }else if (m==8){
//                 uvec[m][1] << -1.0,1.0,0.0;
//                 uvec[m][1] /= uvec[m][1].norm();
//             }else{
//                 Eigen::Vector3d v = uvec[m][0].real();
//                 if (v(0)*v(1)*v(2)>0){
//                     K = sqrt(3.0*pow(Q(0),2)-2.0*Q(0)+3);
//                     if (v(0)>0){
//                         tx = -(Q(0)-1)/K;
//                         ty = (Q(0)+1+K)/2/K;
//                         tz = (Q(0)+1-K)/2/K;
//                     }else if (v(1)>0){
//                         ty = -(Q(0)-1)/K;
//                         tz = (Q(0)+1+K)/2/K;
//                         tx = (Q(0)+1-K)/2/K;
//                     }else{
//                         tz = -(Q(0)-1)/K;
//                         tx = (Q(0)+1+K)/2/K;
//                         ty = (Q(0)+1-K)/2/K;
//                     }
//                     uvec[m][1] << tx,ty,tz;
//                     uvec[m][1] /= uvec[m][1].norm();
//                 }else {
//                     K = sqrt(3.0*pow(Q(0),2)+2.0*Q(0)+3);
//                     if (v(0)<0){
//                         tx = (Q(0)+1)/K;
//                         ty = -(Q(0)-1+K)/2/K;
//                         tz = (-Q(0)+1+K)/2/K;
//                     }else if (v(1)<0){
//                         ty = (Q(0)+1)/K;
//                         tz = -(Q(0)-1+K)/2/K;
//                         tx = (-Q(0)+1+K)/2/K;
//                     }else{
//                         tz = (Q(0)+1)/K;
//                         tx = -(Q(0)-1+K)/2/K;
//                         ty = (-Q(0)+1+K)/2/K;
//                     }
//                     uvec[m][1] << tx,ty,tz;
//                     uvec[m][1] /= uvec[m][1].norm();
//                 }
//             }
//         }
//     }else if (abs(Q(1)-Q(2))<tol && abs(Q(3)<tol)){// [110]
//         for (int m=0; m<NEPS; m++){
//             Eigen::Vector3d v = uvec[m][0].real();
//             if (m==0){
//                 uvec[m][1] << 0.0,0.0,1.0;
//             }else if (m==1 || m==2 || m==7 || m==8){
//                 s = G[m](0)*G[m](1)*G[m](2);
//                 tx = s/2+s*v(2)/2;
//                 ty = -s/2+s*v(2)/2;
//                 tz = -s*v(0);
//                 uvec[m][1] << tx,ty,tz;
//             }else if (m==4 || m==6){
//                 s = G[m](0)*G[m](1)*G[m](2);
//                 K = sqrt(2*pow(Q(0),2)+3);
//                 tx = -s*(4*Q(0)*Q(0)+s*3*Q(0)+3-3*(Q(0)+s)*K);
//                 ty = s*(4*Q(0)*Q(0)-s*3*Q(0)+3-3*(Q(0)-s)*K);
//                 tz = s*2*(Q(0)*Q(0)+3);
//                 uvec[m][1] << tx,ty,tz;
//             }else{
//                 s = G[m](0)*G[m](1)*G[m](2);
//                 K = sqrt(2*pow(Q(0),2)+3);
//                 tx = s*(4*Q(0)*Q(0)-s*3*Q(0)+3+3*(Q(0)+s)*K);
//                 ty = -s*(4*Q(0)*Q(0)+s*3*Q(0)+3+3*(Q(0)-s)*K);
//                 tz = s*2*(Q(0)*Q(0)+3);
//                 uvec[m][1] << tx,ty,tz;
//             }
//         }
//     }else{
//         for (int m=0; m<NEPS; m++){
//             Eigen::Vector3d T;
//             T << 0.0,0.0,1.0;
//             uvec[m][1] = uvec[m][0].cross(T);
//         }
//     }

//     //
//     Eigen::Vector3cd t,p;
//     // std::complex<double> im = std::complex<double>(0.0,1.0);
//     for (int m=0; m<NEPS; m++){
//         uvec[m][1] /= uvec[m][1].norm();
//         uvec[m][2] = uvec[m][0].cross(uvec[m][1]);

//         t = uvec[m][1];
//         p = uvec[m][2];
//         uvec[m][1] = (t+im*p)/sqrt(2);
//         uvec[m][2] = (t-im*p)/sqrt(2);
//     }



//     for (int m=0; m<NEPS; m++){
//         if (abs(uvec[m][0].dot(uvec[m][1]))>tol){
//             std::cout << "Warning::not normal found\n" << Q << "-" << m << std::endl;
//             std::cout << uvec[m][0] << std::endl;
//             std::cout << uvec[m][1] << std::endl;
//             abort();
//         }
//         if (abs(uvec[m][0].dot(uvec[m][2]))>tol){
//             std::cout << "not normal-" << Q << "-" << m << std::endl;
//             std::cout << uvec[m][0] << std::endl;
//             std::cout << uvec[m][2] << std::endl;
//             abort();
//         }
//     }
//     return uvec;
// }

}/* namespace pmx */

















// recycle bin





// // generate plane waves G for fcc in spherical...
// void generate_planewaves (lattice &lat){
//     std::vector<int> g2 {0,3,4,8,11,12,16,19,20,24,27,32,35,36,40,43,44,48,51,52,56,59,64,67,68,72,75,76,80,83,84,88,91,96,99,100};
//     std::vector<int> nnn {1,8,6,12,24,8,6,24,24,24,32,12,48,30,24,24,24,8,48,24,48,72,6,24,48,36,56,24,24,72,48,24,48,24,72,30};
//     std::vector<int>::iterator itr = std::lower_bound(g2.begin(), g2.end(), pow(lat.gcut,2));
//     int N = int(itr - g2.begin());
//     std::vector<int> cnn(N);

//     NPW = 0;
//     for (int i=0; i<N; i++){
//         NPW += nnn[i];
//         if (i==0) cnn[i] = 1;
//         else cnn[i] = cnn[i-1]+nnn[i];
//     }

//     lat.G.reserve(NPW);
//     for (int i=0; i<NPW; i++){
//         lat.G.push_back({0,0,0});
//     }

//     #pragma omp parallel for
//     for (int i=0; i<N; i++){
//         std::vector<Eigen::Vector3d> iG = find_G(g2[i],nnn[i]);
//         int in;
//         if (i==0) in = 0;
//         else in = cnn[i-1];
//         for (int ic=0; ic<nnn[i]; ic++){
//             lat.G[in+ic] = iG[ic];
//         }
//     }
//     #pragma omp barrier

//     // find local-field index
//     find_locGi(lat);

//     return;
// }

// // find n reciprocal lattice vectors with squared magnitude g2
// std::vector<Eigen::Vector3d> find_G (int g2, int n){
//     std::vector<Eigen::Vector3d> G(n);

//     int ia = int(sqrt(double(g2)));
//     int gx2, gy2, gz2, gxy2, gxyz2;
//     int ic = 0;

//     for (int gx=-ia; gx<=ia; gx++){
//         gx2 = gx*gx;
//         for (int gy=-ia; gy<=ia; gy++){
//             gy2 = gy*gy;
//             gxy2 = gx2+gy2;
//             for (int gz=-ia; gz<=ia; gz++){
//                 gz2 = gz*gz;
//                 gxyz2 = gxy2+gz2;
//                 if (gxyz2==g2){
//                     G[ic] << double(gx), double(gy), double(gz);
//                     ic +=1;
//                     if (ic==n) return G;
//                 }
//             }
//         }
//     }
//     return G;
// }

// /*
//     find local-field index i, such that Gi=Gp+Gm

//     idx_offset [NEPS x NPW]

//     index p: planewave summation index
//     Gm: offset wavevector (local-field effect)
// */ 
// void find_locGi (lattice &lat){

//     // initialize
//     lat.loci = new int *[NEPS];
//     for (int m=0; m<NEPS; m++){
//         lat.loci[m] = new int [NPW];
//     }

//     // g(0)
//     for (int p=0; p<NPW; p++){
//         lat.loci[0][p] = p;
//     }
//     #pragma omp parallel for
//     for (int m=1; m<NEPS; m++){

//         // local-field vector Gm
//         Eigen::Vector3d Gm = lat.G[m];
//         for (int p=0; p<NPW; p++){
//             // bra vector Gp
//             Eigen::Vector3d Gp = lat.G[p];

//             // if no such vector found, -1
//             lat.loci[m][p] = -1;
//             for (int i=0; i<NPW; i++){
//                 // ket vector Gi
//                 Eigen::Vector3d Gi = lat.G[i];

//                 // find i, such that Gi = Gp+Gm
//                 if ((Gp-Gi+Gm).norm()<1e-10){
//                     lat.loci[m][p] = i;
//                     break;
//                 }
//             }
//         }
//     }
//     #pragma omp barrier
//     return;
// }