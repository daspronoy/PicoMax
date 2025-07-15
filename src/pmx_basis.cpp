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
    // Use a map to handle unique k-points and avoid duplicates from symmetry operations.
    // The key is a tuple representing the k-vector to handle floating point comparisons.
    std::map<std::tuple<double, double, double>, double> unique_k_weights;
    double tol = 1e-9;

    // 1. Populate the map with original k-points and their weights.
    for (size_t i = 0; i < lat.K.size(); ++i) {
        Eigen::Vector3d k = lat.K[i];
        double w = lat.KW[i];
        std::tuple<double, double, double> k_key = {
            std::round(k.x() / tol) * tol,
            std::round(k.y() / tol) * tol,
            std::round(k.z() / tol) * tol
        };
        unique_k_weights[k_key] += w;
    }

    // 2. Create a new map for the symmetrized points.
    std::map<std::tuple<double, double, double>, double> symmetrized_weights;
    std::set<std::tuple<double, double, double>> processed_keys;

    for (const auto& pair : unique_k_weights) {
        const auto& k_key = pair.first;

        // Skip if this k-point has already been processed as part of a symmetric pair.
        if (processed_keys.count(k_key)) {
            continue;
        }

        Eigen::Vector3d k = {std::get<0>(k_key), std::get<1>(k_key), std::get<2>(k_key)};
        double w = pair.second;

        // Generate the x<->y swapped counterpart.
        Eigen::Vector3d k_swap_vec = {k.y(), k.x(), k.z()};
        std::tuple<double, double, double> k_swap_key = {
            std::round(k_swap_vec.x() / tol) * tol,
            std::round(k_swap_vec.y() / tol) * tol,
            std::round(k_swap_vec.z() / tol) * tol
        };

        // Mark both as processed.
        processed_keys.insert(k_key);
        processed_keys.insert(k_swap_key);

        if (k_key == k_swap_key) {
            // The point is on the symmetry line (x=y), so it's its own pair.
            symmetrized_weights[k_key] += w;
        } else {
            // The point has a distinct symmetric partner.
            // Check if the partner also exists in the original grid.
            double w_swap = 0.0;
            if (unique_k_weights.count(k_swap_key)) {
                w_swap = unique_k_weights[k_swap_key];
            }
            
            // The combined weight is distributed equally between the two points.
            double combined_weight = w + w_swap;
            symmetrized_weights[k_key] += combined_weight / 2.0;
            symmetrized_weights[k_swap_key] += combined_weight / 2.0;
        }
    }

    // 3. Reconstruct the lattice's k-point and weight vectors from the symmetrized map.
    lat.K.clear();
    lat.KW.clear();
    double total_weight = 0.0;

    for (const auto& pair : symmetrized_weights) {
        lat.K.emplace_back(std::get<0>(pair.first), std::get<1>(pair.first), std::get<2>(pair.first));
        lat.KW.push_back(pair.second);
        total_weight += pair.second;
    }

    // 4. Normalize the final weights.
    if (total_weight > 0) {
        for (double& w : lat.KW) {
            w /= total_weight;
        }
    }
    
    // 5. Update the global k-point count.
    NKPT = lat.K.size();
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
    // For Tellurium (hexagonal), define a set of fixed, high-symmetry candidate
    // vectors for the transverse direction. This avoids arbitrary basis construction.
    // These vectors correspond to the primary crystal axes/directions.
    // Note: The coordinates for K and M are based on a common convention for hexagonal lattices.
    const Eigen::Vector3d dir_M(0.5, -0.5 / sqrt(3.0), 0.0); // Direction towards M point
    const Eigen::Vector3d dir_K(2.0 / 3.0, 0.0, 0.0);         // Direction towards K point
    const Eigen::Vector3d dir_A(0.0, 0.0, 1.0);               // Direction towards A point (c-axis)

    // Normalize the input vector to get the longitudinal direction.
    Eigen::Vector3d v_norm = v;
    if (v.norm() > 1e-9) {
        v_norm.normalize();
    } else {
        // If v is the zero vector (Gamma point), default longitudinal to a high-symmetry direction.
        // We choose K here, but any choice is arbitrary; the key is consistency.
        v_norm = dir_K.normalized();
    }

    // Find which candidate vector is most perpendicular to v.
    // The one with the smallest absolute dot product is the most transverse.
    double dot_M = std::abs(v_norm.dot(dir_M.normalized()));
    double dot_K = std::abs(v_norm.dot(dir_K.normalized()));
    double dot_A = std::abs(v_norm.dot(dir_A.normalized()));

    Eigen::Vector3d p; // This will be the chosen transverse vector.

    // Select the high-symmetry direction that is most orthogonal to v
    if (dot_M <= dot_K && dot_M <= dot_A) {
        p = dir_M;
    } else if (dot_K <= dot_M && dot_K <= dot_A) {
        p = dir_K;
    } else {
        p = dir_A;
    }

    // // Project out the longitudinal component from the chosen direction
    // // to make it perfectly orthogonal to v.
    // p = p - v_norm.dot(p) * v_norm;

    // // Normalize the final transverse vector.
    // // This should not fail unless v is perfectly parallel to all three candidate
    // // directions, which is geometrically impossible.
    // if (p.norm() < 1e-9) {
    //     // Fallback for the highly unlikely case that p is a zero vector.
    //     // This can happen if v is perfectly aligned with one of the candidate vectors.
    //     // In this case, we can just pick another candidate.
    //     if (dot_M > 0.999) { // v is along M
    //         p = dir_K - v_norm.dot(dir_K) * v_norm;
    //     } else { // v is along K or A
    //         p = dir_M - v_norm.dot(dir_M) * v_norm;
    //     }
    // }
    
    p.normalize();

    return p;
}

// inline Eigen::Vector3d uvec_T (Eigen::Vector3d v){
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
}

/*
    returns the Cartesian polarization unit vectors (x,y,z)
*/
// std::vector<Eigen::Vector3cd> uvec_xyz(){
//     std::vector<Eigen::Vector3cd> uvec;
//     uvec.reserve(3);
    
// }









