




// // kramers-kronig
// void kramerskronig(env &dat){

//     double realepsL1;
//     double realepsT1;
//     double ds = dat.dfreq*2.0/pi;

    
//     for (int m=0; m<NEPS; m++){ for(int n=0; n<=m; n++){
//         // #pragma omp parallel for
//         for (int f=0; f<NFREQ; f++){
//             realepsL1 = 0.0;
//             realepsT1 = 0.0;
//             for (int s=0; s<NFREQ; s++){
//                 if (f!=s){
//                     realepsL1 += dat.freq[s]*dat.imagepsL[m][n][s]/(dat.freq[s]*dat.freq[s]-dat.freq[f]*dat.freq[f]);
//                     realepsT1 += dat.freq[s]*dat.imagepsT[m][n][s]/(dat.freq[s]*dat.freq[s]-dat.freq[f]*dat.freq[f]);
//                 }
//             }
//             realepsL1 *= ds;
//             realepsT1 *= ds;
//             if (m==n){
//                 realepsL1 += 1.0;
//                 realepsT1 += 1.0;
//             }
//             dat.realepsL[m][n][f] = realepsL1;
//             dat.realepsT[m][n][f] = realepsT1;
            
//         }
        
//     }}
//     // #pragma omp barrier

//     return;
// }





    // // reciprocal lattice vectors b1,b2,b3
    // if (inp.existOption("-b1")){
    //     std::vector<std::string> v;
    //     v = inp.vectorOption("-b1");
    //     for (int i=0; i<v.size(); i++){
    //         dat.b[0](i) = std::stod(v[i]);
    //     }
    // }
    // if (inp.existOption("-b2")){
    //     std::vector<std::string> v;
    //     v = inp.vectorOption("-b2");
    //     for (int i=0; i<v.size(); i++){
    //         dat.b[1](i) = std::stod(v[i]);
    //     }
    // }
    // if (inp.existOption("-b3")){
    //     std::vector<std::string> v;
    //     v = inp.vectorOption("-b3");
    //     for (int i=0; i<v.size(); i++){
    //         dat.b[2](i) = std::stod(v[i]);
    //     }
    // }

    // // number of plane waves
    // if (inp.existOption("-npw")){
    //     NPW = std::stoi(inp.valueOption("-npw"));
    // }




    
//======================================================//
//                 density of states                    //
//======================================================//

// Calculate density of states
// void dos(env &dat){
//     double epsilon = 0.001;

//     std::ofstream dosfile;
// 	dosfile.open("dos.out");
// 	dosfile.setf(std::ios::left); //setf sets format flags
// 	dosfile.setf(std::ios::scientific);
// 	dosfile.setf(std::ios::showpoint);
// 	dosfile.setf(std::ios::showpos);
// 	dosfile.precision(8);
// 	double energy = -5.0;
// 	double dosvalue;
//     int N = dat.numpw;
// 	Eigen::MatrixXcd A(N, N); //Global Matrix

//     std::complex<double> eigenvec[dat.numk][30][N];
// 	Eigen::MatrixXd eigenvalues(dat.numk, 30);
// 	int i, j, k, n; //dummy indices

//     setRefEnergy(dat);

//     Eigen::Vector3d kvector;
//     for (int k=0; k<dat.numk; k++){
//         kvector << dat.kgrid[k][0], dat.kgrid[k][1], dat.kgrid[k][2];

//         for(i=-N/2; i<N/2; i++){for(j=-N/2; j<N/2; j++){
// 			A(i+N/2, j+N/2) = matrixelement(i,j,N,kvector,dat.b,dat.epm,dat.a); 		
// 		}}

//         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(A);
// 		if (eigensolver.info() !=Eigen::Success){ std::cout << "ERROR: Eigensolver failed to solve the matrix." << std::endl; abort(); }
// 		for(n = 0; n < 30; n++){
// 			eigenvalues(k, n) = eigensolver.eigenvalues()(n) - dat.energyoffset;	

// 			for(j=-N/2; j<N/2; j++)	eigenvec[k][n][j+N/2] = eigensolver.eigenvectors().col(n)(j+N/2);
// 		}//end of n
//     }
// }
