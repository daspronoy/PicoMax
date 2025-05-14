Calculate density of states
void dos(env &dat){
    double epsilon = 0.001;

    std::ofstream dosfile;
	dosfile.open("dos.out");
	dosfile.setf(std::ios::left); //setf sets format flags
	dosfile.setf(std::ios::scientific);
	dosfile.setf(std::ios::showpoint);
	dosfile.setf(std::ios::showpos);
	dosfile.precision(8);
	double energy = -5.0;
	double dosvalue;
    int N = dat.numpw;
	Eigen::MatrixXcd A(N, N); //Global Matrix

    std::complex<double> eigenvec[dat.numk][30][N];
	Eigen::MatrixXd eigenvalues(dat.numk, 30);
	int i, j, k, n; //dummy indices

    setRefEnergy(dat);

    Eigen::Vector3d kvector;
    for (int k=0; k<dat.numk; k++){
        kvector << dat.kgrid[k][0], dat.kgrid[k][1], dat.kgrid[k][2];

        for(i=-N/2; i<N/2; i++){for(j=-N/2; j<N/2; j++){
			A(i+N/2, j+N/2) = matrixelement(i,j,N,kvector,dat.b,dat.epm,dat.a); 		
		}}

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(A);
		if (eigensolver.info() !=Eigen::Success){ std::cout << "ERROR: Eigensolver failed to solve the matrix." << std::endl; abort(); }
		for(n = 0; n < 30; n++){
			eigenvalues(k, n) = eigensolver.eigenvalues()(n) - dat.energyoffset;	

			for(j=-N/2; j<N/2; j++)	eigenvec[k][n][j+N/2] = eigensolver.eigenvectors().col(n)(j+N/2);
		}//end of n
    }
}
