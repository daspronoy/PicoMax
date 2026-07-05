/*
    pmx_io.cpp | includes picomax initialization, input, output functionalities

    Author: Jungho Mun
    Date: June 6, 2024
*/

/*
    header files
*/
#include "pmx_io.hpp"
#include "pmx_basis.hpp"


/*
    function definitions
*/
namespace pmx {

/*
    print documentation of picomax
*/
void doc(){
    std::cout   
        << "============================================================\n"
        << "  PicoMax version 0.1\n"
        << "Usage:\n"
        << "  picomax [-option] [value]\n"
        << "============================================================\n"
        << "-h, -help : print help\n"
        << "-v, -version : print PicoMax version\n"
        << "-debug : debug level {0}\n"
        << "-outputfile : output filename\n"
        << "-wdir : working directory {pwd}\n"
        << "-logfile : print output to a logfile {0}\n"
        << "\n"
        << "-switch : calculation switch\n"
        << "    band, epsij, epsLL\n"
        << "\n"
        << "-crystal : crystal structure\n"
        << "    zb, wz, rs, te\n"
        << "-dim : crystal dimension {3}\n"
        << "-a : lattice constant [angstrom]\n"
        << "-f : lattice parameter\n"
        << "-u : atomic basis parameter\n"
        << "-vg : empirical pseudopotential form factors\n"
        << "    ex) ...\n"
        << "-encut : energy cutoff [eV]\n"
        << "-gcut : g-vector cutoff [2*pi/a]\n"
        << "-kpoint : switch for kpoint grids\n"
        << "    0->import kpointfile\n"
        << "    1->full BZ\n"
        << "    2->\n"
        << "-kpointorder : kpoint density NxNxN\n"
        << "-nband : number of electronic bands\n"
        << "-nchi : the size of dielectric matrix NxN\n"
        << "-nvalence : the number of valence electron bands\n"
        << "\n"
        << "-qvec : q-vectors\n"
        << "-qpath : path of q-vectors\n"
        << "         {L,G,X,W,K,G}\n"
        << "-qnum : number of q-vectors in each path\n"
        << "         {20,20,10,10,20}\n"
        << "-refpoint : the reference q-vector for electronic band energy {0,0,0}\n"
        << "-gammazero : the Gamma-point {0,0,0}\n"
        << "    * move away from the zero to avoid singularity for epsLL calculation\n"
        << "\n"
        << "-nfreq : number of frequency points {1}\n"
        << "-dfreq : frequency interval [eV]\n"
        << "-delta : delta function model {0}\n"
        << "-epsilon : delta function smearing [eV] {0.05}\n"
        << "-kk : switch for Kramers-Kronig transform {0}\n"
        << "\n"
        << "============================================================\n"
        << std::endl;
    // std::cout << "\n";
    // std::cout << "\n";
    // std::cout << "\n";
    return;
}

/*
    print version of picomax
*/
void ver(){
    std::cout << "PicoMax version v0.1\n";
    std::cout << "  Date: 2024.06.06\n";
    return;
}

// progress messaging
void progmsg(){
    
}


/*
    initialize lattice parameters
*/
void init_crystal(inputParser inp, lattice &lat){

    inp.fill_int(inp,"-dim",lat.dim);

    if (inp.existOption("-crystal")){
        lat.name = inp.valueOption("-crystal");
    }

    // lattice constant
    if (inp.existOption("-a")){
        a = std::stod(inp.valueOption("-a"));
    }else{
        std::cout << "lattice constant (a) is a mandatory argument!" << std::endl; abort();
    }
    inp.fill_double(inp,"-f",lat.f);// c/a lattice constant for hexagonal
    inp.fill_double(inp,"-u",lat.u);// 

    // default parameters for some crystal structures
    if (lat.name.compare("zb")==0){// zinc-blende, sphalerite
        // lattice vectors in units of [a]
        lat.lattice.reserve(3);
        lat.lattice.push_back({0,0.5,0.5});
        lat.lattice.push_back({0.5,0,0.5});
        lat.lattice.push_back({0.5,0.5,0});

        // reciprocal lattice vectors in units of [2*pi/a]
        lat.reciprocal.reserve(3);
        lat.reciprocal.push_back({-1,1,1});
        lat.reciprocal.push_back({1,-1,1});
        lat.reciprocal.push_back({1,1,-1});
        lat.bz_volume = 4.0;

        // atomic basis vectors in units of [a]
        NATOM = 2;
        lat.atomic.reserve(NATOM);
        lat.atomic.push_back({0,0,0});
        lat.atomic.push_back({0.25,0.25,0.25});
        // lat.atomic.push_back({-0.125,-0.125,-0.125});
        // lat.atomic.push_back({0.125,0.125,0.125});

    }else if (lat.name.compare("rs")==0){// rock-salt, halite
        // lattice vectors in units of [a]
        lat.lattice.reserve(3);
        lat.lattice.push_back({0,0.5,0.5});
        lat.lattice.push_back({0.5,0,0.5});
        lat.lattice.push_back({0.5,0.5,0});

        // reciprocal lattice vectors in units of [2*pi/a]
        lat.reciprocal.reserve(3);
        lat.reciprocal.push_back({-1,1,1});
        lat.reciprocal.push_back({1,-1,1});
        lat.reciprocal.push_back({1,1,-1});
        lat.bz_volume = 4.0;

        // atomic basis vectors in units of [a]
        NATOM = 2;
        lat.atomic.reserve(NATOM);
        lat.atomic.push_back({0,0,0});
        lat.atomic.push_back({0.5,0,0});

    }else if (lat.name.compare("wz")==0){// wurtzite
        // lattice vectors in units of [a]
        lat.lattice.reserve(3);
        lat.lattice.push_back({0.5,-0.5*sqrt(3),0});
        lat.lattice.push_back({0.5,0.5*sqrt(3),0});
        lat.lattice.push_back({0,0,lat.f});

        // reciprocal lattice vectors in units of [2*pi/a]
        lat.reciprocal.reserve(3);
        lat.reciprocal.push_back({1,-1/sqrt(3),0});
        lat.reciprocal.push_back({1,1/sqrt(3),0});
        lat.reciprocal.push_back({0,0,1.0/lat.f});
        lat.bz_volume = abs((lat.reciprocal[0].cross(lat.reciprocal[1])).dot(lat.reciprocal[2]));

        // atomic basis vectors in units of [a]
        NATOM = 4;
        lat.atomic.reserve(NATOM);
        lat.atomic.push_back({0.5,0.5/sqrt(3),0});
        lat.atomic.push_back({0.5,-0.5/sqrt(3),0.5*lat.f});
        lat.atomic.push_back({0.5,0.5/sqrt(3),lat.u*lat.f});
        lat.atomic.push_back({0.5,-0.5/sqrt(3),(0.5+lat.u)*lat.f});

    }else if (lat.name.compare("te")==0){// tellurium
        // lattice vectors in units of [a]
        lat.lattice.reserve(3);
        lat.lattice.push_back({0.5,-0.5*sqrt(3),0});
        lat.lattice.push_back({0.5,0.5*sqrt(3),0});
        lat.lattice.push_back({0,0,lat.f});

        // reciprocal lattice vectors in units of [2*pi/a]
        lat.reciprocal.reserve(3);
        lat.reciprocal.push_back({1,-1/sqrt(3),0});
        lat.reciprocal.push_back({1,1/sqrt(3),0});
        lat.reciprocal.push_back({0,0,1.0/lat.f});
        lat.bz_volume = abs((lat.reciprocal[0].cross(lat.reciprocal[1])).dot(lat.reciprocal[2]));

        // atomic basis vectors in units of [a]
        NATOM = 3;
        lat.atomic.reserve(NATOM);
        // lat.atomic.push_back({-1.5*lat.u,0,0});
        // lat.atomic.push_back({0,0.5*sqrt(3)*lat.u,1.0/3.0*lat.f});
        // lat.atomic.push_back({0,-0.5*sqrt(3)*lat.u,2.0/3.0*lat.f});
        
        lat.atomic.push_back({0,0.5*sqrt(3)*lat.u,-1.0/3.0*lat.f});
        lat.atomic.push_back({1.5*lat.u,0,0});
        lat.atomic.push_back({0,-0.5*sqrt(3)*lat.u,1.0/3.0*lat.f});

    }else if (lat.name.compare("hex")==0){// 2D hexagonal
        lat.dim = 2;
        // lattice vectors in units of [a]
        lat.lattice.reserve(2);
        lat.lattice.push_back({0.5*3,-0.5*sqrt(3),0});
        lat.lattice.push_back({0.5*3,0.5*sqrt(3),0});

        // reciprocal lattice vectors in units of [2*pi/a]
        lat.reciprocal.reserve(2);
        lat.reciprocal.push_back({1.0/3.0,-sqrt(3)/3.0,0});
        lat.reciprocal.push_back({1.0/3.0,sqrt(3)/3.0,0});
        lat.bz_volume = abs(lat.reciprocal[0].dot(lat.reciprocal[1]));

        // atomic basis vectors in units of [a]
        NATOM = 2;
        lat.atomic.reserve(NATOM);
        lat.atomic.push_back({0,0,0});
        lat.atomic.push_back({0.5,0,5*sqrt(3)});
    }

    // manual inputs
    // if (inp.existOption("-lattice_vector"))
    // if (inp.existOption("-reciprocal_vector"))
    // if (inp.existOption("-atomic_vector"))
    inp.fill_int(inp,"-gsym",lat.gsym);//

    // encut: planewave energy cutoff [eV]
    inp.fill_double(inp,"-encut",lat.encut);
    lat.gcut = sqrt(lat.encut*pow(a/2.0/pi,2)/KCON);
    inp.fill_double(inp,"-gcut",lat.gcut);// gcut can be directly provided

    // gammazero
    if (inp.existOption("-gammazero")){
        std::vector<std::string> v;
        v = inp.vectorOption("-gammazero");
        if (v.size()!=3){
            std::cout << "check the number of gammazero!" << std::endl; abort();
        }
        for (int i=0; i<3; i++){
            lat.gammazero[i] = std::stod(v[i]);
        }
    }

    // q-vector
    if (inp.existOption("-qvec")){ // direct qvec input
        std::vector<std::string> v;
        v = inp.vectorOption("-qvec");
        if (v.size()%3!=0){
            std::cout << "check the number of qvec!" << std::endl; abort();
        }
        NQ = v.size()/3;
        Eigen::Vector3d qvec;
        for (int i=0; i<v.size(); i=i+3){
            qvec(0) = std::stod(v[i]);
            qvec(1) = std::stod(v[i+1]);
            qvec(2) = std::stod(v[i+2]);
            lat.Q.push_back(qvec);
        }
    }else if ((inp.existOption("-qpath") && !inp.existOption("-qnum"))
            || (!inp.existOption("-qpath") && inp.existOption("-qnum"))){
        std::cout << "specify both -qpath and -qnum!" << std::endl; abort();
    }else if (inp.existOption("-qpath") && inp.existOption("-qnum")){ // generate qvec based on qpath and qnum
        std::vector<std::string> v;
        v = inp.vectorOption("-qpath");
        for (int i=0; i<v.size(); i++){
            lat.qpath.push_back(v[i]);
        }
        v = inp.vectorOption("-qnum");
        for (int i=0; i<v.size(); i++){
            lat.qnum.push_back(std::stoi(v[i]));
        }

        // construct q-vectors
        NQ = std::accumulate(lat.qnum.begin(),lat.qnum.end(),0)+1;
        const static std::unordered_map<std::string,int> bandind{
            {"G",0},{"X",1},{"L",2},{"W",3},{"U",4},{"K",5}
        };
        Eigen::Vector3d v0,v1,dv;
        for (int i=0; i<lat.qnum.size(); i++){
            if (lat.name.compare("zb")==0 || lat.name.compare("rs")==0){// fcc
                v0 = highsympoint_fcc(lat.qpath[i]);
                v1 = highsympoint_fcc(lat.qpath[i+1]);
            }else if (lat.name.compare("wz")==0 || lat.name.compare("te")==0){// hex
                v0 = highsympoint_hex(lat.qpath[i],lat.f);
                v1 = highsympoint_hex(lat.qpath[i+1],lat.f);
            }else if (lat.name.compare("hex")==0){// hexagonal 2d
                v0 = highsympoint_hex2d(lat.qpath[i]);
                v1 = highsympoint_hex2d(lat.qpath[i+1]);
            }
            dv = v1-v0;
            int nqi = lat.qnum[i];
            for (int j=0; j<nqi; j++){
                lat.Q.push_back(v0 + double(j)/double(nqi)*dv);
                if(lat.Q.back().norm()==0) lat.Q.back() += lat.gammazero;
            }
        }
        lat.Q.push_back(v1);
        if(lat.Q.back().norm()==0) lat.Q.back() += lat.gammazero;
    }// q-vector
}

/*
    initialize material parameters
*/
void init_mater(inputParser inp, mater &mat){
    double tol = 1e-2;

    // empirical pseudopotential form factor (local)
    if (inp.existOption("-vg")){
        std::vector<std::string> v;
        v = inp.vectorOption("-vg");

        if (v.size()%(NATOM+1)!=0){
            std::cout << "check the number of arguments for -vg!" << std::endl; abort();
        }

        mat.NV = v.size()/(NATOM+1);// number of Vg

        mat.g2.reserve(mat.NV);

        mat.vg = new double *[NATOM];
        for (int i=0; i<NATOM; i++){
            mat.vg[i] = new double [mat.NV];
        }
        for (int i=0; i<mat.NV; i++){
            mat.g2.push_back(std::stod(v[i]));
            for (int j=0; j<NATOM; j++){
                mat.vg[j][i] = std::stod(v[(j+1)*mat.NV+i]);
            }
        }

        mat.g2max = mat.g2[mat.NV-1]+tol;
    }

    


}


/*
    initialize and process the input arguments
*/
void inputProcess(inputParser inp, env &dat){

    // debug level
    inp.fill_int(inp,"-debug",dat.debug);

    // outputfile
    if (inp.existOption("-outputfile")){
        const char *filename = inp.valueOption("-outputfile").c_str();
        dat.outputfile = new char[strlen(filename)+1];
        strcpy(dat.outputfile,filename);
    }else{ // default outputfile name
        const std::string s = "picomax_output";;
        dat.outputfile = new char [s.length()+1];
        strcpy(dat.outputfile,s.c_str());
    }

    // working directory
    if (inp.existOption("-wdir")){
        const char *wdir = inp.valueOption("-wdir").c_str();
        dat.wdir = new char[strlen(wdir)+1];
        strcpy(dat.wdir,wdir);
    } else { // insert current path as working directory
        const std::string current_path = std::filesystem::current_path().string();
        dat.wdir = new char[current_path.length()+1];
        strcpy(dat.wdir,current_path.c_str());
    }
    
    // logfile
    inp.fill_int(inp,"-logfile",dat.logfile);

    

    // k-grids generation
    if (inp.existOption("-kpointfile")){
        const char *filename = inp.valueOption("-kpointfile").c_str();
        dat.kpointfile = new char[strlen(filename)+1];
        strcpy(dat.kpointfile,filename);
        dat.kpoint = 0;
    }else{
        inp.fill_int(inp,"-kpoint",dat.kpoint);
        inp.fill_int(inp,"-kpointorder",dat.lat.kpointorder);
    }
    
    // number of electronic bands to consider
    inp.fill_int(inp,"-nband",NBAND);

    // size of dielectric matrix
    inp.fill_int(inp,"-neps",NEPS);

    // number of valence bands
    inp.fill_int(inp,"-nvalence",dat.nvalence);

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
    inp.fill_int(inp,"-kk",dat.kk);
    // delta function epsilon
    inp.fill_double(inp,"-epsilon",EPSILON);
    // delta function model
    inp.fill_int(inp,"-delta",dat.delta);

    if (inp.existOption("-refpoint")){
        std::vector<std::string> v;
        v = inp.vectorOption("-refpoint");
        if (v.size()!=3){
            std::cout << "check the number of refpoint!" << std::endl; abort();
        }
        for (int i=0; i<3; i++){
            dat.refpoint(i) = std::stod(v[i]);
        }
    }

    return;
}

/*
    print input parameters to stdout for debugging purposes
    this function will be invoked if `-debug 1`
    if `-debug 2`, all gpoints and kpoint will be printed
*/
void printInput(env &dat){
    
    std::cout   << "working directory: " << dat.wdir << "\n"
                << "outputfile: " << dat.outputfile << "\n"
                << "logfile: " << dat.logfile << "\n"
                << "number of electronic bands: " << NBAND << "\n"
                << "number of valence electrons: " << dat.nvalence << "\n"
                << "size of dielectric matrix: " << NEPS << "\n"
                << "Material properties" << "\n"
                << "lattice constant [angstrom]: " << a << "\n"
                << "Empirical pseudopotential form factors, Vg:" << "\n"
                << "g2: \t";
                for (int i=0; i<dat.mat.NV; i++){
                    std::cout << dat.mat.g2[i] << "\t";
                }
                std::cout << "\n";
                for (int i=0; i<NATOM; i++){
                    std::cout << "atom " << (i+1) << ": ";
                    for (int j=0; j<dat.mat.NV; j++){
                        std::cout << dat.mat.vg[i][j] << "\t";
                    }
                    std::cout << "\n";
                }
    std::cout << "atomic basis vectors:\n";
    for (int i=0; i<NATOM; i++){
        std::cout << " " << dat.lat.atomic[i].transpose() << "\n";
    }
    std::cout << std::endl;
    
    std::cout << "planewave energy cutoff [eV]: " << dat.lat.encut << std::endl;
    std::cout << "planewave gpoint cutoff [2*pi/a]: " << dat.lat.gcut << std::endl;
    std::cout << "number of plane waves: " << NPW << std::endl;
    std::cout << "Number of gpoint grids: " << NPW << std::endl;
    if (dat.debug>2){
        for(int g=0; g<NPW; g++){
            std::cout << dat.lat.G[g].transpose() << "\n";
        }
        std::cout << std::endl;
    }
    
    
    if (dat.kpoint==0){
        std::cout << "KPOINT imported from kpointfile: " << dat.kpointfile << std::endl;
    }else if (dat.kpoint==1){
        std::cout << "KPOINT generated for full BZ with order: " << dat.lat.kpointorder << std::endl;
    }else if (dat.kpoint==2){
        std::cout << "KPOINT generated for Monkhorst-Pack grids with order: " << dat.lat.kpointorder << std::endl;
    }else if (dat.kpoint==3){
        std::cout << "KPOINT generated for Cohen-Chadi grids with order: " << dat.lat.kpointorder << std::endl;
    }
    std::cout << "Number of K-vector grids: " << NKPT << std::endl;
    if (dat.debug>2){
        for(int k=0; k<NKPT; k++){
            std::cout << dat.lat.K[k].transpose() << "\n";
        }
        std::cout << std::endl;
    }
    
    std::cout << "Number of q-grids: " << NQ << std::endl;
    if (dat.debug>2){
        for(int q=0; q<NQ; q++){
            std::cout << dat.lat.Q[q].transpose() << "\n";
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

/*
    obtain total system memory size
*/
unsigned int getMEMtotal(){
    std::string token;
    std::ifstream file("/proc/meminfo");
    while(file >> token) {
        if(token == "MemTotal:") {
            unsigned int mem;
            if(file >> mem) {
                return mem;
            } else {
                return 0;
            }
        }
        // Ignore the rest of the line
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return 0; // Nothing found
}


/*
    obtain free system memory size
*/
unsigned int getMEMfree(){
    std::string token;
    std::ifstream file("/proc/meminfo");
    while(file >> token) {
        if(token == "MemFree:") {
            unsigned int mem;
            if(file >> mem) {
                return mem;
            } else {
                return 0;
            }
        }
        // Ignore the rest of the line
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return 0; // Nothing found
}

/*
    obtain program memory size
*/
unsigned int getMEMself(){
    std::string token;
    std::ifstream file("/proc/self/status");
    while(file >> token) {
        if(token == "VmSize:") {
            unsigned int mem;
            if(file >> mem) {
                return mem;
            } else {
                return 0;
            }
        }
        // Ignore the rest of the line
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return 0; // Nothing found
}
// VmPeak, VmSize, VmRSS

/*
    obtain cpu model info
*/
std::string getCPUmodel(){
    std::string token, cpumodel;
    std::ifstream file("/proc/cpuinfo");
    while (std::getline(file, cpumodel)) {
        if (cpumodel.substr(0, strlen("Processor")) == "Processor" 
                || cpumodel.substr(0,strlen("model name")) == "model name") {
            cpumodel = std::regex_replace(cpumodel, std::regex("\\model name\t: "), "");
            return cpumodel;
        }
        // Ignore the rest of the line
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return 0; // Nothing found
}


/*
    export electron band structure data
    the outputfile is saved as "${wdir}/${outputfile}.dat"
*/
void export_eband(env &dat){
    std::string outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".dat";
    std::ofstream file(outputfile, std::ios::binary);
    for(int q=0; q<NQ; q++){
        for(int n=0; n<NBAND; n++){
            file.write(reinterpret_cast<char*>(&dat.eband[q][n]), sizeof dat.eband[q][n]);
        }
    }
    file.close();
}


/*
    export susceptibility tensor matrix data
    the outputfile is save as "${wdir}/${outputfile}.dat"
*/
void export_Xij(env &dat){
    std::string outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".dat";// susceptibility outputfile
    std::ofstream file(outputfile, std::ios::binary);
    for (int i=0; i<NTSR; i++){ for (int j=0; j<NTSR; j++){
        for (int q=0; q<NQ; q++){
            for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.ReXij[q][i][j][m][n][f]), sizeof dat.ReXij[q][i][j][m][n][f]);
                }
            }}
        }
    }}
    for (int i=0; i<NTSR; i++){ for (int j=0; j<NTSR; j++){
        for (int q=0; q<NQ; q++){
            for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.ImXij[q][i][j][m][n][f]), sizeof dat.ImXij[q][i][j][m][n][f]);
                }
            }}
        }
    }}
    file.close();


    // // // try hdf5 file export
    // outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".h5";
    // const char* h5filename = outputfile.data();

    // // Create a file
    // const H5std_string FILE_NAME ("test.h5");
    // H5File h5file (FILE_NAME, H5F_ACC_TRUNC);


    // // test
    // int FSPACE_RANK = 2;
    // hsize_t fdim[2];
    // fdim[0] = 2;
    // fdim[1] = 3;
    // int d1[2][3];
    // for (int i=0; i<2; i++){
    //     for (int j=0; j<3; j++){
    //         d1[i][j] = i+j;
    //     }
    // }
    // DataSpace fspace( FSPACE_RANK, fdim );
    // const H5std_string DATASET_NAME ("dset");
    // DataSet dataset = h5file.createDataSet(
    //         DATASET_NAME, PredType::NATIVE_INT, fspace);
    // dataset.write(d1, PredType::NATIVE_INT);


    // // Create property list
    // int fillvalue = 0;
    // H5::DSetCreatPropList plist;
    // plist.setFillValue(H5::PredType::NATIVE_INT, &fillvalue);

    // Create dataspace
    

    // Create a dataset "ReXij"
    // [NQ x NTSR x NTSR x NCHI x NCHI x NFREQ]
    // int FSPACE_RANK = 6;
    // hsize_t fdim[6];
    // fdim[0] = NQ;
    // fdim[1] = NTSR;
    // fdim[2] = NTSR;
    // fdim[3] = NEPS;
    // fdim[4] = NEPS;
    // fdim[5] = NFREQ;
    // //  = {NQ, NTSR, NTSR, NEPS, NEPS, NFREQ};
    // H5::DataSpace fspace( FSPACE_RANK, fdim );
    // H5::DataSet dataset = h5file.createDataSet(
    //         "ReXij", H5::PredType::NATIVE_DOUBLE, fspace);
    // dataset.write(dat.ReXij, H5::PredType::NATIVE_DOUBLE);

    // std::vector<int> d1(50,1);

    // int RANK = 5;


    // HighFive::File test(outputfile, HighFive::File::ReadWrite | HighFive::File::Truncate);

    // 

    // test.createDataSet("/group/ReXij", d1);
    // // test.createDataSet("/group/ImXij", dat.ImXij);

    //https://docs.hdfgroup.org/archive/support/HDF5/doc/cpplus_RM/writedata_8cpp-example.html

    return;
}

/*
    export gpoint data
    the outputfile is save as "${wdir}/${outputfile}.gpt"
*/
void export_gpt(env &dat){
    std::string outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".gpt";// gpoint outputfile
    std::ofstream file(outputfile, std::ios::binary);
    for (int g=0; g<NEPS; g++){
        for (int i=0; i<3; i++){
            file.write(reinterpret_cast<char*>(&dat.lat.G[g](i)), sizeof dat.lat.G[g](i));
        }
    }
    file.close();
    return;
}

/*
    export kpoint data
    the outputfile is save as "${wdir}/${outputfile}.kpt"
*/
void export_kpt(env &dat){
    std::string outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".kpt";// kpoint outputfile
    std::ofstream file(outputfile, std::ios::binary);
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<3; i++){
            file.write(reinterpret_cast<char*>(&dat.lat.K[k](i)), sizeof dat.lat.K[k](i));
        }
        file.write(reinterpret_cast<char*>(&dat.lat.KW[k]), sizeof dat.lat.KW[k]);
    }
    file.close();
    return;
}

}/* namespace pmx */












//===================================================================================================================================================
// recycle bin



// // export susceptibility matrix as a binary file
// void writeEpsilon(env &dat){
//     std::string outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".dat";

//     std::ofstream file(outputfile, std::ios::binary);

//     for (int q=0; q<NQ; q++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     file.write(reinterpret_cast<char*>(&dat.realepsL[q][m][n][f]), sizeof dat.realepsL[q][m][n][f]);
//                 }
//             }
//         }
//     }

//     for (int q=0; q<NQ; q++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     file.write(reinterpret_cast<char*>(&dat.imagepsL[q][m][n][f]), sizeof dat.imagepsL[q][m][n][f]);
//                 }
//             }
//         }
//     }

//     for (int q=0; q<NQ; q++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     file.write(reinterpret_cast<char*>(&dat.realepsT[q][m][n][f]), sizeof dat.realepsT[q][m][n][f]);
//                 }
//             }
//         }
//     }

//     for (int q=0; q<NQ; q++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     file.write(reinterpret_cast<char*>(&dat.imagepsT[q][m][n][f]), sizeof dat.imagepsT[q][m][n][f]);
//                 }
//             }
//         }
//     }

//     file.close();
//     return;
// }


// void writeEpsilon_tensor(env &dat){
//     std::string outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".dat";
//     std::ofstream file(outputfile, std::ios::binary);
//     for (int i=0; i<3; i++){ for (int j=0; j<3; j++){
//         for (int q=0; q<NQ; q++){
//             for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     file.write(reinterpret_cast<char*>(&dat.realepsij[i][j][q][m][n][f]), sizeof dat.realepsij[i][j][q][m][n][f]);
//                 }
//             }}
//         }
//     }}
//     for (int i=0; i<3; i++){ for (int j=0; j<3; j++){
//         for (int q=0; q<NQ; q++){
//             for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     file.write(reinterpret_cast<char*>(&dat.imagepsij[i][j][q][m][n][f]), sizeof dat.imagepsij[i][j][q][m][n][f]);
//                 }
//             }}
//         }
//     }}
//     file.close();
//     return;
// }

// using namespace H5;









// // print electronic bandstructure to stdout
// void printBand(env &dat){
//     for(int q=0; q<NQ; q++){
//         for(int n=0; n<NBAND; n++){
//             std::cout << dat.band[q][n];
//             if(n<NBAND-1) std::cout << ",";
//         }
//         std::cout << std::endl;
//     }
//     return;
// }


// // print susceptibility matrix to stdout
// void printEpsilon(env &dat){

//     for (int q=0; q<NQ; q++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     std::cout << dat.realepsL[q][m][n][f];
//                     if (f<NFREQ-1) std::cout << ",";
//                 }
//                 std::cout << std::endl;
//             }
//         }
//     }

//     for (int q=0; q<NQ; q++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     std::cout << dat.imagepsL[q][m][n][f];
//                     if (f<NFREQ-1) std::cout << ",";
//                 }
//                 std::cout << std::endl;
//             }
//         }
//     }

//     for (int q=0; q<NQ; q++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     std::cout << dat.realepsT[q][m][n][f];
//                     if (f<NFREQ-1) std::cout << ",";
//                 }
//                 std::cout << std::endl;
//             }
//         }
//     }

//     for (int q=0; q<NQ; q++){
//         for (int m=0; m<NEPS; m++){
//             for (int n=0; n<=m; n++){
//                 for (int f=0; f<NFREQ; f++){
//                     std::cout << dat.imagepsT[q][m][n][f];
//                     if (f<NFREQ-1) std::cout << ",";
//                 }
//                 std::cout << std::endl;
//             }
//         }
//     }

//     return;
// }





