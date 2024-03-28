/*
    io.cpp
*/
#include "pmx.hpp"

int NPW; // number of plane waves
int MAXPW;
double a; // lattice constant [angstrom]
int NEPS = 1; // epsilon cutoff
int NBAND = 8; // number of bands
int NKPT; // number of kpoint grids
int NFREQ = 1; // number of frequency
int NQ = 1; // number of q-vectors
int DELTA = 0; // delta function model; 0: lorentzian, 1: gaussian
double EPSILON = 0.1; // delta function epsilon
int ORIGINSHIFT = 0; //
bool EPM_NONLOCAL = false; //
int NVALENCE = 4; // number of valence electronic bands
int KK = 1; // kramers-kronig transform


namespace pmx {

// print documentation
void doc(){
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
void ver(){
    std::cout << "picomax-0.10\n";
    std::cout << "  Date: 2023.10.23\n";
    std::cout << "  Author: Sathwik Bharadwaj\n";
    std::cout << "  Email: sdaralag@purdue.edu\n";
    std::cout << "  https://electrodynamics.org\n";
    return;
}

// progress messaging
void progmsg(){
    
}





// process the input arguments into dat
void inputProcess(inputParser inp, env &dat){

    // debug level
    if (inp.existOption("-debug")){
        dat.debug = std::stoi(inp.valueOption("-debug"));
    }

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
    if (inp.existOption("-logfile")){
        dat.logfile = std::stoi(inp.valueOption("-logfile"));
    }

    // encut: planewave energy cutoff [eV]
    if (inp.existOption("-encut")){
        dat.encut = std::stod(inp.valueOption("-encut"));
    }

    // k-grids generation
    if (inp.existOption("-kpointfile")){
        const char *filename = inp.valueOption("-kpointfile").c_str();
        dat.kpointfile = new char[strlen(filename)+1];
        strcpy(dat.kpointfile,filename);
        dat.kpoint = 0;
    }else{
        if (inp.existOption("-kpoint")){
            dat.kpoint = std::stoi(inp.valueOption("-kpoint"));
        }
        if (inp.existOption("-kpointorder")){
            dat.kpointorder = std::stoi(inp.valueOption("-kpointorder"));
        }
    }
    
    // number of electronic bands to consider
    if (inp.existOption("-nband")){
        NBAND = std::stoi(inp.valueOption("-nband"));
    }

    // size of dielectric matrix
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
    // delta function epsilon
    if (inp.existOption("-epsilon")){
        EPSILON = std::stod(inp.valueOption("-epsilon"));
    }
    // delta function model
    if (inp.existOption("-delta")){
        DELTA = std::stod(inp.valueOption("-delta"));
    }

    // q-vector
    if (inp.existOption("-qvec")){ // direct qvec input
        std::vector<std::string> v;
        v = inp.vectorOption("-qvec");
        if (v.size()%3!=0){
            std::cout << "check the number of qvec!"; abort();
        }
        NQ = v.size()/3;
        Eigen::Vector3d qvec;
        for (int i=0; i<v.size(); i=i+3){
            qvec(0) = std::stod(v[i]);
            qvec(1) = std::stod(v[i+1]);
            qvec(2) = std::stod(v[i+2]);
            dat.Q.push_back(qvec);
        }
    }else if ((inp.existOption("-qpath") && !inp.existOption("-qnum"))
            || (!inp.existOption("-qpath") && inp.existOption("-qnum"))){
        std::cout << "specify both -qpath and -qnum!"; abort();
    }else if (inp.existOption("-qpath") && inp.existOption("-qnum")){ // generate qvec based on qpath and qnum
        std::vector<std::string> v;
        v = inp.vectorOption("-qpath");
        for (int i=0; i<v.size(); i++){
            dat.qpath.push_back(v[i]);
        }
        v = inp.vectorOption("-qnum");
        for (int i=0; i<v.size(); i++){
            dat.qnum.push_back(std::stoi(v[i]));
        }

        // construct q-vectors
        NQ = std::accumulate(dat.qnum.begin(),dat.qnum.end(),0)+1;
        const static std::unordered_map<std::string,int> bandind{
            {"G",0},{"X",1},{"L",2},{"W",3},{"U",4},{"K",5}
        };
        Eigen::Vector3d v0,v1,dv;
        for (int i=0; i<dat.qnum.size(); i++){
            switch (bandind.at(dat.qpath[i])){
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
            switch (bandind.at(dat.qpath[i+1])){
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
            int nqi = dat.qnum[i];
            for (int j=0; j<nqi; j++){
                dat.Q.push_back(v0 + double(j)/double(nqi)*dv);
                if(dat.Q.back().norm()==0) dat.Q.back() += dat.gammazero;
            }
        }
        dat.Q.push_back(v1);
        if(dat.Q.back().norm()==0) dat.Q.back() += dat.gammazero;
    }// q-vector

    // lattice constant
    if (inp.existOption("-a")){
        a = std::stod(inp.valueOption("-a"));
    }else{
        std::cout << "lattice constant (a) should be inserted!"; abort();
    }

    // parameters for empirical pseudopotential method
    if (inp.existOption("-epm")){
        std::vector<std::string> v;
        v = inp.vectorOption("-epm");
        for (int i=0; i<v.size(); i++){
            dat.epm.push_back(std::stod(v[i]));
        }
        if (v.size()==11) // nonlocal EPM
            EPM_NONLOCAL = true;
        else if (v.size()!=6){
            std::cout << "check the number of epm parameters!\n";
            abort();
        }
    }

    // whether to shift lattice origin
    if (inp.existOption("-originshift")){
        ORIGINSHIFT = std::stod(inp.valueOption("-originshift"));
    }
    
    // preprocessing...

    dat.gcut = sqrt(dat.encut*pow(a/2.0/pi,2)/KCON);

    return;
}



// print input to stdout (for debugging purposes)
void printInput(env &dat){
    
    std::cout << "working directory: " << dat.wdir << std::endl;
    std::cout << "outputfile: " << dat.outputfile << std::endl;
    std::cout << "logfile: " << dat.logfile << std::endl;

    
    std::cout << "number of electronic bands: " << NBAND << std::endl;
    std::cout << "size of dielectric matrix: " << NEPS << std::endl;
    std::cout << "Material properties" << std::endl;
    std::cout << "lattice constant [angstrom]: " << a << std::endl;
    std::cout << "Empirical pseudopotential method parameters:" << std::endl;
    std::cout << "vs3\tvs8\tvs11\tva3\tva4\tva11\ta0\tb0\tA2\trl0\trl2" << std::endl;
    for(int i=0; i<dat.epm.size(); i++){
        std::cout << dat.epm[i] << "\t";
    }
    std::cout << std::endl;
    // std::cout << "Lattice vectors" << std::endl;
    // std::cout << "b1\tb2\tb3" << std::endl;
    // for(int i=0; i<3; i++){
    //     std::cout << dat.b[0](i) << "\t" << dat.b[1](i) << "\t" << dat.b[2](i) << std::endl;
    // }
    
    std::cout << "planewave energy cutoff [eV]: " << dat.encut << std::endl;
    std::cout << "planewave g-vector cutoff [2*pi/a]: " << dat.gcut << std::endl;
    std::cout << "number of plane waves: " << NPW << std::endl;
    std::cout << "Number of G-vector grids: " << NPW << std::endl;
    for(int i=0; i<3; i++){
        for(int g=0; g<NPW; g++){
            std::cout << dat.G[g](i) << "\t";
        }
        std::cout << std::endl;
    }
    
    
    if (dat.kpoint==0){
        std::cout << "KPOINT imported from kpointfile: " << dat.kpointfile << std::endl;
    }else if (dat.kpoint==1){
        std::cout << "KPOINT generated for full BZ with order: " << dat.kpointorder << std::endl;
    }else if (dat.kpoint==2){
        std::cout << "KPOINT generated for Monkhorst-Pack grids with order: " << dat.kpointorder << std::endl;
    }else if (dat.kpoint==3){
        std::cout << "KPOINT generated for Cohen-Chadi grids with order: " << dat.kpointorder << std::endl;
    }
    std::cout << "Number of K-vector grids: " << NKPT << std::endl;
    for(int i=0; i<3; i++){
        for(int k=0; k<NKPT; k++){
            std::cout << dat.K[k][i] << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << "Number of q-grids: " << NQ << std::endl;
    for(int i=0; i<3; i++){
        for(int q=0; q<NQ; q++){
            std::cout << dat.Q[q](i) << "\t";
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


// get total system memory size
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

// get free system memory size
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


// get cpu model name
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




// print electronic bandstructure to stdout
void printBand(env &dat){
    for(int q=0; q<NQ; q++){
        for(int n=0; n<NBAND; n++){
            std::cout << dat.band[q][n];
            if(n<NBAND-1) std::cout << ",";
        }
        std::cout << std::endl;
    }
    return;
}


// print susceptibility matrix to stdout
void printEpsilon(env &dat){

    for (int q=0; q<NQ; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    std::cout << dat.realepsL[q][m][n][f];
                    if (f<NFREQ-1) std::cout << ",";
                }
                std::cout << std::endl;
            }
        }
    }

    for (int q=0; q<NQ; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    std::cout << dat.imagepsL[q][m][n][f];
                    if (f<NFREQ-1) std::cout << ",";
                }
                std::cout << std::endl;
            }
        }
    }

    for (int q=0; q<NQ; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    std::cout << dat.realepsT[q][m][n][f];
                    if (f<NFREQ-1) std::cout << ",";
                }
                std::cout << std::endl;
            }
        }
    }

    for (int q=0; q<NQ; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    std::cout << dat.imagepsT[q][m][n][f];
                    if (f<NFREQ-1) std::cout << ",";
                }
                std::cout << std::endl;
            }
        }
    }

    return;
}

// export susceptibility matrix as a binary file
void writeEpsilon(env &dat){
    std::string outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".dat";

    std::ofstream file(outputfile, std::ios::binary);

    for (int q=0; q<NQ; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.realepsL[q][m][n][f]), sizeof dat.realepsL[q][m][n][f]);
                }
            }
        }
    }

    for (int q=0; q<NQ; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.imagepsL[q][m][n][f]), sizeof dat.imagepsL[q][m][n][f]);
                }
            }
        }
    }

    for (int q=0; q<NQ; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.realepsT[q][m][n][f]), sizeof dat.realepsT[q][m][n][f]);
                }
            }
        }
    }

    for (int q=0; q<NQ; q++){
        for (int m=0; m<NEPS; m++){
            for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.imagepsT[q][m][n][f]), sizeof dat.imagepsT[q][m][n][f]);
                }
            }
        }
    }

    file.close();
    return;
}


void writeEpsilon_tensor(env &dat){
    std::string outputfile = std::string(dat.wdir)+"/"+dat.outputfile+".dat";
    std::ofstream file(outputfile, std::ios::binary);
    for (int i=0; i<3; i++){ for (int j=0; j<3; j++){
        for (int q=0; q<NQ; q++){
            for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.realepsij[i][j][q][m][n][f]), sizeof dat.realepsij[i][j][q][m][n][f]);
                }
            }}
        }
    }}
    for (int i=0; i<3; i++){ for (int j=0; j<3; j++){
        for (int q=0; q<NQ; q++){
            for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
                for (int f=0; f<NFREQ; f++){
                    file.write(reinterpret_cast<char*>(&dat.imagepsij[i][j][q][m][n][f]), sizeof dat.imagepsij[i][j][q][m][n][f]);
                }
            }}
        }
    }}
    file.close();
    return;
}

}