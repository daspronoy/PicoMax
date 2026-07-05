/*
    pmx.cpp | main file for picomax

    Author: Jungho Mun
    Date: June 6, 2024
*/


/*
    header files
*/
#include "pmx.hpp"
#include "pmx_io.hpp"
#include "pmx_basis.hpp"
#include "pmx_epm.hpp"
#include "pmx_chi.hpp"


/*
    global variables
*/
int NPW; // number of plane waves
int MAXPW;
double a; // lattice constant [angstrom]
int NEPS = 1; // epsilon cutoff
int NBAND = 8; // number of bands
int NKPT; // number of kpoint grids
int NFREQ = 1; // number of frequency
int NQ = 0; // number of q-vectors
int DELTA = 0; // delta function model; 0: lorentzian, 1: gaussian
double EPSILON = 0.1; // delta function epsilon
int ORIGINSHIFT = 0; //
bool EPM_NONLOCAL = false; //
int KK = 1; // kramers-kronig transform
int NTSR = 3; // tensor dimension (1 or 3)
int NATOM = 0; // number of atoms in the primitive unit cell
std::complex<double> im = std::complex<double>(0.0,1.0);



//======================================================//
//                      main                            //
//======================================================//
int main (int argc, char **argv) {
    std::chrono::time_point<std::chrono::system_clock> time_init, time_0, time_1;
    std::chrono::duration<double> elapsed_time;
    time_init = std::chrono::system_clock::now();
    std::time_t datetime_init = std::chrono::system_clock::to_time_t(time_init);

    // input parser
    pmx::inputParser inp(argc, argv);

    // help
    if (inp.existOption("-h") || inp.existOption("-help")){
        pmx::doc();
        return 0;
    }

    // version
    if (inp.existOption("-v") || inp.existOption("-version")){
        pmx::ver();
        return 0;
    }

    // parameters
    pmx::env dat;
    pmx::inputProcess(inp, dat);
    pmx::init_crystal(inp, dat.lat);
    pmx::init_mater(inp, dat.mat);

    // check parameters before proceeding calculations
    unsigned int memtotal = pmx::getMEMtotal();
    unsigned int memfree = pmx::getMEMfree();
    std::string cpumodel = pmx::getCPUmodel();
    int proc_count = omp_get_num_procs();

    if (dat.logfile){
        std::string logfile = std::string(dat.wdir)+"/"+dat.outputfile+".log";
        std::cout << "Output is redirected to logfile, " << logfile << std::endl;
        if (!std::freopen(logfile.c_str(),"w",stdout)){
            std::cout << "Logfile may have a problem?" << std::endl; abort();
        }
    }

    std::cout   << "*********************************\n"
                << "***PicoMax 0.1 progress output***\n"
                << "*********************************\n"
                << std::ctime(&datetime_init)
                << "PicoMax v" << dat.version << "\n"
                << "  Running on " << cpumodel << "\n"
                << "  Using " << proc_count << " cores\n"
                << "  Total memory: " << (memtotal/1000) << " MB\n"
                << "  Avail memory: " << (memfree/1000) << " MB\n"
                << "  Working dir: " << dat.wdir << "\n"
                << "  Output file: " << dat.outputfile << ".dat\n"
                << "---------------------------------\n"
                << "Parameters:\n"
                << "  encut = " << dat.lat.encut << " eV" << "\n"
                << "  gcut = " << dat.lat.gcut << "\n"
                << "  nchi = " << NEPS << "\n"
                << "  nband = " << NBAND << "\n"
                << "  nfreq = " << NFREQ << "\n"
                << "  dfreq = " << dat.dfreq << "\n"
                << "  epsilon = " << EPSILON << "\n"
                << "---------------------------------" << std::endl;



    //------------------------------------------
    // q-vectors
    std::cout   << "Generated Q-vector\n"
                << "  Path: ";
    for (std::string i: dat.lat.qpath)
        std::cout << i << ",";
    std::cout   << "\n"
                << "  Qnum: ";
    for (int i: dat.lat.qnum)
        std::cout << i << ",";
    std::cout   << "\n"
                << "  Total number: " << NQ << std::endl;


    // generate g-vector grids
    time_0 = std::chrono::system_clock::now();
    if (dat.lat.gsym){
        pmx::generate_gpt(dat.lat);
    }else{
        pmx::generate_gpt_rect(dat.lat);
    }
    

    // pmx::gen_G_fcc_f(dat);
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout   << "Generated G-vectors\n"
                << "  Number of planewaves: " << NPW << "\n"
                << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;


    // construct kpoint (k-vector grids) (not required for bandstructure calculation)
    time_0 = std::chrono::system_clock::now();
    if (dat.kpoint==0){ 
        pmx::import_kpt(dat);
        std::cout << "Imported K-vectors from kpointfile: " << dat.kpointfile << std::endl;
    }else if (dat.kpoint==1){ // generate full BZ
        pmx::generate_kpt_bz(dat.lat);
        std::cout << "Generated K-vectors (Full BZ)\n";
        std::cout << "  Simulation size in kpoint: " << dat.lat.kpointorder << std::endl;
    }else if (dat.kpoint==2){ // generate Monkhorst-Pack grids
        pmx::generate_kpt_monkhorst_pack(dat);
        std::cout << "Generated K-vectors (Monkhorst-Pack grids)\n";
        std::cout << "  Simulation size in kpoint: " << dat.lat.kpointorder << std::endl;
    }else if (dat.kpoint==3){ // generate Cohen-Chadi grids
        pmx::generate_kpt_cohen_chadi(dat.lat);
        std::cout << "Generated K-vectors\n";
        std::cout << "  Using Cohen-Chadi grids with order: " << dat.lat.kpointorder << std::endl;
    }
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout   << "  Number of K-vectors: " << NKPT << "\n"
                << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

    //------------------------------------------
    // print out input parameter (for debugging purposes)
    if (dat.debug>1){
        pmx::printInput(dat);
        return 0;
    }

    //------------------------------------------
    // electronic bandstructure calculation
    if (inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="band")){
        std::cout << "Electronic bandstructure solver" << std::endl;
        time_0 = std::chrono::system_clock::now();
        pmx::empiricalpseudopotentialmethod(dat);
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

        // export bandstructure results
        pmx::export_eband(dat);
    }

    //------------------------------------------
    if (inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="epsij")){
        std::cout << "Susceptibility tensor matrix solver" << std::endl;
        time_0 = std::chrono::system_clock::now();

        // obtain reference energy level
        pmx::setRefEnergy(dat);

        // permittivity calculation
        pmx::chi_tensor(dat);

        // export results
        pmx::export_Xij(dat);
        pmx::export_gpt(dat);

        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;
    }

    if (inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="epsLL")){
        std::cout << "Susceptibility tensor matrix solver" << std::endl;
        time_0 = std::chrono::system_clock::now();

        // obtain reference energy level
        pmx::setRefEnergy(dat);

        // permittivity calculation
        pmx::chi_LL(dat);

        // export results
        pmx::export_Xij(dat);
        pmx::export_gpt(dat);

        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;
    }

    if (inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="epsxyz")){
        std::cout << "Susceptibility tensor matrix solver (xyz-basis)" << std::endl;
        time_0 = std::chrono::system_clock::now();

        // obtain reference energy level
        pmx::setRefEnergy(dat);

        // permittivity calculation
        pmx::chi_xyz(dat);

        // export results
        pmx::export_Xij(dat);
        pmx::export_gpt(dat);

        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;
    }


    return 0;
}