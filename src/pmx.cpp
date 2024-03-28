#include "pmx.hpp"

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

    // check parameters before proceeding calculations
    unsigned int memtotal = pmx::getMEMtotal();
    unsigned int memfree = pmx::getMEMfree();
    std::string cpumodel = pmx::getCPUmodel();
    int proc_count = omp_get_num_procs();

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
                << "  encut = " << dat.encut << " eV" << "\n"
                << "  gcut = " << dat.gcut << "\n"
                << "  nchi = " << NEPS << "\n"
                << "  nband = " << NBAND << "\n"
                << "  nfreq = " << NFREQ << "\n"
                << "  dfreq = " << dat.dfreq << "\n"
                << "  epsilon = " << dat.epsilon << "\n"
                << "---------------------------------" << std::endl;



    //------------------------------------------
    
    // generate g-vector grids
    time_0 = std::chrono::system_clock::now();
    pmx::generate_planewaves(dat);
    // pmx::gen_G_fcc_f(dat);
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout   << "Generated G-vectors\n"
                << "  Number of planewaves: " << NPW << "\n"
                << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;


    // construct kpoint (k-vector grids) (not required for bandstructure calculation)
    time_0 = std::chrono::system_clock::now();
    if (dat.kpoint==0){ 
        pmx::importkpointfile(dat);
        std::cout << "Imported K-vectors from kpointfile: " << dat.kpointfile << std::endl;
    }else if (dat.kpoint==1){ // generate full BZ
        pmx::genkpoint_monkhorst_pack(dat);
        std::cout << "Generated K-vectors (Full BZ)\n";
        std::cout << "  Simulation size in kpoint: " << dat.kpointorder << std::endl;
    }else if (dat.kpoint==2){ // generate Monkhorst-Pack grids
        std::cout << "Generated K-vectors (Monkhorst-Pack grids)\n";
        std::cout << "  Simulation size in kpoint: " << dat.kpointorder << std::endl;
    }else if (dat.kpoint==3){ // generate Cohen-Chadi grids
        pmx::genkpoint_cohen_chadi(dat);
        std::cout << "Generated K-vectors\n";
        std::cout << "  Using Cohen-Chadi grids with order: " << dat.kpointorder << std::endl;
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
        pmx::printBand(dat);
    }

    //------------------------------------------
    // permittivity calculation
    if (inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="eps")){
        std::cout << "Susceptibility matrix solver" << std::endl;
        time_0 = std::chrono::system_clock::now();

        // obtain reference energy level
        pmx::setRefEnergy(dat);

        // permittivity calculation
        pmx::permittivity(dat);

        // pmx::printEpsilon(dat);
        pmx::writeEpsilon(dat);

        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;
    }
    
    if (inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="epsij")){
        std::cout << "Susceptibility tensor matrix solver" << std::endl;
        time_0 = std::chrono::system_clock::now();

        // obtain reference energy level
        pmx::setRefEnergy(dat);

        // permittivity calculation
        pmx::chi_tensor_LF2(dat);

        // pmx::printEpsilon(dat);
        pmx::writeEpsilon_tensor(dat);

        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;
    }

    return 0;
}