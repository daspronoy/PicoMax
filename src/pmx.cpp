#include "pmx.hpp"

//======================================================//
//                      main                            //
//======================================================//
int main (int argc, char **argv) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;

    // input parser
    pmx::inputParser inp(argc, argv);

    // help
    if(inp.existOption("-h") || inp.existOption("-help")){
        pmx::doc();
        return 0;
    }

    // version
    if(inp.existOption("-v") || inp.existOption("-version")){
        pmx::ver();
        return 0;
    }

    // parameters
    pmx::env dat;
    pmx::inputProcess(inp, dat);
    // check parameters before proceeding calculations


    //------------------------------------------
    
    // generate g-vector grids
    std::cout << "generating G-vector grids...\n";
    pmx::generate_planewaves(dat);

    // construct kpoint (k-vector grids) (not required for bandstructure calculation)
    std::cout << "generating K-vector grids...\n";
    if (dat.kpoint==3){ // generate Cohen-Chadi grids
        pmx::genkpoint_cohen_chadi(dat);
    }else if (dat.kpoint==2){ // generate Monkhorst-Pack grids
        pmx::genkpoint_monkhorst_pack(dat);
    }else if (dat.kpoint==1){ // generate full BZ

    }else if (dat.kpoint==0){
        pmx::importkpointfile(dat);
    }

    //------------------------------------------
    // print out input parameter (for debugging purposes)
    if(dat.debug>1){
        pmx::printInput(dat);
        return 0;
    }

    //------------------------------------------
    // electronic bandstructure calculation
    if(inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="band")){
        std::cout << "starting electronic bandstructure calculation...\n";
        start = std::chrono::system_clock::now();

        pmx::empiricalpseudopotentialmethod(dat);

        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";

        // export bandstructure results
        pmx::printBand(dat);
    }

    //------------------------------------------
    // permittivity calculation
    if(inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="eps")){
        std::cout << "starting susceptibility matrix calculation...\n";
        start = std::chrono::system_clock::now();

        // obtain reference energy level
        pmx::setRefEnergy(dat);

        // permittivity calculation
        pmx::permittivity(dat);

        // pmx::printEpsilon(dat);
        pmx::writeEpsilon(dat);

        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";
    }
    
    if(inp.existOption("-switch") 
        && (inp.valueOption("-switch")=="epsij")){
        std::cout << "starting susceptibility tensor matrix calculation...\n";
        start = std::chrono::system_clock::now();

        // obtain reference energy level
        pmx::setRefEnergy(dat);

        // permittivity calculation
        pmx::chi_tensor_LF(dat);

        // pmx::printEpsilon(dat);
        pmx::writeEpsilon_tensor(dat);

        end = std::chrono::system_clock::now();
        elapsed_seconds = end-start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";
    }

    return 0;
}