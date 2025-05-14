#pragma once
#ifndef _PMX_IO_H_
#define _PMX_IO_H_

/*
    pmx_io.cpp | includes picomax initialization, input, output functionalities

    Author: Jungho Mun
    Date: June 6, 2024
*/


/*
    header files
*/
#include "pmx.hpp"


/*
    function declaration
*/
namespace pmx {

    // utility functions
    void doc();
    void ver();
    void progmsg();


    void init_crystal(inputParser inp, lattice &lat);
    void init_mater(inputParser inp, mater &mat);
    void inputProcess(inputParser inp, env &dat);

    void printInput(env &dat);

    unsigned int getMEMtotal();// obtain total memory
    unsigned int getMEMfree();// obtain free memory
    std::string getCPUmodel();// obtain CPU model info

    void export_eband(env &dat);// export electron energy band
    void export_Xij(env &dat);// export susceptibility tensor matrix
    void export_gpt(env &dat);// export gpoint
    void export_kpt(env &dat);// export kpoint

} /* namespace pmx */

#endif /* _PMX_IO_H_ */