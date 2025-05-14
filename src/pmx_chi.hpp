#pragma once
#ifndef _PMX_CHI_HPP_
#define _PMX_CHI_HPP_

/*
    pmx_chi.cpp | includes susceptibility calculation functinalities of picomax

    Author: Jungho Mun
    Date: June 6, 2024
*/


/*
    header files
*/
#include "pmx.hpp"
#include "pmx_chi.hpp"
#include "pmx_basis.hpp"


/*
    function declaration
*/
namespace pmx {

    extern double (*diracdelta) (double value);
    inline double diracdelta_lorentzian(double x);
    inline double diracdelta_gaussian(double x);
    inline double diracdelta_triangle(double x);
    inline double diracdelta_rect(double x);
    void kramerskronigtransform(double *ReX, double *ImX, double *w, double dw);
    void chi_LL(env &dat);
    void chi_tensor(env &dat);
    void chi_xyz(env &dat);

} /* namespace pmx */


#endif /* _PMX_CHI_HPP_ */