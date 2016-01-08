/*
 * =====================================================================================
 *
 *       Filename:  fftPar1D.h
 *
 *    Description:  Header of FFT file
 *
 *        Version:  1.0
 *        Created:  01/07/2016 08:45:49 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Axel Fahy (),
 *   Organization:  HES-SO hepia section ITI
 *
 * =====================================================================================
 */

#ifndef __FFT__PAR__1D
#define __FFT__PAR__1D

#include <vector>
#include <complex>

using namespace std;

void fftPar(vector<complex<double> >& data);

void ifftPar(vector<complex<double> >& data);

#endif
