#ifndef LIBE_IO_H
#define LIBE_IO_H

#include "params.h"


void PrintDoubleTab(const unsigned int N, RR_t const * const Tab);
void PrintDoubleMatrix(const unsigned int N, const RR_t Tab[2*N0][2*N0]);
void PrintIntTab(const unsigned int N, long int const * const Tab);
void PrintFFT(CC_t const * const f_fft);


#endif
