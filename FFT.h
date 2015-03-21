#ifndef LIBE_FFT_H
#define LIBE_FFT_H

#include "params.h"


void FFTStep(CC_t * const f_fft, RR_t const * const f, const unsigned long N, const CC_t w0);
void ReverseFFTStep(CC_t * const f, CC_t const * const f_fft, const unsigned long N, const CC_t w0);
void MyRealReverseFFT(double * const f, CC_t const * const f_fft);
void MyIntReverseFFT(long int * const f, CC_t const * const f_fft);
void MyIntFFT(CC_t * f_FFT, const long int * const f);
void ZZXToFFT(CC_t * f_FFT, const ZZX f);
void FFTToZZX(ZZX& f, CC_t const * const f_FFT);

#endif
