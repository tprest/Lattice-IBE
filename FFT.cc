#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "FFT.h"
#include "params.h"

using namespace std;
using namespace NTL;


void FFTStep(CC_t * const f_fft, RR_t const * const f, const unsigned long N, const CC_t w0)
{
    if(N==1)
    {
        f_fft[0] = f[0];
    }
    else
    {
        if(N==2)
        {
            f_fft[0] = f[0] + ii*f[1];
            f_fft[1] = f[0] - ii*f[1];
        }
        else
        {
            assert(N%2==0);
            RR_t f0[N0/2], f1[N0/2];
            CC_t f0_fft[N0/2], f1_fft[N0/2], wk, w02;
            unsigned int k;
            for(k=0; k<(N/2); k++)
            {
                f0[k] = f[2*k];
                f1[k] = f[2*k+1];
            }

            w02 = w0*w0;
            wk = w0;
            FFTStep(f0_fft, f0, N/2, w02);
            FFTStep(f1_fft, f1, N/2, w02);
            for(k=0; k<N; k++)
            {
                f_fft[k] = f0_fft[k%(N/2)] + wk*f1_fft[k%(N/2)];
                wk *= w02;
            }
        }
    }
}


void ReverseFFTStep(CC_t * const f, CC_t const * const f_fft, const unsigned long N, const CC_t w0)
{
    if(N!=2)
    {
        assert(N%2==0);
        CC_t f0[N0/2], f1[N0/2];
        CC_t f0_fft[N0/2], f1_fft[N0/2], w02, wk;
        unsigned int k;

        w02 = w0*w0;
        wk = w0;

        for(k=0; k<N/2; k++)
        {
            f0_fft[k] = (f_fft[k] + f_fft[k+(N/2)])*0.5l;
            f1_fft[k] = wk*(f_fft[k] - f_fft[k+(N/2)])*0.5l;
            wk *= w02;
        }
        ReverseFFTStep(f0, f0_fft, (N/2), w02);
        ReverseFFTStep(f1, f1_fft, (N/2), w02);

        for(k=0; k<N/2; k++)
        {
            f[2*k] = f0[k];
            f[2*k+1] = f1[k];
        }
    }
    else
    {
        f[0] = (f_fft[0] + f_fft[1])*0.5l;
        f[1] = (f_fft[0] - f_fft[1])*(-0.5l*ii);
    }
}


void MyRealReverseFFT(double * const f, CC_t const * const f_fft)
{
    CC_t fprime[N0];
    unsigned int i;
    ReverseFFTStep(fprime, f_fft, N0, omega_1);

    for(i=0; i<N0; i++)
    {
        f[i] = (fprime[i]).real();
    }
}


void MyIntReverseFFT(long int * const f, CC_t const * const f_fft)
{
    CC_t fprime[N0];
    unsigned int i;

    ReverseFFTStep(fprime, f_fft, N0, omega_1);
    for(i=0; i<N0; i++)
    {
        f[i] = ((long int) round( fprime[i].real() ) );
    }
}


void MyIntFFT(CC_t * f_FFT, const long int * const f)
{
    RR_t f_double[N0];
    unsigned int i;

    for(i=0; i<N0; i++)
    {
        f_double[i] = ( RR_t ( f[i] ) );
    }

    FFTStep(f_FFT, f_double, N0, omega);
}


void ZZXToFFT(CC_t * f_FFT, const ZZX f)
{
    RR_t f_double[N0];
    unsigned int i;

    assert(deg(f)==N0-1);
    assert(MaxBits(f)<900);
    for(i=0; i<N0; i++)
    {
        f_double[i] = ( RR_t ( conv<double>(f[i]) ) );
    }

    FFTStep(f_FFT, f_double, N0, omega);
}


void FFTToZZX(ZZX& f, CC_t const * const f_FFT)
{
    double f_double[N0];
    unsigned int i;
    MyRealReverseFFT(f_double, f_FFT);

    f.SetLength(N0);
    for(i=0; i<N0; i++)
    {
        f[i] = conv<ZZ>(round(f_double[i]));
    }

}
