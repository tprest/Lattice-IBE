#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "io.h"
#include "params.h"

using namespace std;
using namespace NTL;

void PrintDoubleTab(const unsigned int N, RR_t const * const Tab)
{


    cout << "[";
    for(unsigned int j=0; j<N; j++)
    {
        cout << ((long double)Tab[j]) << "    ";
    }
    cout << "]" << endl;
}

void PrintDoubleMatrix(const unsigned int N, const RR_t Tab[2*N0][2*N0])
{

    for(unsigned int j=0; j<N; j++)
    {
        PrintDoubleTab(N,Tab[j]);
    }
    cout << endl;
}


void PrintIntTab(const unsigned int N, long int const * const Tab)
{


    cout << "[";
    for(unsigned int j=0; j<N; j++)
    {
        cout << Tab[j] << " ";
    }
    cout << "]" << endl;
    cout << endl;
}


void PrintFFT(CC_t const * const f_fft)
{
    unsigned int i;
    cout << "[";
    for(i=0; i<N0; i++)
    {
        cout << (f_fft[i]).real() << "+ I*" << (f_fft[i]).imag() << "	";
    }
    cout << "]" << endl;
}
