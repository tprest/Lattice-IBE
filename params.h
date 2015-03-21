
#ifndef LIBE_PARAMS_H
#define LIBE_PARAMS_H

#include <math.h>
#include <complex.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>

using namespace NTL;
using std::complex;



//=====================================================================================
// These are the parameters you need to change
// N0 is the degree of the polynomial ring used. N0 must be a power of 2!
// q0 is the modulus w.r.t. whom the integers are reduced. We suggest to take q0 prime
//=====================================================================================
#define N0 512
#define q0 (1<<20)
//======================================================================================






const ZZ q1 = conv<ZZ>(q0);

//#ifdef  USE_FLOAT128
//    typedef __float128 RR_t;
//#else
    typedef long double RR_t;
    typedef complex<RR_t> CC_t;
//#endif

typedef struct
{
    ZZX PrK[4];
    CC_t PrK_fft[4][N0];
    RR_t GS_Norms[2*N0];
    RR_t sigma;
    RR_t B[2*N0][2*N0];
    RR_t Bstar[2*N0][2*N0];
} MSK_Data;


typedef struct
{
    ZZ_pX h;
    CC_t h_FFT[N0];
} MPK_Data;





//==============================================================================
//          Seed for the RNG    
//==============================================================================

#ifdef __i386
extern __inline__ uint64_t rdtsc(void) {
  uint64_t x;
  __asm__ volatile ("rdtsc" : "=A" (x));
  return x;
}
#elif defined __amd64
extern __inline__ uint64_t rdtsc(void) {
  uint64_t a, d;
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  return (d<<32) | a;
}
#endif



// Useful constants up to ~500 bits
const long double sigma_1= 0.84932180028801904272150283410288961971514109378435394286159953238339383120795466719298223538163406787061691601172910413284884326532697308797136114023L;//sqrt(1/(2*log(2)))
const long double log_2 = 0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516L;
const RR_t Pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081L;
const RR_t PiPrime = 0.39894228040143267793994605993438186847585863116493465766592582967065792589930183850125233390730693643030255886263518268551099195455583724299621273062L; //1/sqrt(2*Pi)
const RR_t LDRMX = ((RR_t)RAND_MAX );
const CC_t ii(1, 1);
const CC_t omega   = exp( ii*(Pi/N0));
const CC_t omega_1 = exp(-ii*(Pi/N0));

#endif
