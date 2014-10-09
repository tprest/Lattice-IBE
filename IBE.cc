/*

Copyright or © or Copr. Thomas Prest.

Thomas.Prest@ENS.fr

This software is a computer program which purpose is to provide to the 
research community a proof-of-concept implementation of the identity-based
encryption scheme over NTRU lattices, described in the paper
"Efffficient Identity-Based Encryption over NTRU Lattices", of
Léo Ducas, Vadim Lyubashevsky and Thomas Prest, available at
http://eprint.iacr.org/2014 or http://www.di.ens.fr/~lyubash/ .

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/




#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
extern "C" {
#include <quadmath.h>
#include <complex.h>
}
#include <time.h>
#include <NTL/mat_ZZ.h>
#include <iomanip>
#include <gmp.h>
#include <cfloat>


using namespace std;
using namespace NTL;

#define Nmax 4096

#ifdef  USE_FLOAT128
    typedef __float128 long_double_t;
#else
    typedef long double long_double_t;
#endif

typedef struct
{
    unsigned int N;
    unsigned long q0;
    ZZ q;
} RingData;

typedef struct
{
    ZZX PrK[4];
    double complex PrK_fft[4][Nmax];
    long_double_t GS_Norms[Nmax];
    long_double_t sigma;
    long_double_t Lookuptable[Nmax][20];
} MSK_Data;


typedef struct
{
    ZZ_pX h;
    double complex h_FFT[Nmax];
} MPK_Data;


long_double_t B[Nmax][Nmax];
long_double_t Bstar[Nmax][Nmax];
long_double_t Bst[Nmax][Nmax];




//*************************************
/*          Seed for the RNG         */
//*************************************

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
const long_double_t Pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081L;
const long_double_t PiPrime = 0.39894228040143267793994605993438186847585863116493465766592582967065792589930183850125233390730693643030255886263518268551099195455583724299621273062L; //1/sqrt(2*Pi)
const long_double_t LDRMX = ((long_double_t)RAND_MAX );





//==============================================================================
//==============================================================================
//                      I/0 AND CONVERSION FUNCTIONS
//==============================================================================
//==============================================================================


void printdoubletab(const unsigned int N, long_double_t * Tab)
{


    cout << "[";
    for(unsigned int j=0; j<N; j++)
    {
        cout << ((long double)Tab[j]) << "    ";
    }
    cout << "]" << endl;
}

void printdoublematrix(const unsigned int N, long_double_t Tab[Nmax][Nmax])
{

    for(unsigned int j=0; j<N; j++)
    {
        printdoubletab(N,Tab[j]);
    }
    cout << endl;
}


void printinttab(const unsigned int N, long int * Tab)
{


    cout << "[";
    for(unsigned int j=0; j<N; j++)
    {
        cout << Tab[j] << " ";
    }
    cout << "]" << endl;
    cout << endl;
}


vec_ZZ Pol_to_Vec(ZZX& source, const unsigned int N)
{
    unsigned int i;
    vec_ZZ rep;
    assert(deg(source)<N);

    rep.SetLength(N);
    for(i=0; i<N; i++)
    {
        rep[i] = source[i];
    }
    return rep;
}


ZZX Vec_to_Pol(vec_ZZ& source)
{
    unsigned int i, n;
    ZZX rep;

    n = source.length();
    rep.SetLength(n);
    for(i=0; i<n; i++)
    {
        rep[i] = source[i];
    }
    return rep;
}


void printFFT(double complex const * const f_fft, const unsigned long N)
{
    unsigned int i;
    cout << "[";
    for(i=0; i<N; i++)
    {
        cout << creal(f_fft[i]) << "+ I*" << cimag(f_fft[i]) << "	";
    }
    cout << "]" << endl;
}


void printFFTnorm(double complex const * const f_fft, const unsigned long N)
{
    unsigned int i;
    cout << "[";
    for(i=0; i<N; i++)
    {
        cout << round(creal(f_fft[i])*creal(f_fft[i]) +cimag(f_fft[i])*cimag(f_fft[i])) << "	";
    }
    cout << "]" << endl;
}


void MyComplexFFT(double complex * const f_fft, double const * const f, const unsigned long N, const double complex w0)
{
    if(N==1)
    {
        f_fft[0] = f[0];
    }
    else
    {
        if(N==2)
        {
            //cout << "f[0] = " << f[0] << endl;
            //cout << "f[1] = " << f[1] << endl;
            f_fft[0] = f[0] + I*f[1];
            f_fft[1] = f[0] - I*f[1];
            //printFFT(f_fft,N);
        }
        else
        {
            assert(N%2==0);
            double f0[N/2], f1[N/2];
            double complex f0_fft[N/2], f1_fft[N/2], wk, w02;
            unsigned int k;
            for(k=0; k<(N/2); k++)
            {
                f0[k] = f[2*k];
                f1[k] = f[2*k+1];
            }

            //w0 = cexp(I*Pi/N);
            w02 = w0*w0;
            wk = w0;
            MyComplexFFT(f0_fft, f0, N/2, w02);
            //printFFT(f0_fft,N/2);
            MyComplexFFT(f1_fft, f1, N/2, w02);
            //printFFT(f1_fft,N/2);
            for(k=0; k<N; k++)
            {
                f_fft[k] = f0_fft[k%(N/2)] + wk*f1_fft[k%(N/2)];
                wk *= w02;
            }
            //printFFT(f_fft,N);
        }
    }
}


void MyReverseFFT(double complex * const f, double complex const * const f_fft, const unsigned long N, const double complex w0)
{
    //cout << "w0 = " << creal(w0) << "+ I*" << cimag(w0) << endl;
    if(N!=2)
    {
        assert(N%2==0);
        double complex f0[N/2], f1[N/2];
        double complex f0_fft[N/2], f1_fft[N/2], w02, wk;
        unsigned int k;

        //w0 = cexp(-I*Pi/N);
        w02 = w0*w0;
        wk = w0;

        for(k=0; k<N/2; k++)
        {
            f0_fft[k] = (f_fft[k] + f_fft[k+(N/2)])/2;
            f1_fft[k] = wk*(f_fft[k] - f_fft[k+(N/2)])/2;
            wk *= w02;
        }
        MyReverseFFT(f0, f0_fft, (N/2), w02);
        MyReverseFFT(f1, f1_fft, (N/2), w02);

        for(k=0; k<N/2; k++)
        {
            f[2*k] = f0[k];
            f[2*k+1] = f1[k];
        }
    }
    else
    {
        f[0] = (f_fft[0] + f_fft[1])/2;
        f[1] = (f_fft[0] - f_fft[1])/(2*I);
    }
}


void MyRealReverseFFT(double * const f, double complex const * const f_fft, const unsigned long N, const double complex w0)
{
    //cout << "w0 = " << creal(w0) << "+ I*" << cimag(w0) << endl;
    double complex fprime[Nmax];
    unsigned int i;
    MyReverseFFT(fprime, f_fft, N, w0);

    //printFFT(fprime, N);
    for(i=0; i<N; i++)
    {
        f[i] = creal(fprime[i]);
    }
}


void MyIntReverseFFT(long int * const f, double complex const * const f_fft, const unsigned long N)
{
    const double complex w0 = cexp(-I*Pi/N);
    double complex fprime[Nmax];
    unsigned int i;

    MyReverseFFT(fprime, f_fft, N, w0);
    //printFFT(fprime, N);
    for(i=0; i<N; i++)
    {
        f[i] = ((long int) round(creal(fprime[i])) );
    }
}


void MyIntFFT(double complex * f_FFT, const long int * const f, const unsigned int N)
{
    const double complex w0 = cexp(I*Pi/N);

    double f_double[Nmax];
    unsigned int i;

    for(i=0; i<N; i++)
    {
        f_double[i] = ( double ( f[i] ) );
    }

    MyComplexFFT(f_FFT, f_double, N, w0);
}


void ZZX_To_FFT(double complex * f_FFT, const ZZX f, const unsigned int N)
{
    const double complex w0 = cexp(I*Pi/N);
    double f_double[Nmax];
    unsigned int i;

    assert(deg(f)==N-1);
    assert(MaxBits(f)<900);
    for(i=0; i<N; i++)
    {
        f_double[i] = conv<double>(f[i]);
    }

    MyComplexFFT(f_FFT, f_double, N, w0);
}


void FFT_To_ZZX(ZZX& f, double complex const * const f_FFT, const unsigned int N)
{
    double f_double[Nmax];
    unsigned int i;
    const double complex w0 = cexp(-I*Pi/N);
    MyRealReverseFFT(f_double, f_FFT, N, w0);

    f.SetLength(N);
    for(i=0; i<N; i++)
    {
        f[i] = conv<ZZ>(round(f_double[i]));
    }

}


vec_ZZ RandomVector(const unsigned int N, const ZZ& q)
{
    vec_ZZ w;
    unsigned int i;
    w.SetLength(N);
    for(i=0; i<N; i++)
    {
        w[i] = conv<ZZ>(rand())%q;
    }
    return w;
}


//==============================================================================
//==============================================================================
//                      GAUSSIAN SAMPLING OVER INTEGERS
//==============================================================================
//==============================================================================


//*****************************************************************************************
/* Takes in input a random value and samples from distribution D_{\sigma_2}^+,           */  
/* Samples an element in Z+ with probability proportionnal to 2^{-x^2}                   */
//*****************************************************************************************
unsigned int Sample_0(unsigned long alea)
{
    if((alea&1UL)==0UL)
    {
        return 0;
    }
    unsigned int i;
    unsigned int k = 1;
    unsigned long mask=0;
    unsigned long aux;
    for(i=1; i<1000;)
    {
        aux = (alea&mask);
        alea = (alea>>k);     
        if(aux)
        {
            return Sample_0(alea);
        }
        else
        {
            if((alea&1UL)==0UL)
            {
                return i;
            }
        }
        i++;
        k += 2;
        mask = (mask<<2)|6UL;
    }
    cout << "ERROR" << endl;
    return 999999;
}



//*****************************************************************************************
/* Samples from distribution D_{k\sigma_2}^+, ie                                         */  
/* Samples an element in Z+ with probability proportionnal to 2^{-(x/k)^2}               */
//*****************************************************************************************
unsigned int Sample_1(const unsigned int k)
{
    unsigned int x, y, z;
    unsigned long alea = rand();

    x = Sample_0(alea);
    y = rand()%k;
    z = k*x + y;
    long double w = y*( (z<<1) - y );
    long double borne =  LDRMX / exp( w*log_2/(k*k) );
    alea = rand();
    if(alea>borne)
    {
        return Sample_1(k);
    }
    else
    {
        return z;
    }
    cout << "ERROR" << endl;
    return 999999;
}




//*****************************************************************************************
/* Samples from distribution D_{k\sigma_2}, ie                                           */  
/* Samples an element in Z with probability proportionnal to 2^{-(x/k)^2}                */
//*****************************************************************************************
signed int Sample_2(const unsigned int k)
{
    signed int signe;
    signed int x;
    unsigned long alea = rand();
    while(1)
    {
        x = Sample_1(k);
        if( (x!=0) || ((alea&1)==1) )
        {
            alea >>= 1;
            signe = 1 - 2*(alea&1);
            x *= signe;
            return x;
        }
        alea >>= 1;
    }
}


//*****************************************************************************************
/* Samples from distribution D_{sigma}, ie                                               */  
/* Samples an element in Z with probability proportionnal to e^{-x^2/2*(sigma^2)}        */
//*****************************************************************************************
signed int Sample_3(const long_double_t sigma128)
{
    signed int x;
    double alea, borne;

    const long double sigma = ((long double) sigma128);
    const unsigned long k = ( (unsigned long) ceil( (long double) sigma/sigma_1 ) );
    while(1)

    {
        x = Sample_2(k);
        alea = ((long double)rand()) / LDRMX;
        borne = exp( -x*x*( 1/(2*sigma*sigma) - 1/(2*k*k*sigma_1*sigma_1) )   );
        //cout << borne << endl;
        assert(borne<=1);
        if(alea<borne)
        {
            //cout << x << endl;
            return x;
        }
    }
}


//*****************************************************************************************
/* Samples from distribution D_{c,sigma}, ie                                             */  
/* Samples an element in Z with probability proportionnal to e^{-(c-x)^2/2*(sigma^2)}    */
//*****************************************************************************************
signed int Sample_4(long_double_t c, long_double_t sigma)
{
    double alea, borne;
    signed int x;
    unsigned int coin;

    const signed int intc = ((signed int)floor((long double)c));
    const double fracc = c-intc;
    coin = rand();
    const double denom = 1/(2*sigma*sigma);

    while(1)
    {
        x = Sample_3(sigma);
        x += (coin&1);
        if(abs(x)>8){cout << x << endl;}
        coin >>= 1;
        borne = exp(-(x-fracc)*(x-fracc)*denom)/ ( exp(-x*x*denom) + exp(-(x-1)*(x-1)*denom) );

        assert(borne<1);
        alea = ((double)rand()) / LDRMX;
        if(alea<borne)
        {
            //cout << (x+intc) << endl;
            return (x+intc);
        }
    }
}



//==============================================================================
//==============================================================================
//                      FUNCTIONS OVER POLYNOMIAL
//                      MAINLY FOR MASTER KEYGEN (eg SETUP)
//==============================================================================
//==============================================================================


/*****************************************************************************************
Computes the squared norm of a polynomial f   
*****************************************************************************************/
ZZ SquaredNorm(const ZZX& f, const unsigned int degree)
{
    unsigned int i;
    ZZ somme;
    for(i=0; i<=degree; i++)
    {
        somme += sqr(f[i]);
    }
    return somme;
}



/*****************************************************************************************
Generates a random polynomial of fixed degree
*****************************************************************************************/
ZZX RandomPoly(const unsigned int degree)
{
    unsigned int i;
    ZZX f;
    f.SetLength(degree+1);
    for(i=0; i<=degree; i++)
    {
        f[i] = rand();
    }
    return f;
}


/*****************************************************************************************
Generates a random polynomial of fixed degree and "approximately" fixed squared norm
*****************************************************************************************/
ZZX RandomPolyFixedSqNorm(const ZZ& SqNorm, const unsigned int degree)
{
    unsigned int i;
    ZZ SqNorm0, Ratio;
    ZZX f;
    f.SetLength(degree+1);

    long_double_t sigma = sqrt( ( (double) conv<double>(SqNorm)/(degree+1) ) );

    for(i=0; i<=degree; i++)
    {
        //f[i] = rand()-(RAND_MAX>>1);
        f[i] = conv<ZZ>(Sample_3(sigma));
    }
    f[degree] |= 1;
    return f;
}


ZZX Cyclo(const unsigned int N)
{
    ZZX phi;
    phi.SetLength(N+1);
    phi[0] = 1;
    phi[N] = 1;
    return phi;
}


ZZX FastMod(const ZZX& f, const unsigned int N)
{
    return (trunc(f,N) - (f>>N));
}


/****************************************************************************************
Verifies that for a parameter N, polynomials f, g are a valid semi-basis for building a NTRU lattice.
If GGCD!=1, then (f,g) isn't a valid pair
****************************************************************************************/
void ValidPair(ZZ& PGCD, ZZ& Alpha, ZZ& Beta, ZZX& rho_f, ZZX& rho_g, const ZZX& f, const ZZX& g, const unsigned int N, const ZZ& q)
{
    ZZX phi, Res_fx, Res_gx, iphi;
    ZZ Res_f, Res_g; 
    phi = Cyclo(N);

    XGCD(Res_f, rho_f, iphi, f, phi, 0);
    //cout << "Res_f = " << Res_f << endl;
    //cout << "GCD(Res_f, q) = " << GCD(Res_f, q) << endl;
    if(GCD(Res_f, q)!=1)
    {
        PGCD = 0;
    }
    else
    {    XGCD(Res_g, rho_g, iphi, g, phi, 0);
         XGCD(PGCD, Alpha, Beta, Res_f, Res_g);
    }
}


/****************************************************************************************
Computes f(1/x) mod (x^N + 1)
If f = a0 + a1*x + ... + a_{N-1}*x^{N-1}, then
Reverse(f) = a0 + a_{N-1}*x + ... + a1*x^{N-1}
****************************************************************************************/
ZZX Reverse(const ZZX& f, const unsigned int N)
{
    assert(deg(f)>=0);
    assert(deg(f)<N);

    ZZX fb;
    unsigned int i;
    fb.SetLength(N);
    fb[0] = f[0];
    fb.SetLength(N);
    for(i=N-deg(f); i<N; i++)
    {
        fb[i] = -f[N-i];
    }
    fb[0] = f[0];
    return fb;
}




/****************************************************************************************
Computes the polynomial k such that (F,G) <-- (F,G) - k*(f,g) minimizes the size of (F,G)
****************************************************************************************/
ZZX ReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, const unsigned int N, unsigned int & mb)
{
    unsigned int i;
    ZZ a;
    ZZX phi, fb, gb, num, den, iden, iphi, k;


    phi = Cyclo(N);
    fb = Reverse(f,N);
    gb = Reverse(g,N);
    num = FastMod( (fb*F+gb*G), N);
    den = FastMod( (f*fb+g*gb), N);
    mb = MaxBits(num);


    XGCD(a, iden, iphi, den, phi);
    k = FastMod(num*iden, N);

    k.SetLength(N);
    for(i=0; i<N; i++)
    {
        k[i] /= a;
    }
    //~fb; ~gb; ~phi; ~num; ~den; ~iden; ~iphi;

    return k;
}


/****************************************************************************************
Computes the polynomial k such that (F,G) <-- (F,G) - k*(f,g) minimizes the size of (F,G)
****************************************************************************************/
ZZX FastReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, const unsigned int N)
{
    //cout << "Fast Reduc" << endl;
    unsigned int i;
    ZZX k;
    double complex f_FFT[Nmax], g_FFT[Nmax], F_FFT[Nmax], G_FFT[Nmax], num_FFT[Nmax], den_FFT[Nmax], k_FFT[Nmax];

    assert(MaxBits(f)<900);
    ZZX_To_FFT(f_FFT, f, N);

    assert(MaxBits(g)<900);
    ZZX_To_FFT(g_FFT, g, N);

    assert(MaxBits(F)<900);
    ZZX_To_FFT(F_FFT, F, N);

    assert(MaxBits(G)<900);
    ZZX_To_FFT(G_FFT, G, N);

    for(i=0; i<N; i++)
    {
        num_FFT[i] = f_FFT[N-1-i]*F_FFT[i] + g_FFT[N-1-i]*G_FFT[i];
        den_FFT[i] = f_FFT[N-1-i]*f_FFT[i] + g_FFT[N-1-i]*g_FFT[i];
        k_FFT[i] = num_FFT[i]/den_FFT[i];
    }

    FFT_To_ZZX(k, k_FFT, N);
    return k;
}



/****************************************************************************************
Returns the anticircular matrix associated to integral polynomial f and integer N
****************************************************************************************/
mat_ZZ AnticircularMatrix(const ZZX& f, const unsigned int N)
{
    unsigned int i,j;
    int df;
    mat_ZZ M;
    M.SetDims(N, N);
    df = deg(f);
    if(df==-1)
    {
        return M;
    }
    unsigned dfu;
    dfu = ((unsigned) df);
    if(dfu>=N)
    {
        cout << "df = " << dfu << endl;
        cout << "f = " << f << endl;
    }
    assert(dfu<N);


    for(i=0; i<N; i++)
    {
        for(j=i; ((j<=dfu+i)&&(j<N)); j++)
        {
            M[i][j] = f[j-i];
        }
        for(j=0; (j+N)<=(dfu+i); j++)
        {
            M[i][j] = -f[j-i+N];
        }
    }
    return M;
}



/****************************************************************************************
Generates a basis from the double pair (f,g), (F,G) and N
This basis has the form :
    |f g|
M = |F G|
****************************************************************************************/
mat_ZZ BasisFromPolynomials(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, const unsigned int N)
{
    unsigned int i,j;
    mat_ZZ A,M;
    M.SetDims(2*N, 2*N);
    A = AnticircularMatrix(f, N);
    for(i=0; i<N; i++){
    for(j=0; j<N; j++){
        M[i][j] = A[i][j];
    }}

    A = AnticircularMatrix(g, N);
    for(i=0; i<N; i++){
    for(j=0; j<N; j++){
        M[i][j+N] = A[i][j];
    }}

    A = AnticircularMatrix(F, N);
    for(i=0; i<N; i++){
    for(j=0; j<N; j++){
        M[i+N][j] = A[i][j];
    }}

    A = AnticircularMatrix(G, N);
    for(i=0; i<N; i++){
    for(j=0; j<N; j++){
        M[i+N][j+N] = A[i][j];
    }}

    return M;
}



/****************************************************************************************
Computes the inverse of f (mod phi) (mod q)
****************************************************************************************/
ZZ_pX inverse(const ZZX& f, const ZZX& phi, const ZZ& q)
{
    ZZ_p::init(q);
    ZZX rho_f, iphi;
    ZZ Res_f;
    ZZ_p Res_f_1;
    XGCD(Res_f, rho_f, iphi, f, phi, 0);    
    inv(Res_f_1, conv<ZZ_p>(Res_f));
    assert(Res_f_1*conv<ZZ_p>(Res_f) == 1);

    return ( Res_f_1 * conv<ZZ_pX>(rho_f) );
}


/****************************************************************************************
Computes h = g/f (mod phi) (mod q)
****************************************************************************************/
ZZ_pX h_from_fg(const ZZX& f, const ZZX& g, const ZZX& phi, const ZZ& q)
{
    ZZ_pX f_1, g0, h0, phi0;
    f_1 = inverse(f,phi,q);
    g0 = conv<ZZ_pX>(g);
    phi0 = conv<ZZ_pX>(phi);
    h0 = (f_1*g0)%phi0;
    return h0;
}



/****************************************************************************************
Computes the Gram-Schmidt norm of the basis B generated from f,g
****************************************************************************************/
void GS_Norm(ZZX fx, ZZX gx, RingData * RD, int& flag)
{
    unsigned int i;
    const unsigned int N = RD->N;
    const unsigned q0 = RD->q0;

    double acc, acc3, Fred[Nmax], Gred[Nmax];
    double complex f[Nmax], g[Nmax], F[Nmax], G[Nmax];


    acc = 0;
    for(i=0; i<N; i++)
    {
        acc += conv<double>(fx[i]*fx[i] + gx[i]*gx[i]);
    }
    acc = sqrt(acc);
    //cout << "|b1| = " << acc << endl;

    ZZX_To_FFT(f, fx, N);
    ZZX_To_FFT(g, gx, N);

    double complex w0 = cexp(-I*Pi/N);

    for(i=0; i<N; i++)
    {
        F[i] = f[i]/(f[i]*f[N-1-i]+g[i]*g[N-1-i]);
        G[i] = g[i]/(f[i]*f[N-1-i]+g[i]*g[N-1-i]);
    }
    MyRealReverseFFT(Fred,F,N,w0);
    MyRealReverseFFT(Gred,G,N,w0);

    acc3 = 0;
    for(i=0; i<N; i++)
    {
        acc3 += Fred[i]*Fred[i] + Gred[i]*Gred[i];
    }
    acc3 = q0*sqrt(acc3);
    //cout << acc << "		"/* << acc2 << "		"*/ << acc3 << endl;
    if(acc3<acc)
    {
        //cout << "match" << endl;
        flag = 1;
    }
}



/****************************************************************************************
Generates a secret basis (f,g),(F,G) from the parameters N,q,Norme
This bases generates a NTRU lattice
****************************************************************************************/
void GenerateBasis(ZZX& f, ZZX& g, ZZX& F, ZZX& G, const ZZ& Norme, RingData * RD)
{
    int i, N;
    N = RD->N;
    ZZX rho_f, rho_g, k, aux, fb, gb, num;
    ZZ PGCD, Alpha, Beta, q;
    q = conv<ZZ>(RD->q0);

    int flag = 0;

    while( (PGCD!=1) || (flag==0) )
    {
        flag = 1;
        //while(flag==0)
        //{
            f = RandomPolyFixedSqNorm(Norme,N-1);
            g = RandomPolyFixedSqNorm(Norme,N-1);
            GS_Norm(f, g, RD, flag);
        //}
        ValidPair(PGCD, Alpha, Beta, rho_f, rho_g, f, g, N, q);
    }
    //cout << f << endl;
    F = -q*Beta*rho_g;
    G = q*Alpha*rho_f;


    f.SetLength(N);
    g.SetLength(N);


    unsigned int mb;
    k = ReductionCoefficient(f, g, F, G, N, mb);
    while(deg(k)>=0)
    {
        i++;

        F = FastMod(F - k*f, N);
        G = FastMod(G - k*g, N);

        fb = Reverse(f,N);
        gb = Reverse(g,N);

        num = FastMod( (fb*F+gb*G), N);
        mb = MaxBits(num);


        k = ReductionCoefficient(f, g, F, G, N, mb);
        k.normalize();
    }

    aux = FastMod(f*G - g*F, N);

    assert(aux[0]==q);
    assert(deg(aux)==0);
    aux.SetLength(N);
}



/****************************************************************************************
Generates from parameters N and q :
 - a public key : polynomial h
 - a private key : polynomials f,g,F,G
****************************************************************************************/
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey, RingData * RD)
{
    const unsigned long q0  = RD->q0;
    const unsigned int N = RD->N;

    ZZ SqNorm, q;
    q = conv<ZZ>(q0);
    ZZX f,g,F,G,phi;
    phi = Cyclo(N);

    SqNorm = conv<ZZ>(1.36*q0/2);

    GenerateBasis(f, g, F, G, SqNorm, RD);
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;

    for(unsigned int i=0; i<4; i++)
    {
            PrivateKey[i].SetLength(N);
    }

    PublicKey = h_from_fg(f, g, phi, q);
}

/****************************************************************************************
Computes the private basis B from private key PrivateKey and parameter N
****************************************************************************************/
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey, const unsigned int N)
{
    ZZX f,g,F,G;
    f = PrivateKey[0];
    g = PrivateKey[1];
    F = PrivateKey[2];
    G = PrivateKey[3];

    f = -f;
    F = -F;

    B = BasisFromPolynomials(g, f, G, F, N);
}



long_double_t DotProduct(const long_double_t * x1, const long_double_t * x2, const unsigned int n)
{
    unsigned int i;
    long_double_t rep = 0;
    for(i=0; i<n; i++)
    {
        rep += x1[i]*x2[i];
    }
    return rep;
}



void rotate(long_double_t * dest, long_double_t * src, unsigned int N)
{
    unsigned int i;
    for(i=0; i<N-1; i++)
    {
        dest[i+1] = src[i];
        dest[N+i+1] = src[N+i];
    }
    dest[0] = -src[N-1];
    dest[N] = -src[2*N-1];
}



void double_MGS(long_double_t Bstar[Nmax][Nmax], const long_double_t B[Nmax][Nmax], const unsigned int n)
{
    long_double_t SquareNorm[Nmax], aux[Nmax];
    unsigned int i,j,k;

    SquareNorm[0] = DotProduct(B[0], B[0], n);
    for(j=0; j<n; j++)
    {

        Bstar[0][j] = B[0][j];
    }

    for(i=1; i<n; i++)
    {
        for(k=0; k<n; k++)
        {
            Bstar[i][k] = B[i][k];
        }
        for(j=0; j<i; j++)
        {
            aux[j]= DotProduct(Bstar[i], Bstar[j],n) / SquareNorm[j];
        }
        for(k=0; k<n; k++)
        {
            for(j=0; j<i; j++)
            {
                Bstar[i][k] -= aux[j]*Bstar[j][k];
            }
        }
        SquareNorm[i] = DotProduct(Bstar[i], Bstar[i], n);
    }
}


void fast_MGS(long_double_t Bst[Nmax][Nmax], const long_double_t B[Nmax][Nmax], const unsigned int N, const unsigned long q0)
{
    long_double_t v[Nmax], v1[Nmax], C[Nmax], D[Nmax], CD[Nmax], C_k, D_k, C_ko, D_ko, aux;
    unsigned int j, k;

    cout << endl;
    //Reducing first vector (obvious)
    for(j=0; j<2*N; j++)
    {    Bst[0][j] = B[0][j];    }


    //Initialising the vector v = b_N - Proj(b_N, (b_1...b_k-2) )
    for(j=0; j<N-1; j++)
    {    v[j] = Bst[0][j+1];
         v[j+N] = Bst[0][j+1+N];    }
    v[N-1] = -Bst[0][0];
    v[2*N-1] = -Bst[0][N];

    for(j=0; j<2*N; j++)
    {    v1[j] = v[j];    }


    //Initialising recurring variables
    C_k = DotProduct(Bst[0], v, 2*N);
    C[0] = C_k;
    D_k = DotProduct(v, v, 2*N);
    D[0] = D_k;
    CD[0] = C[0]/D[0];


    //Reducing b_2 to b_N and updating v at the same time
    //cout << "DOUBLE VERSION" << endl;
    for(k=1; k<N; k++)
    {
        //printdoubletab(2*N, v);
        //b~k <-- r(b~_{k-1}) - <b~_{k-1},b_N>/<v_{k-1},b_N> r(v)
        aux = C_k/D_k;
        //cout << aux << endl;
        //printf ("C[%d]= %f\n", (k-1), C[k-1]);
        //printf ("C[%d]/D[%d]= %f\n", (k-1), (k-1), aux);
        Bst[k][0] = -Bst[k-1][N-1] + aux*v[N-1];
        Bst[k][N] = -Bst[k-1][2*N-1] + aux*v[2*N-1];
        for(j=1; j<N; j++)
        {
            Bst[k][j] = Bst[k-1][j-1] - aux*v[j-1];
            Bst[k][j+N] = Bst[k-1][j+N-1] - aux*v[j+N-1];
        }
        //SquareNorm[k] = SquareNorm[k-1] - aux*aux*sqnorm_v;


        //v <-- v - Proj(v, b~_{k-1} )
        for(j=0; j<2*N; j++)
        {
            v[j] -= aux*Bst[k-1][j];
        }
        //sqnorm_v -= aux*aux*SquareNorm[k-1];

        C_ko = C_k;
        D_ko = D_k;

        C_k = DotProduct(Bst[k], v1, 2*N);
        C[k] = C_k;
        D_k = D_ko - C_ko*C_ko/D_ko;
        D[k] = D_k;
        CD[k] = C[k]/D[k];
        //printf ("D[%d]= %f\n", k, D[k]);
    }



    //Reducing second half!
    //cout << "aux = " << (1<<10)/D[N-1] << endl;
    for(j=0; j<N; j++)
    {    Bst[N][N+j] = Bst[N-1][N-1-j]*q0/D[N-1];
         Bst[N][j] = -Bst[N-1][2*N-1-j]*q0/D[N-1];    }

    //Initialising the vector v = b_N - Proj(b_N, (b_1...b_k-2) )
    for(j=0; j<N-1; j++)
    {    v[j] = Bst[N][j+1];
         v[j+N] = Bst[N][j+1+N];    }
    v[N-1] = -Bst[N][0];
    v[2*N-1] = -Bst[N][N];

    for(j=0; j<2*N; j++)
    {    v1[j] = v[j];    }


    //Initialising recurring variables
    C_k = DotProduct(Bst[N], v1, 2*N);
    C[N] = C_k;
    D_k = DotProduct(Bst[N], Bst[N], 2*N);
    D[N] = D_k;
    CD[N] = C[N]/D[N];


    //Reducing b_2 to b_N and updating v at the same time
    for(k=N+1; k<2*N; k++)
    {
        //b~k <-- r(b~_{k-1}) - <b~_{k-1},b_N>/<v_{k-1},b_N> r(v)
        aux = C_k/D_k;
        Bst[k][0] = -Bst[k-1][N-1] + aux*v[N-1];
        Bst[k][N] = -Bst[k-1][2*N-1] + aux*v[2*N-1];
        for(j=1; j<N; j++)
        {
            Bst[k][j] = Bst[k-1][j-1] - aux*v[j-1];
            Bst[k][j+N] = Bst[k-1][j+N-1] - aux*v[j+N-1];
        }
        //SquareNorm[k] = SquareNorm[k-1] - aux*aux*sqnorm_v;


        //v <-- v - Proj(v, b~_{k-1} )
        for(j=0; j<2*N; j++)
        {
            v[j] -= aux*Bst[k-1][j];
        }
        //sqnorm_v -= aux*aux*SquareNorm[k-1];

        C_ko = C_k;
        D_ko = D_k;

        C_k = DotProduct(Bst[k], v1, 2*N);
        C[k] = C_k;
        D_k = D_ko - C_ko*C_ko/D_ko;
        D[k] = D_k;
        CD[k] = C[k]/D[k];
    }
}



void GPV(long_double_t * v, const long_double_t * const c, const long_double_t s, const unsigned int N, const MSK_Data * const MSKD)
{

    int i;
    unsigned j;
    long_double_t ci[Nmax], zi, cip, sip, aux;

    for(j=0; j<2*N;j++)
    {
        ci[j] = c[j];
    }

    for(j=0; j<2*N; j++)
    {

    }    

    for(i=2*N-1; i>=0; i--)
    {
        aux = (MSKD->GS_Norms)[i];
        cip = DotProduct(ci, Bstar[i], 2*N)/(aux*aux);
        sip = s/aux;
        zi = Sample_4(cip, sip*PiPrime);

        for(j=0; j<2*N; j++)
        {
            //aux = zi*B[i][j];
            //ci[j] -= aux;
            ci[j] -= zi*B[i][j];
        }
    }

    for(j=0; j<2*N; j++)
    {
        v[j] = c[j] - ci[j];
    }

}



//==============================================================================
//==============================================================================
//                            MAIN PROGRAMS
//==============================================================================
//==============================================================================

void IBE_Setup(ZZ_pX& MPK, ZZX * MSK, RingData * RD)
{
    Keygen(MPK, MSK, RD);
}



void CompleteMSK(MSK_Data& MSKD, ZZX * MSK, RingData * RD)
{
    unsigned int i, j, N;
    mat_ZZ B0;
    N = RD->N;

    for(i=0; i<4; i++)
    {
        MSKD.PrK[i] = MSK[i];
        ZZX_To_FFT(MSKD.PrK_fft[i], MSK[i], N);
    }

    CompletePrivateKey(B0, MSK, N);

    for(i=0; i<2*N; i++)
    {
        for(j=0; j<2*N; j++)
        {
            B[i][j] = ( (long_double_t) conv<double>(B0[i][j]) );
        }
    }

    for(i=0; i<1; i++)
    {
        fast_MGS(Bstar, B, N, RD->q0);
    }

    for(i=0; i<2*N; i++)
    {
        MSKD.GS_Norms[i] = sqrtq( DotProduct(Bstar[i], Bstar[i], 2*N) );
    }

    MSKD.sigma = 2*MSKD.GS_Norms[0];

}



void CompleteMPK(MPK_Data& MPKD, ZZ_pX MPK, RingData * RD)
{
    MPKD.h = MPK;
    ZZX_To_FFT(MPKD.h_FFT, conv<ZZX>(MPK), RD->N);
}



void IBE_Extract(ZZX SK_id[2], vec_ZZ id, const MSK_Data * const MSKD, const RingData * const RD)
{
    unsigned int i, N;
    ZZ q;

    q = conv<ZZ>(RD->q0);
    N = RD->N;
    long_double_t c[Nmax], sk[Nmax], sigma;
    ZZX f,g,t,phi,aux;

    phi = Cyclo(N);
    t = conv<ZZX>(id);
    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    sigma = MSKD->sigma;
    //cout << "sigma = " << sigma << endl;
    SK_id[0].SetLength(N);
    SK_id[1].SetLength(N);

    for(i=0;i<N;i++)
    {
        c[i] = ((long_double_t) conv<double>(id[i])) ;
        c[i+N] = 0;
    }

    GPV(sk, c, sigma, N, MSKD);



    for(i=0; i<N; i++)
    {
        sk[i] = c[i] - sk[i];
        sk[i+N] = - sk[i+N];
    }

    //cout << "sk = " << endl;
    //printdoubletab(2*N, sk);

    for(i=0; i<N; i++)
    {
        SK_id[0][i] = sk[i];
        SK_id[1][i] = sk[i+N];
    }
    
}


unsigned long IBE_Verify_Key(const ZZX SK_id[2], const vec_ZZ id, const MSK_Data * const MSKD, const RingData * const RD)
{
    unsigned int i, N;
    ZZ q;
    ZZX f,g,t,phi,aux;

    q = conv<ZZ>(RD->q0);
    N = RD->N;
    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    phi = Cyclo(N);
    
    t = conv<ZZX>(id);
    aux = ((SK_id[0] - t)*f + g*SK_id[1])%phi;

    for(i=0; i<N; i++)
    {
        aux[i] %= q;
    }

    if( IsZero(aux) != 0)
    {
        cout << "The signature (s1,s2) doesn't verify the required equality [ (s1 - t)*f + g*s2 = 0 ] !\nActually, (s1 - t)*f + g*s2 = " << aux << endl << endl;
    }
    return IsZero(aux);
}


void IBE_Encrypt(long C[2][Nmax], const long m[Nmax], const long id0[Nmax], const MPK_Data * const MPKD, const RingData * const RD)
{

    unsigned long i;
    long r[Nmax], e1[Nmax], e2[Nmax];
    double complex r_FFT[Nmax], t_FFT[Nmax], aux1_FFT[Nmax], aux2_FFT[Nmax];
    unsigned long N, q0;

    q0 = RD->q0;
    N = RD->N;

    for(i=0; i<N; i++)
    {
        e1[i] = (rand()%3) - 1;
        e2[i] = (rand()%3) - 1;
        r[i] = (rand()%3) - 1;
    }

    MyIntFFT(r_FFT, r, N);
    MyIntFFT(t_FFT, id0, N);

    for(i=0; i<N; i++)
    {
        aux1_FFT[i] = r_FFT[i]*((MPKD->h_FFT)[i]);
        aux2_FFT[i] = r_FFT[i]*t_FFT[i];
    }


    MyIntReverseFFT(C[0], aux1_FFT, N);
    MyIntReverseFFT(C[1], aux2_FFT, N);

    for(i=0; i<N; i++)
    {
        C[0][i] = (C[0][i] + e1[i]                + q0/2)%q0 - (q0/2);
        C[1][i] = (C[1][i] + e2[i] + (q0/2)*m[i] + q0/2)%q0 - (q0/2);
    }

    

}


void IBE_Decrypt(long w0[Nmax], const long C[2][Nmax], const double complex * const SKid_FFT, const RingData * const RD)
{
    unsigned int i;
    unsigned long N, q0;
    double complex c0_FFT[Nmax], aux_FFT[Nmax];

    N = RD->N;
    q0 = RD->q0;

    MyIntFFT(c0_FFT, C[0], N);

    for(i=0; i<N; i++)
    {
        aux_FFT[i] = c0_FFT[i]*SKid_FFT[i];
    }

    MyIntReverseFFT(w0, aux_FFT, N);

    for(i=0; i<N; i++)
    {
        w0[i] = C[1][i] - w0[i];
        w0[i] = ((unsigned long)(w0[i] ))%q0;
        w0[i] = (w0[i] + (q0>>2) )/(q0>>1);
        w0[i] %= 2;
    }

}



//==============================================================================
//==============================================================================
//                             BENCHES AND TESTS
//                   FOR EXTRACTION AND ENCRYPTION/DECRYPTION
//==============================================================================
//==============================================================================


void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD, RingData * RD)
{
    clock_t t1, t2;
    float diff;
    unsigned int i;
    vec_ZZ id;
    ZZX SK_id[2];

    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        id = RandomVector(RD->N, conv<ZZ>(RD->q0));

        IBE_Extract(SK_id, id, MSKD, RD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "\n\nIt took " << diff << " seconds to extract " << nb_extr << " keys." << endl;
    cout << "That's " << (diff/nb_extr)*1000 << " milliseconds per key." << endl << endl;
}


void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD, RingData * RD)
{
    clock_t t1, t2;
    double diff;
    unsigned int i,j;
    vec_ZZ id;
    ZZX SK_id[2], w;
    double complex SKid_FFT[Nmax];
    long int message[Nmax], decrypted[Nmax];
    long int identity[Nmax], Ciphertext[2][Nmax];


    id = RandomVector(RD->N, conv<ZZ>(RD->q0));
    IBE_Extract(SK_id, id, MSKD, RD);
    IBE_Verify_Key(SK_id, id, MSKD, RD);
    ZZX_To_FFT(SKid_FFT, SK_id[1], RD->N);
    for(i=0; i<(RD->N); i++)
    {
        identity[i] = conv<long int>(id[i]);
    }

    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {

        for(j=0; j < (RD->N); j++)
        {
            message[j] = (rand()%2);
        }

        IBE_Encrypt(Ciphertext, message, identity, MPKD, RD);
        IBE_Decrypt(decrypted, Ciphertext, SKid_FFT, RD);

        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }
    }

    t2 = clock();
    diff = ((double)t2 - (double)t1)/1000000.0l;
    cout << "\n\nIt took " << diff << " seconds to do " << nb_cryp << " encryptions and decryptions." << endl;
    cout << "That's " << (diff/nb_cryp)*1000 << " milliseconds per encryption+decryption." << endl;
    cout << "That's " << (diff/nb_cryp)*1000*1024/(RD->N) << " milliseconds per encryption+decryption per Kilobit." << endl << endl;
}


void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD, RingData * RD)
{
    unsigned int i, rep;
    vec_ZZ id;
    ZZX SK_id[2];

    rep = 0;

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        id = RandomVector(RD->N, conv<ZZ>(RD->q0));

        IBE_Extract(SK_id, id, MSKD, RD);
        rep += IBE_Verify_Key(SK_id, id, MSKD, RD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_extr << " extractions successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_extr << " extractions failed miserabily!" << endl << endl;    }
}


void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD, RingData * RD)
{
    unsigned int i, j, rep;
    vec_ZZ id;
    ZZX SK_id[2], m, w;
    double complex SKid_FFT[Nmax];
    long int id0[Nmax], Ciphertext[2][Nmax];
    long int message[Nmax], decrypted[Nmax];


    id = RandomVector(RD->N, conv<ZZ>(RD->q0));
    IBE_Extract(SK_id, id, MSKD, RD);
    IBE_Verify_Key(SK_id, id, MSKD, RD);
    ZZX_To_FFT(SKid_FFT, SK_id[1], RD->N);

    rep = 0;

    for(i=0; i<(RD->N); i++)
    {
        id0[i] = conv<long int>(id[i]);
    }

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {

        for(j=0; j < (RD->N); j++)
        {
            message[j] = (rand()%2);
        }

        IBE_Encrypt(Ciphertext, message, id0, MPKD, RD);
        IBE_Decrypt(decrypted, Ciphertext, SKid_FFT, RD);
        
        for(j=0; j < (RD->N); j++)
        {
            if(message[j] != decrypted[j])
            {
                cout << "ERROR : Dec(Enc(m)) != m " << endl;
                rep++;
                break;
            }
        }

        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }
    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_cryp << " encryptions+decryptions successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_cryp << " encryptions+decryptions failed miserabily!" << endl << endl;    }
}




//==============================================================================
//==============================================================================
//                                  MAIN
//==============================================================================
//==============================================================================


int main(int argc, char* argv[])
{
    cout << "\n=======================================================================\n";
    cout << "This program is a proof-of concept for efficient IBE over lattices.\n";
    cout << "It generates a NTRU lattice of dimension 2N and associated modulus q,\n";
    cout << "and perform benches and tests, for user key extraction and encryption/decryption.";
    cout << "\n=======================================================================\n\n";
    cout << "The command : ./IBE N0 q0	will run the program for N=N0, q=q0" << endl;
    cout << "The command : ./IBE N0		will run the program for N=N0, q=1048576" << endl;
    cout << "The command : ./IBE		will run the program for N=256, q=1048576" << endl << endl;


    const unsigned int logN = 8;
    unsigned long int N;
    unsigned long q0;
    ZZ q;
    ZZX phi, MSK[4];
    ZZ_pX phiq, MPK;
    unsigned int i;
    float diff;
    RingData RD;
    MSK_Data MSKD;
    MPK_Data MPKD;
    clock_t t1, t2;

    if(argc > 2)
    {
        q0 = atol(argv[2]);
    }
    else
    {
        q0 = (1<<20);
    }
    if(argc > 1)
    {
        N = atoi(argv[1]);
    }
    else
    {
        N = (1<<logN);
    }


    srand(rdtsc()); // initialisation de rand


    q = conv<ZZ>(q0);
    cout << "N = " << N << endl;
    cout << "q = " << q << endl;

    RD.N = N;
    RD.q0 = q0;
    RD.q = q;


    ZZ_p::init(q);
    zz_p::init(q0);


    phi = Cyclo(N);
    phiq = conv<ZZ_pX>(phi);
    ZZ_pXModulus PHI(phiq);


    cout << "\n===================================================================\n KEY GENERATION";
    cout << "\n===================================================================\n";
    t1 = clock();
    for(i=0; i<1; i++)
    {
        Keygen(MPK, MSK, &RD);
    }

    //printdoublematrix(2*N, B);
    //cout << endl;

    CompleteMSK(MSKD, MSK, &RD);
    CompleteMPK(MPKD, MPK, &RD);

    //cout << 1 << endl;
    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "It took " << diff << " seconds to generate the Master Secret Key" << endl;



    /********************************************************/
    //Key extraction bench and encryption/decryption bench
    /********************************************************/
    const unsigned int nb_extrb = 100;
    const unsigned int nb_crypb = 1000;

    cout << "\n===================================================================\n RUNNING EXTRACTION BENCH FOR ";
    cout << nb_extrb << " DIFFERENT IDENTITIES\n===================================================================\n";
    Extract_Bench(nb_extrb, &MSKD, &RD);

    cout << "\n===================================================================\n RUNNING ENCRYPTION BENCH FOR ";
    cout << nb_crypb << " DIFFERENT MESSAGES\n===================================================================\n";
    Encrypt_Bench(nb_crypb, &MPKD, &MSKD, &RD);



    /********************************************************/
    //Key extraction test and encryption/decryption test
    /********************************************************/
    const unsigned int nb_extrt = 100;
    const unsigned int nb_crypt = 100;

    cout << "\n===================================================================\n CHECKING EXTRACTION VALIDIDY FOR ";
    cout << nb_extrt << " DIFFERENT IDENTITIES\n===================================================================\n";
    Extract_Test(nb_extrt, &MSKD, &RD);

    cout << "\n===================================================================\n CHECKING ENCRYPTION VALIDITY FOR ";
    cout << nb_extrt << " DIFFERENT MESSAGES\n===================================================================\n";
    Encrypt_Test(nb_crypt, &MPKD, &MSKD, &RD);



    return 0;
}
