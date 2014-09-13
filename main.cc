#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
# include <stdlib.h>
# include <stdint.h>
#include <stdio.h>
#include <assert.h>
# include <math.h>
# include <time.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <iomanip>
extern "C"{
#include <quadmath.h>
}
#include <complex.h>
#include <gmp.h>



using namespace std;
using namespace NTL;

#define Nmax 2048

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
    double GS_Norms[Nmax];
} MSK_Data;


typedef struct
{
    ZZ_pX h;
    double complex h_FFT[Nmax];
} MPK_Data;


double B[Nmax][Nmax];
double Bstar[Nmax][Nmax];


//************************
/*  Seed for the PRNG   */
//************************
signed int rdtsc()    
{    
    __asm__ __volatile__("rdtsc");    
} 

// Some useful constants
const double sigma_1 = sqrt(1/(2*log(2)));
const double log_2 = log(2);
const double Pi = 3.1415926535897932384626;
const double PiPrime = 0.39894228040143267793994605993438186847; //1/sqrt(2*Pi)





//********************************************
//********************************************
/*  INPUT/OUTPUT AND CONVERSION FUNCTIONS   */
//********************************************
//********************************************

void printdoubletab(const unsigned int N, double * Tab)
{
    cout << "[";
    for(unsigned int j=0; j<N; j++)
    {
        cout << Tab[j] << " ";
    }
    cout << "]" << endl;
    cout << endl;
}

void printdoublematrix(const unsigned int N, double Tab[Nmax][Nmax])
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
        cout << Tab[j] << "	";
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


//******************************************************************
/* Performs FFT of a polynomial over the complex roots of x^N+1   */
//******************************************************************
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
            f_fft[0] = f[0] + I*f[1];
            f_fft[1] = f[0] - I*f[1];
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


//****************
/* Inverse FFT  */
//****************
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


//**********************************************
/* Inverse FFT with output forced to be real  */
//**********************************************
void MyRealReverseFFT(double * const f, double complex const * const f_fft, unsigned long N, const double complex w0)
{
    //cout << "w0 = " << creal(w0) << "+ I*" << cimag(w0) << endl;
    double complex fprime[N];
    MyReverseFFT(fprime, f_fft, N, w0);
    //printFFT(fprime, N);
    for(int i=0; i<N; i++)
    {
        f[i] = creal(fprime[i]);
    }
}


//*************************************************
/* Inverse FFT with output forced to be integer  */
//*************************************************
void MyIntReverseFFT(long int * const f, double complex const * const f_fft, unsigned long N, const double complex w0)
{
    double complex fprime[N];
    MyReverseFFT(fprime, f_fft, N, w0);
    //printFFT(fprime, N);
    for(int i=0; i<N; i++)
    {
        f[i] = ((long int) round(creal(fprime[i])) );
    }
}


//**************************************
/* NTL's integer polynomial ---> FFT  */
//**************************************
void ZZX_To_FFT(double complex * f_FFT, const ZZX f, const unsigned int N)
{
    const double complex w0 = cexp(I*Pi/N);
    double f_double[N];
    unsigned int i;

    assert(deg(f)==N-1);
    assert(MaxBits(f)<900);
    for(i=0; i<N; i++)
    {
        f_double[i] = conv<double>(f[i]);
    }

    MyComplexFFT(f_FFT, f_double, N, w0);
}


//**************************************
/* FFT ---> NTL's integer polynomial  */
//**************************************
void FFT_To_ZZX(ZZX& f, double complex const * const f_FFT, const unsigned int N)
{
    double f_double[N];
    unsigned int i;
    const double complex w0 = cexp(-I*Pi/N);
    MyRealReverseFFT(f_double, f_FFT, N, w0);

    f.SetLength(N);
    for(i=0; i<N; i++)
    {
        f[i] = conv<ZZ>(round(f_double[i]));
    }

}



//*********************************
//*********************************
/*  GAUSSIAN SAMPLING FUNCTIONS  */
//*********************************
//*********************************


//*****************************************************************************************
/* Takes in input a random value and samples from distribution D_{\sigma_2}^+,           */  
/* Samples an element in Z+ with probability proportionnal to 2^{-x^2}                   */
//*****************************************************************************************
unsigned int Sample_0(unsigned long alea)
{
    unsigned long mask, aux;
    unsigned int i,k;
    if((alea&1UL)==0UL)
    {
        return 0;
    }
    k = 1;
    mask=0;
    for(i=1; i<100;)
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
}



//*****************************************************************************************
/* Samples from distribution D_{k\sigma_2}^+, ie                                         */  
/* Samples an element in Z+ with probability proportionnal to 2^{-(x/k)^2}               */
//*****************************************************************************************
unsigned int Sample_1(const unsigned int k)
{
    unsigned int b, x, y, z;
    unsigned long alea = rand();
    //const double sigma1 = k*sigma_1;
    x = Sample_0(alea);
    y = rand()%k;
    z = k*x + y;
    double w = y*(y+2*k*x);
    double borne =  ((double)RAND_MAX)/exp( w*log_2/(k*k) );
    alea = rand();
    if(alea>borne)
    {
        //cout << "redo" << endl;
        Sample_1(k);
    }
    else
    {
        //cout << "continue" << endl;
        return z;
    }
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
signed int Sample_3(const double sigma)
{
    unsigned int k;
    signed int x;
    double alea, borne;
    k = ((unsigned int)ceil(sigma/sigma_1));
    while(1)
    {
        x = Sample_2(k);
        alea = ((double)rand()) / ((double)RAND_MAX);
        borne = exp( -x*x*( 1/(2*sigma*sigma) - 1/(2*k*k*sigma_1*sigma_1) )   );
        //cout << borne << endl;
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
/* It corresponds to the algorithm SampleZ in our article                                */
//*****************************************************************************************
signed int Sample_4(const double c, const double sigma)
{
    double fracc, alea, borne, denom;
    signed int x, intc;
    unsigned int coin;

    intc = ((signed int)floor(c));
    fracc = c-floor(c);
    coin = rand();
    denom = 1/(2*sigma*sigma);

    while(1)
    {
        x = Sample_3(sigma);
        x += (coin&1);
        coin >>= 1;
        borne = exp(-(x-fracc)*(x-fracc)*denom)/( exp(-x*x*denom) + exp(-(x-1)*(x-1)*denom) );
        /*cout << "c = " << c << endl;
        cout << "fracc = " << fracc << endl;
        cout << "x = " << x << endl;
        cout << "num = " << exp(-(x-fracc)*(x-fracc)*denom) << endl;
        cout << "den = " << ( exp(-x*x*denom) + exp(-(x-1)*(x-1)*denom) ) << endl;*/
        /*if(borne>0.8)
        {cout << "borne = " << borne << endl;}*/
        assert(borne<1);
        //cout << "borne = " << borne << endl;
        alea = ((double)rand()) / ((double)RAND_MAX);
        if(alea<borne)
        {
            //cout << (x+intc) << endl;
            return (x+intc);
        }
    }
}




/*****************************************************************************************
Computes the squared norm of a polynomial f   
*****************************************************************************************/
ZZ SquaredNorm(const ZZX& f, const unsigned int degree)
{
    int i = 0;
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
    int i = 0;
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
    int i = 0;
    ZZ SqNorm0, Ratio;
    ZZX f;
    f.SetLength(degree+1);

    double sigma = sqrt(conv<double>(SqNorm)/(degree+1));

    for(i=0; i<=degree; i++)
    {
        //f[i] = rand()-(RAND_MAX>>1);
        f[i] = conv<ZZ>(Sample_3(sigma));
    }
    //cout << f << endl;
    /*Ratio = SqNorm0/SqNorm;  
    SqNorm0 = SquaredNorm(f,degree);
    Ratio = SqrRoot(SqNorm0/SqNorm);
    for(i=0; i<=degree; i++)
    { 
        f[i] /= Ratio;
    }*/
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
    ZZX fb;
    int i;
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
    int i;
    ZZ a;
    RR::SetPrecision(7000);
    RR aux;
    ZZX phi, fb, gb, num, den, iden, iphi, k, k0;


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
        aux = (conv<RR>(k[i]))/(conv<RR>(a));
        k[i] = RoundToZZ(aux);
    }

    return k;
}





/****************************************************************************************
Computes the polynomial k such that (F,G) <-- (F,G) - k*(f,g) minimizes the size of (F,G)
****************************************************************************************/
ZZX FastReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, const unsigned int N)
{
    //cout << "Fast Reduc" << endl;
    int i;
    ZZX k;
    double complex f_FFT[N], g_FFT[N], F_FFT[N], G_FFT[N], num_FFT[N], den_FFT[N], k_FFT[N];

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
    int i,j,df;
    mat_ZZ M;
    M.SetDims(N, N);
    df = deg(f);
    if(df==-1)
    {
        return M;
    }
    if(df>=N)
    {
        cout << "df = " << df << endl;
        cout << "f = " << f << endl;
    }
    assert(df<N);

    for(i=0; i<N; i++)
    {
        for(j=i; ((j<=df+i)&&(j<N)); j++)
        {
            M[i][j] = f[j-i];
        }
        for(j=0; (j+N)<=(df+i); j++)
        {
            M[i][j] = -f[j-i+N];
        }
    }
    return M;
}


/****************************************************************************************
Same but for a real polynomial
****************************************************************************************/
mat_RR AnticircularRealMatrix(const vec_RR& f, const unsigned int N)
{
    int i,j,lf;
    mat_RR M;
    M.SetDims(N, N);
    lf = f.length();
    assert(lf==N);
    for(i=0; i<N; i++)
    {
        for(j=i; j<N; j++)
        {
            M[i][j] = f[j-i];
        }
        for(j=0; j<i; j++)
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
    int i,j;
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
Does modified GSO (more stable than classic GSO) over matrix B
****************************************************************************************/
/*void MGS(mat_RR& C, const mat_ZZ& B)
{
    int i, j, n, p;
    RR r;
    vec_RR w, Cj_Squared;

    n = B.NumRows();
    p = B.NumCols();
    w.SetLength(n);
    Cj_Squared.SetLength(n);
    C.SetDims(n, p);

    for(i=0; i<n; i++)
    {
        //cout << i <<"\n";
        w = conv<vec_RR>(B[i]);
        for(j=0; j<i; j++)
        {
            r = (w*C[j])*Cj_Squared[j];
            w = w - r*C[j];
        }
        C[i] = w;
        Cj_Squared[i] = 1/(C[i]*C[i]);
    }
}*/


/****************************************************************************************
Computes the orthogonalisation B* of matrix B
****************************************************************************************/
/*void Orthogonalise_Basis(mat_RR& Bstar, const mat_ZZ& B)
{
    mat_RR mu;
    vec_RR c,aux;
    unsigned int i,j,n;
    n = B.NumRows();
    Bstar.SetDims(n, n);
    ComputeGS(B, mu, c);
    for(i=0; i<n; i++)
    {
        aux = conv<vec_RR>(B[i]);
        for(j=0; j<i; j++)
        {
            aux -= (mu[i][j])*(Bstar[j]);
        }
        Bstar[i] = aux;
    }
}*/


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
Computes the Gram-Schmidt norm of the basis B generated from f,g, using lemmas 3.4 and 3.6
flag <-- 1 if the GS norm is attained in b1, 0 otherwise
****************************************************************************************/
void GS_Norm(ZZX fx, ZZX gx, RingData * RD, int& flag)
{
    unsigned int i;
    const unsigned int N = RD->N;
    const unsigned q0 = RD->q0;

    double acc, acc2, acc3, Fred[N], Gred[N];
    double complex f[N], g[N], F[N], G[N], k[N];

    flag = 0;
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
        cout << "match" << endl;
        flag = 1;
    }
}



/****************************************************************************************
Generates a secret basis (f,g),(F,G) from the parameters N,q,Norme
This bases generates a NTRU lattice as defined in proposition 3.2
****************************************************************************************/
void GenerateBasis(ZZX& f, ZZX& g, ZZX& F, ZZX& G, const ZZ& Norme, RingData * RD)
{
    RR::SetPrecision(7000);
    int i, j, N;
    N = RD->N;
    ZZX rho_f, rho_g, k, aux, fb, gb, num;
    ZZ PGCD, Alpha, Beta, q;
    q = conv<ZZ>(RD->q0);


    ZZ SqNorm = conv<ZZ>(1.2*(RD->q0)/2);
    int flag = 0;

    while( (PGCD!=1) || (flag==0) )
    {
        flag = 0;
        while(flag==0)
        {
            f = RandomPolyFixedSqNorm(Norme,N-1);
            g = RandomPolyFixedSqNorm(Norme,N-1);
            GS_Norm(f, g, RD, flag);
        }
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
Where N,q,h,f,g,F,G verify the conditions of proposition 3.2
****************************************************************************************/
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey, RingData * RD)
{
    RR::SetPrecision(7000);
    const unsigned int q0  = RD->q0;
    const unsigned int N = RD->N;

    ZZ SqNorm, q;
    q = conv<ZZ>(q0);
    ZZX f,g,F,G,phi;
    phi = Cyclo(N);
    int flag;
    //RoundToZZ(SqNorm, (1/sqrt(120))* conv<RR>(q));

    //cout << SqNorm << endl;
    //RoundToZZ(SqNorm, (1/sqrt(1+N/12.))* conv<RR>(q));
    //cout << SqNorm << endl;

    SqNorm = conv<ZZ>(1.2*q0/2);

    GenerateBasis(f, g, F, G, SqNorm, RD);
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;
    //cout << f << endl;
    for(unsigned int i=0; i<4; i++)
    {
            PrivateKey[i].SetLength(N);
    }

    PublicKey = h_from_fg(f, g, phi, q);
}

/****************************************************************************************
Computes the private basis Lambda^perp from private key PrivateKey and parameter N
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



/*void CompletePublicKey(mat_ZZ& A, const ZZ_pX& PublicKey, const unsigned int N, const ZZ& q)
{
    ZZX One, Zero, Q, h;

    One = conv<ZZX>(1);
    Zero = conv<ZZX>(0);
    Q = conv<ZZX>(q);
    h = conv<ZZX>(PublicKey);
    One.SetLength(N);
    Zero.SetLength(N);
    h.SetLength(N);
    Q.SetLength(N);

    A = BasisFromPolynomials(One, h, Zero, Q, N);  
}*/



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



double DotProduct(double * x1, double * x2, unsigned int n)
{
    unsigned int i;
    double rep = 0;
    for(i=0; i<n; i++)
    {
        rep += x1[i]*x2[i];
    }
    return rep;
}

void printdotproduct(const unsigned int N, double Tab[Nmax][Nmax])
{

    for(unsigned int j=0; j<2; j++)
    {
        cout << sqrt(DotProduct(Tab[N*j/2], Tab[N*j/2], N)) << endl;;
    }
    cout << endl;
}



void double_MGS(double Bstar[Nmax][Nmax], double B[Nmax][Nmax], unsigned int N)
{
    double SquareNorm[2*N], aux[2*N];
    unsigned int i,j,k;

    SquareNorm[0] = DotProduct(B[0], B[0], 2*N);
    for(j=0; j<2*N; j++)
    {

        Bstar[0][j] = B[0][j];
    }

    for(i=1; i<2*N; i++)
    {
        for(k=0; k<2*N; k++)
        {
            Bstar[i][k] = B[i][k];
        }
        for(j=0; j<i; j++)
        {
            //cout << "dot = " << DotProduct(B[i], Bstar[j],2*N) << endl;
            aux[j]= DotProduct(Bstar[i], Bstar[j],2*N) / SquareNorm[j];
            //cout << "sq = " << SquareNorm[j] << endl;
            //cout << "aux = dot/sq = " << aux[j] << endl;
        }
        for(k=0; k<2*N; k++)
        {
            for(j=0; j<i; j++)
            {
                Bstar[i][k] -= aux[j]*Bstar[j][k];
            }
        }
        SquareNorm[i] = DotProduct(Bstar[i], Bstar[i], 2*N);
    }
}



void GPV(double * v, double * c, double s, unsigned int N, MSK_Data * MSKD)
{

    int i;
    unsigned j,k;
    double ci[Nmax], vi[Nmax], zi, cip, sip, aux;

    vec_RR vRR;
    vRR.SetLength(2*N);

    for(j=0; j<2*N;j++)
    {
        vi[j] = 0;
        ci[j] = c[j];
    }

    for(i=2*N-1; i>=0; i--)
    {
        //aux = DotProduct(Bstar[i], Bstar[i], 2*N);
        aux = (MSKD->GS_Norms)[i];
        cip = DotProduct(ci, Bstar[i], 2*N)/(aux*aux);
        sip = s/aux;
        //cout << "sip = " << sip << endl;
        zi = Sample_4(cip, sip*PiPrime);
        //zi = round(cip);
        for(j=0; j<2*N; j++)
        {
            aux = zi*B[i][j];
            ci[j] -= aux;
            vi[j] += aux;
        }
    }

    for(i=0; i<2*N; i++)
    {
        v[i] = vi[i];
    }

}



/*************************************************************************************
                                 ALL THE IBE PART
*************************************************************************************/

/********************************************
/* Algorithm Master_Keygen from appendix A 
*********************************************/
void IBE_Setup(ZZ_pX& MPK, ZZX * MSK, RingData * RD)
{
    Keygen(MPK, MSK, RD);
}


/********************************************
/* Algorithm Extract from appendix A 
*********************************************/
void IBE_Extract(ZZX SK_id[2], vec_ZZ id, MSK_Data * MSKD, RingData * RD)
{
    unsigned int i, j, N;
    N = RD->N;
    double c[2*N], sk[2*N], s;

    s = 5*(MSKD->GS_Norms[0]);
    SK_id[0].SetLength(N);
    SK_id[1].SetLength(N);

    for(i=0;i<N;i++)
    {
        c[i] = conv<double>(id[i]);
        c[i+N] = 0;
    }

    GPV(sk, c, s, N, MSKD);

    for(i=0; i<N; i++)
    {
        sk[i] = c[i] - sk[i];
        sk[i+N] = - sk[i+N];
    }


    for(i=0; i<N; i++)
    {
        SK_id[0][i] = sk[i];
        SK_id[1][i] = sk[i+N];
    }

}


/********************************************
/* Algorithm Encrypt from appendix A 
*********************************************/
void IBE_Encrypt(ZZX C[2], vec_ZZ id, MPK_Data * MPKD, RingData * RD, ZZX SK_id[2])
{
    ZZX t, e, r, h, aux;
    int i, N, borne;
    ZZ q;
    double complex r_FFT[Nmax], t_FFT[Nmax], aux1_FFT[Nmax], aux2_FFT[Nmax];

    q = conv<ZZ>(RD->q0);
    N = RD->N;
    borne = 4;
    h = conv<ZZX>(MPKD->h);

    t = conv<ZZX>(id);
    
    t.SetLength(N);
    e.SetLength(N);
    r.SetLength(N);


    for(i=0; i<N; i++)
    {
        e[i] = (rand()%borne) - (borne>>1);
        r[i] = (rand()%borne) - (borne>>1);
    }
    /*cout << "h = " << h << endl;
    cout << "e = " << e << endl;
    cout << "r = " << r << endl;*/


    ZZX_To_FFT(r_FFT, r, N);
    ZZX_To_FFT(t_FFT, t, N);

    for(i=0; i<N; i++)
    {
        aux1_FFT[i] = r_FFT[i]*((MPKD->h_FFT)[i]);
        aux2_FFT[i] = r_FFT[i]*t_FFT[i];
    }

    FFT_To_ZZX(C[0], aux1_FFT, N);
    FFT_To_ZZX(C[1], aux2_FFT, N);
    C[0] = C[0] + e;

}


/********************************************
/* Algorithm Decrypt from appendix A 
*********************************************/
void IBE_Decrypt(ZZX C[2], ZZX SK_id[2], double complex * SKid_FFT, RingData * RD)
{
    unsigned int i, N;
    ZZ q;
    ZZX w;
    q = conv<ZZ>(RD->q0);
    N = RD->N;
    double complex c0_FFT[Nmax], aux_FFT[Nmax];

    //cout << "c0 = " << C[0] << endl;
    ZZX_To_FFT(c0_FFT, C[0], N);
    /*FFT_To_ZZX(w, c0_FFT, N);
    cout << "c0 = " << w << endl;*/

    for(i=0; i<N; i++)
    {
        aux_FFT[i] = c0_FFT[i]*SKid_FFT[i];
    }
    FFT_To_ZZX(w, aux_FFT, N);
    w = C[1] - w;
    for(i=0; i<N; i++)
    {
        w[i] %= q;
        if( w[i] > (q>>1) )
        {
            w[i] -= q;
        }
        w[i] >>= 21;
    }
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
            //MSKD.B[i][j] = conv<double>(B0[i][j]);
            B[i][j] = conv<double>(B0[i][j]);
        }
    }

    //double_MGS(MSKD.Bstar, MSKD.B, N);
    double_MGS(Bstar, B, N);

    for(i=0; i<2*N; i++)
    {
        //MSKD.GS_Norms[i] = sqrt( DotProduct(MSKD.Bstar[i], MSKD.Bstar[i], 2*N) );
        MSKD.GS_Norms[i] = sqrt( DotProduct(Bstar[i], Bstar[i], 2*N) );
    }
}



void CompleteMPK(MPK_Data& MPKD, ZZ_pX MPK, RingData * RD)
{
    MPKD.h = MPK;
    ZZX_To_FFT(MPKD.h_FFT, conv<ZZX>(MPK), RD->N);
}



int main()
{
    const unsigned int logN = 9;
    const unsigned int N = (1<<logN);

    RR::SetPrecision(7000);
    RR::SetOutputPrecision(200);

    ZZ q, r;
    ZZX f,g,F,G,fb,gb,Fb,Gb,phi, C[2]; 
    ZZ_pX f_1, fq, phiq, aux;
    vec_ZZ w,v,sigma, id;
    mat_ZZ M, A;
    signed int mean, acc, var;
    int i,j;
    unsigned int q0;
    double complex PrK_FFT[4][Nmax], SKid_FFT[Nmax];
    ZZX MSK[4], SK_id[2];
    ZZ_pX MPK, PublicKey;
    time_t t0;
    RingData RD;
    MSK_Data MSKD;
    MPK_Data MPKD;
    clock_t t1, t2;
    float diff;


    q0 = 1<<24;
    q = conv<ZZ>(q0);
    cout << "N = " << N << endl;
    cout << "q = " << q << endl;
    cout << "sqrt(q) = " << sqrt(q0) << endl << endl;

    RD.N = N;
    RD.q0 = q0;
    RD.q = q;


    ZZ_p::init(q);
    zz_p::init(q0);


    phi = Cyclo(N);
    phiq = conv<ZZ_pX>(phi);
    ZZ_pXModulus PHI(phiq);



    srand(rdtsc()); // initialisation de rand

    mean = 0;
    var = 0;

    t1 = clock();
    for(i=0; i<1; i++)
    {
        Keygen(MPK, MSK, &RD);
    }
    CompleteMSK(MSKD, MSK, &RD);
    CompleteMPK(MPKD, MPK, &RD);
    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << endl << "It took " << diff << " seconds to generate the Master Secret Key" << endl;

    /*cout << "f=" << f << endl;
    cout << "g=" << g << endl;
    cout << "F=" << F << endl;
    cout << "G=" << G << endl;*/

    /*cout << "f=" << MSKD.PrK[0] << endl;
    cout << "g=" << MSKD.PrK[1] << endl;
    cout << "F=" << MSKD.PrK[2] << endl;
    cout << "G=" << MSKD.PrK[3] << endl;*/

    printdotproduct(2*N, B);
    cout << endl;
    printdotproduct(2*N, Bstar);

    ZZX t, h;

    vec_ZZ S;
    ZZX x;






    /********************************************************/
    //Key extraction bench
    /********************************************************/
    const unsigned int nb_extr = 100;
    t1 = clock();
    for(i=0; i<nb_extr; i++)
    {
        id = RandomVector(RD.N, conv<ZZ>(RD.q0));

        IBE_Extract(SK_id, id, &MSKD, &RD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << (i+1)/(nb_extr/10) << "0% of extractions completed" << endl;
        }
    }
    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << endl << "It took " << diff << " seconds to extract " << nb_extr << " keys." << endl;
    cout << "That's " << (diff/nb_extr)*1000 << " milliseconds per key." << endl << endl;



    ZZX_To_FFT(SKid_FFT, SK_id[1], RD.N);



    /********************************************************/
    //Encryption and decryption bench
    /********************************************************/
    const unsigned int nb_cryp = 1000;
    t1 = clock();
    for(i=0; i<nb_cryp; i++)
    {
        IBE_Encrypt(C, id, &MPKD, &RD, SK_id);
        IBE_Decrypt(C, SK_id, SKid_FFT, &RD);

        if((i+1)%(nb_cryp/10)==0)
        {
            cout << (i+1)/(nb_cryp/10) << "0% of encryptions and decryptions completed" << endl;
        }
    }
    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << endl << "It took " << diff << " seconds to do " << nb_cryp << " encryptions and decryptions." << endl;
    cout << "That's " << (diff/nb_cryp)*1000 << " milliseconds per encryption+decryption." << endl;
    cout << "That's " << (diff/nb_cryp)*1000*1024/N << " milliseconds per encryption+decryption per Kilobit." << endl << endl;



    return 0;
}
