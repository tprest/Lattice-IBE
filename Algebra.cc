#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "Algebra.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"

using namespace std;
using namespace NTL;


ZZX FastMod(const ZZX& f)
{
    return (trunc(f,N0) - (f>>N0));
}

ZZX Cyclo()
{
    ZZX phi0;
    phi0.SetLength(N0+1);
    phi0[0] = 1;
    phi0[N0] = 1;
    return phi0;
}

const ZZX phi = Cyclo();


//==============================================================================
//Computes the squared norm of a polynomial f   
//==============================================================================
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


//==============================================================================
//Verifies that for a parameter N, polynomials f, g are a valid semi-basis for building a NTRU lattice.
//If PGCD!=1, then (f,g) isn't a valid pair
//==============================================================================
void ValidPair(ZZ& PGCD, ZZ& Alpha, ZZ& Beta, ZZX& rho_f, ZZX& rho_g, const ZZX& f, const ZZX& g)
{
    ZZX Res_fx, Res_gx, iphi;
    ZZ Res_f, Res_g; 

    XGCD(Res_f, rho_f, iphi, f, phi, 0);
    if(GCD(Res_f, q1)!=1)
    {
        PGCD = 0;
    }
    else
    {    XGCD(Res_g, rho_g, iphi, g, phi, 0);
         XGCD(PGCD, Alpha, Beta, Res_f, Res_g);
    }
}


//==============================================================================
//Computes f(1/x) mod (x^N + 1)
//If f = a0 + a1*x + ... + a_{N-1}*x^{N-1}, then
//Reverse(f) = a0 + a_{N-1}*x + ... + a1*x^{N-1}
//==============================================================================
ZZX Reverse(const ZZX& f)
{
    assert(deg(f)>=0);
    assert(deg(f)<N0);

    ZZX fb;
    unsigned int i;
    fb.SetLength(N0);
    fb[0] = f[0];
    fb.SetLength(N0);
    for(i=N0-deg(f); i<N0; i++)
    {
        fb[i] = -f[N0-i];
    }
    fb[0] = f[0];
    return fb;
}




//==============================================================================
//Computes the polynomial k such that (F,G) <-- (F,G) - k*(f,g) minimizes the size of (F,G)
//==============================================================================
ZZX ReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, unsigned int & mb)
{
    unsigned int i;
    ZZ a;
    ZZX fb, gb, num, den, iden, iphi, k;

    fb = Reverse(f);
    gb = Reverse(g);
    num = FastMod(fb*F + gb*G);
    den = FastMod(f*fb + g*gb);
    mb = MaxBits(num);


    XGCD(a, iden, iphi, den, phi);
    k = FastMod(num*iden);

    k.SetLength(N0);
    for(i=0; i<N0; i++)
    {
        k[i] /= a;
    }

    return k;
}


//==============================================================================
//Computes the polynomial k such that (F,G) <-- (F,G) - k*(f,g) minimizes the size of (F,G)
//==============================================================================
ZZX FastReductionCoefficient(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G)
{
    unsigned int i;
    ZZX k;
    CC_t f_FFT[N0], g_FFT[N0], F_FFT[N0], G_FFT[N0], num_FFT[N0], den_FFT[N0], k_FFT[N0];

    assert(MaxBits(f)<900);
    ZZXToFFT(f_FFT, f);

    assert(MaxBits(g)<900);
    ZZXToFFT(g_FFT, g);

    assert(MaxBits(F)<900);
    ZZXToFFT(F_FFT, F);

    assert(MaxBits(G)<900);
    ZZXToFFT(G_FFT, G);

    for(i=0; i<N0; i++)
    {
        num_FFT[i] = f_FFT[N0-1-i]*F_FFT[i] + g_FFT[N0-1-i]*G_FFT[i];
        den_FFT[i] = f_FFT[N0-1-i]*f_FFT[i] + g_FFT[N0-1-i]*g_FFT[i];
        k_FFT[i] = num_FFT[i]/den_FFT[i];
    }

    FFTToZZX(k, k_FFT);
    return k;
}



//==============================================================================
//Returns the anticircular matrix associated to integral polynomial f and integer N
//==============================================================================
mat_ZZ AnticircularMatrix(const ZZX& f)
{
    unsigned int i,j;
    int df;
    mat_ZZ M;
    M.SetDims(N0, N0);
    df = deg(f);
    if(df==-1)
    {
        return M;
    }
    unsigned dfu;
    dfu = ((unsigned) df);
    if(dfu>=N0)
    {
        cout << "df = " << dfu << endl;
        cout << "f = " << f << endl;
    }
    assert(dfu<N0);


    for(i=0; i<N0; i++)
    {
        for(j=i; ((j<=dfu+i)&&(j<N0)); j++)
        {
            M[i][j] = f[j-i];
        }
        for(j=0; (j+N0)<=(dfu+i); j++)
        {
            M[i][j] = -f[j-i+N0];
        }
    }
    return M;
}



//==============================================================================
//Generates a basis from the double pair (f,g), (F,G) and N
//This basis has the form :
//    |f g|
//M = |F G|
//==============================================================================
mat_ZZ BasisFromPolynomials(const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G)
{
    unsigned int i,j;
    mat_ZZ A,M;
    M.SetDims(2*N0, 2*N0);
    A = AnticircularMatrix(f);
    for(i=0; i<N0; i++){
    for(j=0; j<N0; j++){
        M[i][j] = A[i][j];
    }}

    A = AnticircularMatrix(g);
    for(i=0; i<N0; i++){
    for(j=0; j<N0; j++){
        M[i][j+N0] = A[i][j];
    }}

    A = AnticircularMatrix(F);
    for(i=0; i<N0; i++){
    for(j=0; j<N0; j++){
        M[i+N0][j] = A[i][j];
    }}

    A = AnticircularMatrix(G);
    for(i=0; i<N0; i++){
    for(j=0; j<N0; j++){
        M[i+N0][j+N0] = A[i][j];
    }}

    return M;
}



//==============================================================================
//Computes the Inverse of f (mod phi) (mod q)
//==============================================================================
ZZ_pX Inverse(const ZZX& f)
{
    ZZ_p::init(q1);
    ZZX rho_f, iphi;
    ZZ Res_f;
    ZZ_p Res_f_1;
    XGCD(Res_f, rho_f, iphi, f, phi, 0);    
    inv(Res_f_1, conv<ZZ_p>(Res_f));
    assert(Res_f_1*conv<ZZ_p>(Res_f) == 1);

    return ( Res_f_1 * conv<ZZ_pX>(rho_f) );
}


//==============================================================================
//Computes h = g/f (mod phi) (mod q)
//==============================================================================
ZZ_pX Quotient(const ZZX& f, const ZZX& g)
{
    ZZ_pX f_1, g0, h0, phi0;
    f_1 = Inverse(f);
    g0 = conv<ZZ_pX>(g);
    phi0 = conv<ZZ_pX>(phi);
    h0 = (f_1*g0)%phi0;
    return h0;
}



//==============================================================================
//Computes the Gram-Schmidt norm of the basis B generated from f,g
//==============================================================================
void GS_Norm(const ZZX fx, const ZZX gx, int& flag)
{
    unsigned int i;

    double acc, acc3, Fred[N0], Gred[N0];
    CC_t f[N0], g[N0], F[N0], G[N0];

    acc = 0;
    for(i=0; i<N0; i++)
    {
        acc += conv<double>(fx[i]*fx[i] + gx[i]*gx[i]);
    }
    acc = sqrt(acc);

    ZZXToFFT(f, fx);
    ZZXToFFT(g, gx);

    for(i=0; i<N0; i++)
    {
        F[i] = f[i]/(f[i]*f[N0-1-i]+g[i]*g[N0-1-i]);
        G[i] = g[i]/(f[i]*f[N0-1-i]+g[i]*g[N0-1-i]);
    }
    MyRealReverseFFT(Fred, F);
    MyRealReverseFFT(Gred, G);

    acc3 = 0;
    for(i=0; i<N0; i++)
    {
        acc3 += Fred[i]*Fred[i] + Gred[i]*Gred[i];
    }
    acc3 = q0*sqrt(acc3);
    if(acc3<acc)
    {
        flag = 1;
    }
}



//==============================================================================
//Generates a secret basis (f,g),(F,G) from the parameters N,q,Norme
//This bases generates a NTRU lattice
//==============================================================================
void GenerateBasis(ZZX& f, ZZX& g, ZZX& F, ZZX& G, const ZZ& Norme)
{
    int i;
    ZZX rho_f, rho_g, k, aux, fb, gb, num;
    ZZ PGCD, Alpha, Beta;

    int flag = 0;

    while( (PGCD!=1) || (flag==0) )
    {
        flag = 1;
        f = RandomPolyFixedSqNorm(Norme,N0-1);
        g = RandomPolyFixedSqNorm(Norme,N0-1);
        GS_Norm(f, g, flag);
        ValidPair(PGCD, Alpha, Beta, rho_f, rho_g, f, g);
    }
    F = -q1*Beta*rho_g;
    G = q1*Alpha*rho_f;

    f.SetLength(N0);
    g.SetLength(N0);

    unsigned int mb;
    k = ReductionCoefficient(f, g, F, G, mb);
    while(deg(k)>=0)
    {
        i++;

        F = FastMod(F - k*f);
        G = FastMod(G - k*g);

        fb = Reverse(f);
        gb = Reverse(g);

        num = FastMod(fb*F + gb*G);
        mb = MaxBits(num);


        k = ReductionCoefficient(f, g, F, G, mb);
        k.normalize();
    }

    aux = FastMod(f*G - g*F);

    assert(aux[0]==q1);
    assert(deg(aux)==0);
    aux.SetLength(N0);
}


RR_t DotProduct(const RR_t * x1, const RR_t * x2)
{
    unsigned int i;
    RR_t rep = 0;
    for(i=0; i<2*N0; i++)
    {
        rep += x1[i]*x2[i];
    }
    return rep;
}


void Rotate(RR_t * const dest, RR_t const * const src)
{
    unsigned int i;
    for(i=0; i<N0-1; i++)
    {
        dest[i+1] = src[i];
        dest[N0+i+1] = src[N0+i];
    }
    dest[0] = -src[N0-1];
    dest[N0] = -src[2*N0-1];
}



void ClassicMGS(RR_t Bstar[2*N0][2*N0], const RR_t B[2*N0][2*N0])
{
    RR_t SquareNorm[2*N0], aux[2*N0];
    unsigned int i,j,k;

    SquareNorm[0] = DotProduct(B[0], B[0]);
    for(j=0; j<2*N0; j++)
    {

        Bstar[0][j] = B[0][j];
    }

    for(i=1; i<2*N0; i++)
    {
        for(k=0; k<2*N0; k++)
        {
            Bstar[i][k] = B[i][k];
        }
        for(j=0; j<i; j++)
        {
            aux[j]= DotProduct(Bstar[i], Bstar[j]) / SquareNorm[j];
        }
        for(k=0; k<2*N0; k++)
        {
            for(j=0; j<i; j++)
            {
                Bstar[i][k] -= aux[j]*Bstar[j][k];
            }
        }
        SquareNorm[i] = DotProduct(Bstar[i], Bstar[i]);
    }
}


void FastMGS(RR_t Bst[2*N0][2*N0], const RR_t B[2*N0][2*N0])
{
    RR_t v[2*N0], v1[2*N0], C_k, D_k, C_ko, D_ko, aux;
    //RR_t C[2*N0], D[2*N0];
    unsigned int j, k;

    cout << endl;
    //Reducing first vector (obvious)
    for(j=0; j<2*N0; j++)
    {    Bst[0][j] = B[0][j];    }


    //Initialising the vector v = b_N - Proj(b_N, (b_1...b_k-2) )
    for(j=0; j<N0-1; j++)
    {    v[j] = Bst[0][j+1];
         v[j+N0] = Bst[0][j+1+N0];    }
    v[N0-1] = -Bst[0][0];
    v[2*N0-1] = -Bst[0][N0];

    for(j=0; j<2*N0; j++)
    {    v1[j] = v[j];    }


    //Initialising recurring variables
    C_k = DotProduct(Bst[0], v);
    D_k = DotProduct(v, v);

    //C[0] = C_k;
    //D[0] = D_k;
    //CD[0] = C[0]/D[0];


    //Reducing b_2 to b_N and updating v at the same time
    for(k=1; k<N0; k++)
    {
        //b~k <-- r(b~_{k-1}) - <b~_{k-1},b_N>/<v_{k-1},b_N> r(v)
        aux = C_k/D_k;
        Bst[k][0] = -Bst[k-1][N0-1] + aux*v[N0-1];
        Bst[k][N0] = -Bst[k-1][2*N0-1] + aux*v[2*N0-1];
        for(j=1; j<N0; j++)
        {
            Bst[k][j] = Bst[k-1][j-1] - aux*v[j-1];
            Bst[k][j+N0] = Bst[k-1][j+N0-1] - aux*v[j+N0-1];
        }

        //v <-- v - Proj(v, b~_{k-1} )
        for(j=0; j<2*N0; j++)
        {
            v[j] -= aux*Bst[k-1][j];
        }
        //sqnorm_v -= aux*aux*SquareNorm[k-1];

        C_ko = C_k;
        D_ko = D_k;

        C_k = DotProduct(Bst[k], v1);
        D_k = D_ko - C_ko*C_ko/D_ko;

        //C[k] = C_k;
        //D[k] = D_k;
        //CD[k] = C[k]/D[k];
        //printf ("C[%d]= %Lf		", k, C_k);
        //printf ("D[%d]= %Lf\n", k, D_k);
    }



    //Reducing second half!
    //cout << "aux = " << (1<<10)/D[N0-1] << endl;
    for(j=0; j<N0; j++)
    {    Bst[N0][N0+j] = Bst[N0-1][N0-1-j]*q0/D_k;
         Bst[N0][j] = -Bst[N0-1][2*N0-1-j]*q0/D_k;    }

    //Initialising the vector v = b_N - Proj(b_N, (b_1...b_k-2) )
    for(j=0; j<N0-1; j++)
    {    v[j] = Bst[N0][j+1];
         v[j+N0] = Bst[N0][j+1+N0];    }
    v[N0-1] = -Bst[N0][0];
    v[2*N0-1] = -Bst[N0][N0];

    for(j=0; j<2*N0; j++)
    {    v1[j] = v[j];    }


    //Initialising recursive variables
    C_k = DotProduct(Bst[N0], v1);
    D_k = DotProduct(Bst[N0], Bst[N0]);

    //C[N0] = C_k;
    //D[N0] = D_k;
    //CD[N0] = C[N0]/D[N0];


    //Reducing b_2 to b_N and updating v at the same time
    for(k=N0+1; k<2*N0; k++)
    {
        //b~k <-- r(b~_{k-1}) - <b~_{k-1},b_N>/<v_{k-1},b_N> r(v)
        aux = C_k/D_k;
        Bst[k][0] = -Bst[k-1][N0-1] + aux*v[N0-1];
        Bst[k][N0] = -Bst[k-1][2*N0-1] + aux*v[2*N0-1];
        for(j=1; j<N0; j++)
        {
            Bst[k][j] = Bst[k-1][j-1] - aux*v[j-1];
            Bst[k][j+N0] = Bst[k-1][j+N0-1] - aux*v[j+N0-1];
        }
        //SquareNorm[k] = SquareNorm[k-1] - aux*aux*sqnorm_v;


        //v <-- v - Proj(v, b~_{k-1} )
        for(j=0; j<2*N0; j++)
        {
            v[j] -= aux*Bst[k-1][j];
        }
        //sqnorm_v -= aux*aux*SquareNorm[k-1];

        C_ko = C_k;
        D_ko = D_k;

        C_k = DotProduct(Bst[k], v1);
        D_k = D_ko - C_ko*C_ko/D_ko;

        //C[k] = C_k;
        //D[k] = D_k;
        //CD[k] = C[k]/D[k];
    }
}
