#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "Sampling.h"
#include "params.h"

using namespace std;
using namespace NTL;



//==============================================================================
// Takes in input a random value and samples from distribution D_{\sigma_2}^+,   
// Samples an element in Z+ with probability proportionnal to 2^{-x^2}       
//==============================================================================
unsigned int Sample0(unsigned long alea)
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
            return Sample0(alea);
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



//==============================================================================
// Samples from distribution D_{k\sigma_2}^+, ie 
// Samples an element in Z+ with probability proportionnal to 2^{-(x/k)^2} 
//==============================================================================
unsigned int Sample1(const unsigned int k)
{
    unsigned int x, y, z;
    unsigned long alea = rand();

    x = Sample0(alea);
    y = rand()%k;
    z = k*x + y;
    RR_t w = y*( (z<<1) - y );
    RR_t borne =  LDRMX / exp( w*log_2/(k*k) );
    alea = rand();
    if(alea>borne)
    {
        return Sample1(k);
    }
    else
    {
        return z;
    }
    cout << "ERROR" << endl;
    return 999999;
}




//==============================================================================
// Samples from distribution D_{k\sigma_2}, ie                
// Samples an element in Z with probability proportionnal to 2^{-(x/k)^2} 
//==============================================================================
signed int Sample2(const unsigned int k)
{
    signed int signe;
    signed int x;
    unsigned long alea = rand();
    while(1)
    {
        x = Sample1(k);
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


//==============================================================================
// Samples from distribution D_{sigma}, ie                                       
// Samples an element in Z with probability proportionnal to e^{-x^2/2*(sigma^2)}
//==============================================================================
signed int Sample3(const RR_t sigma128)
{
    signed int x;
    double alea, borne;

    const RR_t sigma = sigma128;
    const unsigned long k = ( (unsigned long) ceil( (RR_t) sigma/sigma_1 ) );
    while(1)

    {
        x = Sample2(k);
        alea = ((RR_t)rand()) / LDRMX;
        borne = exp( -x*x*( 1/(2*sigma*sigma) - 1/(2*k*k*sigma_1*sigma_1) )   );
        assert(borne<=1);
        if(alea<borne)
        {
            return x;
        }
    }
}


//==============================================================================
// Samples from distribution D_{c,sigma}, ie                                              
// Samples an element in Z with probability proportionnal to e^{-(c-x)^2/2*(sigma^2)}    
//==============================================================================
signed int Sample4(RR_t c, RR_t sigma)
{
    RR_t alea, borne;
    signed int x;
    unsigned int coin;

    const signed int intc = ( (signed int) floor(c) );
    const RR_t fracc = c-intc;
    coin = rand();
    const RR_t denom = 1/(2*sigma*sigma);

    while(1)
    {
        x = Sample3(sigma);
        x += (coin&1);
        if(abs(x)>8){cout << x << endl;}
        coin >>= 1;
        borne = exp(-(x-fracc)*(x-fracc)*denom)/ ( exp(-x*x*denom) + exp(-(x-1)*(x-1)*denom) );

        assert(borne<1);
        alea = ( (RR_t)rand() ) / LDRMX;
        if(alea<borne)
        {
            return (x+intc);
        }
    }
}
