/*

Copyright or © or Copr. Thomas Prest.

Thomas.Prest@ENS.fr

This software is a computer program which purpose is to provide to the 
research community a proof-of-concept implementation of the identity-based
encryption scheme over NTRU lattices, described in the paper
"Efficient Identity-Based Encryption over NTRU Lattices", of
Léo Ducas, Vadim Lyubashevsky and Thomas Prest, available at
homepages.cwi.nl/~ducas/ , www.di.ens.fr/~lyubash/
and www.di.ens.fr/~prest/ .

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



#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>


#include "params.h"
#include "io.h"
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"
#include "Algebra.h"
#include "Scheme.h"


using namespace std;
using namespace NTL;




//==============================================================================
//==============================================================================
//                                  MAIN
//==============================================================================
//==============================================================================


int main()
{
    cout << "\n=======================================================================\n";
    cout << "This program is a proof-of concept for efficient IBE over lattices.\n";
    cout << "It generates a NTRU lattice of dimension 2N and associated modulus q,\n";
    cout << "and perform benches and tests, for user key extraction and encryption/decryption.";
    cout << "\n=======================================================================\n\n";

    ZZX MSK[4];
    ZZ_pX phiq, MPK;
    unsigned int i;
    float diff;
    MSK_Data * MSKD = new MSK_Data;
    MPK_Data * MPKD = new MPK_Data;
    clock_t t1, t2;
    const ZZX phi = Cyclo();

    srand(rdtsc()); // initialisation of rand

    cout << "N = " << N0 << endl;
    cout << "q = " << q0 << endl;

    ZZ_p::init(q1);
    zz_p::init(q0);

    phiq = conv<ZZ_pX>(phi);
    ZZ_pXModulus PHI(phiq);


    cout << "\n===================================================================\n KEY GENERATION";
    cout << "\n===================================================================\n";
    t1 = clock();
    for(i=0; i<1; i++)
    {
        Keygen(MPK, MSK);
    }

    CompleteMSK(MSKD, MSK);
    CompleteMPK(MPKD, MPK);

    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "It took " << diff << " seconds to generate the Master Secret Key" << endl;



    //==============================================================================
    //Key extraction bench and encryption/decryption bench
    //==============================================================================
    const unsigned int nb_extrb = 100;
    const unsigned int nb_crypb = 1000;

    cout << "\n===================================================================\n RUNNING EXTRACTION BENCH FOR ";
    cout << nb_extrb << " DIFFERENT IDENTITIES\n===================================================================\n";
    Extract_Bench(nb_extrb, MSKD);

    cout << "\n===================================================================\n RUNNING ENCRYPTION BENCH FOR ";
    cout << nb_crypb << " DIFFERENT MESSAGES\n===================================================================\n";
    Encrypt_Bench(nb_crypb, MPKD, MSKD);



    ///==============================================================================
    //Key extraction test and encryption/decryption test
    //==============================================================================
    const unsigned int nb_extrt = 100;
    const unsigned int nb_crypt = 100;

    cout << "\n===================================================================\n CHECKING EXTRACTION VALIDITY FOR ";
    cout << nb_extrt << " DIFFERENT IDENTITIES\n===================================================================\n";
    Extract_Test(nb_extrt, MSKD);

    cout << "\n===================================================================\n CHECKING ENCRYPTION VALIDITY FOR ";
    cout << nb_extrt << " DIFFERENT MESSAGES\n===================================================================\n";
    Encrypt_Test(nb_crypt, MPKD, MSKD);

    free(MSKD);
    free(MPKD);
    return 0;
}
