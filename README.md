Identity-Based Encryption over NTRU Lattices
===========

This software is a proof-of-concept implementation of an identity-based encryption scheme over NTRU lattices, described in the paper "Efficient Identity-Based Encryption over NTRU Lattices", of Léo Ducas, Vadim Lyubashevsky and Thomas Prest, available at http://eprint.iacr.org/2014 or http://www.di.ens.fr/~lyubash/ .

How to use?
===========

Compile the file using a C++ compiler, and linking to the GMP, NTL and quadmath libraries, and then run the executable.
Example on an Unix machine with gcc:

g++ -Ofast IBE.cc -o IBE -lntl -lgmp -lquadmath
./IBE


If GMP, NTL and quadmath are not in a standard directory, you have to indicate where they are upon compilation.
Example:

g++ -Ofast -I$HOME/sw/include IBE.cc -o IBE -L$HOME/sw/lib -lntl -lgmp -lquadmath

./IBE
