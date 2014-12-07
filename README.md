Identity-Based Encryption over NTRU Lattices
===========

This software is a proof-of-concept implementation of an identity-based encryption scheme over NTRU lattices, described in the paper "Efficient Identity-Based Encryption over NTRU Lattices", of Léo Ducas, Vadim Lyubashevsky and Thomas Prest, available at http://eprint.iacr.org/2014/794 , http://www.di.ens.fr/~lyubash/ and http://www.di.ens.fr/~prest/ .

Warning
=======
This code is not to be considered nice, secure or efficient. Its purpose is not to be used for actual encryption, but to provide the research community a tool to verify, analyze and reproduce the statements made in our paper.

How to use?
===========

Compile the file using a C++ compiler with the C++11 norm, and linking to the GMP and NTL libraries, and then run the executable.
Example on an Unix machine with gcc:
```
$ g++ -Ofast IBE.cc -o IBE -std=c++11 -lntl -lgmp
$ ./IBE
```

If GMP, NTL and quadmath are not in a standard directory, you have to indicate where they are upon compilation.
Example:
```
$ g++ -Ofast -isystem $HOME/sw/include IBE.cc -o IBE -std=c++11 -L$HOME/sw/lib -lntl -lgmp
$ ./IBE
```

A note on efficiency
====================
Since the publication of "Efficient Identity-Based Encryption over NTRU Lattices", this code has been updated several times, so do the timings claimed. Currently, the timings for encryption/decryption are ~10 times faster than what was claimed in "Efficient Identity-Based Encryption over NTRU Lattices".
