Identity-Based Encryption over NTRU Lattices
===========

This software is a proof-of-concept implementation of an identity-based encryption scheme over NTRU lattices, described in the paper "Efficient Identity-Based Encryption over NTRU Lattices", of Léo Ducas, Vadim Lyubashevsky and Thomas Prest, available at http://homepages.cwi.nl/~ducas/ , http://www.di.ens.fr/~lyubash/ and http://www.di.ens.fr/~prest/ .

Warning
=======
This code is not to be considered secure, efficient or fully portable. Its purpose is not to be used for actual encryption, but to provide the research community a tool to verify, analyze and reproduce the statements made in our paper.

How to use?
===========

To modify the parameters, edit the values N0 and q0 in params.h.

To run on an Unix machine with g++:
```
$ make
$ ./IBE
```

If GMP and NTL are not in a standard directory, you have to modify the CCFLAGS and LDFLAGS in the Makefile to indicate where they are.


A note on efficiency
====================
Since the publication of "Efficient Identity-Based Encryption over NTRU Lattices", this code has been updated several times, so do the timings claimed. Currently, the timings for encryption/decryption are ~10 times faster than what was claimed in the Proceedings version of "Efficient Identity-Based Encryption over NTRU Lattices".
