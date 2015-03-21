#ifndef LIBE_SAMPLING_H
#define LIBE_SAMPLING_H

#include "params.h"


unsigned int Sample0(unsigned long alea);
unsigned int Sample1(const unsigned int k);
signed int Sample2(const unsigned int k);
signed int Sample3(const RR_t sigma128);
signed int Sample4(RR_t c, RR_t sigma);


#endif
