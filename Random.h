#ifndef LIBE_RANDOM_H
#define LIBE_RANDOM_H

#include "params.h"
#include "Sampling.h"

vec_ZZ RandomVector();
ZZX RandomPoly(const unsigned int degree);
ZZX RandomPolyFixedSqNorm(const ZZ& SqNorm, const unsigned int degree);

#endif
