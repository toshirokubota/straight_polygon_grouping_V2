#ifndef _Svd_h_
#define _Svd_h_

#include "mytype.h"
#include "Misc.h"
#include <math.h>
#include <stdlib.h>

static real at, bt, ct;

/* PHTHAG computes sqrt(a*a+b*b) w/o destructive overflow
or underflow */
 
#define PYTHAG(a,b) ((at=Abs(a) ) > (bt=Abs(b)) ? \
                     (ct=bt/at,at*sqrt(1.0+ct*ct)) : \
                     (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)) : 0))

#define SIGN(a,b) ((b) >= 0 ? Abs(a) : -Abs(a))

void
SVD_Image(RealImage& A, RealImage& W, RealImage& V);

 
RealImage
SVD_Approximate(const RealImage& U, const RealImage& W, 
                  const RealImage& V, int order);

void
SVOSD_Image(const RealImage& A, RealImage& hf, RealImage& vf, 
            int order);

RealImage
SVOSD_ApproximationMatrix(const RealImage& fil);

void
SVOSD_ShuffleMatrix(RealImage& U, RealImage& W, RealImage& V);

RealImage
SVOSD_Approximate(const RealImage& cmmn, const RealImage& indp, int order);

#endif /* Svd_h */
