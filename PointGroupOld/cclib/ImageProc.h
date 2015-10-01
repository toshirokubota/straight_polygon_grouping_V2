#ifndef _ImageProc_h_
#define _ImageProc_h_

#include <assert.h>
#include <math.h>

#include "mytype.h"
#include "Complex.h"

RealImage Complex2RealImage(const ComplexImage& image, Complex2Real_Op type);

ComplexImage Real2ComplexImage(const RealImage& image, Complex2Real_Op type);

void NormalizeImage(RealImage& image);

void NormalizeImage(RealImage& image, real& mean, real& var);

void UnnormalizeImage(RealImage& image, real mean, real std);

RealImage UnnormalizeImage(const RealImage& image, real mean, real std);

void NormalizeImageRange(RealImage& image, real high, real low);

real ComputeMeanImage(const RealImage& image); 

RealImage AdjustMeanImage(const RealImage& image, real mean); 

RealImage NitzbergGradient(const RealImage& image);

RealImage LGNGradient(const RealImage& image);

RealImage AndoGradient(const RealImage& image);

RealImage GaussianSmoothing(const RealImage& image, real sgm, real epsilon);

RealImage RotateImage(const RealImage& image, ByteImage& valid, real angle);

#endif /* ImageProc_h */
