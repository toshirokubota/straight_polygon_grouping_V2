#ifndef _Feature_h_
#define _Feature_h_

#include "mytype.h"
#include "Complex.h"
#include "Convolve.h"
#include "ImageProc.h"
#include "Misc.h"
#include "Fourier.h"
#include <math.h>
#include <mytype.h>

RealImage
LawTextureFilter(const RealImage& image, 
                  bool L5, bool E5, bool S5, bool W5, bool R5);

ComplexImage
FanFilter(int bands, int rows, int cols, 
          real freq_low, real freq_high, real gauss_sgm);

RealImage
GaussianFilter2D(real ysgm, real xsgm, real theta, real epsilon);

RealImage
GaussianFilter1D(real sgm, real epsilon);

ComplexImage
GaborFilter(real ysgm, real xsgm, real alpha, real theta, real epsilon);

ComplexImage
GaborFilters(int numfilter, real ysgm, real xsgm, real alpha, real epsilon);


RealImage
FanFilterFeature(const RealImage& image, int num_orient, real sgm, 
                 real freq_low, real freq_high);

RealImage
LawTextureFeature(const RealImage& image, real gauss_sgm,
                  bool L5, bool E5, bool S5, bool W5, bool R5);

#endif /* Feature_h */
