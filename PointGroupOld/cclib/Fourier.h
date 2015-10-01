#ifndef _Fourier_h_
#define _Fourier_h_

#include "Complex.h"
#include "Misc.h"
#include <math.h>

ComplexImage
FourierTransform(const ComplexImage& image, int direction);

ComplexImage
FourierTransform(const RealImage& image, int direction);

#endif /* Fourier_h */
