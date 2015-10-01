#ifndef _GRADIENT_H_
#define _GRADIENT_H_
#include <mytype.h>
#include <Misc.h>
#include <ImageProc.h>
#include <Feature.h>
//#include <canny.h>

RealImage
NitzbergGradientPolar(const RealImage& image);

RealImage
CannyGradientPolar(const RealImage& image, 
              real sigma, real low, real high);

RealImage
CannyNitzbergGradientPolar(const RealImage& image, 
              real sigma, real low, real high);

RealImage
MarrGradientPolar(const RealImage& image, real sigma);

#endif /* _GRADIENT_H_ */