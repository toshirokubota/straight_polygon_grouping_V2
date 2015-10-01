#ifndef _Misc_h_
#define _Misc_h_
#include "mytype.h"
#include "Complex.h"
//#include <kim.h>
#include <math.h>

#define LET_USE_THIS_ONE_INSTEAD_OF_KIM_ONES
#ifdef LET_USE_THIS_ONE_INSTEAD_OF_KIM_ONES
double rndm(int);
double gasdev(long jd);
#endif

RealImage AddNoiseImage(const RealImage& image, real sgm1, real sgm2);

RealImage ResizeImagePowerTwo(const RealImage& image);

ComplexImage ResizeImagePowerTwo(const ComplexImage& image);

template<class Item>
inline Item Max(Item a, Item b) {return ((a>b) ? a: b);}
template<class Item>
inline Item Min(Item a, Item b) {return ((a<b) ? a: b);}
template<class Item>
inline Item Abs(Item a) {return ((a>0) ? a: -a);}

inline int Round(real a) {return (int) ((a) > 0 ? ((a) + 0.5):((a) - 0.5));}

//const real PI=3.1415265358979323844;
#undef PI
const real PI=atan2(0,-1.0);
const real TwoPi=2.0*PI;
const real HalfPi=.5*PI;
const real QtrPi=.25*PI;
const real Qtr3Pi=.75*PI;

#endif /* Misc_h */
