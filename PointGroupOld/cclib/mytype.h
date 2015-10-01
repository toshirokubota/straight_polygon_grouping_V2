#ifndef _mytype_h_
#define _mytype_h_

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

#define _BOOL
#ifndef _BOOL
enum bool { false, true};
#define _BOOL
#endif

typedef double real;
typedef unsigned char uchar;

#ifdef AIMAGE

#include "Aimage.h"

typedef Aimage<float> RealImage;
typedef Aimage<unsigned char> ByteImage;

#define Image Aimage

#else

#include "Image.h"

typedef Image<real> RealImage;
typedef Image<uchar> ByteImage;
typedef Image<int> IntImage;

typedef vector<float> vReal;
typedef vector<int> vInt;

#endif /* AIMAGE*/


#endif /* mytype_h */
