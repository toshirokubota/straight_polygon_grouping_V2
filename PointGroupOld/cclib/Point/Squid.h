#ifndef _Squid_h_
#define _Squid_h_

#include <iostream>
#include <fstream>
using namespace std;

#include <stdlib.h>
#include <assert.h>
#include <mytype.h>
#include <Cpoint.h>
#include <PointLink.h>

Cpoint*
ReadSquidPoints(char* filename, int& h, int& w, int& n);

Cpoint*
ReadPetrakisPoints(char* filename, int& h, int& w, int& n, int scale);

Cpoint*
ReadBooksteinPoints(char* filename, int& h, int& w, real scale, int& n);

#endif /* Squid_h */
