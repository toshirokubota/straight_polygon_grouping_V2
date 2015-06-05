#ifndef _GaussJordan_H_
#define _GaussJordan_H_

#include <iostream>
#include <vector>
using namespace std;

#include <myDataType.h>

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

/*
solve Ax = v by Gauss-Jordan elimination.
the result is returned in v.
*/
bool
GaussJordan(const vector<double>& A, vector<double>& v, int n);

#endif /* _GaussJordan_H_ */