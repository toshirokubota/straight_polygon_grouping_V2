#ifndef _Cluster_h_
#define _Cluster_h_

#include <assert.h>
#include <math.h>

#include "mytype.h"

real 
ComputeClusterDistance(const RealImage& image, 
                       const RealImage& centroid,
                       int index, int y, int x);

IntImage 
KmeanImage(const RealImage& image, RealImage& cent, real thr, int iter);

RealImage 
FuzzyCmeanImage(const RealImage& image, RealImage& cent, 
                real thr, int m, int iter);

#endif /* Cluster_h */
