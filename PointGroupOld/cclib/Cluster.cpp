/*
This routine performs various clustering algorithms.
*/

#include <iostream>
#include <fstream>
using namespace std;

#include "mytype.h"
#include "Feature.h"
#include "Misc.h"
#include "Cluster.h"

#include <stdlib.h>
#include <math.h>

/*
compute the Euclidean distance SQUARED.
*/
real 
ComputeClusterDistance(const RealImage& image, 
             const RealImage& centroid,
             int index, int y, int x){
  int i;
  real d, dist;

  dist = 0;
  for(i = 0; i < image.NumBands(); ++i) {
    d = image.GetPixel(i,y,x)-centroid.GetPixel(index,i,0);
    dist += d * d;
  }
  return(dist);
}

/*
This routine performs kmean clustering.
im: input image with multiple planes where across the planes forms
    a feature vector.
nclust: # of clusters to be obtained
thr: convergence criterion.  If the change from the previous iteration
     is smaller than 'thr', the iteration stops.
iter: the maximum # of iterations.
cent: the routine returns the center of each cluster.
      if cent is NULL, then the routine doesn't return the centroids.
ret: if 1, the clusterd image (plus distance ratio) are returned,
     o.w. the result is replaced in 'im'.
init: if 1 the content of cent is used for the initial clustering.
*/

IntImage
KmeanImage(const RealImage& image,
           RealImage& cent,
           real thr, int iter) {
  IntImage res(1,image.NumRows(),image.NumCols());
  int n, done;
  int i,j,k,c;
  for(n=0,done=0; n<iter && !done; ++n) {
    /* compute the distance from every centroid */
    /* and determine which centroid is the closest */
    real dist;
    for(i = 0; i < image.NumRows(); ++i) {
      for(j = 0; j < image.NumCols(); ++j) {
        int mink=0;
        real min = ComputeClusterDistance(image, cent, 0, i, j);
        for(c=1; c<cent.NumBands(); ++c) {
          dist = ComputeClusterDistance(image, cent, c, i, j);
          if(dist < min) {
            min = dist;
            mink=c;
          }
        }
        res.SetPixel(0,i,j,mink);
      }
    }

    /* update the centroids */
    RealImage previous=cent;
    for(k=0; k<cent.NumBands(); ++k) {
      for(i=0; i<cent.NumRows(); ++i) {
        for(j=0; j<cent.NumCols(); ++j) {
          cent.SetPixel(k,i,j,0);
        }
      }
    }
    for(i=0; i<image.NumRows(); ++i) {
      for(j=0; j<image.NumCols(); ++j) {
        for(k = 0; k < image.NumBands(); ++k) {
          real val=image.GetPixel(k,i,j);
          c=res.GetPixel(0,i,j);
          cent.AddPixel(c,k,0,val);
          cent.AddPixel(c,k,1,val*val);
          cent.AddPixel(c,k,2,1.0);
        }
      }
    }
    for(i=0; i<cent.NumBands(); ++i) {
      for(j=0; j<cent.NumRows(); ++j) {
        real nump=cent.GetPixel(i,j,2);
        if(nump) {
          real sum=cent.GetPixel(i,j,0);
          real var=cent.GetPixel(i,j,1);
          real average=sum/nump;
          real std=sqrt(Max(.0,var/nump-average*average));
          cent.SetPixel(i,j,0,average);
          cent.SetPixel(i,j,1,std);
        }
      }
    }

    previous-=cent;
    /* compute the new total distance */
    dist = 0;
    for(i=0; i<previous.NumBands(); ++i) {
      for(j=0; j<previous.NumRows(); ++j) {
        real dd=previous.GetPixel(i,j,0);
        dist+=dd*dd;
      }
    }
    printf("Iteration %d: dist = %f\n", n, dist);

    if(dist < thr)
    done = 1;
  }
  return res;
}

/*
this routine computes the Maximum Absolute Column Sum Norm (MACSN)
for an image.
here, we consider a band of the image as a column vector
and data across the band at the same pixel location as a
row vector.

the formula for the MACSN is
||A||=max_k sum_{x,y} |a_kxy|
*/
//real
//MaxAbsColSumNorm(const RealImage& A) {
//}

void
UpdateUMatrix(RealImage& U,
              const RealImage& image,
              const RealImage& cent,
              int m) {

  int i,j,k;
  real* dist=new real[cent.NumBands()];
  for(i=0; i<image.NumRows(); ++i) {
    for(j=0; j<image.NumCols(); ++j) {
      real sum=0;
      for(k=0; k<cent.NumBands(); ++k) {
        dist[k]=ComputeClusterDistance(image,cent,k,i,j);
        sum+=1.0/dist[k]; //assuming m=2
      }
      for(k=0; k<cent.NumBands(); ++k) {
        U.SetPixel(k,i,j,1.0/(dist[k]*sum));
      }
    }
  }
}

void
UpdateFuzzyCentroid(RealImage& cent,
                    const RealImage& image,
                    const RealImage& U,
                    int m) {
  int i,j,k,l;
  for(k=0; k<cent.NumBands(); ++k) {
    for(l=0; l<cent.NumRows(); ++l) {
      real deno=0;
      real nume=0;
      for(i=0; i<image.NumRows(); ++i) {
        for(j=0; j<image.NumCols(); ++j) {
          real uval=U.GetPixel(k,i,j);
          deno+=uval*uval;
          nume+=uval*uval*image.GetPixel(l,i,j);
        }
      }
      cent.SetPixel(k,l,0,nume/deno);
    }
  }
}
  
/*
This routine performs kmean clustering.
im: input image with multiple planes where across the planes forms
    a feature vector.
nclust: # of clusters to be obtained
thr: convergence criterion.  If the change from the previous iteration
     is smaller than 'thr', the iteration stops.
iter: the maximum # of iterations.
m: a measure of fuzzyness (1.1<m<5)
cent: the routine returns the center of each cluster.
      if cent is NULL, then the routine doesn't return the centroids.
*/

RealImage
FuzzyCmeanImage(const RealImage& image,
                RealImage& cent,
                real thr, int m, int iter) {
  RealImage U(cent.NumBands(),image.NumRows(),image.NumCols());
  int n;
  for(n=0; n<iter; ++n) {
    //RealImage oldU=U;
    UpdateUMatrix(U,image,cent,m);
    //real eps=MaxAbsColSumNorm(U-oldU);
    //if(eps<thr)
    //  break;
    UpdateFuzzyCentroid(cent,image,U,m);
  }
  return U;
}
