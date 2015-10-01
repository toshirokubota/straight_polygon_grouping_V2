#ifndef _Eigen_h_
#define _Eigen_h_
/*
Eigen.cpp implements various eigenvalue/eigenvector related
routines, in particular Principle Component Analysis.
*/

#include "mytype.h"
#include "Jacobi.h"

RealImage
TransposeImage(const RealImage image);

RealImage
LinearTransform(const RealImage& image, const RealImage& tr);

RealImage
CovarianceMatrixImageBand(const RealImage& image);

RealImage
ComputeEigenVectors(const RealImage& covariance, bool whitening);

RealImage
PrincipleComponentAnalysis(const RealImage& image, RealImage& eigen);

RealImage
WhiteningFilter(const RealImage& image, RealImage& eigen);

RealImage
MaximumNoiseFraction(const RealImage& image);

#endif /* Eigen_h */
