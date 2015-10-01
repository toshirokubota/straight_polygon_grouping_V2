/*
Eigen.cpp implements various eigenvalue/eigenvector related
routines, in particular Principle Component Analysis.
*/

#include "Eigen.h"

#include <assert.h>

RealImage
CovarianceMatrixImageBand(const RealImage& image){
	RealImage res(1,image.NumBands(),image.NumBands(),.0);
	int nump=image.NumRows()*image.NumCols();
	real wgt=1.0/sqrt((real)nump);
	int i,j;
	for(i=0; i<image.NumRows(); ++i) {
		for(j=0; j<image.NumRows(); ++j) {
			real v1=wgt*image.GetPixel(0,i,j);
			real v2=wgt*image.GetPixel(1,i,j);
			real v3=wgt*image.GetPixel(2,i,j);
			res.AddPixel(0,0,0,v1*v1);
			res.AddPixel(0,0,1,v1*v2);
			res.AddPixel(0,0,2,v1*v3);
			res.AddPixel(0,1,0,v2*v1);
			res.AddPixel(0,1,1,v2*v2);
			res.AddPixel(0,1,2,v2*v3);
			res.AddPixel(0,2,0,v3*v1);
			res.AddPixel(0,2,1,v3*v2);
			res.AddPixel(0,2,2,v3*v3);
		}
	}
	return res;
}

RealImage
ComputeEigenVectors(const RealImage& cov, bool whitening) {
	int i,j;

	RealImage res(1,cov.NumCols(),cov.NumCols());

	int nrot;
	float* a=new float[cov.NumCols()*cov.NumRows()];
	for(i=0; i<cov.NumRows()*cov.NumCols(); ++i)
		a[i]=cov.GetPixel(i);
	float* v=new float[cov.NumCols()*cov.NumRows()];
	float* d=new float[cov.NumRows()];
	jacobi(a,cov.NumCols(),d,v,&nrot);
	
	if(whitening) {
		for(i=0; i<cov.NumRows()*cov.NumCols(); ++i)
			res.SetPixel(i,v[i]/sqrt(d[i%cov.NumRows()]));
	}
	else {
		for(i=0; i<cov.NumRows()*cov.NumCols(); ++i)
			res.SetPixel(i,v[i]);
	}
		
	delete [] a;
	delete [] v;
	delete [] d;
	return res;
}

RealImage
TransposeImage(const RealImage image) {
	RealImage res(image.NumBands(),image.NumCols(),image.NumRows());
	int i,j,k;
	for(k=0; k<image.NumBands(); ++k) {
    for(i=0; i<image.NumRows(); ++i) {
      for(j=0; j<image.NumCols(); ++j) {
      	res.SetPixel(k,j,i,image.GetPixel(k,i,j));
      }
    }
  }
  return res;
}

RealImage
LinearTransform(const RealImage& image, const RealImage& tr) {

	RealImage res(tr.NumCols(),image.NumRows(),image.NumCols());
	real* mat=new real[tr.NumRows()*tr.NumCols()];
	real* vec=new real[tr.NumCols()];
	//real mat[9];
	//real vec[3];
	int i,j,k,n;
	for(i=0; i<tr.NumRows(); ++i) 
		for(j=0; j<tr.NumCols(); ++j) 
			mat[i*tr.NumCols()+j]=tr.GetPixel(0,i,j);
			
	for(i=0; i<image.NumRows(); ++i) {
		for(j=0; j<image.NumCols(); ++j) {
			for(k=0; k<image.NumBands(); ++k) {
				vec[k]=image.GetPixel(k,i,j);
			}
			for(k=0; k<res.NumBands(); ++k) {
			  real sum=0;
			  for(n=0; n<tr.NumCols(); ++n) {
			  	sum+=mat[n*tr.NumCols()+k]*vec[n];
			  }
			  res.SetPixel(k,i,j,sum);
			}
		}
	}
	delete [] mat;
	delete [] vec;
	
	return res;
}

RealImage
PrincipleComponentAnalysis(const RealImage& image, RealImage& eigen){
	int i,j;
	RealImage res;
	RealImage cov=CovarianceMatrixImageBand(image);
	eigen=ComputeEigenVectors(cov,false);
	res=LinearTransform(image,eigen);
	cout << cov;
	cout << eigen;
	
	return res;
}	

RealImage
WhiteningFilter(const RealImage& image, RealImage& eigen){
	int i,j;
	RealImage res;
	RealImage cov=CovarianceMatrixImageBand(image);
	eigen=ComputeEigenVectors(cov,true);
	res=LinearTransform(image,eigen);
	
	return res;
}	

RealImage
TakeDifference(const RealImage& image) {
	RealImage res(image.NumBands(),image.NumRows(),image.NumCols());
	for(int k=0; k<image.NumBands(); ++k) {
		for(int i=0; i<image.NumRows(); ++i) {
			for(int j=0; j<image.NumRows(); ++j) {
				real v1=image.GetPixel(k,i,j);
				real v2=image.GetPixelRepeat(k,i,j+1);
				res.AddPixel(k,i,j,v1-v2);
			}
		}
	}
	return res;
}

RealImage
MaximumNoiseFraction(const RealImage& image){
	int i,j;
	RealImage res,res2;
	RealImage dif=TakeDifference(image);
	RealImage cov=CovarianceMatrixImageBand(dif);
	RealImage eigen1=ComputeEigenVectors(cov,true);
	res=LinearTransform(image,eigen1);

	RealImage cov2=CovarianceMatrixImageBand(res);
	RealImage eigen2=ComputeEigenVectors(cov,false);
	res2=LinearTransform(image,eigen1);
	res2=LinearTransform(image,eigen2);
	
	return res2;
}	

