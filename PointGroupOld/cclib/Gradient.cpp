
#include "Gradient.h"

RealImage
NitzbergGradientPolar(const RealImage& image) {
  RealImage result(image.NumBands()*2,image.NumRows(),image.NumCols());

  real coeff=(sqrt(2.0)-1)/(2-sqrt(2.0));
  real wgt=1.0/(2.0+4*coeff);
  int i,j,k;
  for(k=0; k<image.NumBands(); ++k) {
    for(i=0; i<image.NumRows(); ++i) {
      for(j=0; j<image.NumCols(); ++j) {
        real s=image.GetPixel(k,i,j);
        real a=image.GetPixelDefault(k,i-1,j-1,s);
        real b=image.GetPixelDefault(k,i-1,j,s);
        real c=image.GetPixelDefault(k,i-1,j+1,s);
        real d=image.GetPixelDefault(k,i,j-1,s);
        real e=image.GetPixelDefault(k,i,j+1,s);
        real f=image.GetPixelDefault(k,i+1,j-1,s);
        real g=image.GetPixelDefault(k,i+1,j,s);
        real h=image.GetPixelDefault(k,i+1,j+1,s);
        real gx= wgt*((e-d)+coeff*((c-a)+(h-f)));
        real gy= wgt*((b-g)+coeff*((a-f)+(c-h)));
        real mag=sqrt(gx*gx+gy*gy);
		real theta=atan(gx/gy);
        result.SetPixel(2*k,i,j,mag);
        if(mag==0)
          result.SetPixel(2*k+1,i,j,.0);
        else 
          result.SetPixel(2*k+1,i,j,theta);
      }
    }
  }
  return result;
}

RealImage
CannyGradientPolar(const RealImage& image, 
              real sigma, real low, real high) {
  double* data=new double[image.NumRows()*image.NumCols()];
  double* deriv=new double[image.NumRows()*image.NumCols()];
  double* mag=new double[image.NumRows()*image.NumCols()];
  double* ori=new double[image.NumRows()*image.NumCols()];
  int i;
  for(i=0; i<image.NumRows()*image.NumCols();++i) {
    data[i]=(double)image.GetPixel(i);
    deriv[i]=0;
    mag[i]=0;
    ori[i]=0;
  }    
  
  CannyEdgeDetect(data,deriv,mag,ori,sigma, 
                  image.NumRows(),image.NumCols(),low,high);
 
  RealImage edge(2,image.NumRows(),image.NumCols());
  for(i=0; i<image.NumRows()*image.NumCols();++i) {
    edge.SetPixel(i,data[i]);
    //ori[i]=-ori[i]+HalfPi;
	ori[i]=atan(cos(ori[i])/sin(ori[i]));
    //if(ori[i]>HalfPi) ori[i]-=PI;
    //else if(ori[i]<-HalfPi) ori[i]+=PI;
    edge.SetPixel(i+image.NumRows()*image.NumCols(),ori[i]);
  }
  delete [] data;
  delete [] deriv;
  delete [] mag;
  delete [] ori;

  return edge;
}

RealImage
CannyNitzbergGradientPolar(const RealImage& image, 
              real sigma, real low, real high) {
	RealImage grad1=NitzbergGradientPolar(image);
  RealImage grad2=CannyGradientPolar(image,sigma,low,high);
  grad1.extractROI(1,1,0,grad1.NumRows()-1,0,grad1.NumCols()-1);
  grad2.insert(grad1,1,0,0);
  return grad2;
}

RealImage
MarrFilter(real sgm, real thr) {
	int radius=Round(-log(thr)*2*sgm*sgm);
	RealImage filter(1,2*radius+1,2*radius+1,0);
	for(int i=0; i<filter.NumRows(); ++i) {
		real y=(real)(i-radius);
		for(int j=0; j<filter.NumCols(); ++j) {
			real x=(real)(j-radius);
			real r2=y*y+x*x;
			real val=(1-r2/(sgm*sgm))*exp(-r2/(2.*sgm*sgm));
			filter.SetPixel(0,i,j,val);
		}
	}
	return filter;
}

RealImage
MarrGradientPolar(const RealImage& image, real sigma) {
	real thres=0.001;
	RealImage marr=MarrFilter(sigma,thres);
	RealImage gauss=GaussianFilter1D(sigma,thres);
	RealImage lap=ImageConvolution(image,marr);
	InteractiveImageWrite(lap,"temp2.pgm",1);
	RealImage grad=NitzbergGradientPolar(image);
	RealImage res(2,image.NumRows(),image.NumCols(),0);
	
	for(int i=0; i<res.NumRows(); ++i) {
		for(int j=0; j<res.NumCols(); ++j) {
			real val=lap.GetPixel(0,i,j);
			if(val*lap.GetPixelRepeat(0,i+1,j)<=0 ||
				val*lap.GetPixelRepeat(0,i,j+1)<=0) {
				res.SetPixel(0,i,j,Max(Abs(val),sigma));
				res.SetPixel(1,i,j,grad.GetPixel(0,i,j));
			}
		}
	}
	return res;
}

