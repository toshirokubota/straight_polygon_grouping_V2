
#include "ImageProc.h"
#include "Feature.h"
#include "Convolve.h"

#define DEBUG_ImageProc

//
// This routine converts a complex image into a real image.
//
RealImage
Complex2RealImage(const ComplexImage& image, Complex2Real_Op type) {
  RealImage result(image.NumBands(), image.NumRows(), image.NumCols());

  Complex cval;
  real rval;
  int num = image.NumBands()*image.NumRows()*image.NumCols();
  int i;
  switch(type) {
  case Complex2Real_Real:
    for(i=0; i<num; ++i) {
      cval = image.GetPixel(i);
      rval = cval.GetReal();
      result.SetPixel(i, rval);
    }
    break;
  case Complex2Real_Imag:
    for(i=0; i<num; ++i) {
      cval = image.GetPixel(i);
      rval = cval.GetImag();
      result.SetPixel(i, rval);
    }
    break;
  case Complex2Real_Mag:
    for(i=0; i<num; ++i) {
      cval = image.GetPixel(i);
      rval = cval.Magnitude();
      result.SetPixel(i, rval);
    }
    break;
  case Complex2Real_Phase:
    for(i=0; i<num; ++i) {
      cval = image.GetPixel(i);
      rval = cval.Phase();
      result.SetPixel(i, rval);
    }
    break;
  default:
    cerr << "Invalid Complex to Real operator type." << endl;
    assert(1);
  }

  return result;
}

//
// This routine converts a real image into a complex image.
//
ComplexImage
Real2ComplexImage(const RealImage& image, Complex2Real_Op type) {
  ComplexImage result(image.NumBands(), image.NumRows(), image.NumCols());

  Complex cval;
  int num = image.NumBands()*image.NumRows()*image.NumCols();
  int i;
  switch(type) {
  case Complex2Real_Real:
    for(i=0; i<num; ++i) {
      cval.SetReal(image.GetPixel(i));
      cval.SetImag(0);
      result.SetPixel(i, cval);
    }
    break;
  case Complex2Real_Imag:
    for(i=0; i<num; ++i) {
      cval.SetReal(0);
      cval.SetImag(image.GetPixel(i));
      result.SetPixel(i, cval);
    }
    break;
  default:
    cerr << "Invalid Real to Complex operator type." << endl;
    assert(1);
  }

  return result;
}

/*
This routine adjust the given image so that it has zero mean and unit
variance.
*/
void
NormalizeImage(RealImage& image) {
  real av = .0;
  real var = .0;
  int i;
  int num_pixels = image.NumBands()*image.NumRows()*image.NumCols();
  for(i=0; i<num_pixels; ++i) {
    av += image.GetPixel(i);
    var += image.GetPixel(i)*image.GetPixel(i);
  }
  av /= (real)(num_pixels);
  var = var/(real)(num_pixels)- av*av;
  cout << "average=" << av << " variance=" << var << endl;

  var = 1.0 / sqrt(var);

  for(i=0; i<num_pixels; ++i) {
    image.SetPixel(i, (image.GetPixel(i)-av)*var);
  }
}

/*
This routine adjust the given image so that it has zero mean and unit
variance.
*/
void
NormalizeImage(RealImage& image, real& mean, real& std) {
  mean = .0;
  real var = .0;
  int i;
  int num_pixels = image.NumBands()*image.NumRows()*image.NumCols();
  for(i=0; i<num_pixels; ++i) {
    mean += image.GetPixel(i);
    var += image.GetPixel(i)*image.GetPixel(i);
  }
  mean /= (real)(num_pixels);
  var = var/(real)(num_pixels)- mean*mean;
  std=sqrt(var);
  cout << "average=" << mean << " std=" << std << endl;

  real ivar = 1.0 / std;

  for(i=0; i<num_pixels; ++i) {
    image.SetPixel(i, (image.GetPixel(i)-mean)*ivar);
  }
}

/*
Reverse the operation done by NormalizeImage().
*/
void
UnnormalizeImage(RealImage& image, real mean, real std) {
  int i;
  int num_pixels = image.NumBands()*image.NumRows()*image.NumCols();
  for(i=0; i<num_pixels; ++i) {
    real val=image.GetPixel(i);
    val=val*std+mean;
    image.SetPixel(i,val);
  }
}

/*
Reverse the operation done by NormalizeImage().
*/
RealImage
UnnormalizeImage(const RealImage& image, real mean, real std) {
  int i;
  int num_pixels = image.NumBands()*image.NumRows()*image.NumCols();
  RealImage res(image.NumBands(),image.NumRows(),image.NumCols());
  for(i=0; i<num_pixels; ++i) {
    real val=image.GetPixel(i);
    val=val*std+mean;
    res.SetPixel(i,val);
  }
  return res;
}

/*
This routine adjust the given image so that the maximum = high
and the minimum=low.
*/
void
NormalizeImageRange(RealImage& image, real low, real high) {
  real max,min,val,wgt;
  int i;
  int num_pixels = image.NumBands()*image.NumRows()*image.NumCols();
  max=min=image.GetPixel(0);
  for(i=1; i<num_pixels; ++i) {
    val = image.GetPixel(i);
    max = (val > max) ? val: max;
    min = (val < min) ? val: min;
  }

  if(max - min > 0)
    wgt = (high-low)/(max-min);
  else
    wgt = 0;
  for(i=0; i<num_pixels; ++i) {
    val = image.GetPixel(i);
    val = (val - min) * wgt + low;
    image.SetPixel(i, val);
  }
}

real
ComputeMeanImage(const RealImage& image) {
  int numpixel=image.NumBands()*image.NumRows()*image.NumCols();
  real sum=0;
  for(int i=0; i<numpixel; ++i){
  	sum+=image.GetPixel(i);
  }
  return sum/numpixel;
}

RealImage
AdjustMeanImage(const RealImage& image, real mean) {
  real mean2=ComputeMeanImage(image);
  RealImage res(image.NumBands(),image.NumRows(),image.NumCols());
  int numpixel=image.NumBands()*image.NumRows()*image.NumCols();
  for(int i=0; i<numpixel; ++i){
  	res.SetPixel(i,image.GetPixel(i)+mean-mean2);
  }
  return res;
}

  

RealImage
NitzbergGradient(const RealImage& image) {
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
        result.SetPixel(2*k,i,j,gx);
        result.SetPixel(2*k+1,i,j,gy);
      }
    }
  }
  return result;
}

RealImage
LGNGradient(const RealImage& image) {
  RealImage result(image.NumBands()*2,image.NumRows(),image.NumCols());

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
        real gy=Abs((2.*s+d+e)-.5*(a+2.*b+c+f+2.*g+h));
        real gx=Abs((2.*s+b+g)-.5*(a+c+2.*d+2.*e+f+h));
        result.SetPixel(2*k,i,j,gx);
        result.SetPixel(2*k+1,i,j,gy);
      }
    }
  }
  return result;
}

RealImage
AsanoGradient(const RealImage& image) {
  //cout << "Asano Gradient:" << endl;
  RealImage result(image.NumBands()*2,image.NumRows(),image.NumCols());

  real wgt1=0.112727;
  real wgt2=0.274526;
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
        real gx= wgt1*(c+h-a-f)+wgt2*(e-d);
        real gy= wgt1*(f+h-a-c)+wgt2*(g-b);
        result.SetPixel(2*k,i,j,gx);
        result.SetPixel(2*k+1,i,j,gy);
      }
    }
  }
  return result;
}


/*
This routine performs Gaussian smoothing.
*/
RealImage 
GaussianSmoothing(const RealImage& image, real sgm, real epsilon) {

  RealImage gauss=GaussianFilter1D(sgm,epsilon);
  RealImage result=SeparableImageConvolution(image,gauss,gauss);

  return result;
}
 

/*
This routine rotate an image by 'angle' (in radian) using
binary interpolation.
The auxiliary file, 'valid', record which region has
pixels computed through interpolation of the original image
(valid=1), and which region has invalid data where
zeros are padded in (valid=0).
The routine also sets valid=0 to at 1 pixel width boarder
between valid and invalid regions so that gradient operations
can be performed totally within the valid region by consulting
'valid'.
*/

RealImage
RotateImage(const RealImage& image, ByteImage& valid, real angle) {
  assert(image.NumBands()==valid.NumBands() &&
         image.NumRows()==valid.NumRows() &&
         image.NumCols()==valid.NumCols());
  RealImage result=RealImage(image.NumBands(),image.NumRows(),
                             image.NumCols(), .0);

  real cy=(real)result.NumRows()*.5;
  real cx=(real)result.NumCols()*.5;
  int k;
  for(k=0; k<result.NumBands(); ++k) {
    for(int i=0; i<result.NumRows(); ++i) {
      real y=(real)i-cy;
      for(int j=0; j<result.NumCols(); ++j) {
        real x=(real)j-cx;
        real ux=x*cos(angle)-y*sin(angle)+cx;
        real uy=x*sin(angle)+y*cos(angle)+cy;
        int bx=(int)floor(ux);
        int by=(int)floor(uy);
        int ex=(int)ceil(ux);
        int ey=(int)ceil(uy);
        if(bx>=0 && by>=0 && ex<image.NumCols() && ey<image.NumRows()) {
          valid.SetPixel(k,i,j,1);
          real a=image.GetPixel(k,by,bx);
          real b=image.GetPixel(k,by,ex);
          real c=image.GetPixel(k,ey,bx);
          real d=image.GetPixel(k,ey,ex);
          real ty=uy-floor(uy);
          real tx=ux-floor(ux);
          real e=tx*b+(1.0-tx)*a;
          real f=tx*d+(1.0-tx)*c;
          real val=ty*f+(1.0-ty)*e;
          result.SetPixel(k,i,j,val);
        }
        else
          valid.SetPixel(k,i,j,0);
      }
    }
  }
  ByteImage valid2=valid;
  for(k=0; k<valid.NumBands(); ++k) {
    for(int i=0; i<valid.NumRows(); ++i) {
      for(int j=0; j<valid.NumCols(); ++j) {
        if(valid2.GetPixel(k,i,j)) { //erosion of valid
          if(!valid2.GetPixelZero(k,i-1,j-1)||
             !valid2.GetPixelZero(k,i-1,j)||
             !valid2.GetPixelZero(k,i-1,j+1)||
             !valid2.GetPixelZero(k,i,j-1)||
             !valid2.GetPixelZero(k,i,j+1)||
             !valid2.GetPixelZero(k,i+1,j-1)||
             !valid2.GetPixelZero(k,i+1,j)||
             !valid2.GetPixelZero(k,i+1,j+1))
            valid.SetPixel(k,i,j,0);
        }
      }
    }
  }
  return result;
}
