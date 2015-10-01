
#include "Feature.h"

#include <assert.h>

const real Epsilon=1.0/256.0; //insignificant below this value in filter coeff.

/*
This routine computes Law's texture metrics.
L5, E5, S5, W5 and R5 is a boolean to switch on/off the cooresponding
filters.  They are:
L5 = [1,4,6,4,1]
E5 = [-1,-2,0,2,1]
S5 = [-1,0,2,0,-1]
W5 = [-1,2,0,-2,1]
R5 = [1,-4,6,-4,1]
*/

RealImage
LawTextureFilter(const RealImage& image, 
                  bool L5, bool E5, bool S5, bool W5, bool R5) {

  int bands=0;
  if(L5)  bands++;
  if(E5)  bands++;
  if(S5)  bands++;
  if(W5)  bands++;
  if(R5)  bands++;
  
  assert(bands!=0);
  RealImage result(bands*bands, image.NumRows(),image.NumCols());
  RealImage *filter;
  filter = new RealImage[bands];

  //
  // setting up the filter
  //
  int i;
  int k=0;

  if(L5) {
    real fil_l5[]={1.0/16, 4.0/16, 6.0/16, 4.0/16, 1.0/16};
    filter[k] = RealImage(1,1,5);
    for(i=0; i<filter[k].NumCols(); ++i)
      filter[k].SetPixel(0,0,i,fil_l5[i]);
    k++;
  }
  if(E5) {
    real fil_e5[]={-1.0/6, -2.0/6, .0, 2.0/6, 1.0/6};
    filter[k] = RealImage(1,1,5);
    for(i=0; i<filter[k].NumCols(); ++i)
      filter[k].SetPixel(0,0,i,fil_e5[i]);
    k++;
  }
  if(S5) {
    real fil_s5[]={-1.0/4, 0, 2.0/4, 0, -1.0/4};
    filter[k] = RealImage(1,1,5);
    for(i=0; i<filter[k].NumCols(); ++i)
      filter[k].SetPixel(0,0,i,fil_s5[i]);
    k++;
  }
  if(W5) {
    real fil_w5[]={-1.0/6, 2.0/6, 0, -2.0/6, 1.0/6};
    filter[k] = RealImage(1,1,5);
    for(i=0; i<filter[k].NumCols(); ++i)
      filter[k].SetPixel(0,0,i,fil_w5[i]);
    k++;
  }
  if(R5) {
    real fil_r5[]={1.0/16, -4.0/16, 6.0/16, -4.0/16, 1.0/16};
    filter[k] = RealImage(1,1,5);
    for(i=0; i<filter[k].NumCols(); ++i)
      filter[k].SetPixel(0,0,i,fil_r5[i]);
    k++;
  }
  //
  // Perform filtering
  //
  RealImage temp_im;
  for(k=0; k<bands; ++k) {
    for(int k2=0; k2<bands; ++k2) {
      temp_im = SeparableImageConvolution(image,filter[k],filter[k2]);
      result.insert(temp_im, k*bands+k2,0,0);
    }
  }

  delete [] filter;

  return result;
}

RealImage
GaussianFilter2D(real ysgm, real xsgm, real theta, real epsilon) {

  assert(epsilon < 1.0);
  int ysize, xsize;
  if(ysgm > xsgm)
    xsize=ysize=(int)(sqrt(-2.0*ysgm*ysgm*log(epsilon)) + .5); 
  else
    xsize=ysize=(int)(sqrt(-2.0*xsgm*xsgm*log(epsilon)) + .5); 

  cout << "Gauss Filter: size=(" << ysize << "," << xsize << ")\n";
  RealImage result(1, ysize, xsize);

  int i, j;
  int xcenter = (xsize-1)/2;
  int ycenter = (ysize-1)/2;
  ysgm = 1.0/(2.0*ysgm*ysgm);
  xsgm = 1.0/(2.0*xsgm*xsgm);
  real val;
  theta=-theta;  // the other direction for the coordinate mapping 
  for(i=0; i<ysize; ++i) {
    for(j=0; j<xsize; ++j) {
      real ux = cos(theta)*(real)(j-xcenter)+sin(theta)*(real)(i-ycenter);
      real uy = -sin(theta)*(real)(j-xcenter)+cos(theta)*(real)(i-ycenter);
      val = uy*uy*ysgm+ux*ux*xsgm;
      val = exp(-val);
      result.SetPixel(0,i,j, val);
    }
  }

  return result;
}

RealImage
GaussianFilter1D(real gauss_sgm, real epsilon) {

  assert(epsilon < 1.0);

  int fil_size=(int)(sqrt(-2.0*gauss_sgm*gauss_sgm*log(epsilon)) + .5);

  fil_size = 2*(fil_size-1)+1;
  cout << "Gaussian filter size= " << fil_size << endl;

  RealImage result(1, 1, fil_size);

  int i;
  int center = (fil_size-1)/2;
  real inv_sgm = 1.0/(2.0*gauss_sgm*gauss_sgm);

  real val;
  real sum=0.0;
  for(i=0; i<fil_size; ++i) {
    val = (i-center)*(i-center)*inv_sgm;
    val = exp(-val);
    result.SetPixel(i, val);
    sum+=val;
  }
  result *=(1.0/sum);
  return result;
}

RealImage
GaussianDerivative1D(real gauss_sgm, real epsilon) {

  assert(epsilon < 1.0);

  int fil_size=(int)(sqrt(-2.0*gauss_sgm*gauss_sgm*log(epsilon)) + .5);

  fil_size = 2*(fil_size-1)+1;
  cout << "Gaussian filter size= " << fil_size << endl;

  RealImage result(1, 1, fil_size);

  int i;
  int center = (fil_size-1)/2;
  real inv_sgm = 1.0/(2.0*gauss_sgm*gauss_sgm);

  real val;
  for(i=0; i<fil_size; ++i) {
    real x=(i-center);
    val = x*x*inv_sgm;
    val = -x*exp(-val)/gauss_sgm;
    result.SetPixel(i, val);
  }

  return result;
}

ComplexImage
FanFilter(int bands, int rows, int cols, 
          real freq_low, real freq_high, real gauss_sgm) {
  assert(bands > 3 && bands <= 8);

  real angle = PI / (real) bands;
  int center_y = rows / 2;
  int center_x = cols / 2;
  real low_angle = -angle*.5;
  real high_angle = low_angle + angle;
  Complex zero(0, 0), one(1.0, 0);
  ComplexImage result(bands+2, rows, cols, zero); ///

  int i, j, k;
  for(k=0; k<bands; ++k) {
    for(i=0; i<rows; ++i) {
      real y=(real)(i-center_y);
      for(j=0; j<cols; ++j) {
        real x=(real)(j-center_x);
        real ang=atan2((real)(i-center_y), (real)(j-center_x));

        real dist=sqrt(x*x+y*y);
        if(dist > freq_low && dist < freq_high) {
          if(ang < -angle*.5) ang += PI;
          else if(ang > PI-angle*.5) ang -= PI;
          if(ang >= low_angle && ang < high_angle)
            result.SetPixel(k,i,j, one);
        }
        else if(dist <= freq_low)
          result.SetPixel(result.NumBands()-2,i,j, one);
        else
          result.SetPixel(result.NumBands()-1,i,j, one);
      }
    }
    low_angle += angle;
    high_angle += angle;
  }

  // Gaussian smoothing of the filters

  RealImage filter=GaussianFilter1D(gauss_sgm, Epsilon);
  //filter = UnitSumNormalizeImage(filter);
  result = SeparableImageConvolution(result, filter, filter);

  return result;
}

ComplexImage
GaborFilter(real ysgm, real xsgm, real alpha, real theta, real epsilon) {

  assert(epsilon < 1.0);

  RealImage gauss = GaussianFilter2D(ysgm, xsgm, theta, epsilon);
  ComplexImage result(1, gauss.NumRows(), gauss.NumCols());

  int xcenter = (result.NumRows()-1)/2;
  int ycenter = (result.NumCols()-1)/2;
  for(int i=0; i<result.NumRows(); ++i) {
    for(int j=0; j<result.NumCols(); ++j) {
      real ux = sin(theta)*(real)(j-xcenter)+cos(theta)*(real)(i-ycenter);
      real r=cos(ux*alpha)*gauss.GetPixel(0,i,j);
      real im=sin(ux*alpha)*gauss.GetPixel(0,i,j);
      Complex cval(r, im);
      result.SetPixel(0,i,j, cval);
    }
  }

  return result;
}

ComplexImage
GaborFilters(int numfilter, real ysgm, real xsgm, real alpha, real epsilon) {

  assert(epsilon < 1.0);

  ComplexImage result;

  real theta=.0;
  real dtheta = PI / numfilter;

  for(int k=0; k<numfilter; ++k) {
    ComplexImage gabor=GaborFilter(ysgm,xsgm,alpha,theta,epsilon);
    if(k==0)
      result=ComplexImage(numfilter,gabor.NumRows(),gabor.NumCols());
    result.Insert(gabor, k, 0, 0);
    
    theta += dtheta;
  }

  return result;
}

RealImage
FanFilterFeature(const RealImage& image, int num_orient, real sgm, 
                 real freq_low, real freq_high) {

  RealImage resized = ResizeImagePowerTwo(image);
  ComplexImage filter = FanFilter(num_orient, 
                                  resized.NumRows(), resized.NumCols(),
                                  freq_low, freq_high, sgm);

  filter.circularShift(0, filter.NumRows()/2, filter.NumCols()/2);

  ComplexImage filtered = FourierTransform(resized, 0);

  filtered *= filter;

  filtered = FourierTransform(filtered, 1);

  RealImage result = Complex2RealImage(filtered, Complex2Real_Real);

  result = resized - result;
  //for(int i=0; i<result.NumBands()*result.NumRows()*result.NumCols(); ++i)
  //  result.SetPixel(i, abs(result.GetPixel(i)));
  result.extractROI(0, result.NumBands()-1, 0, image.NumRows()-1, 
                    0, image.NumCols()-1);
  return result;
}

RealImage
LawTextureFeature(const RealImage& image, real gauss_sgm,
                 bool L5, bool E5, bool S5, bool W5, bool R5) {

  RealImage result = LawTextureFilter(image, L5, E5, S5, W5, R5);

  for(int i=0; i<result.NumBands()*result.NumRows()*result.NumCols(); ++i){
    real val= result.GetPixel(i);
    val = (val < 0) ? -val: val;
    result.SetPixel(i, val);
  }
  result *= 2.0;

  RealImage filter = GaussianFilter1D(gauss_sgm, Epsilon);
  //filter = UnitSumNormalizeImage(filter);
  result = SeparableImageConvolution(result, filter, filter);
  if(L5==true)
    result.Insert(image, 0, 0, 0);

  return result;
}
