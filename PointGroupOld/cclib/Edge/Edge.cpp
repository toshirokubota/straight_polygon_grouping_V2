#include "Dipole.h"
//#include "Table.h"
//#define ALVAREZ

const real FieldBeta=10.0;

/*
A field definition.
*/
/*Vector Dipole::ComputeField(real y, real x, real beta, real sigma) {

  real yr = y-Y;
  real xr = x-X;
  real cs=cos(Theta);
  real sn=sin(Theta);
  real u=cs*xr+sn*yr;
  real v=-sn*xr+cs*yr;
  //real val=Mag*exp(-(u*u)/(2.*beta*beta)-(v*v)/(2.*sigma*sigma));
  real val=Mag*exp(-Abs(u)/(2.*beta*beta)-Abs(v)/(2.*sigma*sigma));
  return Vector(Theta, val);
}*/

/*
A field definition.
*/
/*Vector Dipole::ComputeField(real y, real x) {

  real yr = y-Y;
  real xr = x-X;
  real alpha=0;
  if(yr!=0 || xr!=0)
    alpha=atan2(yr,xr);
  real diff=alpha-Theta;
  real sn=sin(diff);
  //real snn=Abs(sn);
  real snn=sn*sn;
  real val;

#ifdef ALVAREZ
  val = Mag*exp(-(yr*yr+xr*xr)/32.0); //Alvarez diffusion
#else
  //val = Mag*exp(-FieldBeta*snn-(yr*yr+xr*xr)/32.0);
  val = Mag*exp(-FieldBeta*snn-(yr*yr+xr*xr)/32.0);
#endif
  Vector v(Theta, val);
  return v;
}*/

/*
More general field definition.  Used in the default UpdateField().
*/
Vector Dipole::ComputeField(real y, real x, real beta, real sigma) {
  real yr = y-Y;
  real xr = x-X;
  real alpha=0;
  if(yr!=0 || xr!=0)
    alpha=atan2(yr,xr);
  real diff=alpha-Theta;
  real sn=sin(diff);
  real snn=sn*sn;
  real val;
  val = Mag*exp(-beta*snn-(yr*yr+xr*xr)/(2*sigma*sigma));

  Vector v(Theta, val);
  return v;
}

/*
Table based field generation.
*/
Vector Dipole::ComputeFieldTable(real y, real x,
                                 real** radial_table,
                                 real* angular_table,
                                 real** atan_table) {
  int yi=(int)(Y-y);
  int xi=(int)(X-x);
  real val=Mag*radial_table[Abs(yi)][Abs(xi)];
  //real val=Mag*exp(-((real)yi*yi+xi*xi)/32.0);
  //val*=GetAngularTerm(Theta,yi,xi,angular_table,atan_table);

  Vector v(Theta, val);
  return v;
}

Vector Dipole::ComputeFieldTable(real y, real x,
                                 real** field_table,
                                 real** angle_table) {
  int range=3;
  int yi=(int)(Y-y)+range;
  int xi=(int)(X-x)+range;
  int k=yi*(2*range+1)+xi;
  int num_bin=37;
  real wgt=(real)(num_bin-1)/TwoPi;
  int theta=Round(wgt*(Theta+PI));
  real val=Mag*field_table[k][theta];
  //real phi=angle_table[k][theta];
  //Vector v(phi, val);
  Vector v(Theta, val);
  return v;
}


RealImage
FieldImage2RealImage(const FieldImage& field) {
  RealImage im=RealImage(2*field.NumBands(),field.NumRows(),field.NumCols());
  for(int k=0; k<field.NumBands(); ++k) {
    for(int i=0; i<im.NumRows(); ++i) {
      for(int j=0; j<im.NumCols(); ++j) {
        Field f=field.GetPixel(k,i,j);
        im.SetPixel(2*k,i,j,f.GetMagnitude());
        im.SetPixel(2*k+1,i,j,f.GetTheta());
      }
    }
  }

  return im;
}

RealImage
DipoleImage2RealImage(const DipoleImage& dipole) {
  RealImage im=RealImage(4,dipole.NumRows(),dipole.NumCols());
  for(int i=0; i<im.NumRows(); ++i) {
    for(int j=0; j<im.NumCols(); ++j) {
      Dipole dp=dipole.GetPixel(0,i,j);
      im.SetPixel(0,i,j,dp.GetMagnitude());
      im.SetPixel(1,i,j,dp.GetTheta());
      im.SetPixel(2,i,j,dp.GetY());
      im.SetPixel(3,i,j,dp.GetX());
    }
  }
  return im;
}

VectorImage
DipoleImage2VectorImage(const DipoleImage& dipole) {
  VectorImage im(1,dipole.NumRows(),dipole.NumCols());
  for(int i=0; i<im.NumRows(); ++i) {
    for(int j=0; j<im.NumCols(); ++j) {
      Vector v=dipole.GetPixel(0,i,j);
      im.SetPixel(0,i,j,v);
     }
  }
  return im;
}

VectorImage
FieldImage2VectorImage(const FieldImage& field) {
  VectorImage im(1,field.NumRows(),field.NumCols());
  for(int i=0; i<im.NumRows(); ++i) {
    for(int j=0; j<im.NumCols(); ++j) {
      Vector v=field.GetPixel(0,i,j);
      im.SetPixel(0,i,j,v);
     }
  }
  return im;
}
