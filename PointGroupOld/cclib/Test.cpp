/*

*/

#include <iostream.h>
#include <fstream.h>
#include <vector.h>

#include <mytype.h>
#include <ImageProc.h>
#include <Feature.h>
#include <Convolve.h>
#include <Misc.h>

#include "Dipole.h"
#include "Update.h"
#include "Postp.h"
#include "Diffuse.h"
#include "Initialize.h"

//#define GL_DISPLAY
#ifdef GL_DISPLAY
#include "glr.h"
#endif

#include <stdlib.h>
#include <math.h>

typedef vector<Dipole> DipoleVector;
typedef Image<DipoleVector> DipoleVectorImage;

int MaxIter=20;
int FieldRange=3;
real FieldThres=.0;
real FieldSigma=.5;
real FieldBeta=0.0;
real IncDelta=0.25;
real CornerThres=5.0;
int ScaleParam=4;

char infile[256];
char outfile[256];

void
Usage(char* program) {
  cerr << program << ": usage" << endl;
  cerr << "Required arguments:" << endl;
  cerr << "\t -i <input file name>" << endl;
  cerr << "\t -o <output file name>" << endl;

  cerr << "Optional arguments:" << endl;
  cerr << "\t -m (int) set the number of update iterations" << endl;
  cerr << "\t \t default value is " << MaxIter << endl;

  cerr << "\t -r (int) set the range of the field" << endl;
  cerr << "\t \t default value is " << FieldRange << endl;

  cerr << "\t -b (real) set the beta value for the field" << endl;
  cerr << "\t \t default value is " << FieldBeta << endl;

  cerr << "\t -g (real) set the sigma value for the field" << endl;
  cerr << "\t \t default value is " << FieldSigma << endl;

  cerr << "\t -t (real) set the threshold value for the field" << endl;
  cerr << "\t \t default value is " << FieldThres << endl;

  cerr << "\t -c (real) set the threshold value for the corner" << endl;
  cerr << "\t \t default value is " << FieldThres << endl;

  cerr << "\t -d (real) set the incremental delta" << endl;
  cerr << "\t \t default value is " << IncDelta << endl;

  cerr << "\t -sf (int) set the scale factor"<<endl;
  cerr << "\t \t default value is " << ScaleParam << endl;
}

void
PrintParameters() {
  cout << "input file: " << infile << endl;
  cout << "output file: " << outfile << endl;
  cout << "# of update iterations: " << MaxIter << endl;
  cout << "field beta parameter: " << FieldBeta << endl;
  cout << "field range parameter: " << FieldRange << endl;
  cout << "field sigma parameter: " << FieldSigma << endl;
  cout << "field threshold parameter: " << FieldThres << endl;
  cout << "incremental delta value: " << IncDelta << endl;
  cout << "corner threshold parameter: " << CornerThres << endl;
  cout << "scaling factor: " << ScaleParam << endl;
}
  

int
ReadArguments(int argc, char* argv[]) {
  int in_found=0;
  int out_found=0;
  for(int i=1; i<argc; ++i) {
    if(strcmp(argv[i],"-i")==0) {
      if(argc>i+1) {
        in_found=1;
        strcpy(infile,argv[++i]);
      }
      else
        return 1; // no file specified
    }
    else if(strcmp(argv[i],"-o")==0) {
      if(argc>i+1) {
        out_found=1;
        strcpy(outfile,argv[++i]);
      }
      else
        return 2; // no file specified
    }
    else if(strcmp(argv[i],"-m")==0) {
      if(argc>i+1) {
        MaxIter=atoi(argv[++i]);
      }
      else
        return 3; // no value specified
    }
    else if(strcmp(argv[i],"-b")==0) {
      if(argc>i+1) {
        FieldBeta=strtod(argv[++i],NULL);
      }
      else
        return 4; // no value specified
    }
    else if(strcmp(argv[i],"-g")==0) {
      if(argc>i+1) {
        FieldSigma=strtod(argv[++i],NULL);
      }
      else
        return 5; // no value specified
    }
    else if(strcmp(argv[i],"-t")==0) {
      if(argc>i+1) {
        FieldThres=strtod(argv[++i],NULL);
      }
      else
        return 6; // no value specified
    }
    else if(strcmp(argv[i],"-r")==0) {
      if(argc>i+1) {
        FieldRange=atoi(argv[++i]);
      }
      else
        return 7; // no value specified
    }
    else if(strcmp(argv[i],"-d")==0) {
      if(argc>i+1) {
        IncDelta=strtod(argv[++i],NULL);
      }
      else
        return 8; // no value specified
    }
    else if(strcmp(argv[i],"-c")==0) {
      if(argc>i+1) {
        CornerThres=strtod(argv[++i],NULL);
      }
      else
        return 8; // no value specified
    }
    else if(strcmp(argv[i],"-sf")==0) {
      if(argc>i+1) {
        ScaleParam=atoi(argv[++i]);
      }
      else
        return 9; // no value specified
    }
  }
  if(!in_found)
    return 11;
  if(!out_found)
    return 12;
  else
    return 0;
}
void
InitializeDipole(DipoleVectorImage& dipole, 
                 const RealImage& image, real thr) {

  RealImage grd=NitzbergGradient(image);
  int i,j;
  for(i=0; i<dipole.NumRows(); ++i) {
    for(j=0; j<dipole.NumCols(); ++j) {
      Dipole d;
      d.SetX((real) j);
      d.SetY((real)(dipole.NumRows()-i-1));
      real gh=grd.GetPixel(0,i,j);
      real gv=grd.GetPixel(1,i,j);
      real mag=sqrt(gv*gv+gh*gh);
      if(mag>thr) {
        real phi=atan2(gh,-gv);
        d.SetTheta(phi);
        d.SetMagnitude(mag);
      }
      else {
        d.SetMagnitude(0);
      }
      DipoleVector vd=dipole.GetPixel(0,i,j);
      vd.push_back(d);
      dipole.SetPixel(0,i,j,vd);
    }
  }
}


/*
This routine computes the edge field at given location (y,x)
*/
Vector
ComputeField(const DipoleVectorImage& dipole,
             real y, real x, int h, int w, 
             real beta, real sigma, int range) {
  Vector f(0,0);
  for(int i=h-range; i<h+range; i++) {
    if(i>=0 && i<dipole.NumRows()) {
      for(int j=w-range; j<w+range; j++) {
        if(j>=0 && j<dipole.NumCols()) {
          DipoleVector dv=dipole.GetPixel(0,i,j);
          DipoleVector::iterator p;
          for(p=dv.begin(); p<dv.end(); ++p) {
            Vector v=p->ComputeField(y,x,beta,sigma);
            f+=v;
          }
        }
      }
    }
  }

  return f;
}

/*
This routine computes the cornerness meature at given location (y,x)
*/
real
ComputeCorner(const DipoleVectorImage& dipole,
             real y, real x, int h, int w, 
             real beta, real sigma, int range) {
  Vector f(0,0);
  real exx=.0;
  real eyy=.0;
  real exy=.0;
  int count=0;
  for(int i=h-range; i<h+range; i++) {
    if(i>=0 && i<dipole.NumRows()) {
      for(int j=w-range; j<w+range; j++) {
        if(j>=0 && j<dipole.NumCols()) {
          DipoleVector dv=dipole.GetPixel(0,i,j);
          DipoleVector::iterator p;
          for(p=dv.begin(); p<dv.end(); ++p) {
            Vector v=p->ComputeField(y,x,beta,sigma);
            real ex=v.GetXComp();
            real ey=v.GetYComp();
            exx+=ex*ex;
            eyy+=ey*ey;
            exy+=ex*ey;
            count++;
          }
        }
      }
    }
  }

  real corner=(exx+eyy-sqrt((exx-eyy)*(exx-eyy)+4.*exy*exy))/count;
  return corner;
}

void
hereUpdateField(FieldImage& field, const DipoleVectorImage& dipole, 
                real thres, real beta, real sigma, int range) {
  Field f;
  int i,j,n,m;
  // Initialize the field to 0
  for(i=0; i<field.NumRows(); ++i) {
    for(j=0; j<field.NumCols(); ++j) {
      f = field.GetPixel(0,i,j);
      real yf=f.GetY();
      real xf=f.GetX();
      real cnr;
      Vector v=ComputeField(dipole,yf,xf,i,j,beta,sigma,range);
      f.SetMagnitude(v.GetMagnitude());
      f.SetTheta(v.GetTheta());
      field.SetPixel(0,i,j,f);
    }
  }
}

void
UpdateSeparableField(FieldImage& field, const DipoleVectorImage& dipole, 
                     real sigma, real epsilon) {
  RealImage fcomp(2,dipole.NumRows(),dipole.NumCols(),0);
  for(int i=0; i<dipole.NumRows(); ++i) {
    for(int j=0; j<dipole.NumCols(); ++j) {
      DipoleVector dv = dipole.GetPixel(0,i,j);
      DipoleVector::iterator p;
      p=dv.begin();  // assume there is only one dipole each site
      fcomp.SetPixel(0,i,j,p->GetXComp());
      fcomp.SetPixel(1,i,j,p->GetYComp());
    }
  }
  RealImage gauss=GaussianFilter1D(sigma,epsilon);
  RealImage res=SeparableImageConvolution(fcomp,gauss,gauss);

  for(int i=0; i<field.NumRows(); ++i) {
    for(int j=0; j<field.NumCols(); ++j) {
      real vx=res.GetPixel(0,i,j);
      real vy=res.GetPixel(1,i,j);
      real mag=sqrt(vx*vx+vy*vy);
      real theta=atan2(vy,vx);
      
      Field f = field.GetPixel(0,i,j);
      real yf=f.GetY();
      real xf=f.GetX();
      f.SetMagnitude(mag);
      f.SetTheta(theta);
      field.SetPixel(0,i,j,f);
    }
  }
}

int
FindMaximumField(Vector* v, int num) {
  real max=v[0].GetMagnitude();
  int index=0;
  for(int i=0; i<num; ++i) {
    real val=v[i].GetMagnitude();
    if(val>max) {
      max=val;
      index=i;
    }
  }
  return index;
}

/* 
This routine adjust not only the theta of the dipole but 
the location, thus performs sub-pixel edge detection.
*/

UpdateDipoleVector(DipoleVectorImage& dipole,
                  real beta, real sigma, real thres, 
                  real delta, real cthres, int range) {
  int count=0;
  int count2=0;
  //DipoleVectorImage res(1,dipole.NumRows(),dipole.NumCols());
  RealImage corner(1,dipole.NumRows(),dipole.NumCols(),.0);
  for(int i=0; i<dipole.NumRows(); ++i) {
    for(int j=0; j<dipole.NumCols(); ++j) { 
      DipoleVector dv = dipole.GetPixel(0,i,j);
      DipoleVector::iterator p;
      for(p=dv.begin(); p<dv.end(); ++p) {
        if(p->GetMagnitude()>thres){
          real x[3], y[3], cnr;
          Vector f[3];
          y[0]=p->GetY();
          x[0]=p->GetX();
          cnr=ComputeCorner(dipole,y[0],x[0],i,j,beta,sigma,range);
          corner.AddPixel(0,i,j,cnr);
          if(cnr>cthres) {
            Dipole dp1=*p;
            real theta=p->GetTheta();
            dp1.SetY(y[0]+delta*sin(theta));
            dp1.SetX(x[0]+delta*cos(theta));
            dv.push_back(dp1);
            Dipole dp2=*p;
            dp2.SetY(y[0]-delta*sin(theta));
            dp2.SetX(x[0]-delta*cos(theta));
            dv.push_back(dp2);
            count2+=2;
          }
          int index;
          f[0]=ComputeField(dipole,y[0],x[0],i,j,beta,sigma,range);
          real theta=p->GetTheta();
          y[1]=y[0]+delta*sin(theta+HalfPi);
          x[1]=x[0]+delta*cos(theta+HalfPi);
          f[1]=ComputeField(dipole,y[1],x[1],i,j,beta,sigma,range);
          y[2]=y[0]-delta*sin(theta+HalfPi);
          x[2]=x[0]-delta*cos(theta+HalfPi);
          f[2]=ComputeField(dipole,y[2],x[2],i,j,beta,sigma,range);
          index=FindMaximumField(f,3);
          if(index!=0) count++;
          p->SetY(y[index]);
          p->SetX(x[index]);
          p->SetTheta(f[index].GetTheta());
          //compute new array index 
          real yy=p->GetY();
          real xx=p->GetX();
          int m=Round((real)dipole.NumRows()-yy-1.0);
          m=Max(0,Min(m,dipole.NumRows()-1));
          int n=Round(xx);
          n=Max(0,Min(n,dipole.NumCols()-1));
          if(i!=m || j!=n) {
            DipoleVector dv2=dipole.GetPixel(0,m,n);
            dv2.push_back(*p);
            dipole.SetPixel(0,m,n,dv2);
            dv.erase(p);
            p--;  //so that for loop doesn't skip the next item
          }
        }
      }
      dipole.SetPixel(0,i,j,dv);
    }
  }
  cout << "Total number of moves= " << count << endl;
  cout << "Newly created edges= " << count2 << endl;
  corner.WritePnmFile("corner.pgm",IO_PGM,1);
}

DipoleVectorImage
DipoleExpandVector(const DipoleVectorImage& dipole, int scale) {
  
  DipoleVectorImage res(1,dipole.NumRows()*scale,dipole.NumCols()*scale);

  for(int i=0; i<dipole.NumRows(); ++i) {
    for(int j=0; j<dipole.NumCols(); ++j) {
      DipoleVector dv = dipole.GetPixel(0,i,j);
      DipoleVector::iterator p;
      for(p=dv.begin(); p<dv.end(); ++p) {
        // get the new coordinate
        real y=p->GetY()*scale;
        real x=p->GetX()*scale;
        int m=Round(res.NumRows()-y-1);
        m=Max(0,Min(m,res.NumRows()-1));
        int n=Round(x);
        n=Max(0,Min(n,res.NumCols()-1));
        
        //check if there is something already in it.
        DipoleVector dv2=res.GetPixel(0,m,n);
        if(dv2.size()>0) {
        // compare with the existing one
          DipoleVector::iterator p2;
          p2=dv2.begin();  // there is at most one
          if(p->GetMagnitude()>p2->GetMagnitude()) {
            Dipole dp(p->GetTheta(),p->GetMagnitude(),y,x);
            dv2.clear();  // replace with the old one
            dv2.push_back(dp);
          }
        }
        else {
          Dipole dp(p->GetTheta(),p->GetMagnitude(),y,x);
          dv2.push_back(dp);
        }
        res.SetPixel(0,m,n,dv2);
      }
    }
  }
  
  return res;
}

void main(int argc,char *argv[])
{
  if(ReadArguments(argc,argv)!=0) {
    Usage(argv[0]);
    abort();
  }
  else
    PrintParameters();
  RealImage image;
  image.ReadPnmFile(infile);

  DipoleVectorImage dipole(1,image.NumRows(),image.NumCols());
  RealImage dim;

  InitializeDipole(dipole, image, FieldThres);

  // Localizing edges
  for(int i=0; i<MaxIter; ++i) {
    UpdateDipoleVector(dipole,FieldBeta,FieldSigma,FieldThres,
                       IncDelta,CornerThres,FieldRange);
  }        

#ifdef GL_DISPLAY
  SquareImage square= InitializeSquare(image);
  FieldImage field=InitializeField(dipole);
  hereUpdateField(field,dipole,FieldThres,FieldBeta,FieldSigma,FieldRange);
  DisplayField(argc, argv, image, square, dipole, field);
#else
  dipole=DipoleExpandVector(dipole,ScaleParam);

  //dim=DipoleVectorImage2RealImage(dipole);
  //cout << "Expanded dipole config:" << endl;
  //InteractiveImageWrite(dim,outfile,1);

  FieldImage field=InitializeField(dipole.NumRows(),dipole.NumCols());
  UpdateSeparableField(field,dipole,FieldSigma*ScaleParam,0.001);
  dim=FieldImage2RealImage(field);
  cout << "Interpolation Field:" << endl;
  InteractiveImageWrite(dim,outfile,1);
  
  field=NonMaximumFieldSuppression(field,FieldThres);
  dim=FieldImage2RealImage(field);
  cout << "Final Interpolation Field:" << endl;
  InteractiveImageWrite(dim,outfile,1);

  dim.extractROI(0,0,0,dim.NumRows()-1,0,dim.NumCols()-1);
  dim.WritePnmFile(outfile,IO_PGM,1);
#endif
}
 