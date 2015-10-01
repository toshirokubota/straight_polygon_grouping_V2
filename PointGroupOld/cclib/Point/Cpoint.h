#ifndef _Cpoint_h_
#define _Cpoint_h_

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include <stdlib.h>
#include <assert.h>
#include <mytype.h>
#include <Complex.h>
#include <Misc.h>
#include <Draw.h>

class Cpoint{
public:
	Cpoint(real x=0, real y=0, real a=0, real b=0, real c=0, real d=0) {
		X=x;
		Y=y;
		A=a;
		B=b;
		C=c;
		D=d;
	}
	
	Cpoint(const Cpoint& cp) {
		X=cp.GetX();
		Y=cp.GetY();
		A=cp.GetA();
		B=cp.GetB();
		C=cp.GetC();
		D=cp.GetD();
	}
	
	const Cpoint& operator=(const Cpoint& cp) {
		X=cp.GetX();
		Y=cp.GetY();
		A=cp.GetA();
		B=cp.GetB();
		C=cp.GetC();
		D=cp.GetD();
		
		return *this;
	}
	
	inline real GetX() const {return X;}
	inline real GetY() const {return Y;}
	inline real GetA() const {return A;}
	inline real GetB() const {return B;}
	inline real GetC() const {return C;}
	inline real GetD() const {return D;}
	inline void SetX(real x) {
		X=x;
	}
	inline void SetY(real y) {
		Y=y;
	}
	inline void SetA(real a) {
		A=a;
	}
	inline void SetB(real b) {
		B=b;
	}
	inline void SetC(real c) {
		C=c;
	}
	inline void SetD(real d) {
		D=d;
	}
	
	real EvaluateX(real t) const {
		return X+A*t*t+B*t;
	}
	
	real EvaluateY(real t) const {
		return Y+C*t*t+D*t;
	}
	
	Complex Tangent(real t) const {
		return Complex(2.*A*t+B,2.*C*t+D);
	}
	
	real PseudoCurvature() const {
		return B*C-A*D;
	}
	
	real Curvature(real t) const {
        Complex tt=Tangent(t);
        real numerator=2.*(tt.GetReal()*C-A*tt.GetImag());
        real dd=tt.Power();
        dd=Max(dd,1.0e-30);
        real ee=numerator/dd;
        ee/=sqrt(dd);
        return ee;
	}
	
	real ArcLength(real t1, real t2, int order) const {
		int n=2*order;
		real inc=(t2-t1)/n;
		real len=.0;
		int i;
		for(i=0; i<=n; ++i) {
			real t=t1+i*inc;
			real tx=2.*A*t+B;
			real ty=2.*C*t+D;
			real f=sqrt(tx*tx+ty*ty);
			if(i==0 || i==n)
				len+=f;
			else if(i%2==0)
				len+=2.*f;
			else
				len+=4.*f;
		}
		len*=(inc/3.0);
		return len;
	}
	
	// Artithmetic overloaded operators
	Cpoint operator+(const Cpoint& s) const {
		Cpoint p(X+s.GetX(),Y+s.GetY(),\
			A+s.GetA(),B+s.GetB(),C+s.GetC(),D+s.GetD());
		return p;
	}
	Cpoint operator-(const Cpoint& s) const {
		Cpoint p(X-s.GetX(),Y-s.GetY(),\
			A-s.GetA(),B-s.GetB(),C-s.GetC(),D-s.GetD());
		return p;
	}
	
	void operator+=(const Cpoint& s) {
		X+=s.GetX(); Y+=s.GetY(); 
		A+=s.GetA(); B+=s.GetB(); C+=s.GetC(); D+=s.GetD();
	}
	void operator-=(const Cpoint& s) {
		X-=s.GetX(); Y-=s.GetY(); 
		A-=s.GetA(); B-=s.GetB(); C-=s.GetC(); D-=s.GetD();
	}
	
	//real-complex operators
	Cpoint operator*(real v) const {
		return Cpoint(v*X,v*Y,v*A,v*B,v*C,v*D);
	}
	void operator*=(real v) {
		X*=v; Y*=v; A*=v; B*=v; C*=v; D*=v;
	}
	
 protected:
	 real X;
	 real Y;
	 real A;
	 real B;
	 real C;
	 real D;
};

// Standard Output
ostream& operator<<(ostream& s, const Cpoint& v);
istream& operator>>(istream& s, Cpoint& v);

//
// Interface to Image class so that we can use a variety of File IO routines.
//
typedef Image<Cpoint> CpointImage;
typedef vector<Cpoint> vPoint;

CpointImage
Cpoint2CpointImage(Cpoint* p, int n);

Cpoint*
CpointImage2Cpoint(const CpointImage& pim);

void
OutputVectorCpoints(vector<vPoint>& plink, char* filename);

vector<vPoint>
InputVectorCpoints(char* filename);

//
// Miscellaneous routines
//
Cpoint* CopyPoints(Cpoint* pnt, int num);

real ComputeTurnAngle(const Cpoint& p1, const Cpoint& p2, const Cpoint& p3);

real ComputeArea(const Cpoint& p1, const Cpoint& p2, const Cpoint& p3);

real ComputeDistance(const Cpoint& p1, const Cpoint& p2);

bool IsClockwise(Cpoint* pnts, int num);

Cpoint* ReorientCpoints(Cpoint* pnts, int num);

Cpoint* RotateCpoints(Cpoint* pnts, int num, int sh);

Cpoint* ReflectCpoints(Cpoint* pnts, int num);

Cpoint* NormalizeCpoints(Cpoint* pnts, int num);

real ComputeBestScale(int height, int width, Cpoint* pnts, int nump);

void DelineateContour(ByteImage& flag, Cpoint* pnts, int nump, int band, 
					  real scale, unsigned char val);

void DelineateContour(ByteImage& flag, vPoint& pnts, int nump, int band, 
					  real scale, unsigned char val);

void DelineateApproximateContour(ByteImage& flag, Cpoint* pnts, int nump, 
								 int band1, int band2, real scale, unsigned char val);

Cpoint* MakeNoisyCircle(real rad, int num, real nvar);

Cpoint* MakeNoisyTwoEllipses(real rad, int num, real nvar); 

Cpoint* MakeNoisySquare(real rad, int num, real nvar); 

#endif /* Cpoint_h */
