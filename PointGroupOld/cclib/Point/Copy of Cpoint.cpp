#include "Cpoint.h"

// Standard Output
ostream& operator<<(ostream& s, const Cpoint& v) {
    s << v.GetX();
	s << " " << v.GetY();
    s << " " << v.GetA();
    s << " " << v.GetB();
    s << " " << v.GetC();
    s << " " << v.GetD() << endl;
    return s;
}

istream& operator>>(istream& s, Cpoint& v) {
	real val;
	s>>val; v.SetX(val);
	s>>val; v.SetY(val);
	s>>val; v.SetA(val);
	s>>val; v.SetB(val);
	s>>val; v.SetC(val);
	s>>val; v.SetD(val);
    return s;
}

CpointImage
Cpoint2CpointImage(Cpoint* p, int n) {
	CpointImage im(1,1,n);
	int i;
	for(i=0; i<n; ++i)
		im.SetPixel(i,p[i]);
	return im;
}

Cpoint*
CpointImage2Cpoint(const CpointImage& pim) {
	Cpoint* pnt=new Cpoint[pim.NumRows()*pim.NumCols()];
	int i;
	for(i=0; i<pim.NumRows()*pim.NumCols(); ++i)
		pnt[i]=pim.GetPixel(i);
	return pnt;
}

real 
ComputeDistance(const Cpoint& p1, const Cpoint& p2) {
	real dx=p1.GetX()-p2.GetX();
	real dy=p1.GetY()-p2.GetY();
	return sqrt(dx*dx+dy*dy);
}

real
ComputeTurnAngle(const Cpoint& p1, const Cpoint& p2, const Cpoint& p3) {
	real x1=p1.GetX();
	real y1=p1.GetY();
	real x2=p2.GetX();
	real y2=p2.GetY();
	real x3=p3.GetX();
	real y3=p3.GetY();
	real phi1=atan2(y2-y1,x2-x1);
	real phi2=atan2(y3-y2,x3-x2);
	return phi2-phi1;
}

real
ComputeArea(const Cpoint& p1, const Cpoint& p2, const Cpoint& p3) {
	real x1=p1.GetX();
	real y1=p1.GetY();
	real x2=p2.GetX();
	real y2=p2.GetY();
	real x3=p3.GetX();
	real y3=p3.GetY();
	return (x2-x1)*(y3-y2)-(y2-y1)*(x3-x2);
}

bool 
IsClockwise(Cpoint* pnts, int num) {
	int i;
	real sum=.0;
	for(i=0; i<num; ++i) {
		Cpoint p1=pnts[i];
		Cpoint p2=pnts[(i+1)%num];
		sum+=(p1.GetX()*p2.GetY()-p2.GetX()*p1.GetY());
	}
	if(sum>0)
		return false;
	else
		return true;
}

Cpoint*
RotateCpoints(Cpoint* pnts, int num, int sh) {
	Cpoint* res=new Cpoint[num];
	int i;
	for(i=0; i<num; ++i) {
		res[i]=pnts[(i+sh+num)%num];
	}
	return res;
}

Cpoint*
ReflectCpoints(Cpoint* pnts, int num) {
	Cpoint* res=new Cpoint[num];
	int i;
	for(i=0; i<num; ++i) {
		res[i]=pnts[i];
		res[i].SetY(-res[i].GetY());
	}
	return res;
}

Cpoint* 
ReorientCpoints(Cpoint* pnts, int num) {
	Cpoint* res=new Cpoint[num];
	int i;
	for(i=0; i<num; ++i) {
		res[i]=pnts[num-1-i];
	}
	return res;
}

Cpoint* 
NormalizeCpoints(Cpoint* pnts, int num) {
	Cpoint* res=new Cpoint[num];
	int i;
	// Compute the arclength and centroid
	real cy=.0, cx=.0, arclen=.0;
	for(i=0; i<num; ++i) {
		cy+=pnts[i].GetY();
		cx+=pnts[i].GetX();
		arclen+=ComputeDistance(pnts[i],pnts[(i+1)%num]);
	}
	cy/=num;
	cx/=num;
	arclen/=TwoPi;
	//Normalize it
	for(i=0; i<num; ++i) {
		res[i]=pnts[i];
		real y=pnts[i].GetY();
		real x=pnts[i].GetX();
		y=(y-cy)/arclen;
		x=(x-cx)/arclen;
		res[i].SetY(y); 
		res[i].SetX(x);
	}
	return res;
}

real
ComputeBestScale(int height, int width, Cpoint* pnts, int nump) {
	real maxx, minx, maxy, miny;
	int i;
	for(i=0; i<nump; ++i) {
		real y=pnts[i].GetY();
		real x=pnts[i].GetX();
		if(i==0) {
			maxx=minx=x;
			maxy=miny=y;
		}
		else {
			if(y>maxy) maxy=y;
			else if(y<miny) miny=y;
			if(x>maxx) maxx=x;
			else if(x<minx) minx=x;
		}
	}
	real scalex=Min(width/(2.*Abs(minx)),width/(2.*Abs(maxx)));
	real scaley=Min(height/(2.*Abs(miny)),height/(2.*Abs(maxy)));
	real scale=0.975*Min(scalex,scaley);
	
	return scale;
}

/*
Brute force implementation of interpolating two adjcent edge points.
*/
void
DelineateContour(ByteImage& flag, Cpoint* pnts, int nump, int band, real scale) {
	int i;
	for(i=0; i<nump; ++i) {
		real y=pnts[i].GetY();
		real x=pnts[i].GetX();
		real a=pnts[i].GetA();
		real b=pnts[i].GetB();
		real c=pnts[i].GetC();
		real d=pnts[i].GetD();
		int np=15; //# of sample points
		int n;
		real t, inc;
		inc=1.0/(np-1);
		for(n=0, t=-0.5; n<=np; ++n, t+=inc) {
			real xx=(x+a*t*t+b*t);
			real yy=(y+c*t*t+d*t);
			int iy=Round(scale*yy)+flag.NumRows()/2;
			int ix=Round(scale*xx)+flag.NumCols()/2;
			if(iy>=0 && iy<flag.NumRows() && ix>=0 && ix<flag.NumCols() && 
				!flag.GetPixel(band,iy,ix)) {
				flag.SetPixel(band,iy,ix,1);
			}
		}
	}
}

void
DelineateApproximateContour(ByteImage& flag, Cpoint* pnts, int nump, 
							int band1, int band2, real scale) {
	for(int i=0; i<nump; ++i) {
		real y1=pnts[i].GetY();
		real x1=pnts[i].GetX();
		real y2=pnts[(i+1)%nump].GetY();
		real x2=pnts[(i+1)%nump].GetX();
		int iy1=Round(y1*scale)+flag.NumRows()/2;
		int ix1=Round(x1*scale)+flag.NumCols()/2;
		int iy2=Round(y2*scale)+flag.NumRows()/2;
		int ix2=Round(x2*scale)+flag.NumCols()/2;
		DrawLine(flag,band1,iy1,ix1,iy2,ix2,1);
		DrawCross(flag,band2,iy1,ix1,5,1);
	}
}
