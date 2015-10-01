/*
This routine implements grouping of moving points.
*/

#include <iostream>
#include <fstream>
using namespace std;

#include <mytype.h>
#include <ImageProc.h>
#include <Feature.h>
#include <Convolve.h>
#include <Gradient.h>
#include <Misc.h>

#include <Draw.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "MCpoint.h"

int Factor=4;
int Width=256;
int Height=256;
real FrameRate=1.0;
real Rate=0.25;
int NumIter=20;
real NoiseVar=1.0;
int NumGroup=3;
int PointsPerGroup=6;
int NoisePoints=20;
int NumPoints=NumGroup*PointsPerGroup+NoisePoints;
real Epsilon=0.00001;
real Threshold1=.01;
real Threshold2=.01;
int Seed=54;
real Sigma1=sqrt(2.0)*NoiseVar;
real Sigma2=1.0;
real Sigma3=1.0;

const int MaxNumGroups=10;
//const int MaxCandidates=10;
real MaxRange=25.0;

real VelocityY[]={20.0, 20.0, -20.0};
real VelocityX[]={-20.0, 20.0, 20.0};

char infile[256];
char outfile[256];

//
// Usage - argument parsing routines
//
void
Usage(char* program) {
	cerr << program << ": usage" << endl;
	cerr << "Required arguments:" << endl;
	//cerr << "\t -i <input file name>" << endl;
	cerr << "\t -o <output file name>" << endl;
	
	cerr << "Optional arguments:" << endl;
	cerr << "\t -m (int) set the # of iterations" << endl;
	cerr << "\t \t default value is " << NumIter << endl;
	cerr << "\t -rt (int) set the update rate (>0)" << endl;
	cerr << "\t \t default value is " << Rate << endl;
	cerr << "\t -nv (int) set the noise variance" << endl;
	cerr << "\t \t default value is " << NoiseVar << endl;
	cerr << "\t -np (int) set the # of sample points" << endl;
	cerr << "\t \t default value is " << NumPoints << endl;
	cerr << "\t -seed (int) set the seed value for random number generator" << endl;
	cerr << "\t \t default value is " << Seed << endl;
}

void
PrintParameters() {
	cout << "input file: " << infile << endl;
	cout << "output file: " << outfile << endl;
	cout << "# of iterations: " << NumIter << endl;
	cout << "the update rate: " << Rate << endl;
	cout << "the noise variance: " << NoiseVar << endl;
	cout << "# of sample points: " << NumPoints << endl;
	cout << "Seed number: " << Seed << endl;
}


int
ReadArguments(int argc, char* argv[]) {
	int in_found=0;
	int out_found=0;
	real tempval;
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
				NumIter=atoi(argv[++i]);
			}
			else
				return 6; // no value specified
		}
		else if(strcmp(argv[i],"-seed")==0) {
			if(argc>i+1) {
				Seed=atoi(argv[++i]);
			}
			else
				return 6; // no value specified
		}
		else if(strcmp(argv[i],"-np")==0) {
			if(argc>i+1) {
				NumPoints=atoi(argv[++i]);
			}
			else
				return 6; // no value specified
		}
		else if(strcmp(argv[i],"-rt")==0) {
			if(argc>i+1) {
				tempval=strtod(argv[++i],NULL);
				if(tempval<0)
					return 8;
				else
					Rate=tempval;
			}
			else
				return 6; // no value specified
		}
		else if(strcmp(argv[i],"-nv")==0) {
			if(argc>i+1) {
				NoiseVar=strtod(argv[++i],NULL);
				//MaxRange=NoiseVar*NoiseVar;
				Sigma1=sqrt(2.)*NoiseVar;
			}
			else
				return 6; // no value specified
		}
	}
	
	//if(!in_found)
	//  return 11;
	rndm(Seed);
	if(!out_found)
		return 12;
	else
		return 0;
}

RealImage
ShowPoints(const vMCpoint& points, int h, int w, int nump) {
	RealImage im(1,h,w,0);
	int i;
	int length=5;
	for(i=0; i<nump; ++i) {
		real x=points[i].GetX();
		real y=points[i].GetY();
		int ix=Round(x);
		int iy=Round(y);
		if(iy>=0 && iy<h && ix>=0 && ix<w)
			im.SetPixel(0,iy,ix,1.0);
		//DrawCross(im,0,iy,ix,0);
	}
	return im;
}

void
TraceLink(const vMCpoint& points, ByteImage& flag, RealImage& im, 
          int i, real color) {
   int nump=points.size();
	for(int j=0; j<nump; ++j) {
		if(flag.GetPixel(0,i,j))
			continue;
		for(int m=0; m<points[i].GetNumCandidates(); ++m) {
			for(int n=0; n<points[j].GetNumCandidates(); ++n) {
				if(points[i].GetLink(j,m,n)>=.5) {
					flag.SetPixel(0,i,j,1);
					flag.SetPixel(0,j,i,1);
					int y1=Round(points[i].GetY());
					int x1=Round(points[i].GetX());
					int y2=Round(points[j].GetY());
					int x2=Round(points[j].GetX());
					DrawLine(im,0,y1,x1,y2,x2,color);
					TraceLink(points,flag,im,j,color);
				}
			}
		}
		flag.SetPixel(0,i,j,1);
		flag.SetPixel(0,j,i,1);
	}
}

RealImage
ShowLinks(const vMCpoint points, int h, int w) {
   int nump=points.size();
	RealImage im(1,h,w,0);
	ByteImage flag(1,nump,nump,0);
	int i;
	for(i=0; i<nump; ++i)
		flag.SetPixel(0,i,i,1);
	real color=1.0;
	for(i=0; i<nump; ++i) {
			TraceLink(points,flag,im,i,color);
			color+=1.0;
	}
	return im;
}

RealImage
ShowOnlyPoints(const vMCpoint& points, int h, int w) {
   int nump=points.size();
	RealImage im(1,h,w,0);
	int i;
	real color=1.0;
	for(i=0; i<nump; ++i) {
		int y1=Round(points[i].GetY());
		int x1=Round(points[i].GetX());
		if(y1>=0 && y1<h && x1>=0 && x1<w)
			DrawCross(im,0,y1,x1,5);
	}
	return im;
}

RealImage
ShowGroundTruthLinks(const vMCpoint& points, int h, int w) {
   int nump=points.size();
	RealImage im(NumGroup,h,w,0);
	int i,j,j2;
	for(i=0; i<NumGroup; ++i) {
		for(j=0; j<PointsPerGroup; ++j) {
			real y1=points[i*PointsPerGroup+j].GetY();
			real x1=points[i*PointsPerGroup+j].GetX();
			for(j2=0; j2<PointsPerGroup; ++j2) {
				if(j!=j2) {
					real y2=points[i*PointsPerGroup+j2].GetY();
					real x2=points[i*PointsPerGroup+j2].GetX();
					DrawLine(im,i,Round(y1),Round(x1),Round(y2),Round(x2),1);
				}
			}
		}
	}
	for(i=NumGroup*PointsPerGroup; i<NumPoints; ++i) {
		real y1=points[i].GetY();
		real x1=points[i].GetX();
		DrawCross(im,0,Round(y1),Round(x1),5);
		DrawCross(im,1,Round(y1),Round(x1),5);
		DrawCross(im,2,Round(y1),Round(x1),5);
	}
	return im;
}

vMCpoint
MovePoints(const vMCpoint& points, real rate, real var) {
   int nump=points.size();
	vMCpoint res=points;
	
	int i;
	for(i=0; i<nump; ++i) { 
		res[i].Move(rate,var);
	}
	return res;
}

vMCpoint
ShufflePoints(const vMCpoint& points) {
   int nump=points.size();
	vMCpoint res(nump);
	
	int i;
	for(i=0; i<nump; ++i) { 
		res[nump-i-1]=points[i];
	}
	return res;
}

real
Dist(real y1, real x1, real y2, real x2) {
	return (y1-y2)*(y1-y2)+(x1-x2)*(x1-x2);
}

void
PlacePoints(vMCpoint& points, int height, int width) {
   int nump=points.size();
	int i,j,k,m;
	int offy=30;
	int offx=30;
	for(i=0,m=0; i<NumGroup; ++i) {
		for(j=0; j<PointsPerGroup; ++j,m++) {
			real y;
			real x;
			points[m].SetNumPoints(nump);
			y=(height-2*offy)*rndm(0)+offy;
			x=(width-2*offx)*rndm(0)+offx;
			points[m].SetY(y);
			points[m].SetX(x);
			points[m].SetVelocityY(VelocityY[i]);
			points[m].SetVelocityX(VelocityX[i]);
		}
	}
	for(i=NumGroup*PointsPerGroup; i<NumPoints; ++i){
		points[i].SetNumPoints(nump);
		points[i].SetY(rndm(0)*height);
		points[i].SetX(rndm(0)*width);
		points[i].SetVelocityY(rndm(0)*50);
		points[i].SetVelocityX(rndm(0)*50);
	}
}

/*
Sort y, x and index according to the dist using bubble sort.
*/
void
SortCandidates(real* dist, real* y, real* x, int* index, int count) {
	for(int i=0; i<count; ++i) {
		real max=dist[i];
		int maxi=i;
		for(int j=i; j<count; ++j) {
			if(max<dist[j]) {
				maxi=j;
				max=dist[i];
			}
		}
		if(maxi!=i) {
			real temp=dist[i];
			dist[i]=dist[maxi];
			dist[maxi]=temp;
			temp=y[i];
			y[i]=y[maxi];
			y[maxi]=temp;
			temp=x[i];
			x[i]=x[maxi];
			x[maxi]=temp;
			int tempi=index[i];
			index[i]=index[maxi];
			index[maxi]=tempi;
		}
	}
}

const int MaxNumPoints=256;

void
EstimateTransform(vMCpoint& p1, const vMCpoint& p2) {
	int nump=p1.size();
	real dist[MaxNumPoints];
	real est_vy[MaxNumPoints];
	real est_vx[MaxNumPoints];
	int index[MaxNumPoints];
	int i,j,k;
	for(i=0; i<nump; ++i) {
		real x0=p1[i].GetX();
		real y0=p1[i].GetY();
		int count=0;
		for(j=0; j<nump; ++j) {
			real x1=p2[j].GetX();
			real y1=p2[j].GetY();
			real dd=(x0-x1)*(x0-x1)+(y0-y1)*(y0-y1);
			if(dd<2.*MaxRange*MaxRange) {
				est_vy[count]=y1-y0;
				est_vx[count]=x1-x0;
				index[count]=j;
				count++;
			}
		}
		if(count==0) 
			cout << "Count=0 @ " << i << endl;
		cout << i << ":" << count << endl;
		assert(count>0);
		SortCandidates(dist,est_vy,est_vx,index,count);
		count=Min(count,MCpointMaxCandidates);
		p1[i].SetNumCandidates(count);
		for(k=0; k<count; ++k) {
			p1[i].SetEstimateVelocityY(k,est_vy[k]);
			p1[i].SetEstimateVelocityX(k,est_vx[k]);
			p1[i].SetIndex(k,index[k]);
			p1[i].SetProb(k,1.0/(count+1));
		}
		p1[i].SetProb(count,1.0/(count+1));
	}
   cout << "Done estimate transform.\n";
}

real
CompatibilityTransform(const MCpoint& p1, const MCpoint& p2, int m, int n, real sgm) {
	real vy1=p1.GetEstimateVelocityY(m);
	real vx1=p1.GetEstimateVelocityX(m);
	real vy2=p2.GetEstimateVelocityY(n);
	real vx2=p2.GetEstimateVelocityX(n);
	real dd=Dist(vy1,vx1,vy2,vx2);
	if(dd>2.*sgm*sgm)
		return .0;
	else
		return exp(-dd/(2.*sgm*sgm));
}

real
InterframeCompatibilityTransform(const MCpoint& p1, const MCpoint& p2, int m, int n, real sgm){
	real vy1=p1.GetEstimateVelocityY(m);
	real vx1=p1.GetEstimateVelocityX(m);
	real vy2=-p2.GetEstimateVelocityY(n);
	real vx2=-p2.GetEstimateVelocityX(n);
	real dd=Dist(vy1,vx1,vy2,vx2);
	if(dd>2.*sgm*sgm)
		return .0;
	else
		return exp(-dd/(2.*sgm*sgm));
}

vMCpoint
UpdateTransformEstimate(const vMCpoint& points, real sigma, real rate) {
   int nump=points.size();
	vMCpoint res=points;
	real sum[MCpointMaxCandidates];
	int i,m,n;
	real nvy=0;
	real nvx=0;
	real alpha=0.5;
	for(i=0; i<nump; ++i) {
		real total_sum=0;
		MCpoint pnt1=points[i];
		for(m=0; m<pnt1.GetNumCandidates(); ++m) {
			sum[m]=.0;
			real evy1=pnt1.GetEstimateVelocityY(m);
			real evx1=pnt1.GetEstimateVelocityX(m);
			for(int j=0; j<nump; ++j) {
				if(i!=j) {
					MCpoint pnt2=points[j];
					for(n=0; n<pnt2.GetNumCandidates(); ++n) {
						real evy2=pnt2.GetEstimateVelocityY(n);
						real evx2=pnt2.GetEstimateVelocityX(n);
						real lnk1=pnt1.GetLink(j,m,n);
						real lnk2=pnt2.GetLink(i,n,m);
						real dd=CompatibilityTransform(pnt1,pnt2,m,n,sigma);
						real dy=evy2-evy1;
						real dx=evx2-evx1;
						real pp=pnt2.GetProb(n);
						nvy+=dd*pp*lnk1*lnk2*dy;
						nvx+=dd*pp*lnk1*lnk2*dx;
					}
				}
			}
			res[i].SetEstimateVelocityY(m,evy1+rate*nvy);
			res[i].SetEstimateVelocityX(m,evx1+rate*nvx);
		}
	}
	return res;
}

vMCpoint
UpdateLinkWeights(const vMCpoint& points, real thres, real sigma) {
   int nump=points.size();
	vMCpoint res=points;
	for(int i=0; i<nump; ++i) {
		MCpoint p1=points[i];
		int numc1=p1.GetNumCandidates();
		for(int m=0; m<numc1; ++m) {
			for(int j=0; j<nump; ++j) {
				MCpoint p2=points[j];
				int numc2=p2.GetNumCandidates();
				for(int n=0; n<numc2; ++n) {
					if(i!=j) {
						real dd=CompatibilityTransform(p1,p2,m,n,sigma);
						real lnk1=p1.GetLink(j,m,n);
						real lnk2=p2.GetLink(i,n,m);
						real pp=p2.GetProb(n);
#ifdef BAUM_EAGON
						//real ee=pp*lnk1*lnk2;
						//real pp2=ee*dd/(ee*dd+(1.-ee)*thres);
#else
						real ee=pp*dd*lnk2; // for EM
						real pp2=ee/(ee+thres); //for EM
#endif
						res[i].SetLink(j,m,n,pp2);
					}
				}
			}
		}
	}
	return res;
}

vMCpoint
UpdateProbMeasure(const vMCpoint& frame1, const vMCpoint& frame2, 
                  real thres, real sigma) {
   int nump=frame1.size();
	vMCpoint res=frame1;
	real sum[MCpointMaxCandidates];
	for(int i=0; i<nump; ++i) {
		MCpoint p1=frame1[i];
		int numc1=p1.GetNumCandidates();
		real total_sum=0;
      for(int m=0; m<numc1; ++m) { //compute fitness for each candidate
			sum[m]=.0;
			MCpoint p3=frame2[p1.GetIndex(m)];
			real match_prob=.0;
         //check if the candidate also consider it as a candidate
			for(int m2=0; m2<p3.GetNumCandidates(); ++m2) {
				if(p3.GetIndex(m2)==i) {
					match_prob=p3.GetProb(m2);
					break;
				}
			}
         if(match_prob>0) { //mutual candidancy
				for(int j=0; j<nump; ++j) {
               //compute the support from other INTRA-frame points
               //measured by 1. transformation similarity
               //            2. strength of links
					if(i!=j) {
						MCpoint p2=frame1[j];
						int numc2=p2.GetNumCandidates();
						for(int n=0; n<numc2; ++n) {
							real dd=CompatibilityTransform(p1,p2,m,n,sigma);
							real lnk1=p1.GetLink(j,m,n);
							real lnk2=p2.GetLink(i,n,m);
							real pp=p2.GetProb(n);
							real ss=pp*lnk1*lnk2;
							sum[m]+=ss*match_prob;
						}
					}
				}
#ifdef BAUM_EAGON
				sum[m]*=p1.GetProb(m); //activate this for Baum-Eagon
#endif
				total_sum+=sum[m];
			}
		}
#ifdef BAUM_EAGON
		for(m=0; m<numc1; ++m)
			res[i].SetProb(m,sum[m]/(total_sum+p1.GetProb(numc1)*thres));
		res[i].SetProb(numc1,p1.GetProb(numc1)*thres/(total_sum+p1.GetProb(numc1)*thres));
#else
		for(m=0; m<numc1; ++m)
			res[i].SetProb(m,sum[m]/(total_sum+thres));
		res[i].SetProb(numc1,thres/(total_sum+thres));
#endif
	}
	return res;
}

void main(int argc,char *argv[]) {
	
	if(ReadArguments(argc,argv)!=0) {
		Usage(argv[0]);
		abort();
	}
	else
		PrintParameters();
	
	vMCpoint frame1(NumPoints);
	vMCpoint frame2;
	vMCpoint tmp_frm;
	PlacePoints(frame1,Height,Width);
	frame2=MovePoints(frame1,FrameRate,NoiseVar);
	frame2=ShufflePoints(frame2); //to make sure the array index is not helping
	EstimateTransform(frame1,frame2);
	EstimateTransform(frame2,frame1);
	RealImage map;
	
   map=ShowGroundTruthLinks(frame1,Height,Width);
	map.WritePnmFile("truth.ppm",IO_PPM,1);
	map=ShowGroundTruthLinks(frame2,Height,Width);
	map.WritePnmFile("truth2.ppm",IO_PPM,1);
	map=ShowOnlyPoints(frame1,Height,Width);
	map.WritePnmFile("points1.pgm",IO_PGM,1);
	map=ShowOnlyPoints(frame2,Height,Width);
	map.WritePnmFile("points2.pgm",IO_PGM,1);
	//InteractiveImageWrite(map,"truth.pgm",1);
	
	for(int i=0; i<NumIter; ++i) {	
		cout << i << endl;
		frame1=UpdateLinkWeights(frame1,Threshold1,Sigma1);
		frame2=UpdateLinkWeights(frame2,Threshold1,Sigma1);
		frame1=UpdateTransformEstimate(frame1,Sigma1,Rate);
		frame2=UpdateTransformEstimate(frame2,Sigma1,Rate);
		tmp_frm=UpdateProbMeasure(frame1,frame2,Threshold2,Sigma1);
		frame2=UpdateProbMeasure(frame2,tmp_frm,Threshold2,Sigma1);
		frame1=tmp_frm;
		map=ShowLinks(frame1,Height,Width);
		//sprintf(outfile,"out%2.2d.pgm",i);
		map.WritePnmFile(outfile,IO_PGM,1);
	}

   for(i=0; i<frame1.size(); ++i)
      cout << "Point #" << i << endl << frame1[i];
}