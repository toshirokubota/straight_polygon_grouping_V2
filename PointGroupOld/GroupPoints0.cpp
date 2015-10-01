/*
This routine implements grouping of moving points.
*/

#include <iostream.h>
//#include <fstream.h>
//using namespace std;

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
#include <Mpoint.h>

int Factor=4;
int Width=256;
int Height=256;
real FrameRate=1.0;
real Rate=0.1;
int NumIter=20;
real NoiseVar=.0; 
int NumPoints=100;
real Epsilon=0.00001;
real Threshold=.5;
int Seed=0;
real Sigma=1.0;

const int NumGroups=3;
const int MaxCandidates=10;

real VelocityY[]={0.0, 2.0, -2.0};
real VelocityX[]={0.0, 2.0, 2.0};

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

void
CopyPoints(Mpoint* dest, Mpoint* src, int nump) {
	for(int i=0; i<nump; ++i) 
		dest[i]=src[i];
}

RealImage
ShowPoints(Mpoint* points, int h, int w, int nump) {
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

RealImage
ShowLinks(Mpoint* points, int h, int w, int nump) {
	RealImage im(1,h,w,0);
	int i,j;
	int length=5;
	for(i=0; i<nump; ++i) {
		real x=points[i].GetX();
		real y=points[i].GetY();
		int ix=Round(x);
		int iy=Round(y);
		for(j=0; j<nump; ++j) {
			if(i!=j) {
				real x2=points[j].GetX();
				real y2=points[j].GetY();
				int ix2=Round(x);
				int iy2=Round(y);
				DrawLine(im,0,iy,ix,iy2,ix2,points[i].GetLink(j));
			}
		}
	}
	return im;
}

void
MovePoints(Mpoint* points, real rate, int nump) {
	
	int i;
	for(i=0; i<nump; ++i) {
		points[i].Move(rate);
	}
}

real
Dist(real y1, real x1, real y2, real x2) {
	return (y1-y2)*(y1-y2)+(x1-x2)*(x1-x2);
}

void
InitializePoints(Mpoint* points, Mpoint* oldpoints, int height, int width, int nump, real rate, real var) {
	
	int i,j,k;
	for(i=0; i<nump; ++i) {
		real y=height*rndm(0);
		real x=width*rndm(0);
		int gid=(int)(rndm(0)*NumGroups);
		points[i].SetNumPoints(nump,1.0);
		points[i].SetY(y);
		points[i].SetX(x);
		points[i].SetVelocityY(VelocityY[gid]);
		points[i].SetVelocityX(VelocityX[gid]);
	}
	CopyPoints(oldpoints,points,nump);
	MovePoints(points,rate,nump);
	
	real* est_vy=new real[MaxCandidates];
	real* est_vx=new real[MaxCandidates];
	for(i=0; i<nump; ++i) {
		real x0=oldpoints[i].GetX();
		real y0=oldpoints[i].GetY();
		for(j=0; j<nump; ++j) {
			real x1=points[j].GetX();
			real y1=points[j].GetY();
			int count=0;
			if(Dist(x0,y0,x1,y1)<9.0 && count<MaxCandidates) {
				est_vy[count]=y1-y0;
				est_vx[count]=x1-x0;
				count++;
			}
			assert(count>0);
			points[i].SetNumCandidates(count);
			for(k=0; k<count; ++k) {
				points[i].SetEstimateVelocityY(k,est_vy[k]);
				points[i].SetEstimateVelocityX(k,est_vx[k]);
				points[i].SetProb(k,1.0/count);
			}
		}
	}
}

real
CompatibilityVelocity(const Mpoint& p1, const Mpoint& p2, int m, int n, real sgm) {
	real vy1=p1.GetEstimateVelocityY(m);
	real vx1=p1.GetEstimateVelocityX(m);
	real vy2=p2.GetEstimateVelocityY(n);
	real vx2=p2.GetEstimateVelocityX(n);
	real dd=Dist(vy1,vx1,vy2,vx2);
	if(dd>sgm)
		return .0;
	else
		return exp(-dd/(2.*sgm*sgm));
}

real
CompatibilityPoints(const Mpoint& point1, const Mpoint& point2, real thres, real sgm) {
	real maxc=.0;
	for(int i=0; i<point1.GetNumCandidates(); ++i) {
		real vy=point1.GetEstimateVelocityY(i);
		real vx=point1.GetEstimateVelocityX(i);
		real p1=point1.GetProb(i);
		if(p1<thres)
			continue;
		for(int j=0; j<point2.GetNumCandidates(); ++j) {
			real vy2=point2.GetEstimateVelocityY(j);
			real vx2=point2.GetEstimateVelocityX(j);
			real dd=Dist(vy,vx,vy2,vx2);
			real p2=point2.GetProb(j);
			if(p2<thres)
				continue;
			real cc=exp(-dd/(2.*sgm*sgm));
			maxc=Max(cc,maxc);
		}
	}
	return maxc;
}


void
UpdateVelocityEstimate(Mpoint* points, Mpoint* prev, int nump, real sigma) {
	real* sum=new real[MaxCandidates];
	int i,m,n;
	for(i=0; i<nump; ++i) {
		real total_sum=0;
		for(m=0; m<points[i].GetNumCandidates(); ++m) {
			sum[m]=.0;
			for(int j=0; j<nump; ++j) {
				if(i!=j) {
					for(n=0; n<points[j].GetNumCandidates(); ++n) {
						real cc=CompatibilityVelocity(prev[i],prev[j],m,n,sigma);
						real pp=prev[i].GetLink(j);
						sum[m]+=pp*cc;
					}
				}
			}
			sum[m]*=prev[i].GetProb(m);
			total_sum+=sum[m];
		}
		for(m=0; m<points[i].GetNumCandidates(); ++m) {
			points[i].SetProb(m,sum[m]/total_sum);
		}
	}
}

void
UpdateLinkWeights(Mpoint* points, Mpoint* prev, int nump, real thres) {
	for(int i=0; i<nump; ++i) {
		for(int j=0; j<nump; ++j) {
			if(i!=j) {
				real cc=CompatibilityPoints(prev[i],prev[j],0.1,1.0);
				real pp=prev[i].GetLink(j);
				real p2=pp*cc/(pp*cc+(1.-pp)*thres);
				points[i].SetLink(j,p2);
			}
		}
	}
}

void main(int argc,char *argv[]) {
	
	if(ReadArguments(argc,argv)!=0) {
		Usage(argv[0]);
		abort();
	}
	else
		PrintParameters();
	
	Mpoint* points=new Mpoint[NumPoints];
	Mpoint* previous=new Mpoint[NumPoints];
	InitializePoints(points,previous,Height,Width,NumPoints,FrameRate,NoiseVar);
	//CopyPoints(orig,points,NumPoints);
	RealImage map=ShowLinks(points,Height,Width,NumPoints);
	InteractiveImageWrite(map,"temp.pgm",1);
	
	for(int i=0; i<NumIter; ++i) {
		CopyPoints(previous,points,NumPoints);
		MovePoints(points,FrameRate,NumPoints);
		UpdateVelocityEstimate(points,previous,NumPoints,Sigma);
		UpdateLinkWeights(points,previous,NumPoints,Threshold);
		map=ShowLinks(points,Height,Width,NumPoints);
		sprintf(outfile,"out%2.2d.pgm",i);
		map.WritePnmFile(outfile,IO_PGM,1);
	}
}