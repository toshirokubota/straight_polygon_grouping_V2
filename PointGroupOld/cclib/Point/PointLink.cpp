#include "PointLink.h"

void
CopyPoints(Cpoint* dest, Cpoint* src, int nump) {
	for(int i=0; i<nump; ++i) 
		dest[i]=src[i];
}

void
ShowPoints(RealImage& map, Cpoint* points, int nump, bool closed) {
	int h=map.NumRows();
	int w=map.NumCols();
	int* ix=new int[nump];
	int* iy=new int[nump];
	int i;
	real xmax,xmin,ymax,ymin;
	for(i=0; i<nump; ++i) {
		real x=points[i].GetX();
		real y=points[i].GetY();
		if(i==0) {
			xmax=x; xmin=x;
			ymax=y; ymin=y;
		}
		else {
			if(y>ymax) ymax=y;
			else if(y<ymin) ymin=y;
			if(x>xmax) xmax=x;
			else if(x<xmin) xmin=x;
		}
	}
	real facx=0.9*w/(xmax-xmin);
	real facy=0.9*h/(ymax-ymin);
	real fac=Min(facx,facy);
	for(i=0; i<nump; ++i) {
		real x=points[i].GetX();
		real y=points[i].GetY();
		//ix[i]=Round(fac*(x-xmin)+.05*w);
		//iy[i]=Round(fac*(y-ymin)+.05*h); //reverse it for display
		ix[i]=Round(4*x);
		iy[i]=Round(4*y); //reverse it for display
	}
	for(i=0; i<nump-1; ++i) {
		DrawLine(map,0,iy[i],ix[i],iy[i+1],ix[i+1]);
	}
	if(closed)
		DrawLine(map,0,iy[nump-1],ix[nump-1],iy[0],ix[0]);
	delete [] iy;
	delete [] ix;
}

bool
StartTrack(const RealImage& edge, const ByteImage& flag, int y, int x) {
	if(flag.GetPixel(0,y,x)==0 && //has not been tracked
		edge.GetPixel(0,y,x)>0){ // and edge
		return true;
	}
	else
		return false;
}

bool
FindNextPoint(const RealImage& edge, const ByteImage& flag, 
							int y, int x, int& y2, int& x2) {
	real max=0;
	for(int i=y-1; i<=y+1; ++i) {
		for(int j=x-1; j<=x+1; ++j) {
			if(!(i==y && j==x)) {
				if(flag.GetPixelDefault(0,i,j,1)==0) {
					real theta=edge.GetPixel(1,i,j);
					real cs=cos(theta+HalfPi);
					if(max<edge.GetPixel(0,i,j)*cs*cs) {
						max=edge.GetPixel(0,i,j);
						y2=i;
						x2=j;
					}
				}
			}
		}
	}
	if(max>0)
		return true;
	else
		return false;				
}

bool
FindNextPointNoCanny(const RealImage& edge, const ByteImage& flag, 
							int y, int x, int& y2, int& x2) {
	real max=0;
	int done=0;
	for(int i=y-1; !done && i<=y+1; ++i) {
		for(int j=x-1; !done && j<=x+1; ++j) {
			if(!(i==y && j==x)) {
				if(flag.GetPixelDefault(0,i,j,1)==0) {
					if(max<edge.GetPixel(0,i,j)) {
						max=edge.GetPixel(0,i,j);
						y2=i;
						x2=j;
						done=1;
					}
				}
			}
		}
	}
	if(max>0)
		return true;
	else
		return false;				
}

bool
FindNextPointOld(const RealImage& edge, const ByteImage& flag, 
							int y, int x, int& y2, int& x2) {
	real max=0;
	for(int i=y-1; i<=y+1; ++i) {
		for(int j=x-1; j<=x+1; ++j) {
			if(!(i==y && j==x)) {
				if(flag.GetPixelDefault(0,i,j,1)==0) {
					if(max<edge.GetPixel(0,i,j)) {
						max=edge.GetPixel(0,i,j);
						y2=i;
						x2=j;
					}
				}
			}
		}
	}
	if(max>0)
		return true;
	else
		return false;				
}

Cpoint*
Combine(Cpoint* p1, Cpoint* p2, int c1, int c2) {
	Cpoint* res=new Cpoint[c1+c2];
	CopyPoints(res,p1,c1);
	CopyPoints(res+c1,p2,c2);
	return res;
}

Cpoint*
Track(const RealImage& edge, ByteImage& flag, 
			int y, int x, int& cnt) {
	flag.SetPixel(0,y,x,1);
	int y2,x2;
	Cpoint p2[1];
	p2[0]=Cpoint((real)x,(real)y,0,0,0,0);
	int n=0;
	Cpoint* pnt;
	if(FindNextPoint(edge,flag,y,x,y2,x2)) {
		pnt=Track(edge,flag,y2,x2,n);
	}
	cnt=n+1;
	Cpoint* res=Combine(p2,pnt,1,n);
	if(n>0)
  		delete [] pnt;
	return res;
}

Cpoint*
TrackNoCanny(const RealImage& edge, ByteImage& flag, 
			int y, int x, int& cnt) {
	flag.SetPixel(0,y,x,1);
	int y2,x2;
	Cpoint p2[1];
	p2[0]=Cpoint((real)x,(real)y,0,0,0,0);
	int n=0;
	Cpoint* pnt;
	if(FindNextPointNoCanny(edge,flag,y,x,y2,x2)) {
		pnt=TrackNoCanny(edge,flag,y2,x2,n);
	}
	cnt=n+1;
	Cpoint* res=Combine(p2,pnt,1,n);
	if(n>0)
  		delete [] pnt;
	return res;
}

void
WritePointLinks(struct_link& seg, char* filename) {		
  fstream out;
  //out.open (filename, ios::ascii| ios::out);
  out.open (filename, ios::out);
  vector<int>::iterator np;
  vector<Cpoint>::iterator lp;
  out << seg.num.size() << endl;
  for(np=seg.num.begin(); np<seg.num.end(); ++np)
  	out << *np << " ";
  out << endl;
  for(lp=seg.link.begin(); lp<seg.link.end(); ++lp)
		out << *lp;
	out.close();
}

struct_link
ReadPointLinks(char* filename) {		
  fstream in;
  struct_link seg;
  in.open (filename, ios::in);
  int num_seg, i, j, m;
  in >> num_seg;
  //in >> tmp_char;
  for(i=0; i<num_seg; ++i) {
  	in >> m;
  	seg.num.push_back(m);
  }
  //in >> tmp_char;
  for(i=0; i<num_seg; ++i) {
  	for(j=0; j<seg.num[i]; ++j) {
  		Cpoint pnt;
  		in>>pnt;
  		seg.link.push_back(pnt);
  	}
  }
	in.close();
	return seg;
}