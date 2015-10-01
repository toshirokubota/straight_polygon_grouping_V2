#include "EdgeProc.h"

real
EdgeDistance(const Edge& e1, const Edge& e2) {
	real x1=e1.GetX();
	real y1=e1.GetY();
	real x2=e2.GetX();
	real y2=e2.GetY();
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//#define NOT_WTA

EdgeImage
ProbabilityProjection(const EdgeImage& hyper) {

  EdgeImage res=hyper;
  for(int i=0; i<res.NumRows(); ++i) {
    for(int j=0; j<res.NumCols(); ++j) {
      real sum=.0;
      for(int k=0; k<res.NumBands(); ++k) {
        Edge ei=res.GetPixel(k,i,j);
#ifdef NOT_WTA
        sum+=ei.GetMagnitude();
#else
        sum+=ei.GetMagnitude()*ei.GetProb();
#endif
      }
      if(sum>0) {
        for(int k=0; k<res.NumBands(); ++k) {
          Edge ei=res.GetPixel(k,i,j);
#ifdef NOT_WTA
          real val=ei.GetMagnitude();
#else
          real val=ei.GetMagnitude()*ei.GetProb();
#endif
          ei.SetProb(val/sum);
          res.SetPixel(k,i,j,ei);
        }
      }
      else {
        for(int k=0; k<res.NumBands(); ++k) {
          Edge ei=res.GetPixel(k,i,j);
          ei.SetProb(1.0/hyper.NumBands());
          res.SetPixel(k,i,j,ei);
        }
      }
    }
  }
  return res;
}

EdgeImage
newProbabilityProjection(EdgeImage& hyper) {

  EdgeImage res=hyper;
  for(int i=0; i<hyper.NumRows(); ++i) {
    for(int j=0; j<hyper.NumCols(); ++j) {
      real max=.0;
      int maxk=-1;
      for(int k=0; k<hyper.NumBands(); ++k) {
        Edge ei=hyper.GetPixel(k,i,j);
        if(max<ei.GetMagnitude()) {
        	max=ei.GetMagnitude();
        	maxk=k;
        }
      }
      if(maxk>=0) {
        for(int k=0; k<hyper.NumBands(); ++k) {
          Edge ei=hyper.GetPixel(k,i,j);
          if(k==maxk)
	          ei.SetProb(1.0);
	        else
	        	ei.SetProb(0.0);
          res.SetPixel(k,i,j,ei);
        }
      }
      else {
        for(int k=0; k<hyper.NumBands(); ++k) {
          Edge ei=hyper.GetPixel(k,i,j);
          ei.SetProb(1.0/hyper.NumBands());
          res.SetPixel(k,i,j,ei);
        }
      }
    }
  }
  return res;
}

ByteImage
ExpandEdges(const EdgeImage& hyper, int factor) {
  ByteImage res(1,factor*hyper.NumRows(),factor*hyper.NumCols(),0);
  
  for(int i=0; i<hyper.NumRows(); ++i) {
    for(int j=0; j<hyper.NumCols(); ++j) {
    	real yy=0;
    	real xx=0;
    	bool is_edge=false;
      for(int k=0; k<hyper.NumBands(); ++k) {
      	Edge ed=hyper.GetPixel(k,i,j);
	      if(ed.IsEdge()) {
	      	is_edge=true;
	      }
        yy+=ed.GetY()*factor*ed.GetProb();
        xx+=ed.GetX()*factor*ed.GetProb();
      }
      if(is_edge) {
      	int row=Round(yy);
      	int col=Round(xx);
        if(row>=0 && row<res.NumRows() && col>=0 && col<res.NumCols()) {        
          res.AddPixel(0,row,col,1);
        }
      }
    }
  }
  
  return res;
}

EdgeImage
SetEdges(const EdgeImage& hyper, const RealImage& grd, real thres) {
  EdgeImage res=hyper;
 
  for(int k=0; k<res.NumBands(); ++k) {
    for(int i=0; i<res.NumRows(); ++i) {
      for(int j=0; j<res.NumCols(); ++j) {
        real mag=grd.GetPixel(0,i,j);
        Edge ed=res.GetPixel(k,i,j);
        if(mag>thres) {
          ed.SetFlag(true);////
        }
        else {
          ed.SetFlag(false);////
        }
        res.SetPixel(k,i,j,ed);
      }
    }
  }
  return res;
}

/*
the routine check if multiple edges reside in a small area (specified
by argument 'dist', and remove them except the one with the largest
magnitude.
The purpose of this routine is to simplify the edge grouping process.
*/
EdgeImage 
RearrangeEdgeImage(const EdgeImage& edge, int range, real dist) {
	EdgeImage res=edge;
  for(int i=0; i<edge.NumRows(); ++i) {
    for(int j=0; j<edge.NumCols(); ++j) {
    	real mag=.0;
    	Edge e1,e2,e3;
    	int maxk;
    	for(int k=0; k<edge.NumBands(); ++k) {
		    e1=edge.GetPixel(k,i,j);
		    if(mag<e1.GetMagnitude()) {
		    	mag=e1.GetMagnitude();
		    	maxk=k;
		    }
		  }
		  if(mag>0) {
			  e1=res.GetPixel(maxk,i,j);
				bool done=false;
				for(int i2=i-range; !done && i2<=i+range; ++i2) {
					if(i2>=0 && i2<edge.NumRows()) {
						for(int j2=j-range; !done && j2<=j+range; ++j2) {
							if(j2>=0 && j2<edge.NumCols() && (i!=i2 || j!=j2)) {
					    	for(int k2=0; !done && k2<edge.NumBands(); ++k2) {
							    e2=res.GetPixel(k2,i2,j2);
									if(EdgeDistance(e1,e2)<dist) {
								    if(mag<=e2.GetMagnitude()) {
											for(int k3=0; k3<edge.NumBands(); ++k3) {
												e3=res.GetPixel(k3,i,j);
												e3.SetMagnitude(.0);
												e3.SetProb(.0);
												e3.SetFlag(false);
												res.SetPixel(k3,i,j,e3);
												done=true;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return res;
}

void
InteractiveEdgeImageWrite(const EdgeImage& edge, 
                          char* filename, int normalize) {
  RealImage image(1,edge.NumBands()*edge.NumRows(),edge.NumCols());
  int nump=image.NumBands()*image.NumRows()*image.NumCols();
  int choice;
  int i,j,k,n;
  for(int done=0; !done; ) {
    cout << "Which info to view:\n";
    cout << "(1): X coordinate\n";
    cout << "(2): Y coordinate\n";
    cout << "(3): Theta\n";
    cout << "(4): Magnitude\n";
    cout << "(5): Probability\n";
    cout << "(6): Flag\n";
    cin >> choice;
    switch (choice) {
    case 1:
      n=0;
      for(k=0; k<edge.NumBands(); ++k) {
        for(i=0; i<edge.NumRows(); ++i,++n) {
          for(j=0; j<edge.NumCols(); ++j) {
            Edge ed=edge.GetPixel(k,i,j);
            image.SetPixel(0,n,j,ed.GetX()-j);
          }
        }
      }
      image.WritePnmFile(filename,IO_PGM,normalize);
      break;
    case 2:
      n=0;
      for(k=0; k<edge.NumBands(); ++k) {
        for(i=0; i<edge.NumRows(); ++i,++n) {
          for(j=0; j<edge.NumCols(); ++j) {
            Edge ed=edge.GetPixel(k,i,j);
            image.SetPixel(0,n,j,ed.GetY()-i);
          }
        }
      }
      image.WritePnmFile(filename,IO_PGM,normalize);
      break;
    case 3:
      for(i=0; i<nump; ++i) {
        Edge ed=edge.GetPixel(i);
        image.SetPixel(i,ed.GetTheta());
      }
      image.WritePnmFile(filename,IO_PGM,normalize);
      break;
    case 4:
      for(i=0; i<nump; ++i) {
        Edge ed=edge.GetPixel(i);
        image.SetPixel(i,ed.GetMagnitude());
      }
      image.WritePnmFile(filename,IO_PGM,normalize);
      break;
    case 5:
      for(i=0; i<nump; ++i) {
        Edge ed=edge.GetPixel(i);
        image.SetPixel(i,ed.GetProb());
      }
      image.WritePnmFile(filename,IO_PGM,normalize);
      break;
    case 6:
      for(i=0; i<nump; ++i) {
        Edge ed=edge.GetPixel(i);
        image.SetPixel(i,(real)ed.IsEdge());
      }
      image.WritePnmFile(filename,IO_PGM,normalize);
      break;
    default:
      if(choice<0)
        done=1;
      break;
    }
  }
}
