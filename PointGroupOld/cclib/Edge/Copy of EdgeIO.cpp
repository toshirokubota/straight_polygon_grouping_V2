#include "EdgeIO.h"

EdgeStructImage
EdgeImage2EdgeStructImage(const EdgeImage& hyper) {
	EdgeStructImage edge(hyper.NumBands(),hyper.NumRows(),hyper.NumCols());
	for(int k=0; k<hyper.NumBands(); ++k) {
	  for(int i=0; i<hyper.NumRows(); ++i) {
	    for(int j=0; j<hyper.NumCols(); ++j) {
	    	Edge ed=hyper.GetPixel(k,i,j);
	    	EdgeStruct ee=edge.GetPixel(k,i,j);
	    	ee.X=ed.GetX();
	    	ee.Y=ed.GetY(); //hyper.NumRows()-1-ed.GetY();
	    	ee.Theta=ed.GetTheta();
	    	ee.Prob=ed.GetProb();
	    	ee.Mag=ed.GetMagnitude();
	    	ee.On=ed.IsEdge()? 1: 0;
	    	edge.SetPixel(k,i,j,ee);
	    }
	  }
	}
	return edge;
}

EdgeImage
EdgeStructImage2EdgeImage(const EdgeStructImage& sedge) {
	EdgeImage edge(sedge.NumBands(),sedge.NumRows(),sedge.NumCols());
	for(int k=0; k<edge.NumBands(); ++k) {
		for(int i=0; i<edge.NumRows(); ++i) {
			for(int j=0; j<edge.NumCols(); ++j) {
				Edge ed=edge.GetPixel(k,i,j);
				EdgeStruct ee=sedge.GetPixel(k,i,j);
				ed.SetX(ee.X);
				ed.SetY(ee.Y);
				ed.SetTheta(ee.Theta);
				ed.SetMagnitude(ee.Mag);
				ed.SetProb(ee.Prob);
				ed.SetFlag(ee.Flag);
				edge.SetPixel(k,i,j,ed);
			}
		}
	}
	return edge;
}

/*EdgeStructImage
ReadEdgeImageRaw(char* filename) {
	printf("Reading data from %s\n", filename);
	int ne=6;
	RealImage buf;
	buf.ReadRawFile(filename);
	EdgeStructImage edge(buf.NumBands()/ne,buf.NumRows(),buf.NumCols());
	for(int k=0; k<edge.NumRows(); ++k) {
		for(int i=0; i<edge.NumRows(); ++i) {
			for(int j=0; j<edge.NumCols(); ++j) {
				Edge ed=edge.GetPixel(k,i,j);
				ed.X=buf.GetPixel(ne*k,i,j);
				ed.Y=buf.GetPixel(ne*k+1,i,j);
				ed.Theta=buf.GetPixel(ne*k+2,i,j);
				ed.Mag=buf.GetPixel(ne*k+3,i,j);
				ed.Prob=buf.GetPixel(ne*k+4,i,j);
				ed.On=(unsigned char)buf.GetPixel(ne*k+5,i,j);
				edge.SetPixel(k,i,j,ed);
			}
		}
	}
	return edge;
}*/

EdgeImage
ReadEdgeImageRaw(char* filename) {
	printf("Reading data from %s\n", filename);
	int ne=6;
	RealImage buf;
	buf.ReadRawFile(filename);
	EdgeImage edge(buf.NumBands()/ne,buf.NumRows(),buf.NumCols());
	for(int k=0; k<edge.NumBands(); ++k) {
		for(int i=0; i<edge.NumRows(); ++i) {
			for(int j=0; j<edge.NumCols(); ++j) {
				Edge ed=edge.GetPixel(k,i,j);
				ed.SetX(buf.GetPixel(ne*k,i,j));
				ed.SetY(buf.GetPixel(ne*k+1,i,j));
				ed.SetTheta(buf.GetPixel(ne*k+2,i,j));
				ed.SetMagnitude(buf.GetPixel(ne*k+3,i,j));
				ed.SetProb(buf.GetPixel(ne*k+4,i,j));
				ed.SetFlag((unsigned char)buf.GetPixel(ne*k+5,i,j));
				edge.SetPixel(k,i,j,ed);
			}
		}
	}
	//EdgeImage edge;
	//edge.ReadRawFile(filename);
	return edge;
}

WriteEdgeRaw(const EdgeImage& hyper, char* filename) {

	EdgeStructImage edge(hyper.NumBands(),hyper.NumRows(),hyper.NumCols());
	for(int k=0; k<hyper.NumBands(); ++k) {
	  for(int i=0; i<hyper.NumRows(); ++i) {
	    for(int j=0; j<hyper.NumCols(); ++j) {
	    	Edge ed=hyper.GetPixel(k,i,j);
	    	EdgeStruct ee=edge.GetPixel(k,i,j);
	    	ee.X=ed.GetX();
	    	ee.Y=ed.GetY(); //hyper.NumRows()-1-ed.GetY();
	    	ee.Theta=ed.GetTheta();
	    	ee.Prob=ed.GetProb();
	    	ee.Mag=ed.GetMagnitude();
	    	ee.On=ed.IsEdge()? 1: 0;
	    	edge.SetPixel(k,i,j,ee);
	    }
	  }
	}
	edge.WriteRawFile(filename);
}

/*void
WriteEdgeImageRaw(const EdgeImage& edge, char* filename) {
	printf("Writing data to %s\n", filename);
	int ne=6;
	RealImage buf(ne*edge.NumBands(),edge.NumRows(),edge.NumCols());
	for(int k=0; k<edge.NumBands(); ++k) {
		for(int i=0; i<edge.NumRows(); ++i) {
			for(int j=0; j<edge.NumCols(); ++j) {
				Edge ed=edge.GetPixel(k,i,j);
				buf.SetPixel(ne*k,i,j,ed.GetX());
				buf.SetPixel(ne*k+1,i,j,ed.GetY());
				buf.SetPixel(ne*k+2,i,j,ed.GetTheta());
				buf.SetPixel(ne*k+3,i,j,ed.GetMagnitude());
				buf.SetPixel(ne*k+4,i,j,ed.GetProb());
				buf.SetPixel(ne*k+5,i,j,(real)ed.IsEdge());
			}
		}
	}
	buf.WriteRawFile(filename);
	//edge.WriteRawFile(filename);
}*/
