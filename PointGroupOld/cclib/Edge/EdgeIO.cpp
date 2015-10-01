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
				ed.SetFlag(ee.On);
				edge.SetPixel(k,i,j,ed);
			}
		}
	}
	return edge;
}


EdgeImage
ReadEdgeImageRaw(char* filename) {
	printf("Reading data from %s\n", filename);
	EdgeStructImage buf;
	buf.ReadRawFile(filename);
	EdgeImage edge=EdgeStructImage2EdgeImage(buf);
	return edge;
}

void
WriteEdgeImageRaw(const EdgeImage& hyper, char* filename) {

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