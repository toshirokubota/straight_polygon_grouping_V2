#include "Draw.h"

void
DrawLine(RealImage& image, int b, int y1, int x1, int y2, int x2, real val) {
  int inc;
  int i,n;
  if(Abs(y1-y2)>Abs(x1-x2)) {
    if(y1>y2)
      inc=-1;
    else
      inc=1;
    real dx=(real)x1;
    real delta=(real)(x2-x1)/Abs(y2-y1);
    for(i=y1,n=0; n<=Abs(y1-y2); i+=inc,n++) {
      if(i>=0 && i<image.NumRows()) {
	      int j=Round(dx);
	      if(j>=0 && j<image.NumCols())
	        image.SetPixel(b,i,j,val); //Max(val,image.GetPixel(b,i,j)));
	      dx+=delta;
      }
    }
  }
  else {
    if(x1>x2)
      inc=-1;
    else
      inc=1;
    real dy=(real)y1;
    real delta=(real)(y2-y1)/Abs(x2-x1);
    for(i=x1,n=0; n<=Abs(x1-x2); i+=inc,n++) {
      if(i>=0 && i<image.NumCols()) {
	      int j=Round(dy);
	      if(j>=0 && j<image.NumRows())
	        image.SetPixel(b,j,i,val); //Max(val,image.GetPixel(b,j,i)));
	      dy+=delta;
	    }
    }
  }
}

void
DrawCross(RealImage& image, int b, int y, int x, int length, real val) {
  for(int i=y-length/2; i<=y+length/2; ++i) {
    if(i>=0 && i<image.NumRows() && x>=0 && x<image.NumCols())
      image.SetPixel(b,i,x,val); 
  }
  for(int j=x-length/2; j<=x+length/2; ++j) {
    if(j>=0 && j<image.NumCols() && y>=0 && y<image.NumRows())
      image.SetPixel(b,y,j,val); 
  }
}

void
DrawLine(RealImage& image, int b, int y1, int x1, int y2, int x2) {
  int inc;
  int i,n;
  if(Abs(y1-y2)>Abs(x1-x2)) {
    if(y1>y2)
      inc=-1;
    else
      inc=1;
    real dx=(real)x1;
    real delta=(real)(x2-x1)/Abs(y2-y1);
    for(i=y1,n=0; n<=Abs(y1-y2); i+=inc,n++) {
      if(i>=0 && i<image.NumRows()) {
	      int j=Round(dx);
	      if(j>=0 && j<image.NumCols())
	        image.SetPixel(b,i,j,LineValue);
	      dx+=delta;
      }
    }
  }
  else {
    if(x1>x2)
      inc=-1;
    else
      inc=1;
    real dy=(real)y1;
    real delta=(real)(y2-y1)/Abs(x2-x1);
    for(i=x1,n=0; n<=Abs(x1-x2); i+=inc,n++) {
      if(i>=0 && i<image.NumCols()) {
	      int j=Round(dy);
	      if(j>=0 && j<image.NumRows())
	        image.SetPixel(b,j,i,LineValue);
	      dy+=delta;
      }
    }
  }
}

void
DrawCross(RealImage& image, int b, int y, int x, int length) {
  for(int i=y-length/2; i<=y+length/2; ++i) {
    if(i>=0 && i<image.NumRows() && x>=0 && x<image.NumCols())
      image.SetPixel(b,i,x,CornerValue); 
  }
  for(int j=x-length/2; j<=x+length/2; ++j) {
    if(j>=0 && j<image.NumCols() && y>=0 && y<image.NumRows())
      image.SetPixel(b,y,j,CornerValue); 
  }
}

void
DrawLine(ByteImage& image, int b, int y1, int x1, int y2, int x2, unsigned char val) {
  int inc;
  int i,n;
  if(Abs(y1-y2)>Abs(x1-x2)) {
    if(y1>y2)
      inc=-1;
    else
      inc=1;
    real dx=(real)x1;
    real delta=(real)(x2-x1)/Abs(y2-y1);
    for(i=y1,n=0; n<=Abs(y1-y2); i+=inc,n++) {
      if(i>=0 && i<image.NumRows()) {
	      int j=Round(dx);
	      if(j>=0 && j<image.NumCols())
	        image.SetPixel(b,i,j,val); //Max(val,image.GetPixel(b,i,j)));
	      dx+=delta;
      }
    }
  }
  else {
    if(x1>x2)
      inc=-1;
    else
      inc=1;
    real dy=(real)y1;
    real delta=(real)(y2-y1)/Abs(x2-x1);
    for(i=x1,n=0; n<=Abs(x1-x2); i+=inc,n++) {
      if(i>=0 && i<image.NumCols()) {
	      int j=Round(dy);
	      if(j>=0 && j<image.NumRows())
	        image.SetPixel(b,j,i,val); //Max(val,image.GetPixel(b,j,i)));
	      dy+=delta;
	    }
    }
  }
}

void
DrawCross(ByteImage& image, int b, int y, int x, int length, unsigned char val) {
  for(int i=y-length/2; i<=y+length/2; ++i) {
    if(i>=0 && i<image.NumRows() && x>=0 && x<image.NumCols())
      image.SetPixel(b,i,x,val); 
  }
  for(int j=x-length/2; j<=x+length/2; ++j) {
    if(j>=0 && j<image.NumCols() && y>=0 && y<image.NumRows())
      image.SetPixel(b,y,j,val); 
  }
}
