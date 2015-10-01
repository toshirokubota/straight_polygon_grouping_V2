#ifndef _EDGEIO_H_
#define _EDGEIO_H_
#include <mytype.h>
#include "Edge.h"

typedef struct struct_edge {
  real X;
  real Y;
  real Mag;
  real Theta;
  real Prob;
  unsigned char On;
} EdgeStruct;

typedef Image<EdgeStruct> EdgeStructImage;

EdgeStructImage
EdgeImage2EdgeStructImage(const EdgeImage& edges);

EdgeImage
EdgeStructImage2EdgeImage(const EdgeStructImage& edges);

EdgeImage
ReadEdgeImageRaw(char* filename);

void
WriteEdgeImageRaw(const EdgeImage& edge, char* filename);

void
InteractiveEdgeImageWrite(const EdgeImage& edge, 
                          char* filename, int normalize);

ostream& operator<<(ostream& s, const Edge& ed);

#endif /* _EDGEIO_H_ */ 