#ifndef _EDGE_PROC_
#define _EDGE_PROC_

#include "Edge.h"

EdgeImage
ProbabilityProjection(const EdgeImage& hyper);

EdgeImage
newProbabilityProjection(EdgeImage& hyper);

ByteImage
ExpandEdges(const EdgeImage& hyper, int factor);

EdgeImage
SetEdges(const EdgeImage& hyper, const RealImage& grd, real thres);

EdgeImage 
RearrangeEdgeImage(const EdgeImage& edge, int range, real dist);

void
InteractiveEdgeImageWrite(const EdgeImage& edge, 
                          char* filename, int normalize);

#endif /* _EDGE_PROC_ */