#ifndef _NMS_h_
#define _NMS_h_

#include <Image.h>
#include <mytype.h>
#include <Edge.h>

RealImage
NonMaximumSuppressionQuality(const RealImage& qq, 
														const EdgeImage& edge);

#endif /* NMS_h */