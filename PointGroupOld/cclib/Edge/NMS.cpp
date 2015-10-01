#include "NMS.h"

RealImage
NonMaximumSuppressionQuality(const RealImage& qq, const EdgeImage& edge) {
  
  RealImage q2(qq.NumBands(),qq.NumRows(),qq.NumCols(),0);
  for(int k=0; k<qq.NumBands(); ++k) {
    for(int i=0; i<qq.NumRows(); ++i) {
      for(int j=0; j<qq.NumCols(); ++j) {
        real gc=qq.GetPixel(k,i,j);
			  real gn=qq.GetPixelZero(k,i-1,j);
			  real gs=qq.GetPixelZero(k,i+1,j);
			  real gw=qq.GetPixelZero(k,i,j-1);
			  real ge=qq.GetPixelZero(k,i,j+1);
			  real gne=qq.GetPixelZero(k,i-1,j+1);
			  real gse=qq.GetPixelZero(k,i+1,j+1);
			  real gsw=qq.GetPixelZero(k,i+1,j-1);
			  real gnw=qq.GetPixelZero(k,i-1,j-1);
                
        real tt=(edge.GetPixel(k,i,j)).GetTheta();
        real uy=-cos(tt);
        real ux=sin(tt);
        real g0;
			  if(ux*uy>0) {
			    if(Abs(ux)<Abs(uy)) {
				    if((g0=Abs(uy*gc))<Abs(ux*gse+(uy-ux)*gs) ||
		            g0<=Abs(ux*gnw+(uy-ux)*gn))
		          continue;
		      }
			    else {
				    if((g0=Abs(ux*gc))<Abs(uy*gse+(ux-uy)*ge) ||
		            g0<=Abs(uy*gnw+(ux-uy)*gw))
		          continue;
		      }
		    }
			  else {
			    if(Abs(ux)<Abs(uy)) {
				    if((g0=Abs(uy*gc))<Abs(ux*gne-(uy+ux)*gn) ||
		            g0<=Abs(ux*gsw-(uy+ux)*gs))
				      continue;
		      }
			    else {
				    if((g0=Abs(ux*gc))<Abs(uy*gne-(ux+uy)*ge) ||
		            g0<=Abs(uy*gsw-(ux+uy)*gw))
				      continue;
		      }
		    }
			  q2.SetPixel(k,i,j,qq.GetPixel(k,i,j));
      }
    }
  }
  for(int i=0; i<qq.NumRows(); ++i) {
    for(int j=0; j<qq.NumCols(); ++j) {
      real max=0;
      int maxk=0;
      for(int k=0; k<qq.NumBands(); ++k) {
        if(q2.GetPixel(k,i,j)>max) {
          max=q2.GetPixel(k,i,j);
          maxk=k;
        }
      }
      if(max>0) {
        for(int k=0; k<qq.NumBands(); ++k) {
	        if(k!=maxk) {
	          q2.SetPixel(k,i,j,.0);
	        }
	      }
	    }
	  }
	}
  return q2;
}
