
#include "Misc.h"


#ifdef LET_USE_THIS_ONE_INSTEAD_OF_KIM_ONES
double 
rndm(int jd)
{
   /*static int mdig = 16;*/	/* # digits used for integers */
/*
   M1= 2**(MDIG-2) + (2**(MDIG-2)-1),  M2 = 2**(MDIG/2)
*/

   static int m[17] =
   {30788, 23052, 2053, 19346, 10646, 19427, 23975,
   19049, 10949, 19693, 29746, 26748, 2796, 23890, 29168, 31924, 16499};

   static int i = 4, j = 16, m1 = 32767, m2 = 256;
   static double invm1 = .00003051850947599719;  /* 1.0 / m1 */

   int k, jseed, k0, k1, j0, j1;

   if (jd != 0)
   {
      jseed = (jd < 0) ? -jd : jd;
      jseed = (jseed < m1) ? jseed : m1;
      if (jseed % 2 == 0)
	 --jseed;
      k0 = 9069 % m2;
      k1 = 9069 / m2;
      j0 = jseed % m2;
      j1 = jseed / m2;
#if 0
      printf("jseed= %d, k0= %d, k1= %d, j0= %d, j1= %d\n", jseed, k0, k1, j0, j1);
#endif
      for (i = 0; i < 17; ++i)
      {
	 jseed = j0 * k0;
	 j1 = ((jseed / m2) + j0 * k1 + j1 * k0) % (m2 / 2);
	 j0 = jseed % m2;
	 m[i] = j0 + m2 * j1;
#if 0
	 printf("i= %d, m[i]= %d\n", i, m[i]);
#endif
      }
      i = 4;
      j = 16;
   }
   k = m[i] - m[j];
   if (k < 0)
      k += m1;
   m[j] = k;
   if (--i < 0)
      i = 16;
   if (--j < 0)
      j = 16;
   return (double) k *  invm1;
}


double gasdev(long jd)
{
	static int iset=0;
	static real gset;
	double fac,r,v1,v2;
        double ran1();

	if  (iset == 0) {
		do {
			v1=2.0*rndm(jd)-1.0;
			v2=2.0*rndm(jd)-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
		/*return .2132 * v2*fac;*/
	} else {
		iset=0;
		return gset;
		/*return .2132 * gset;*/
	}
}
#endif

/*
this routine add both multiplicative and additive noise to a given
image.
image: the input image
sgm1: the level of the multiplicative noise
sgm2: the level of the additive noise
*/

RealImage
AddNoiseImage(const RealImage& image, real sgm1, real sgm2) {

  RealImage result(image.NumBands(),image.NumRows(),image.NumCols());

  for(int i=0; i<result.NumBands()*result.NumRows()*result.NumCols(); ++i) {
    real val = image.GetPixel(i);
    //val = (1.0+sgm1*gasdev(0))*val + sgm2*gasdev(0);
    val = (1.0+sgm1*rndm(0))*val + sgm2*rndm(0);
    result.SetPixel(i, val);
  }
  return result;
}

RealImage
ResizeImagePowerTwo(const RealImage& image) {
  int rows=image.NumRows();
  int cols=image.NumCols();

  if(rows==1 || rows==2)
	  ;
  if(rows >= 3 && rows <= 4)
    rows = 4;
  else if(rows <= 8)
    rows = 8;
  else if(rows <= 16)
    rows = 16;
  else if(rows <= 32)
    rows = 32;
  else if(rows <= 64)
    rows = 64;
  else if(rows <= 128)
    rows = 128;
  else if(rows <= 256)
    rows = 256;
  else if(rows <= 512)
    rows = 512;
  else if(rows <= 1024)
    rows = 1024;
  else if(rows>1024){
    cerr << "ResizeImagePowerTwo: Image size is " << rows << endl;
    cerr << "Has to be no more than 1024." << endl;
    cerr << "The number of rows truncated to 1024." << endl;
    rows = 1024;
  }

  if(cols==1 || cols==2)
	  ;
  if(cols >= 3 && cols <= 4)
    cols = 4;
  else if(cols <= 8)
    cols = 8;
  else if(cols <= 16)
    cols = 16;
  else if(cols <= 32)
    cols = 32;
  else if(cols <= 64)
    cols = 64;
  else if(cols <= 128)
    cols = 128;
  else if(cols <= 256)
    cols = 256;
  else if(cols <= 512)
    cols = 512;
  else if(cols <= 1024)
    cols = 1024;
  else if(cols>1024){
    cerr << "ResizeImagePowerTwo: Image size is " << cols << endl;
    cerr << "Has to be no more than 1024." << endl;
    cerr << "The number of cols truncated to 1024." << endl;
    cols = 1024;
  }
  
  RealImage result(image.NumBands(), rows, cols, 0);
  result.insert(image, 0, 0, 0);

  return result;
}


ComplexImage
ResizeImagePowerTwo(const ComplexImage& image) {
  int rows=image.NumRows();
  int cols=image.NumCols();

  if(rows==1 || rows==2)
	  ;
  else if(rows > 3 && rows <= 4)
    rows = 4;
  else if(rows <= 8)
    rows = 8;
  else if(rows <= 16)
    rows = 16;
  else if(rows <= 32)
    rows = 32;
  else if(rows <= 64)
    rows = 64;
  else if(rows <= 128)
    rows = 128;
  else if(rows <= 256)
    rows = 256;
  else if(rows <= 512)
    rows = 512;
  else if(rows <= 1024)
    rows = 1024;
  else {
    cerr << "ResizeImagePowerTwo: Image size is " << rows << endl;
    cerr << "Has to be no more than 1024." << endl;
    cerr << "The number of rows truncated to 1024." << endl;
    rows = 1024;
  }

  if(cols==1 || cols==2)
	  ;
  else if(cols > 3 && cols <= 4)
    cols = 4;
  else if(cols <= 8)
    cols = 8;
  else if(cols <= 16)
    cols = 16;
  else if(cols <= 32)
    cols = 32;
  else if(cols <= 64)
    cols = 64;
  else if(cols <= 128)
    cols = 128;
  else if(cols <= 256)
    cols = 256;
  else if(cols <= 512)
    cols = 512;
  else if(cols <= 1024)
    cols = 1024;
  else {
    cerr << "ResizeImagePowerTwo: Image size is " << cols << endl;
    cerr << "Has to be no more than 1024." << endl;
    cerr << "The number of cols truncated to 1024." << endl;
    cols = 1024;
  }
  
  Complex zero(0, 0);
  ComplexImage result(image.NumBands(), rows, cols, zero);
  result.insert(image, 0, 0, 0);

  return result;
}
