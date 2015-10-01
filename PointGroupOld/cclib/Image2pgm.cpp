
#include <iostream.h>
#include <fstream.h>

#include "mytype.h"
#include "Feature.h"
#include "Misc.h"

#include <stdlib.h>
#include <math.h>


int
main(int argc, char* argv[])
{
  ifstream in(argv[1], ios::in);
  if(!in) {
    cerr << argv[0] << " :File open failure." << endl;
    cerr << "Cannot open " << argv[1] << endl;
  }
  in.close();

  RealImage image;

  in >> image;

  InteractiveImageWrite(image, "temp", 1);
}
