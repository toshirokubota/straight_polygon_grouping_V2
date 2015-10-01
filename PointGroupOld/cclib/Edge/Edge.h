#ifndef _Edge_h_
#define _Edge_h_

#include "Image.h"
#include "mytype.h"
#include "Misc.h"
#include <stdlib.h>
#include <assert.h>

const int NumProbLevels=8;

class Edge{
 public:
  Edge(real theta=0, real mag=0, real x=0, real y=0, real prob=0,
       bool on=false) {
    Theta=theta;
    Mag=mag;
    X=x;
    Y=y;
    SetProb(prob);
    On=on;
  }

  Edge(const Edge& edge) {
    Theta=edge.GetTheta();
    Mag=edge.GetMagnitude();
    X=edge.GetX();
    Y=edge.GetY();
    SetProb(edge.GetProb());
    On=edge.IsEdge();
  }

  const Edge& operator=(const Edge& v) {
    Theta=v.GetTheta();
    Mag=v.GetMagnitude();
    X=v.GetX();
    Y=v.GetY();
    SetProb(v.GetProb());
    On=v.IsEdge();
        
    return *this;
  }

  inline real GetTheta() const {return Theta;}
  inline real GetMagnitude() const {return Mag;}
  inline real GetX() const {return X;}
  inline real GetY() const {return Y;}
  inline real GetProb() const {return Prob;}
  inline bool IsEdge() const {return On;}
  inline void SetTheta(real theta) {
    Theta=theta;
  }
  inline void SetMagnitude(real mag) {
    Mag=mag;
  }
  inline void SetX(real x) {
    X=x;
  }
  inline void SetY(real y) {
    Y=y;
  }
  inline void SetProb(real prob) {
    //Prob=(real)Round(prob*NumProbLevels)/NumProbLevels;
    Prob=prob;
  }
  inline void SetFlag(bool flag) {
    On=flag;
  }
  inline real GetXComp() const {return Mag*cos(Theta);}
  inline real GetYComp() const {return Mag*sin(Theta);}

  // Artithmetic overloaded operators

  /*void operator +=(const Edge& v) {
    real theta=v.GetTheta();
    real y = Mag*sin(Theta);
    real x = Mag*cos(Theta);
    y += v.GetMagnitude()*sin(theta);
    x += v.GetMagnitude()*cos(theta);
    Mag=sqrt(y*y+x*x);
    Theta=atan2(y,x);
  }*/

  // Standard Output
  friend ostream& operator<<(ostream& s, const Edge& v) {
    s <<"("<<v.Mag<<","<<v.Theta;
    s << "," << v.X << "," << v.Y;
    s << "," << v.Prob << ")\n";
    return s;
  }

 protected:
  real Theta;
  real Mag;
  real X;
  real Y;
  real Prob;
  bool On;
};

typedef Image<Edge> EdgeImage;

#endif /* Edge_h */
