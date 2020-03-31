#pragma once

#include "MyVec.h"

#include "Ray.h"

namespace quad
{
  // For evd
  struct BPPt2D
  {
    BPPt2D(double _x, double _z, int _view, double _energy, int _idx) : x(_x), z(_z), view(_view), energy(_energy), hitIdx(_idx), priv(0) {}
    
    bool operator<(const BPPt2D& p) const {return z < p.z;}

    double x, z;
    int view;
    double energy;
    int hitIdx;
    int priv; // for use by algorithms
  };

  struct BPPt
  {
    BPPt(MyVec _origin, MyVec _dir, int _view, int _idx) : ray(_origin, _dir), view(_view), hitIdx(_idx), priv(0) {}
    BPPt(MyVec _origin, UnitVec _dir, int _view, int _idx) : ray(_origin, _dir), view(_view), hitIdx(_idx), priv(0) {}

    Ray ray;
    int view;
    float z;

    int hitIdx;
    int priv; // for use by algorithms
  };
}
