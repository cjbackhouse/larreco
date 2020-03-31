#pragma once

#include "lardataobj/RecoBase/Hit.h"

#include "BPPt.h"

#include "TMatrixD.h"

#include <vector>

namespace geo{class GeometryCore;}
namespace detinfo{class DetectorProperties;}

namespace quad
{
  std::vector<BPPt> GetPts3D(const std::vector<recob::Hit>& hits,
                             const geo::GeometryCore* geom,
                             const detinfo::DetectorProperties* detprop,
                             std::vector<UnitVec>* views = 0,
                             std::vector<UnitVec>* perps = 0);

  /// TMatrixD::operator() does various sanity checks and shows up in profiles
  class TMatrixAccessor
  {
  public:
    TMatrixAccessor(TMatrixD& m)
      : fN(m.GetNrows()),
        fArray(m.GetMatrixArray())
    {
    }
    inline double operator()(unsigned int i, unsigned int j) const
    {
      return fArray[fN*i+j];
    }
    inline double& operator()(unsigned int i, unsigned int j)
    {
      return fArray[fN*i+j];
    }
  protected:
    unsigned int fN;
    double* fArray;
  };
}

