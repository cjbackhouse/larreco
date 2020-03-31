#include "Util.h"

#include "larcorealg/Geometry/GeometryCore.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace{template<class T> T sqr(T x){return x*x;}}

namespace quad
{
  std::vector<BPPt> GetPts3D(const std::vector<recob::Hit>& hits,
                             const geo::GeometryCore* geom,
                             const detinfo::DetectorProperties* detprop,
                             std::vector<UnitVec>* views,
                             std::vector<UnitVec>* perps)
  {
    std::vector<UnitVec> localViews, localPerps;
    if(!views) views = &localViews;
    if(!perps) perps = &localPerps;

    std::vector<BPPt> ret;

    views->clear();
    perps->clear();

    unsigned int hitIdx = 0;
    for(const recob::Hit& hit: hits){
      const geo::WireID wire = hit.WireID();

      const double xpos = detprop->ConvertTicksToX(hit.PeakTime(), wire);

      const MyVec r0(geom->WireEndPoints(wire).start());
      const MyVec r1(geom->WireEndPoints(wire).end());

      MyVec vd = (r1-r0);
      if(vd.y < 0) vd *= -1;
      const UnitVec dir(vd);

      int view = -1;
      for(unsigned int vi = 0; vi < views->size(); ++vi){
        if(dir.Dot((*views)[vi]) > 0.99){
          view = vi;
          break;
        }
      }
      if(view == -1){
        view = views->size();
        if(view == 3){
          std::cout << "TOO MANY VIEWS!" << std::endl;
          for(const UnitVec& v: *views) v.Print();
          dir.Print();
          abort();
        }
        views->push_back(dir);

        // Compute the direction perpendicular to the wires
        // We want to ultimately have a positive z component in "perp"
        if(dir.Y() > 0){
          perps->emplace_back(MyVec(0, -dir.Z(), +dir.Y()));
        }
        else{
          perps->emplace_back(MyVec(0, +dir.Z(), -dir.Y()));
        }
      }

      ret.emplace_back(MyVec(xpos, r0.y, r0.z), dir, view, hitIdx++);
      ret.back().z = (*perps)[view].Dot(ret.back().ray.Origin());
    } // end for hits

    std::cout << ret.size() << " hits in " << views->size() << " views" << std::endl;

    return ret;
  }

} // namespace
