// Christopher Backhouse - c.backhouse@ucl.ac.uk - Aug 2019

// Infinitesimally faster, maybe
//#pragma GCC optimize("O3", "fast-math")

// C/C++ standard libraries
#include <string>
#include <iostream>

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"

// LArSoft libraries
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

//#include "Breakpoint.h"
#include "Util.h"
#include "MyVec.h"
#include "Ray.h"
#include "BPPt.h"

#include "TDecompLU.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "FastRand.h"

template<class T> inline T sqr(T x){return x*x;}
template<class T> inline T cube(T x){return x*x*x;}


const float maxdL = 2; // 2cm in 3D
const float maxd = .25; // 2.5mm laterally
const float maxdsq = sqr(maxd);

const int kMinMatch = 6;//10;//20; // Don't make tracks smaller than this


namespace quad
{

class QuadPts : public art::EDProducer
{
public:
  explicit QuadPts(const fhicl::ParameterSet& pset);

private:
  void produce(art::Event& evt) override;
  void beginJob() override;
  void endJob() override;

  std::string fHitLabel;

  const detinfo::DetectorProperties* detprop;
  const geo::GeometryCore* geom;
};

DEFINE_ART_MODULE(QuadPts)

// ---------------------------------------------------------------------------
QuadPts::QuadPts(const fhicl::ParameterSet& pset) :
  EDProducer(pset),
  fHitLabel(pset.get<std::string>("HitLabel"))
{
  produces<std::vector<recob::Track>>();
  produces<art::Assns<recob::Track, recob::Hit>>();
}

// ---------------------------------------------------------------------------
void QuadPts::beginJob()
{
  detprop = art::ServiceHandle<const detinfo::DetectorPropertiesService>()->provider();
  geom = art::ServiceHandle<const geo::Geometry>()->provider();
}

// ---------------------------------------------------------------------------
void QuadPts::endJob()
{
}

// ---------------------------------------------------------------------------
Ray PointsToRay(const std::array<BPPt, 4>& pts)
{
  TMatrixD M(4, 4); // my, mz, cy, cz
  TVectorD v(4);


  for(int i = 0; i < 4; ++i){
    const UnitVec d = pts[i].ray.Dir();
    const MyVec o   = pts[i].ray.Origin();
    // https://en.wikipedia.org/wiki/Skew_lines#Distance
    //
    // (mx,my,mz).Cross(d).Dot(o - (cx,cy,cz)) == 0
    //
    // (my*d.z - mz*d.y)*(o.x-cx) +
    // (mz*d.x - mx*d.z)*(o.y-cy) +
    // (mx*d.y - my*d.x)*(o.z-cz) == 0
    //
    // we know d.x = 0, and we can set mx = 1 and cx = 0
    //
    // (my*d.z - mz*d.y)*o.x + (-d.z)*(o.y-cy) + (d.y)*(o.z-cz) == 0

    M(i, 0) =  d.Z()*o.x; // my coeff
    M(i, 1) = -d.Y()*o.x; // mz coeff
    M(i, 2) =  d.Z();     // cy coeff
    M(i, 3) = -d.Y();     // cz coeff

    v(i) = d.Z()*o.y - d.Y()*o.z; // const part
  }

  TDecompLU d(M);
/*bool ok = */d.Solve(v);

//  std::cout << "OK? " << ok << std::endl;

  const Ray ret(MyVec(0, v[2], v[3]), MyVec(1, v[0], v[1]));

  for(int i = 0; i < 4; ++i){
    const UnitVec d = pts[i].ray.Dir();
    const MyVec o   = pts[i].ray.Origin();
    std::cout << ret.Dir().Cross(d).Dot(o-ret.Origin()) << std::endl;
  }

  return Ray(MyVec(0, v[2], v[3]), MyVec(1, v[0], v[1]));

/*
  try{
    M.Invert();
  }
  catch(...){
    std::cout << "Singular matrix in QuadPts::PointsToRay" << std::endl;
    M.Print();
    // Matrix will probably be garbage, but it's OK because we'll just get a
    // bad score for the line.
  }

  const TVectorD x = M*v;

  return Ray(MyVec(0, x[2], x[3]), MyVec(1, x[0], x[1]));
*/
}

struct MinimalPt
{
  MinimalPt(float _z, float _x) : z(_z), x(_x) {}

  float DistSq(const MinimalPt& p) const {return sqr(z-p.z) + sqr(x-p.x);}

  float z, x;
};

struct MinimalLine
{
  MinimalLine(MinimalPt a, MinimalPt b)
    : dzdx((b.z-a.z)/(b.x-a.x)),
      z0(a.z-dzdx*a.x)
  {
  }

  MinimalLine(float _dzdx, float _z0) : dzdx(_dzdx), z0(_z0) {}

  float dzdx, z0;
};

class VPLeaf
{
public:
  VPLeaf(MinimalPt _origin, float _radius,
         std::vector<MinimalPt>::iterator begin,
         std::vector<MinimalPt>::iterator end)
    : origin(_origin), radius(_radius), pts(begin, end)
  {
  }

  MinimalPt origin;
  float radius;
  std::vector<MinimalPt> pts;
};

class VPTree
{
public:
  VPTree(std::vector<MinimalPt>::iterator begin,
         std::vector<MinimalPt>::iterator end)
  {
    // Bigger hack
    kids.emplace_back(*begin, std::numeric_limits<float>::infinity(), begin, end);
    return;

    // Adjustable parameter. This one gets ~10 pts per leaf
    const float kRadius = 5;

    while(begin != end){
      MinimalPt origin = *begin;
      // HACK HACK HACK
      //      auto half = begin;
      //      for(int i = 0; i < 10; ++i) if(half != end) ++half;
      auto half = std::partition(begin, end, [origin, kRadius](const MinimalPt& p){return p.DistSq(origin) < sqr(kRadius);});
      kids.emplace_back(origin, kRadius, begin, half);
      begin = half;
    }
  }

  VPTree(const VPTree&) = delete;

  ~VPTree()
  {
  }

  std::vector<VPLeaf> kids;
};

class TreeLayer
{
public:
  std::bitset<16> kids;
};

class VPTree2
{
public:
  VPTree2(const std::vector<MinimalPt>& pts)
  {
    fLayers[std::make_tuple(0, 0, 0)] = TreeLayer(); // top level

    zoff = std::numeric_limits<float>::infinity();
    xoff = std::numeric_limits<float>::infinity();

    for(MinimalPt pt: pts){
      zoff = std::min(zoff, pt.z);
      xoff = std::min(xoff, pt.x);
    }

    xoff -= .5;
    zoff -= .5;

    zoff *= -1;
    xoff *= -1;

    for(MinimalPt pt: pts) Add(pt);
  }

  void Add(MinimalPt pt)
  {
    pt.z += zoff;
    pt.x += xoff;

    int stride = 1;
    const float cellSize = 1; // cm

    // True global position
    const int i = int(pt.x/cellSize);
    const int j = int(pt.z/cellSize);
    fPts[std::make_pair(i, j)].push_back(pt);

    for(int level = 9; level >= 0; --level){
      fSets[level].emplace(i/stride, j/stride);
      stride *= 2;
    }

    return;
    /*
    std::cout << pt.x << " " << pt.z << " -> " << i << " " << j << std::endl;

    TreeLayer* layer = &fLayers[std::make_tuple(0, 0, 0)];
    for(int level = 0; level < 10; ++level){
      const int i0 = (i/stride)*stride;
      const int j0 = (j/stride)*stride;

      const int di = i-i0;
      const int dj = j-j0;

      std::cout << " = " << i0 << "+" << di << ", " << j0 << "+" << dj << std::endl;

      const auto global_key = std::make_tuple(level, i0, j0);

      if(layer->kids[dj*4+di]){
        layer = &fLayers[global_key];
      }
      else{
        layer->kids[dj*4+di] = 1;
        if(level != 9){
          fLayers[global_key] = TreeLayer();
          layer = &fLayers[global_key];
        }
      }

      stride /= 2;
    } // end for level

    abort();
    */
  }

  int Count(MinimalLine line) const
  {
    line.z0 += zoff - xoff*line.dzdx; // correct for point offset

    int stride = 512; // TODO 2^9
    const float cellSize = 1; // cm

    const float angleFactor = sqrt(1+sqr(line.dzdx));
    const float maxdz = maxd * angleFactor;

    std::vector<std::pair<int, int>> todo, next;
    todo.emplace_back(0, 0);

    for(int level = 0; level < 9/*10*/; ++level){
      const float Rz = stride*cellSize/sqrt(2) * angleFactor; // cover the corners

      for(std::pair<int, int> c: todo){
        const float x0 = (c.first+.5)*stride*cellSize;
        const float z0 = (c.second+.5)*stride*cellSize;

        const float dz = fabs(line.dzdx * x0 + line.z0 - z0);
        //        std::cout << d << " " << R << std::endl;
        if(dz < Rz + maxdz){
          for(int di = 0; di <= 1; ++di){
            for(int dj = 0; dj <= 1; ++dj){
              const auto key = std::make_pair(c.first*2+di, c.second*2+dj);
              if(fSets[level+1].count(key)) next.push_back(key);
            } // end for dj
          } // end for di
        } // end if
      } // end for c (todo)

      //      std::cout << "Level " << level << " todo " << todo.size() << " next " << next.size() << std::endl;

      todo.swap(next);
      //      std::swap(todo, next);
      next.clear();

      stride /= 2;
    } // end for level

    int ret = 0;
    for(std::pair<int, int> c: todo){
      for(MinimalPt p: fPts.find(c)->second){
        const float dz = fabs(line.dzdx * p.x + line.z0 - p.z);
        if(dz < maxdz) ++ret;
      }
    }

    // TODO actually inspect the points
    return ret;
    //    const int stride = 1024; // TODO

    //    return Count(line, 0, 0, 0, stride);
  }

protected:
  int Count(MinimalLine line,
            int level, int i0, int j0,
            int stride) const
  {
    const double cellSize = 1; // cm

    const double R = stride*cellSize/sqrt(2); // cover the corners

    const double mfactor = 1/sqrt(1+sqr(line.dzdx)); // TODO hoist?

    // We know this exists because the caller checked the kids mask
    const TreeLayer layer = fLayers.find(std::make_tuple(level, i0, j0))->second;

    int ret = 0;

    for(int i = 0; i < 4; ++i){ // TODO loop to 16?
      for(int j = 0; j < 4; ++j){
        if(!layer.kids[j*4+i]) continue;

        const double z = (i0+(i+.5)*stride)*cellSize;
        const double x = (j0+(j+.5)*stride)*cellSize;

        const double d = fabs(line.dzdx * x + line.z0 - z)*mfactor;

        if(d < R){ // line hits the circle - TODO padding
          if(level == 10){
            // TODO check them individually
            ret += fPts.find(std::make_pair(i0*i*stride, j0+j*stride))->second.size();
          }
          else{
            // Recurse
            ret += Count(line, level+1, i0+i*stride, j0*j+stride, stride/2);
          }
        }
      }
    }

    return ret;
  }

  float zoff, xoff;

  std::array<std::set<std::pair<int, int>>, 10> fSets;

  std::map<std::tuple<int, int, int>, TreeLayer> fLayers;

  // Final layer
  std::map<std::pair<int, int>, std::vector<MinimalPt>> fPts;
};

class View
{
public:
  View(std::vector<MinimalPt>& _pts,
       const UnitVec& _dir,
       const UnitVec& _perp) :
    npts(_pts.size()),
    dir(_dir), perp(_perp),
    vptree(_pts.begin(), _pts.end())
  {
  }

  unsigned int npts;

  UnitVec dir, perp;

  VPTree vptree;
};







// ---------------------------------------------------------------------------
int CountClosePoints(const View& view, const MinimalLine& line) throw()
{
  int ret = 0;

  // Inititialize a lot of constants we can hoist out of the loop
  //  const UnitVec bd = line.Dir().Cross(view.dir);

  //  MyVec n2 = view.dir.Cross(bd);
  //  n2 *= 1./line.Dir().Dot(n2);

  //  // dz/dx
  //  const float m = line.Dir().Dot(view.perp) / line.Dir().X();

  // z(x=0)
  //  const float c = line.Origin().Dot(view.perp) - line.Origin().x / line.Dir().X() * line.Dir().Dot(view.perp);

  const float angleFactor = sqrt(1+sqr(line.dzdx));
  const float maxdz = maxd * angleFactor;

  //  const float lambda0 = line.Origin().Dot(n2);

  //  const float lambdax = n2.x;
  //  const float lambdaz = view.perp.Y()*n2.y + view.perp.Z()*n2.z;

  // This is all extremely hot and worth closely optimizing
  for(const VPLeaf& leaf: view.vptree.kids){
    const float dzo = fabs(line.dzdx*leaf.origin.x + line.z0 - leaf.origin.z);

    // If the line misses the circle plus maxd padding then we can dismiss
    // all the children.
    if(dzo > (leaf.radius+maxd) * angleFactor) continue;

    for(const MinimalPt& p: leaf.pts){
      // Surprisingly fabs() here is substantially faster than sqr()
      const float dz = fabs(line.dzdx*p.x + line.z0 - p.z);

      // Hmm, seems not really
      //      const float dz2 = sqr(line.dzdx*p.x + line.z0 - p.z);

      if(dz < maxdz){
        ++ret;
        //        const float lambda = p.x * lambdax + p.z * lambdaz - lambda0;
        //        Ls.emplace_back(lambda, v);
      }
    } // end for p
  } // end for leaf

  return ret;
}


// ---------------------------------------------------------------------------
int CountClosePoints(const View& view, const Ray& line) throw()
{
  //  int ret = 0;

  // Inititialize a lot of constants we can hoist out of the loop
  //  const UnitVec bd = line.Dir().Cross(view.dir);

  //  MyVec n2 = view.dir.Cross(bd);
  //  n2 *= 1./line.Dir().Dot(n2);

  const MyVec r0 = line.Origin();
  const MyVec r1 = line.Origin() + 10*line.Dir();

  const double x0 = r0.x;
  const double x1 = r1.x;
  const double z0 = r0.Dot(view.perp);
  const double z1 = r1.Dot(view.perp);

  const double m2 = (z1-z0)/(x1-x0);
  const double c2 = z0-m2*x0;

  // dz/dx
  const float m = line.Dir().Dot(view.perp) / line.Dir().X();

  // z(x=0)
  const float c = line.Origin().Dot(view.perp) - m * line.Origin().x;

  std::cout << "  m1 c1 " << m << " " << c << " " << m2 << " " << c2 << std::endl;

  MinimalLine ml(m, c);
  return CountClosePoints(view, ml);

  /*
  const float angleFactor = sqrt(1+sqr(m));
  const float maxdz = maxd * angleFactor;

  //  const float lambda0 = line.Origin().Dot(n2);

  //  const float lambdax = n2.x;
  //  const float lambdaz = view.perp.Y()*n2.y + view.perp.Z()*n2.z;

  // This is all extremely hot and worth closely optimizing
  for(const VPLeaf& leaf: view.vptree.kids){
    const float dzo = fabs(m*leaf.origin.x + c - leaf.origin.z);

    // If the line misses the circle plus maxd padding then we can dismiss
    // all the children.
    if(dzo > (leaf.radius+maxd) * angleFactor) continue;

    for(const MinimalPt& p: leaf.pts){
      // Surprisingly fabs() here is substantially faster than sqr()
      const float dz = fabs(m*p.x + c - p.z);
      if(dz < maxdz){
        ++ret;
        //        const float lambda = p.x * lambdax + p.z * lambdaz - lambda0;
        //        Ls.emplace_back(lambda, v);
      }
    } // end for p
  } // end for leaf

  return ret;
  */
}


// ---------------------------------------------------------------------------
int LongestChainCount(const std::array<View, 3>& views,
                      const Ray& line,
                      unsigned int target,
                      std::vector<int> viewIdxs = {0, 1, 2}) throw()
{
  static std::vector<std::pair<float, int>> Ls;
  Ls.clear(); // try to avoid frequent re-allocation

  for(int v: viewIdxs){
    const View& view = views[v];

    // Inititialize a lot of constants we can hoist out of the loop
    const UnitVec bd = line.Dir().Cross(view.dir);

    MyVec n2 = view.dir.Cross(bd);
    n2 *= 1./line.Dir().Dot(n2);

    // dz/dx
    const float m = line.Dir().Dot(view.perp) / line.Dir().X();

    // z(x=0)
    const float c = line.Origin().Dot(view.perp) - line.Origin().x / line.Dir().X() * line.Dir().Dot(view.perp);

    const float angleFactor = sqrt(1+sqr(m));
    const float maxdz = maxd * angleFactor;

    const float lambda0 = line.Origin().Dot(n2);

    const float lambdax = n2.x;
    const float lambdaz = view.perp.Y()*n2.y + view.perp.Z()*n2.z;

    // This is all extremely hot and worth closely optimizing
    for(const VPLeaf& leaf: view.vptree.kids){
      const float dzo = fabs(m*leaf.origin.x + c - leaf.origin.z);

      // If the line misses the circle plus maxd padding then we can dismiss
      // all the children.
      if(dzo > (leaf.radius+maxd) * angleFactor) continue;

      for(const MinimalPt& p: leaf.pts){
        // Surprisingly fabs() here is substantially faster than sqr()
        const float dz = fabs(m*p.x + c - p.z);
        if(dz < maxdz){
          const float lambda = p.x * lambdax + p.z * lambdaz - lambda0;
          Ls.emplace_back(lambda, v);
        }
      } // end for p
    } // end for leaf
  } // end for v

  if(Ls.size() <= target) return -1;

  // Must have at least two hits in each view
  int nv[3] = {0, 0, 0};
  for(auto L: Ls) ++nv[L.second];
  for(int v: viewIdxs) if(nv[v] < 2) return -1;

  std::array<int, 3> nvcur = {0, 0, 0};
  std::array<int, 3> nvlong = {0, 0, 0};

  std::sort(Ls.begin(), Ls.end());

  float prevL[3] = {-1e10, -1e10, -1e10};

  // The longest run found so far
  int longest = 0;
  int current = 0;

  for(auto it = Ls.begin(); it != Ls.end(); ++it){
    const float L = it->first;

    bool gap = false;
    for(int v: viewIdxs) gap = gap || (L - prevL[v] > maxdL);

    // Gap in any one view is big enough to end segment
    if(gap){
      // Begin a new segment
      current = 0;
      nvcur = {0, 0, 0};
      for(int v = 0; v < 3; ++v) prevL[v] = L;

      // The remaining hits will be insufficient to beat the record
      if(Ls.end()-it < longest) break;
    }

    ++current;
    ++nvcur[it->second];
    if(current > longest){
      longest = current;
      nvlong = nvcur;
    }
    prevL[it->second] = L;
  } // end for (L, view)

  // Must still have at least two hits in each view
  for(int v: viewIdxs) if(nvlong[v] < 2) return -1;

  return longest;
}

// ---------------------------------------------------------------------------
// Reorders the chain to have the longest chain found at the end,
// and updates 'begin' to point at the start of that chain
int LongestChainApply(std::vector<BPPt>::iterator& begin,
                      const std::vector<BPPt>::iterator end,
                      const Ray& line, unsigned int target)
{
  // The longest run found so far
  static std::vector<std::vector<BPPt>::iterator> longest, current;
  static std::vector<std::pair<float, std::vector<BPPt>::iterator>> Ls;
  longest.clear();
  current.clear();
  Ls.clear(); // try to avoid frequent re-allocation

  for(auto it = begin; it != end; ++it){
    if(line.SqrDistanceToRay(it->ray) < maxdsq){
      const float distAlong = line.ClosestLambdaToRay(it->ray);
      Ls.emplace_back(distAlong, it);
    }
  }

  if(Ls.size() < target) return -1;

  std::sort(Ls.begin(), Ls.end());

  float prevL[3] = {-std::numeric_limits<float>::max(),
                    -std::numeric_limits<float>::max(),
                    -std::numeric_limits<float>::max()};

  for(auto it: Ls){
    const float L = it.first;

    bool gap = false;
    for(int v = 0; v < 3; ++v) gap = gap || (L - prevL[v] > maxdL);

    // Gap in any one view is big enough to end segment
    if(gap){
      // Begin a new segment
      current.clear();
      for(int v = 0; v < 3; ++v) prevL[v] = L;
    }

    current.push_back(it.second);
    if(current.size() > longest.size()) longest = current;
    prevL[it.second->view] = L;
  } // end for (L, view)

  if(current.size() > longest.size()) longest = current;

  if(longest.size() < target) return -1;

  // Mark all the points on the longest run we found
  for(auto it: longest) it->priv = 1;

  // Move all the ones we've marked to the end
  begin = std::partition(begin, end, [](const BPPt& p){return p.priv == 0;});

  // Reset the flag
  for(auto it = begin; it != end; ++it) it->priv = 0;

  return longest.size();
}

// ---------------------------------------------------------------------------
void RandomPairs(std::vector<BPPt>::iterator begin,
                 std::vector<BPPt>::iterator end,
                 FastRand& r,
                 std::vector<std::pair<BPPt, BPPt>>& ret)
{
  const int kMaxPairs = 10000;

  ret.clear();
  ret.reserve(kMaxPairs);

  if(((end-begin)*(end-begin-1))/2 <= kMaxPairs){
    //    std::cout << "Deterministic pairs!" << std::endl;
    for(auto it = begin; it != end; ++it){
      for(auto it2 = it+1; it2 != end; ++it2){
        ret.emplace_back(*it, *it2);
      }
    }
  }
  else{
    for(int trial = 0; trial < kMaxPairs; ++trial){
      ret.emplace_back(*(begin + r.xorshift32()%(end-begin)),
                       *(begin + r.xorshift32()%(end-begin)));
    }
  }
}

// ---------------------------------------------------------------------------
void RandomQuads(int N,
                 std::vector<BPPt>::iterator begin,
                 std::vector<BPPt>::iterator end,
                 FastRand& r,
                 std::vector<std::array<BPPt, 4>>& ret)
{
  ret.clear();
  ret.reserve(N);

  if((double(end-begin)*double(end-begin-1)*double(end-begin-2)*double(end-begin-3))/24 <= N){
    std::cout << "Deterministic quads!!!" << std::endl;
    for(auto it = begin; it != end; ++it){
      for(auto it2 = it+1; it2 != end; ++it2){
        for(auto it3 = it2+1; it3 != end; ++it3){
          for(auto it4 = it3+1; it4 != end; ++it4){
            ret.push_back({*it, *it2, *it3, *it4});
          }
        }
      }
    }
  }
  else{
    for(int trial = 0; trial < N; ++trial){
      ret.push_back({*(begin + r.xorshift32()%(end-begin)),
                     *(begin + r.xorshift32()%(end-begin)),
                     *(begin + r.xorshift32()%(end-begin)),
                     *(begin + r.xorshift32()%(end-begin))});
    }
  }
}

// ---------------------------------------------------------------------------
void CandidateLinesRandom(int N,
                          std::vector<BPPt>::iterator begin,
                          std::vector<BPPt>::iterator end,
                          const std::array<View, 3>& views,
                          FastRand& r,
                          std::vector<std::pair<Ray, int>>& lines)
{
  std::vector<std::array<BPPt, 4>> quads;
  RandomQuads(N, begin, end, r, quads);

  for(auto seed: quads){
    // Avoid cases that give a degenerate solution
    bool ok = true;
    for(int i = 0; i < 3; ++i)
      for(int j = i+1; j < 4; ++j)
        if(seed[i].ray.Origin().x == seed[j].ray.Origin().x &&
           seed[i].view == seed[j].view) ok = false;
    if(!ok) continue;

    // 2+2+0 or 2+1+1 is OK, 3/4+anything is not
    int nperview[3] = {0, 0, 0};
    for(const BPPt& s: seed) ++nperview[s.view];
    for(int v = 0; v < 3; ++v) if(nperview[v] >= 3) ok = false;
    if(!ok) continue;

    const Ray line = PointsToRay(seed);

    // Very horizontal lines can cause problems
    if(sqr(line.Dir().X()) < 1e-10*(sqr(line.Dir().Y() + line.Dir().Z()))) continue;


    const int nmatch = LongestChainCount(views, line, kMinMatch-1);

    if(nmatch >= kMinMatch) lines.emplace_back(line, nmatch);
  } // end for trial
}

// ---------------------------------------------------------------------------
void CandidateLinesGreedy(std::vector<BPPt>::iterator begin,
                          std::vector<BPPt>::iterator end,
                          const std::array<View, 3>& views,
                          FastRand& r,
                          std::vector<std::pair<Ray, int>>& lines)
{
  std::array<BPPt, 4> bests = {*begin, *begin, *begin, *begin};

  int bestScore3D = kMinMatch-1;

  for(int view = 0; view < 3; ++view){
    auto view_end = std::partition(begin, end, [view](const BPPt& p){return p.view == view;});

    int bestScore2D = 0;
    BPPt bestA = *begin;
    BPPt bestB = *begin;

    std::vector<std::pair<BPPt, BPPt>> pairs;
    RandomPairs(begin, view_end, r, pairs);

    for(auto& it: pairs){
      const BPPt a = it.first;
      const BPPt b = it.second;

      if(a.ray.Origin().x == b.ray.Origin().x) continue;

      const Ray line2D(a.ray.Origin(), b.ray.Origin() - a.ray.Origin());

      const int score = LongestChainCount(views, line2D, bestScore2D, {view});
      if(score > bestScore2D){
        bestScore2D = score;
        bestA = a;
        bestB = b;
      }
    }

    // TODO - possible/worthwhile to scan optimize the endpoints in 2D first?

    RandomPairs(view_end, end, r, pairs);

    BPPt bestC = *view_end;
    BPPt bestD = *view_end;

    for(auto& it: pairs){
      const BPPt c = it.first;
      const BPPt d = it.second;
      if(c.view == d.view && c.ray.Origin().x == d.ray.Origin().x) continue;

      const Ray line3D = PointsToRay({bestA, bestB, c, d});

      // Very horizontal lines can cause problems
      if(sqr(line3D.Dir().X()) < 1e-10*(sqr(line3D.Dir().Y() + line3D.Dir().Z()))) continue;

      const int score = LongestChainCount(views, line3D, bestScore3D);
      if(score > bestScore3D){
        bestScore3D = score;
        bestC = c;
        bestD = d;
        bests = {bestA, bestB, bestC, bestD};
      }
    }
  } // end for view

  /*
  std::cout << bestScore3D << " -> ";

  while(true){
    bool any = false;
    for(int i = 0; i < 4; ++i){
      for(auto it = begin; it != end; ++it){
        auto maybe = bests;
        maybe[i] = *it;

        // Avoid cases that give a degenerate solution
        bool ok = true;
        for(int i = 0; i < 3; ++i)
          for(int j = i+1; j < 4; ++j)
            if(maybe[i].ray.Origin().x == maybe[j].ray.Origin().x &&
               maybe[i].view == maybe[j].view) ok = false;
        if(!ok) continue;

        // 2+2+0 or 2+1+1 is OK, 3/4+anything is not
        int nperview[3] = {0, 0, 0};
        for(BPPt& s: maybe) ++nperview[s.view];
        for(int v = 0; v < 3; ++v) if(nperview[v] >= 3) ok = false;
        if(!ok) continue;

        const Ray line3D = PointsToRay(maybe);

        // Very horizontal lines can cause problems
        if(sqr(line3D.Dir().X()) < 1e-10*(sqr(line3D.Dir().Y() + line3D.Dir().Z()))) continue;

        const int score = LongestChainCount(views, line3D, bestScore3D);
        if(score > bestScore3D){
          any = true;
          bestScore3D = score;
          bests = maybe;
        }
      }
    }
    if(!any) break;
  }

  std::cout << bestScore3D << " after opt" << std::endl;
  */

  if(bestScore3D >= kMinMatch)
    lines.emplace_back(PointsToRay(bests), bestScore3D);
}

// ---------------------------------------------------------------------------
void CandidateLinesSeedsIncremental(const std::vector<MyVec>& seeds,
                                    const std::array<View, 3>& views,
                                    std::vector<std::pair<Ray, int>>& lines)
{
  // Must use one of the two new endpoints, but not both
  for(unsigned int i = seeds.size()-2; i < seeds.size(); ++i){
    for(unsigned int j = 0; j < seeds.size()-2; ++j){
      const Ray line3D(seeds[i], seeds[j]-seeds[i]);

      const int score = LongestChainCount(views, line3D, kMinMatch);

      if(score >= kMinMatch){
        lines.emplace_back(line3D, score);
      }
    }
  }
}

// ---------------------------------------------------------------------------
void AddArtTrack(const Ray& line,
                 std::vector<BPPt>::iterator begin,
                 std::vector<BPPt>::iterator end,
                 MyVec r0, MyVec r1,
                 std::vector<recob::Track>* trkcol,
                 art::Assns<recob::Track, recob::Hit>* assns,
                 const art::Event& evt,
                 const art::Handle<std::vector<recob::Hit>>& hits)
{
  // Magic up a pointer to the track we're about to make
  const art::ProductID id = evt.getProductID<std::vector<recob::Track>>("");
  const art::EDProductGetter* pg = evt.productGetter(id);
  const art::Ptr<recob::Track> ptrk(id, trkcol->size(), pg);

  for(auto it = begin; it != end; ++it)
    assns->addSingle(ptrk, art::Ptr<recob::Hit>(hits, it->hitIdx));

  std::vector<geo::Point_t> tps;
  tps.emplace_back(r0.x, r0.y, r0.z);
  tps.emplace_back(r1.x, r1.y, r1.z);

  const recob::TrackTrajectory traj(std::move(tps),
                                    std::vector<geo::Vector_t>(tps.size()), // momenta
                                    std::vector<recob::TrajectoryPointFlags>(tps.size()),
                                    false);

  trkcol->emplace_back(traj, 0, 0, 0, recob::Track::SMatrixSym55(), recob::Track::SMatrixSym55(), trkcol->size()+1);
}

// ---------------------------------------------------------------------------
void AddArtTrack2(const Ray& line,
                  //                 std::vector<BPPt>::iterator begin,
                  //                 std::vector<BPPt>::iterator end,
                  //                 MyVec r0, MyVec r1,
                 std::vector<recob::Track>* trkcol,
                 art::Assns<recob::Track, recob::Hit>* assns,
                 const art::Event& evt,
                 const art::Handle<std::vector<recob::Hit>>& hits)
{
  // Magic up a pointer to the track we're about to make
  const art::ProductID id = evt.getProductID<std::vector<recob::Track>>("");
  const art::EDProductGetter* pg = evt.productGetter(id);
  const art::Ptr<recob::Track> ptrk(id, trkcol->size(), pg);

  //  for(auto it = begin; it != end; ++it)
    //    assns->addSingle(ptrk, art::Ptr<recob::Hit>(hits, it->hitIdx));

  const MyVec r0 = line.Origin() - 10000*line.Dir();
  const MyVec r1 = line.Origin() + 10000*line.Dir();

  std::vector<geo::Point_t> tps;
  tps.emplace_back(r0.x, r0.y, r0.z);
  tps.emplace_back(r1.x, r1.y, r1.z);

  const recob::TrackTrajectory traj(std::move(tps),
                                    std::vector<geo::Vector_t>(tps.size()), // momenta
                                    std::vector<recob::TrajectoryPointFlags>(tps.size()),
                                    false);

  trkcol->emplace_back(traj, 0, 0, 0, recob::Track::SMatrixSym55(), recob::Track::SMatrixSym55(), trkcol->size()+1);
}

// ---------------------------------------------------------------------------
void QuadPts::produce(art::Event& evt)
{
  auto trkcol = std::make_unique<std::vector<recob::Track>>();
  auto assns = std::make_unique<art::Assns<recob::Track, recob::Hit>>();

  art::Handle<std::vector<recob::Hit>> hits;
  evt.getByLabel(fHitLabel, hits);

  std::vector<UnitVec> dirs, perps;
  std::vector<BPPt> pts3d = GetPts3D(*hits, geom, detprop, &dirs, &perps);

  for(int i = 0; i < 3; ++i){
    std::cout << dirs[i] << " " << perps[i] << std::endl;
  }

  std::array<std::vector<BPPt>, 3> pts_by_view;

  std::array<std::vector<MinimalPt>, 3> mpts;

  for(const BPPt& pt: pts3d){
    // TODO ideally GetPts3D would already do this
    pts_by_view[pt.view].push_back(pt);

    mpts[pt.view].emplace_back(pt.z, pt.ray.Origin().x);
  }

  const std::array<View, 3> views = {View(mpts[0], dirs[0], perps[0]),
                                     View(mpts[1], dirs[1], perps[1]),
                                     View(mpts[2], dirs[2], perps[2])};


  // TODO skip any event with a very small view

  //  std::vector<Ray> lines;

  std::vector<BPPt> allBestPts;

  for(int mainView = 0; mainView < 3; ++mainView){

    VPTree2 vp2(mpts[mainView]);

    int bestScore = 0;
    //    Ray bestLine(MyVec(0, 0, 0), MyVec(0, 0, 0));

    //    const int otherView = (mainView+1)%3;
    //    const BPPt other_a = pts_by_view[otherView][0];
    //    const BPPt other_b = pts_by_view[otherView][1];

    const View& view = views[mainView]; // TODO confusing naming

    std::pair<int, int> bestPts;

    const std::vector<MinimalPt>& mptsv = mpts[mainView];
    const unsigned int N = mptsv.size();
    for(unsigned int i = 0; i < N; ++i){
      for(unsigned int j = i+1; j < N; ++j){
        const MinimalLine line(mptsv[i], mptsv[j]);

        const int score = CountClosePoints(view, line);
        const int score2 = vp2.Count(line);
        if(score != score2){
          std::cout << score << " != " << score2 << std::endl;
          //          abort();
        }

        if(score > bestScore){
          bestScore = score;
          bestPts = std::make_pair(i, j);
        }
      }
    }

    const MinimalLine line(mptsv[bestPts.first], mptsv[bestPts.second]);
    std::cout << "m c " << line.dzdx << " " << line.z0 << std::endl;
    std::cout << " -> " << CountClosePoints(view, line) << std::endl;

    std::cout << "vp: " << vp2.Count(line) << std::endl;


    allBestPts.push_back(pts_by_view[mainView][bestPts.first]);
    allBestPts.push_back(pts_by_view[mainView][bestPts.second]);

    /*

    const unsigned int N = main_pts.size();
    for(unsigned int i = 0; i < N; ++i){
      for(unsigned int j = i+1; j < N; ++j){

        const Ray line = PointsToRay({main_pts[i], main_pts[j], other_a, other_b});
        const int score = CountClosePoints(view, line);
        if(score > bestScore){
          bestScore = score;
          bestLine = line;
        }
        //        bestScore = std::max(bestScore, score);

        //        lines.push_back(line);
      } // end for j
    } // end for i
    */

    std::cout << "Best score " << bestScore << std::endl;
    //    std::cout << "  " << CountClosePoints(view, bestLine) << std::endl;
  } // end for mainView


  const Ray bestLine = PointsToRay({allBestPts[0], allBestPts[1], allBestPts[2], allBestPts[3]});

  AddArtTrack2(bestLine,// begin, end, r0, r1,
               trkcol.get(), assns.get(), evt, hits);

  //  std::cout << pts3d.size() << " points -> " << lines.size() << " lines" << std::endl;

  evt.put(std::move(trkcol));
  evt.put(std::move(assns));

  return;

#if 0
  FastRand r;

  std::vector<BPPt>::iterator begin = pts3d.begin();
  std::vector<BPPt>::iterator end = pts3d.end();

  std::vector<MyVec> seeds;

  std::vector<std::pair<Ray, int>> lines;

  while(end-begin >= kMinMatch){ // find as many tracks as possible

    // TODO - keep these up-to-date as we go, somehow
    static std::array<std::vector<MinimalPt>, 3> mpts;
    for(int view = 0; view < 3; ++view) mpts[view].clear();

    for(auto it = begin; it != end; ++it)
      mpts[it->view].emplace_back(it->z, it->ray.Origin().x);


    // Remove points that are so far away from others that they can't possibly
    // be in a valid track.
    static std::array<std::vector<MinimalPt>, 3> mpts_dense;

    for(int view = 0; view < 3; ++view){
      mpts_dense[view].clear();

      for(const MinimalPt& p: mpts[view]){
        bool ok = false;
        for(const MinimalPt& q: mpts[view]){
          if(&q != &p && p.DistSq(q) <= sqr(maxdL)){
            ok = true;
            break;
          }
        } // end for q
        if(ok) mpts_dense[view].push_back(p);
      } // end for p
    } // end for view


    const std::array<View, 3> views = {View(mpts[0], dirs[0], perps[0]),
                                       View(mpts[1], dirs[1], perps[1]),
                                       View(mpts[2], dirs[2], perps[2])};

    CandidateLinesRandom(10000, begin, end, views, r, lines);
    CandidateLinesGreedy(begin, end, views, r, lines);
    CandidateLinesSeedsIncremental(seeds, views, lines);

    std::cout << lines.size() << " candidate lines" << std::endl;
    if(lines.empty()) break;

    while(!lines.empty()){
      // Put the best candidate last and the second-best candidate penultimate
      if(lines.size() >= 2){
        std::nth_element(lines.begin(), lines.end()-2, lines.end(),
                         [](const std::pair<Ray, int>& a,
                            const std::pair<Ray, int>& b)
                         {
                           return a.second < b.second; // Longest last
                         });
      }

      // Update the score of the line that's top of the pile
      lines.back().second = LongestChainCount(views, lines.back().first, kMinMatch);
      // If it's now terrible, discard it
      if(lines.back().second < kMinMatch){lines.pop_back(); continue;}
      // If there's no competition, go with it
      if(lines.size() == 1) break;
      // If it beats the closest competition (even though that one may not yet
      // be updated (can only get smaller)) then we go with it. Otherwise go
      // around the loop and let this candidate get shuffled lower.
      if(lines.back().second >= (&lines.back()-1)->second) break;
    }

    if(lines.empty()) break;

    const int expect = lines.back().second;

    if(expect < kMinMatch) break;

    const Ray line = lines.back().first;
    lines.pop_back();

    const int xcheck = LongestChainApply(begin, end, line, kMinMatch);

    if(xcheck != expect){
      std::cout <<  "** Mismatch between Count: " << expect
                << " and Apply: " << xcheck << " **" << std::endl;
    }

    float L0 = +std::numeric_limits<float>::max();
    float L1 = -std::numeric_limits<float>::max();

    for(auto it = begin; it != end; ++it){
      L0 = std::min(L0, line.ClosestLambdaToRay(it->ray));
      L1 = std::max(    line.ClosestLambdaToRay(it->ray), L1);
    }

    const MyVec r0 = line.Origin() + L0 * line.Dir();
    const MyVec r1 = line.Origin() + L1 * line.Dir();

    AddArtTrack(line, begin, end, r0, r1,
                trkcol.get(), assns.get(), evt, hits);

    seeds.push_back(r0);
    seeds.push_back(r1);

    // Now focus on the subset of points that aren't in this track
    end = begin;
    begin = pts3d.begin();

    std::cout << "Made track with " << xcheck << " hits"
              << ", " << (end-begin) << " remain"
              << std::endl;
  } // end while

  // IDEA - consider leaving the points at the extreme ends of this track
  // to be used in future tracks - ie they could share one point at the
  // ends. Note that the track points currently aren't sorted.

  // IDEA - try extending existing lines by adding just two more points to
  // them.

  // IDEA - scan all 2-point combinations in each view, then try combinations
  // of the best. Can set a bound on the score of the 3D track from the sum of
  // the two views and the best possible from the third view

  evt.put(std::move(trkcol));
  evt.put(std::move(assns));
#endif
}

} // end namespace quad
