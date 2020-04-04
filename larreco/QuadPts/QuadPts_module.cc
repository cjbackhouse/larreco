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
#include "lardataobj/RecoBase/SpacePoint.h"
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

  produces<std::vector<recob::SpacePoint>>();
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

  try{
    TDecompLU d(M);
    /*bool ok = */d.Solve(v);
  }
  catch(...){
    std::cout << "Singular matrix in QuadPts::PointsToRay" << std::endl;
    //    M.Print();
    // garbage return is OK, will just score badly
  }

//  std::cout << "OK? " << ok << std::endl;

/*
  const Ray ret(MyVec(0, v[2], v[3]), MyVec(1, v[0], v[1]));

  for(int i = 0; i < 4; ++i){
    const UnitVec d = pts[i].ray.Dir();
    const MyVec o   = pts[i].ray.Origin();
    std::cout << ret.Dir().Cross(d).Dot(o-ret.Origin()) << std::endl;
  }
*/

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

Ray PointsToRay(BPPt a, BPPt b, BPPt c, BPPt d)
{
  return PointsToRay({a, b, c, d});
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

class QuadTree
{
public:
  const float cellSize = 64;//1; // cm

  QuadTree(const std::vector<MinimalPt>& pts)
  {
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

  int Count(MinimalLine line) const throw()
  {
    line.z0 += zoff - xoff*line.dzdx; // correct for point offset

    int stride = 512; // TODO 2^9

    const float angleFactor = sqrt(1+sqr(line.dzdx));
    const float maxdz = maxd * angleFactor;

    static std::vector<int> todo, next;
    //    todo.reserve(1024*1024);
    //    next.reserve(1024*1024);
    todo.clear();
    next.clear();
    todo.push_back(0);

    for(int level = 0; level < 9/*10*/; ++level){
      //      const float Rz = stride*cellSize/sqrt(2) * angleFactor; // cover the corners

      for(int k: todo){
        //        const float x0 = (k/1024+.5)*stride*cellSize;
        //        const float z0 = (k%1024+.5)*stride*cellSize;

        //        const float dz = fabs(line.dzdx * x0 + line.z0 - z0);

        const float xlo = (k/1024  )*stride*cellSize - maxd;
        const float xhi = (k/1024+1)*stride*cellSize + maxd;
        const float zlo = (k%1024  )*stride*cellSize - maxd;
        const float zhi = (k%1024+1)*stride*cellSize + maxd;
        const float zin  = line.dzdx*xlo + line.z0;
        const float zout = line.dzdx*xhi + line.z0;

        // Intersects the (padded) box
        if((zin > zlo || zout > zlo) && (zin < zhi || zout < zhi)){

          //        if(dz < Rz + maxdz){
          for(int dij = 0; dij < 4; ++dij){
            const int di = dij/2;
            const int dj = dij%2;
            //          for(int di = 0; di <= 1; ++di){
            //            for(int dj = 0; dj <= 1; ++dj){
            const int key = 1024*((k/1024)*2+di) + ((k%1024)*2+dj);
            if(fSets[level+1][key]/*.count(key)*/) next.push_back(key);
          } // end for dij
        } // end if
      } // end for k (todo)

      todo.swap(next);
      next.clear();

      stride /= 2;
    } // end for level

    /*
    // Putting all the math together like this is a little faster
    static std::vector<MinimalPt> pts;
    pts.clear();
    for(int k: todo){
      const auto& fpk = fPts[k];
      pts.insert(pts.end(), fpk.begin(), fpk.end());
    }

    int ret = 0;
    for(MinimalPt p: pts){
      if(fabs(line.dzdx * p.x + line.z0 - p.z) < maxdz) ++ret;
    }
    */

    int ret = 0;
    for(int k: todo){
      for(MinimalPt p: fPts[k]){//fPts.find(k)->second){
        const float dz = fabs(line.dzdx * p.x + line.z0 - p.z);
        if(dz < maxdz) ++ret;
      }
    }

    return ret;
  }

protected:
  void Add(MinimalPt pt)
  {
    pt.z += zoff;
    pt.x += xoff;

    int stride = 1;

    // True global position
    const int i = int(pt.x/cellSize);
    const int j = int(pt.z/cellSize);
    fPts[1024*i+j].push_back(pt);

    for(int level = 9; level >= 0; --level){
      const int key = 1024*(i/stride) + j/stride;
      fSets[level][key] = true;//.insert(key);
      stride *= 2;
    }
  }

  float zoff, xoff;

  //  std::array<std::array<bool, 1024*1024>, 10> fSets;
  std::array<std::bitset<1024*1024>, 10> fSets;
  //  std::array<std::unordered_set<int>, 10> fSets;

  //  std::unordered_map<int, std::vector<MinimalPt>> fPts;

  std::array<std::vector<MinimalPt>, 1024*1024> fPts;
};

class View
{
public:
  View(std::vector<MinimalPt>& _pts,
       const UnitVec& _dir,
       const UnitVec& _perp) :
    pts(_pts),
    npts(_pts.size()),
    dir(_dir), perp(_perp),
    vptree(_pts.begin(), _pts.end())
  {
  }

  std::vector<MinimalPt> pts;
  unsigned int npts;

  UnitVec dir, perp;

  VPTree vptree;
};







// ---------------------------------------------------------------------------
int CountClosePoints(const std::vector<MinimalPt>& pts, const MinimalLine& line) throw()
{
  int ret = 0;

  const float angleFactor = sqrt(1+sqr(line.dzdx));
  const float maxdz = maxd * angleFactor;

  for(const MinimalPt& p: pts){
    // Surprisingly fabs() here is substantially faster than sqr()
    const float dz = fabs(line.dzdx*p.x + line.z0 - p.z);

    // Hmm, seems not really
    //      const float dz2 = sqr(line.dzdx*p.x + line.z0 - p.z);

    if(dz < maxdz) ++ret;
  } // end for p

  return ret;
}


// ---------------------------------------------------------------------------
int CountClosePointsBulk(const std::vector<MinimalPt>& pts,
                         const std::vector<MinimalLine>& lines) throw()
{
  int ret = 0;

  for(const MinimalLine& line: lines){
    const float angleFactor = 1+sqr(line.dzdx);
    const float maxdzsq = maxdsq * angleFactor;

    int score = 0;
    for(const MinimalPt& p: pts){
      const float dzsq = sqr(line.dzdx*p.x + line.z0 - p.z);
      if(dzsq < maxdzsq) ++score;
    } // end for p
    ret = std::max(ret, score);
  }

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

  /*
  const MyVec r0 = line.Origin();
  const MyVec r1 = line.Origin() + 10*line.Dir();

  const double x0 = r0.x;
  const double x1 = r1.x;
  const double z0 = r0.Dot(view.perp);
  const double z1 = r1.Dot(view.perp);

  const double m2 = (z1-z0)/(x1-x0);
  const double c2 = z0-m2*x0;
  */

  // dz/dx
  const float m = line.Dir().Dot(view.perp) / line.Dir().X();

  // z(x=0)
  const float c = line.Origin().Dot(view.perp) - m * line.Origin().x;

  //  std::cout << "  m1 c1 " << m << " " << c << " " << m2 << " " << c2 << std::endl;

  MinimalLine ml(m, c);
  return CountClosePoints(view.pts, ml);

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
                  double lambda0, double lambda1,
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

  const MyVec r0 = line.Origin() + lambda0 * line.Dir();
  const MyVec r1 = line.Origin() + lambda1 * line.Dir();

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
void AddSpacePoints(const std::vector<MyVec>& sps,
                    std::vector<recob::SpacePoint>* spscol)
{
  for(MyVec sp: sps){
    const double xyz[3] = {sp.x, sp.y, sp.z};
    const double exyz[6] = {0, 0, 0, 0, 0, 0};
    spscol->emplace_back(xyz, exyz, 0);
  }
}

// ---------------------------------------------------------------------------
struct Result
{
  Result(int s, int a, int b) : score(s), hitA(a), hitB(b) {}
  int score, hitA, hitB;

  bool operator<(const Result& r) const {return score < r.score;}
};

// ---------------------------------------------------------------------------
void PopulateResultsTable(const std::array<std::vector<MinimalPt>, 3>& mpts,
                          std::array<std::vector<Result>, 3>& results)
{
  for(int view = 0; view < 3; ++view){
    const std::vector<MinimalPt>& mptsv = mpts[view];
    const unsigned int N = mptsv.size();
    std::vector<Result>& local_results = results[view];
    local_results.reserve((N*(N+1))/2);

    for(unsigned int i = 0; i < N; ++i){
      for(unsigned int j = i+1; j < N; ++j){
        const MinimalLine line(mptsv[i], mptsv[j]);
        const int score = CountClosePoints(mptsv, line);

        local_results.emplace_back(score, i, j);
      } // end for j
    } // end for i

    // Needs to be sorted so we can safely reject those that can't reach target
    std::sort(local_results.rbegin(), local_results.rend());
  } // end for view
}


// ---------------------------------------------------------------------------
std::vector<BPPt> ExtractClosePointsCopy(const Ray& ray,
                                         const std::array<std::vector<BPPt>, 3>& pts)
{
  // TODO this can be done with 2D math...
  //
  // or for two of the views possibly with a giant bitmask?

  static std::vector<BPPt> ret;
  ret.clear();

  for(int view = 0; view < 3; ++view){
    for(const BPPt& pt: pts[view]){
      // https://en.wikipedia.org/wiki/Skew_lines#Distance
      const MyVec n = ray.Dir().Cross(pt.ray.Dir());
      const MyVec dos = ray.Origin()-pt.ray.Origin();
      const double dsq = sqr(n.Dot(dos))/n.Mag2();

      if(dsq < maxdsq){
        ret.push_back(pt);
      }
    }
  }

  return ret;
}


double ProjectPointToTrackLambda(const Ray& ray, const BPPt& pt)
{
  // https://en.wikipedia.org/wiki/Skew_lines#Distance
  const MyVec dos = pt.ray.Origin()-ray.Origin();
  const MyVec n = ray.Dir().Cross(pt.ray.Dir());
  const MyVec n1 = pt.ray.Dir().Cross(n);
  return dos.Dot(n1) / ray.Dir().Dot(n1);
}

// TODO also score a chisq so that an identical count which is closer to the points is preferred
class Score
{
public:
  Score() {Reset();}

  int GetMin() const {return std::min(fVal[0], std::min(fVal[1], fVal[2]));}
  int GetTot() const {return fVal[0]+fVal[1]+fVal[2];}

  void Increment(int idx) {++fVal[idx];}
  void Reset() {fVal[0] = fVal[1] = fVal[2] = 0;}


  bool operator<(const Score& s) const
  {
    const int m = GetMin(), sm = s.GetMin();
    if(m != sm) return m < sm;
    return GetTot() < s.GetTot();
  }

  bool operator>(const Score& s) const
  {
    const int m = GetMin(), sm = s.GetMin();
    if(m != sm) return m > sm;
    return GetTot() > s.GetTot();
  }
protected:
  std::array<int, 3> fVal;
};

std::ostream& operator<<(std::ostream& os, const Score& s)
{
  os << s.GetMin() << " (" << s.GetTot() << ")";
  return os;
}

// TODO allow one additional hit of each view at the end of track - think about
// symmetry of the algorithm
Score LongestGoodSubsetHelper(const Ray& ray, const std::vector<BPPt>& pts,
                              std::vector<BPPt>* trk_pts,
                              std::vector<BPPt>* reject_pts)
{
  // TODO consider giving BPPt a float priv field and using that

  struct PtInfo
  {
    PtInfo(double l, const BPPt* p) : lambda(l), pt(p) {}

    bool operator<(const PtInfo& i) const {return lambda < i.lambda;}

    // This thing gets sorted, make it as small as possible
    float lambda;
    const BPPt* pt;
  };

  static std::vector<PtInfo> infos;
  //  infos.reserve(pts.size());
  infos.clear();

  for(const BPPt& pt: pts){
    infos.emplace_back(ProjectPointToTrackLambda(ray, pt), &pt);
  }
  std::sort(infos.begin(), infos.end());

  bool inTrack = false;
  auto start = infos.end();
  Score score;
  auto bestStart = infos.end();
  auto bestEnd = infos.end();
  Score bestScore;

  float prevLambda[3] = {infos[0].lambda, infos[0].lambda, infos[0].lambda};

  // Deliberately encounter infos.end() in the loop
  for(auto it = infos.begin(); /*it != infos.end()*/; ++it){
    bool ok = (it != infos.end());

    // Safe because if it is invalid then OK is already false
    for(int view = 0; view < 3; ++view) ok = ok && (it->lambda < prevLambda[view] + maxdL);

    if(ok && !inTrack){
      inTrack = true;
      start = it;
      //      std::cout << "Starting new track candidate at " << it-infos.begin() << std::endl;
    }
    else if(!ok && inTrack){
      inTrack = false;
      //      std::cout << "Track candidate ends at " << it-infos.begin() << " for score " << score << std::endl;
      if(score > bestScore){
        bestScore = score;
        bestStart = start;
        bestEnd = it;
      }
      score.Reset();
    }

    if(it == infos.end()) break;

    prevLambda[it->pt->view] = it->lambda;
    if(inTrack) score.Increment(it->pt->view);
  } // end for it

  if(trk_pts && reject_pts){
    for(auto it = infos.begin(); it != bestStart; ++it) reject_pts->push_back(*it->pt);
    for(auto it = bestStart; it < bestEnd; ++it) trk_pts->push_back(*it->pt);
    for(auto it = bestEnd; it < infos.end(); ++it) reject_pts->push_back(*it->pt);
  }

  return bestScore;
}


std::vector<BPPt> LongestGoodSubset(const Ray& ray, std::vector<BPPt>& pts)
{
  std::vector<BPPt> trk_pts, reject_pts;
  trk_pts.reserve(pts.size());
  reject_pts.reserve(pts.size());
  const Score score = LongestGoodSubsetHelper(ray, pts, &trk_pts, &reject_pts);
  if(score.GetMin() < kMinMatch) return {};
  pts.swap(reject_pts);
  return trk_pts;
}

Score CountLongestGoodSubset(const Ray& ray, const std::vector<BPPt>& pts)
{
  return LongestGoodSubsetHelper(ray, pts, 0, 0);
}

// ---------------------------------------------------------------------------
// NB mutates results
Ray BestRay(std::array<std::vector<Result>, 3>& results,
            const std::array<View, 3>& views,
            const std::array<std::vector<BPPt>, 3>& pts_by_view)
{
  Score bestScore;
  Ray bestRay(MyVec(0, 0, 0), MyVec(0, 0, 0));

  FastRand r;

  // All the ways to have 2 views (order doesn't matter) and then the third (special) view
  const int viewIdxs[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};

  int nNoProg = 0;

  while(nNoProg < 1e5){//1e6){
    for(int vi = 0; vi < 3; ++vi){
      const int viewA = viewIdxs[vi][0];
      const int viewB = viewIdxs[vi][1];
      const int viewC = viewIdxs[vi][2];

      const std::vector<BPPt>& ptsA = pts_by_view[viewA];
      const std::vector<BPPt>& ptsB = pts_by_view[viewB];

      if(results[viewA].empty() || results[viewB].empty()) return bestRay;

      const Result resA = results[viewA][r.xorshift32()%results[viewA].size()];
      const Result resB = results[viewB][r.xorshift32()%results[viewB].size()];

      const Ray ray = PointsToRay(ptsA[resA.hitA], ptsA[resA.hitB],
                                  ptsB[resB.hitA], ptsB[resB.hitB]);

      const Score score3d = CountLongestGoodSubset(ray, ExtractClosePointsCopy(ray, pts_by_view));
      /*
      const int scoreC = CountClosePoints(views[viewC], ray);
      const int totScore = resA.score + resB.score + scoreC;
      const int minScore = std::min(std::min(resA.score, resB.score), scoreC);


      if(minScore > bestMinScore ||
         (minScore == bestMinScore && totScore > bestTotScore)){


        std::cout << "*** new best score " << minScore << std::endl;
        std::cout << "  " << resA.score << " + " << resB.score << " + " << scoreC << std::endl;
        bestMinScore = minScore;
        bestTotScore = totScore;
      */

      if(score3d > bestScore){//bestMinScore){
        std::cout << "*** new best score " << score3d << std::endl;

        bestScore = score3d;

        bestRay = ray;
        nNoProg = 0;

        // TODO figure out how to do this with lower_bound et al
        for(int view = 0; view < 3; ++view){
          while(!results[view].empty() && results[view].back().score < bestScore.GetMin()) results[view].pop_back();
          if(results[view].empty()) return bestRay;
        }

        std::cout << " now " << results[viewA].size() << " " << results[viewB].size() << " " << results[viewC].size() << std::endl;
        // TODO potentially exhaustive search once the product of the two smallest is minimized
      }
      else{
        ++nNoProg;
        if(nNoProg%10000 == 0) std::cout << "No prog: " << nNoProg << std::endl;
      }
    } // end for vi
  } // end while

  return bestRay;
}

// ---------------------------------------------------------------------------
std::vector<BPPt> ExtractClosePoints(const Ray& ray, std::vector<BPPt>& pts)
{
  std::vector<BPPt> ret;

  std::vector<BPPt> rejects;
  rejects.reserve(pts.size());

  for(const BPPt& pt: pts){
    // https://en.wikipedia.org/wiki/Skew_lines#Distance
    const MyVec n = ray.Dir().Cross(pt.ray.Dir());
    const MyVec dos = ray.Origin()-pt.ray.Origin();
    const double dsq = sqr(n.Dot(dos))/n.Mag2();

    if(dsq < maxdsq){
      ret.push_back(pt);
    }
    else{
      rejects.push_back(pt);
    }
  }

  pts.swap(rejects);

  return ret;
}

// ---------------------------------------------------------------------------
std::vector<MyVec> ProjectPoints(const Ray& ray, const std::vector<BPPt>& pts,
                                 double& lambda0, double& lambda1)
{
  std::vector<MyVec> ret;

  lambda0 = +std::numeric_limits<double>::infinity();
  lambda1 = -std::numeric_limits<double>::infinity();

  std::vector<BPPt> rejects;
  rejects.reserve(pts.size());

  for(const BPPt& pt: pts){
    // https://en.wikipedia.org/wiki/Skew_lines#Distance
    const UnitVec n = ray.Dir().Cross(pt.ray.Dir());
    const MyVec dos = ray.Origin()-pt.ray.Origin();

    const MyVec n2 = ray.Dir().Cross(n);
    const double lambda_pt = dos.Dot(n2) / pt.ray.Dir().Dot(n2);
    const MyVec closest = pt.ray.Origin() + lambda_pt * pt.ray.Dir();
    ret.push_back(closest);

    const double lambda_trk = ProjectPointToTrackLambda(ray, pt);
    lambda0 = std::min(lambda0, lambda_trk);
    lambda1 = std::max(lambda1, lambda_trk);
  }

  return ret;
}

// ---------------------------------------------------------------------------
void QuadPts::produce(art::Event& evt)
{
  auto trkcol = std::make_unique<std::vector<recob::Track>>();
  auto assns = std::make_unique<art::Assns<recob::Track, recob::Hit>>();

  auto spscol = std::make_unique<std::vector<recob::SpacePoint>>();

  art::Handle<std::vector<recob::Hit>> hits;
  evt.getByLabel(fHitLabel, hits);

  std::vector<UnitVec> dirs, perps;
  std::vector<BPPt> pts3d = GetPts3D(*hits, geom, detprop, &dirs, &perps);

  for(int i = 0; i < 3; ++i){
    std::cout << dirs[i] << " " << perps[i] << std::endl;
  }

  while(true){
    std::array<std::vector<BPPt>, 3> pts_by_view;

    std::array<std::vector<MinimalPt>, 3> mpts;

    for(const BPPt& pt: pts3d){
      // TODO ideally GetPts3D would already do this
      pts_by_view[pt.view].push_back(pt);
      mpts[pt.view].emplace_back(pt.z, pt.ray.Origin().x);
    }

    // Skip any event with a very small view
    if(mpts[0].size() < kMinMatch ||
       mpts[1].size() < kMinMatch ||
       mpts[2].size() < kMinMatch) break;

    const std::array<View, 3> views = {View(mpts[0], dirs[0], perps[0]),
                                       View(mpts[1], dirs[1], perps[1]),
                                       View(mpts[2], dirs[2], perps[2])};

    std::array<std::vector<Result>, 3> results;
    PopulateResultsTable(mpts, results);
    if(results[0].empty() || results[1].empty() || results[2].empty()) break;

    const Ray bestRay = BestRay(results, views, pts_by_view);

    std::vector<BPPt> trkpts = ExtractClosePoints(bestRay, pts3d);
    std::cout << "Found " << trkpts.size() << " 3D points close to track" << std::endl;
    const unsigned int nClose = trkpts.size();
    const std::vector<BPPt> goodtrkpts = LongestGoodSubset(bestRay, trkpts);

    std::cout << goodtrkpts.size() << " points good of " << nClose << std::endl;
    if(goodtrkpts.empty()){
      std::cout << "No run in track exceeded " << kMinMatch << std::endl;
      break;

      /*
      // Show the rejected ray
      AddArtTrack2(bestRay, -1e6, +1e6,
                   trkcol.get(), assns.get(), evt, hits);

      break;
      */
    }

    // Put the rejects back in the heap
    pts3d.insert(pts3d.end(), trkpts.begin(), trkpts.end());

    double lambda0, lambda1;
    const std::vector<MyVec> sps = ProjectPoints(bestRay, goodtrkpts,
                                                 lambda0, lambda1);

    AddArtTrack2(bestRay, lambda0, lambda1,
                 trkcol.get(), assns.get(), evt, hits);

    AddSpacePoints(sps, spscol.get());

  } // end for itrk

  evt.put(std::move(trkcol));
  evt.put(std::move(assns));

  evt.put(std::move(spscol));

  // TODO could already require no big 2D gaps when constructing the score
  // table

  // TODO require tracks to overlap at least a little in x. But one stray hit
  // can trick that.

  // TODO vastly reduce number of 2D tracks by requiring that none are subsets
  // of each others' hits.
}

} // end namespace quad
