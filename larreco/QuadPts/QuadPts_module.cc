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
const float maxd = .1;//25; // 1mm laterally
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

  return Ray(MyVec(0, v[2], v[3]), MyVec(1, v[0], v[1]));
}

Ray PointsToRay(BPPt a, BPPt b, BPPt c, BPPt d)
{
  return PointsToRay({a, b, c, d});
}

struct MinimalPt
{
  MinimalPt(float _z, float _x, int _n) : z(_z), x(_x), nTrk(_n) {}

  float DistSq(const MinimalPt& p) const {return sqr(z-p.z) + sqr(x-p.x);}

  float z, x;
  int nTrk;
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
    dir(_dir), perp(_perp)
  {
  }

  std::vector<MinimalPt> pts;
  unsigned int npts;

  UnitVec dir, perp;
};


// ---------------------------------------------------------------------------
int CountClosePoints(const std::vector<MinimalPt>& pts,
                     const MinimalLine& line) throw()
{
  int ret = 0;

  const float angleFactor = 1+sqr(line.dzdx);
  const float maxdzsq = maxdsq * angleFactor;

  for(const MinimalPt& p: pts){
    const float dzsq = sqr(line.dzdx*p.x + line.z0 - p.z);

    if(dzsq < maxdzsq) ret += 2-p.nTrk;
  } // end for p

  return ret;
}

typedef float v4f __attribute__ ((vector_size (16)));
typedef int   v4i __attribute__ ((vector_size (16)));

// ---------------------------------------------------------------------------
int CountClosePoints_vec(const std::vector<MinimalPt>& pts,
                         const MinimalLine& line) throw()
{
  // NB putting a const on any of these variables seems to mis-compile(!)

  v4i rets = {0, 0, 0, 0};

  const float angleFactor = 1+sqr(line.dzdx);
  const float maxdzsq = maxdsq * angleFactor;

  const unsigned int N = pts.size();
  //  std::cout << N << std::endl;
  for(unsigned int i = 0; i+3 < N; i += 4){
    v4f zs = {pts[i].z, pts[i+1].z, pts[i+2].z, pts[i+3].z};
    v4f xs = {pts[i].x, pts[i+1].x, pts[i+2].x, pts[i+3].x};
    v4i nTrks = {pts[i].nTrk, pts[i+1].nTrk, pts[i+2].nTrk, pts[i+3].nTrk};

    v4f dzsqs = sqr(line.dzdx*xs + line.z0 - zs);

    rets += (dzsqs < maxdzsq) ? 2-nTrks : 0;
  } // end for p

  // TODO extra loop to mop up the last up-to-three points

  const int sumret = rets[0]+rets[1]+rets[2]+rets[3];

  /*
  const int truth = CountClosePoints(pts, line);
  if(sumret > truth || sumret+6 < truth){
    std::cout << sumret << " vs " << truth << std::endl;
    abort();
  }
  */

  return sumret;
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
      if(dzsq < maxdzsq) score += 2-p.nTrk;
    } // end for p
    ret = std::max(ret, score);
  }

  return ret;
}


// ---------------------------------------------------------------------------
int CountClosePoints(const View& view, const Ray& line) throw()
{
  // dz/dx
  const float m = line.Dir().Dot(view.perp) / line.Dir().X();

  // z(x=0)
  const float c = line.Origin().Dot(view.perp) - m * line.Origin().x;

  MinimalLine ml(m, c);
  return CountClosePoints_vec(view.pts, ml);
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
    local_results.reserve((N*(N-1))/2);

    for(unsigned int i = 0; i < N; ++i){
      for(unsigned int j = i+1; j < N; ++j){
        const MinimalLine line(mptsv[i], mptsv[j]);
        const int score = CountClosePoints(mptsv, line);

        local_results.emplace_back(score, i, j);
      } // end for j
    } // end for i

    std::cout << view << " " << mpts[view].size() << " " << local_results.size() << std::endl;

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

  void Increment(int idx, int delta = 1) {fVal[idx] += delta;}
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
  if(pts.empty()) return Score();

  static std::vector<const BPPt*> ppts;
  ppts.clear();
  for(const BPPt& pt: pts){
    pt.privF = ProjectPointToTrackLambda(ray, pt);
    ppts.push_back(&pt);
  }

  // Oddly stable_sort is quite a bit faster than regular sort
  std::stable_sort(ppts.begin(), ppts.end(),
                   [](const BPPt* a, const BPPt* b){return a->privF < b->privF;});

  bool inTrack = false;
  auto start = ppts.end();
  Score score;
  auto bestStart = ppts.end();
  auto bestEnd = ppts.end();
  Score bestScore;

  float prevLambda[3] = {ppts[0]->privF, ppts[0]->privF, ppts[0]->privF};

  // Deliberately encounter ppts.end() in the loop
  for(auto it = ppts.begin(); /*it != ppts.end()*/; ++it){
    bool ok = (it != ppts.end());

    const BPPt& pt = ok ? **it : **ppts.begin();

    // Safe because if 'it' is invalid then OK is already false
    for(int view = 0; view < 3; ++view) ok = ok && (pt.privF < prevLambda[view] + maxdL);

    if(ok && !inTrack){
      inTrack = true;
      start = it;
      //      std::cout << "Starting new track candidate at " << it-ppts.begin() << std::endl;
    }
    else if(!ok && inTrack){
      inTrack = false;
      //      std::cout << "Track candidate ends at " << it-ppts.begin() << " for score " << score << std::endl;
      if(score > bestScore){
        bestScore = score;
        bestStart = start;
        bestEnd = it;
      }
      score.Reset();
    }

    if(it == ppts.end()) break;

    prevLambda[pt.view] = pt.privF;
    if(inTrack) score.Increment(pt.view, 2-pt.nTrk);
  } // end for it

  if(trk_pts && reject_pts){
    for(auto it = ppts.begin(); it != bestStart; ++it) reject_pts->push_back(**it);
    for(auto it = bestStart; it < bestEnd; ++it) trk_pts->push_back(**it);
    for(auto it = bestEnd; it < ppts.end(); ++it) reject_pts->push_back(**it);
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
void RandomIndices(FastRand& r, int N, int& i1, int& i2)
{
  do{
    i1 = r.xorshift32()%N;
    i2 = r.xorshift32()%N;
  } while(i1 == i2);

  if(i1 > i2) std::swap(i1, i2);
}

// ---------------------------------------------------------------------------
Ray BestRay(const std::array<View, 3>& views,
            const std::array<std::vector<BPPt>, 3>& pts_by_view)
{
  Score bestScore;
  Ray bestRay(MyVec(0, 0, 0), MyVec(0, 0, 0));

  FastRand r;

  // All the ways to have 2 views (order doesn't matter) and then the third (special) view
  const int viewIdxs[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};

  // TODO this could be a more efficient structure, but doesn't seem to cost
  // much
  std::array<std::vector<std::vector<int>>, 3> cache;
  for(int view = 0; view < 3; ++view){
    const unsigned int N = pts_by_view[view].size();
    cache[view].resize(N);
    for(unsigned int i = 0; i < N; ++i) cache[view][i].resize(N);
  }

  int nAttempts = 0;
  int lastProg = -1;

  while(nAttempts < std::max(10*1000, 3*lastProg)){
    for(int vi = 0; vi < 3; ++vi){
      ++nAttempts;

      const int viewA = viewIdxs[vi][0];
      const int viewB = viewIdxs[vi][1];
      const int viewC = viewIdxs[vi][2];

      const std::vector<BPPt>& ptsA = pts_by_view[viewA];
      const std::vector<BPPt>& ptsB = pts_by_view[viewB];

      if(ptsA.empty() || ptsB.empty()) return bestRay;

      int iA1, iA2, iB1, iB2;
      RandomIndices(r, ptsA.size(), iA1, iA2);
      RandomIndices(r, ptsB.size(), iB1, iB2);

      int nA = cache[viewA][iA1][iA2];
      int nB = cache[viewB][iB1][iB2];
      // already evaluated this 3D line
      if(nA > 0 && nB > 0) continue;

      // Guaranteed to be worse than current best
      if(nA > 0 && nA < bestScore.GetMin()) continue;
      if(nB > 0 && nB < bestScore.GetMin()) continue;

      // Evaluate and compare to current best
      if(nA == 0){
        const MinimalPt mptA1 = views[viewA].pts[iA1];
        const MinimalPt mptA2 = views[viewA].pts[iA2];

        int nA = CountClosePoints_vec(views[viewA].pts, MinimalLine(mptA1, mptA2));
        cache[viewA][iA1][iA2] = nA;
        if(nA < bestScore.GetMin()) continue;
      }

      if(nB == 0){
        const MinimalPt mptB1 = views[viewB].pts[iB1];
        const MinimalPt mptB2 = views[viewB].pts[iB2];

        int nB = CountClosePoints_vec(views[viewB].pts, MinimalLine(mptB1, mptB2));
        cache[viewB][iB1][iB2] = nB;
        if(nB < bestScore.GetMin()) continue;
      }

      const Ray ray = PointsToRay(ptsA[iA1], ptsA[iA2], ptsB[iB1], ptsB[iB2]);

      // Would third view fall below limit?
      const int nC = CountClosePoints(views[viewC], ray);
      if(nC < bestScore.GetMin()) continue;
      // TODO is it worthwhile reusing any of the info from these calls for
      // ExtractClosePoints?

      const Score score3d = CountLongestGoodSubset(ray, ExtractClosePointsCopy(ray, pts_by_view));

      if(score3d > bestScore){
        lastProg = nAttempts;

        std::cout << "*** new best score " << score3d << std::endl;

        bestScore = score3d;

        bestRay = ray;
      }
      else{
        if(nAttempts%10000 == 0) std::cout << "Attempts / last: " << nAttempts << " " << lastProg << std::endl;
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
struct FinalTrack
{
  FinalTrack(const Ray& r, float l0, float l1) : ray(r), lambda0(l0), lambda1(l1) {}
  Ray ray;
  float lambda0, lambda1;
};

// ---------------------------------------------------------------------------
FinalTrack ToFinalTrack(const Ray& ray, const std::vector<BPPt>& pts)
{
  float lambda0 = +std::numeric_limits<double>::infinity();
  float lambda1 = -std::numeric_limits<double>::infinity();

  for(const BPPt& pt: pts){
    const float lambda_trk = ProjectPointToTrackLambda(ray, pt);
    lambda0 = std::min(lambda0, lambda_trk);
    lambda1 = std::max(lambda1, lambda_trk);
  }

  return FinalTrack(ray, lambda0, lambda1);
}

// ---------------------------------------------------------------------------
bool ProjectPoint(const std::vector<FinalTrack>& trks, const BPPt& pt,
                  MyVec& ret)
{
  ret = MyVec(0, 0, 0);
  float totW = 0;

  for(const FinalTrack& trk: trks){
    const float lambda_trk = ProjectPointToTrackLambda(trk.ray, pt);
    if(lambda_trk < trk.lambda0 || lambda_trk > trk.lambda1) continue;

    // https://en.wikipedia.org/wiki/Skew_lines#Distance
    const UnitVec n = trk.ray.Dir().Cross(pt.ray.Dir());
    const MyVec dos = trk.ray.Origin()-pt.ray.Origin();

    const double d = fabs(n.Dot(dos));
    if(d > maxd) continue;

    const MyVec n2 = trk.ray.Dir().Cross(n);
    const double lambda_pt = dos.Dot(n2) / pt.ray.Dir().Dot(n2);
    const MyVec closest = pt.ray.Origin() + lambda_pt * pt.ray.Dir();

    const float w = maxd-d; // TODO what should the weighting function be?
    ret += w*closest;
    totW += w;
  }

  if(totW == 0) return false;

  ret *= 1./totW;
  return true;
}

// ---------------------------------------------------------------------------
void AddArtTrack(const FinalTrack& trk,
                 //                 std::vector<BPPt>::iterator begin,
                 //                 std::vector<BPPt>::iterator end,
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

  const MyVec r0 = trk.ray.Origin() + trk.lambda0 * trk.ray.Dir();
  const MyVec r1 = trk.ray.Origin() + trk.lambda1 * trk.ray.Dir();

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
void AddSpacePoint(const MyVec& sp,
                   std::vector<recob::SpacePoint>* spscol)
{
  const double xyz[3] = {sp.x, sp.y, sp.z};
  const double exyz[6] = {0, 0, 0, 0, 0, 0};
  spscol->emplace_back(xyz, exyz, 0);
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
  const std::vector<BPPt> pts3d_all = pts3d;

  for(int i = 0; i < 3; ++i){
    std::cout << dirs[i] << " " << perps[i] << std::endl;
  }

  std::vector<FinalTrack> finals;

  while(true){
    std::array<std::vector<BPPt>, 3> pts_by_view;

    std::array<std::vector<MinimalPt>, 3> mpts;

    for(const BPPt& pt: pts3d){
      // TODO ideally GetPts3D would already do this
      pts_by_view[pt.view].push_back(pt);
      mpts[pt.view].emplace_back(pt.z, pt.ray.Origin().x, pt.nTrk);
    }

    // Skip any event with a very small view
    if(mpts[0].size() < kMinMatch ||
       mpts[1].size() < kMinMatch ||
       mpts[2].size() < kMinMatch) break;

    const std::array<View, 3> views = {View(mpts[0], dirs[0], perps[0]),
                                       View(mpts[1], dirs[1], perps[1]),
                                       View(mpts[2], dirs[2], perps[2])};

    const Ray bestRay = BestRay(views, pts_by_view);

    std::vector<BPPt> trkpts = ExtractClosePoints(bestRay, pts3d);
    std::cout << "Found " << trkpts.size() << " 3D points close to track" << std::endl;
    const unsigned int nClose = trkpts.size();
    std::vector<BPPt> goodtrkpts = LongestGoodSubset(bestRay, trkpts);

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

    // And selected points that aren't in enough tracks yet
    for(BPPt& pt: goodtrkpts){
      ++pt.nTrk;
      if(pt.nTrk < 2) pts3d.push_back(pt);
    }

    finals.push_back(ToFinalTrack(bestRay, goodtrkpts));
  } // end for itrk

  for(const FinalTrack& t: finals){
    AddArtTrack(t, trkcol.get(), assns.get(), evt, hits);
  }

  for(const BPPt& pt: pts3d_all){
    MyVec sp(0, 0, 0);
    if(ProjectPoint(finals, pt, sp)) AddSpacePoint(sp, spscol.get());
  }


  evt.put(std::move(trkcol));
  evt.put(std::move(assns));

  evt.put(std::move(spscol));

  // TODO could already require no big 2D gaps when constructing the score
  // table

  // TODO vastly reduce number of 2D tracks by requiring that none are subsets
  // of each others' hits / have similar gradient and intercept

  // TODO write out rejected hits
}

} // end namespace quad
