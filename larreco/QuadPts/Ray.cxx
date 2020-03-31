#include "Ray.h"

#if 0
/*
double Ray::DistanceToRay(const Ray& r) const
{
  const MyVec bd = dir.Cross(r.dir); // perp to both

  if(bd.Mag2() == 0){
    std::cout << "BD = 0" << std::endl;
    return spine.SqrDistPt(pt); // TODO!
  }

  const MyVec n2 = pt.dir.Cross(bd);
  const double lambda = (pt.origin - spine.Start()).Dot(n2) / spine.Dir().Dot(n2);

    const MyVec bdmu = pt.dir.Cross(spine.Dir()); // perp to both
    const MyVec n2mu = spine.Dir().Cross(bdmu);
    const double mu = (spine.Start() - pt.origin).Dot(n2mu) / pt.dir.Dot(n2mu);

    const MyVec trkpt = Eval(lambda);
    const MyVec closest = pt.origin + mu * pt.dir;

    if(isnan((trkpt-closest).Mag2())){
      std::cout << "NAN? " << lambda << " " << mu << std::endl;
      return 0; // ....

      std::cout << "NAN? " << std::endl;
      trkpt.Print();
      std::cout << mu << std::endl;
      closest.Print();
      abort();
    }

    return (trkpt-closest).Mag2();

  return 0;
}
*/




double Ray::SqrDistanceToRay(const Ray& r) const
{
  const MyVec bd = dir.Cross(r.dir);

  // https://en.wikipedia.org/wiki/Skew_lines#Distance
  if(bd.Mag2() == 0){
    // Special case for parallel lines
    const MyVec d = r.origin - origin;
    const MyVec bd = d - dir.Dot(d)*dir;
    if(bd.Mag2() == 0) abort();
    return sqr(bd.Dot(d)) / bd.Mag2();
  }

  return sqr(bd.Dot(r.origin - origin)) / bd.Mag2();
}

double Ray::ClosestLambdaToRay(const Ray& r) const
{
  // https://en.wikipedia.org/wiki/Skew_lines#Nearest_Points
  const MyVec bd = dir.Cross(r.dir); // perp to both
  const MyVec n2 = r.dir.Cross(bd);
  if(dir.Dot(n2) == 0){
    // TODO - figure out what we should do in this case
    return 0;
  }
  else{
    return (r.origin - origin).Dot(n2) / dir.Dot(n2);
  }
}

MyVec Ray::ClosestPtToRay(const Ray& r) const
{
  return origin + ClosestLambdaToRay(r) * dir;
}

#endif
