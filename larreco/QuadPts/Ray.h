#pragma once

#include "MyVec.h"

class Ray
{
public:
  Ray(const MyVec& o, const MyVec& d) : origin(o), dir(d) {}
  Ray(const MyVec& o, const UnitVec& d) : origin(o), dir(d) {}

  float SqrDistanceToRay(const Ray& r) const;

  float ClosestLambdaToRay(const Ray& r) const;
  MyVec ClosestPtToRay(const Ray& r) const;

  const MyVec& Origin() const {return origin;}
  const UnitVec& Dir() const {return dir;}
protected:
  template<class T> static inline T sqr(T x){return x*x;}

  MyVec origin;
  UnitVec dir;
};

// This all needs to be inlined for speed

inline float Ray::SqrDistanceToRay(const Ray& r) const
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

inline float Ray::ClosestLambdaToRay(const Ray& r) const
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

inline MyVec Ray::ClosestPtToRay(const Ray& r) const
{
  return origin + ClosestLambdaToRay(r) * dir;
}
