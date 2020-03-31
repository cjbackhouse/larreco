#pragma once

#include "TVector3.h"

#include <cassert>
#include <iostream>

class MyVec;
std::ostream& operator<<(std::ostream& os, const MyVec& v);

class UnitVec;
std::ostream& operator<<(std::ostream& os, const UnitVec& v);

// TVector3 really sucks
class MyVec
{
public:
  MyVec(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}

  MyVec(const TVector3& v) : x(v.X()), y(v.Y()), z(v.Z()) {}

  inline float Mag2() const {return x*x + y*y + z*z;}
  inline float Mag() const {return sqrt(Mag2());}

  inline float Dot(MyVec v) const {return x*v.x + y*v.y + z*v.z;}
  inline float Dot(UnitVec v) const;

  inline MyVec Cross(MyVec v) const {return MyVec(y * v.z - z * v.y,
                                                  z * v.x - x * v.z,
                                                  x * v.y - y * v.x);}

  inline float operator[](int i) const
  {
    if(i == 0) return x;
    if(i == 1) return y;
    if(i == 2) return z;
    abort();
  }

  void Print() const {std::cout << *this << std::endl;}

  inline MyVec operator+(MyVec v) const {return MyVec(x+v.x, y+v.y, z+v.z);}
  inline MyVec& operator+=(MyVec v) {x += v.x; y += v.y; z += v.z; return *this;}
  inline MyVec operator-(MyVec v) const {return MyVec(x-v.x, y-v.y, z-v.z);}
  inline MyVec& operator-=(MyVec v) {x -= v.x; y -= v.y; z -= v.z; return *this;}

  inline MyVec operator*(float m) const {return MyVec(m*x, m*y, m*z);}
  inline MyVec& operator*=(float m) {x *= m; y *= m; z *= m; return *this;}

  float x, y, z;
};

inline MyVec operator*(float m, MyVec v){return v*m;}

inline std::ostream& operator<<(std::ostream& os, const MyVec& v)
{
  os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
  return os;
}

class UnitVec
{
public:
  UnitVec(MyVec v)
  {
    const float m = v.Mag();
    if(m == 0){
      // This is the old behaviour
      x = y = z = 0;
    }
    else{
      x = v.x/m;
      y = v.y/m;
      z = v.z/m;
    }
  }

  inline float X() const {return x;}
  inline float Y() const {return y;}
  inline float Z() const {return z;}

  inline float operator[](int i) const
  {
    if(i == 0) return x;
    if(i == 1) return y;
    if(i == 2) return z;
    abort();
  }

  inline float Dot(UnitVec v) const {return x*v.x + y*v.y + z*v.z;}
  inline float Dot(MyVec v) const {return x*v.x + y*v.y + z*v.z;}

  inline MyVec Cross(MyVec v) const {return MyVec(y * v.z - z * v.y,
                                                  z * v.x - x * v.z,
                                                  x * v.y - y * v.x);}
  inline MyVec Cross(UnitVec v) const {return MyVec(y * v.z - z * v.y,
                                                    z * v.x - x * v.z,
                                                    x * v.y - y * v.x);}

  void Print() const {std::cout << *this << std::endl;}

  inline MyVec operator*(float m) const {return MyVec(m*x, m*y, m*z);}
protected:
  float x, y, z;
};

inline MyVec operator*(float m, UnitVec v){return v*m;}

inline float MyVec::Dot(UnitVec v) const {return x*v.X() + y*v.Y() + z*v.Z();}

inline std::ostream& operator<<(std::ostream& os, const UnitVec& v)
{
  os << "{" << v.X() << ", " << v.Y() << ", " << v.Z() << "}";
  return os;
}
