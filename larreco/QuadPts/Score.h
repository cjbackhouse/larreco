#pragma once

#include <array>
#include <ostream>

namespace quad
{
  class Score
  {
  public:
    Score() {Reset();}

    int GetMin() const {return std::min(fVal[0], std::min(fVal[1], fVal[2]));}
    int GetTot() const {return fVal[0]+fVal[1]+fVal[2];}
    float GetTotDist() const {return fDist[0]+fDist[1]+fDist[2];}

    void Increment(int idx, int delta = 1) {fVal[idx] += delta;}
    void IncrementDist(int idx, float delta) {fDist[idx] += delta;}
    void Reset() {fVal[0] = fVal[1] = fVal[2] = 0; fDist[0] = fDist[1] = fDist[2] = 0;}

    bool operator<(const Score& s) const
    {
      const int m = GetMin(), sm = s.GetMin();
      if(m != sm) return m < sm;
      const int t = GetTot(), st = s.GetTot();
      if(t != st) return t < st;
      return GetTotDist() > s.GetTotDist(); // NB large dist is bad
    }

    bool operator>(const Score& s) const
    {
      const int m = GetMin(), sm = s.GetMin();
      if(m != sm) return m > sm;
      const int t = GetTot(), st = s.GetTot();
      if(t != st) return t > st;
      return GetTotDist() < s.GetTotDist(); // NB large dist is bad
    }
  protected:
    std::array<int, 3> fVal;
    std::array<float, 3> fDist;
  };

  std::ostream& operator<<(std::ostream& os, const Score& s)
  {
    os << s.GetMin() << " (" << s.GetTot() << ", " << s.GetTotDist() << ")";
    return os;
  }
}
