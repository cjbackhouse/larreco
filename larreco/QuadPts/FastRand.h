#ifndef MYPMA_FASTRAND_H
#define MYPMA_FASTRAND_H

// https://en.wikipedia.org/wiki/Xorshift
class FastRand
{
public:
  FastRand(uint32_t seed = 123456)
  {
    a = seed;
  }

  uint32_t xorshift32()
  {
    /* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
    uint32_t x = a;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return a = x;
  }

protected:
  uint32_t a;
};

#endif
