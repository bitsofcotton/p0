#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <assert.h>

#if !defined(_FLOAT_BITS_)
  #include <complex>
  #include <cmath>
  using namespace std;
  typedef long double num_t;
#else
  #include "ifloat.hh"
  template <typename T> using complex = Complex<T>;
# if _FLOAT_BITS_ == 8
    typedef uint8_t myuint;
    typedef int8_t  myint;
    typedef SimpleFloat<myuint, uint16_t, 8, int64_t> num_t;
    #define mybits 8
# elif _FLOAT_BITS_ == 16
    typedef uint16_t myuint;
    typedef int16_t  myint;
    typedef SimpleFloat<myuint, uint32_t, 16, int64_t> num_t;
    #define mybits 16
# elif _FLOAT_BITS_ == 32
    typedef uint32_t myuint;
    typedef int32_t  myint;
    typedef SimpleFloat<myuint, uint64_t, 32, int64_t> num_t;
    #define mybits 32
# elif _FLOAT_BITS_ == 64
    typedef uint64_t myuint;
    typedef int64_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 64>, 64, int64_t> num_t;
    #define mybits 64
# elif _FLOAT_BITS_ == 128
    typedef DUInt<uint64_t, 64> uint128_t;
    typedef Signed<uint128_t, 128> int128_t;
    typedef uint128_t myuint;
    typedef int128_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 128>, 128, int64_t> num_t;
    #define mybits 128
# elif _FLOAT_BITS_ == 256
    typedef DUInt<uint64_t, 64> uint128_t;
    typedef DUInt<uint128_t, 128> uint256_t;
    typedef Signed<uint256_t, 128> int256_t;
    typedef uint256_t myuint;
    typedef int256_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 256>, 256, int64_t> num_t;
    #define mybits 256
# else
#   error cannot handle float
# endif
#endif

#include "simplelin.hh"
#include "p0.hh"

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int range(20);
  if(1 < argc)
    range = std::atoi(argv[1]);
  P0B<num_t> p(range);
  num_t d(0);
  auto  d0(d);
  auto  d1(d);
  auto  d2(d);
  auto  d3(d);
  auto  d4(d);
  auto  bbd(d);
  auto  M(d);
  auto  M0(d);
  int   t(0);
  auto  bet0(t);
  auto  bet1(t);
  while(std::getline(std::cin, s, '\n')) {
    const auto bd(d);
    std::stringstream ins(s);
    ins >> d;
    d0 += (d - bd) * M;
    if(d != bd) {
      const auto dd(d + bd);
      d1 += (d - bd)  * M0 * num_t(bet0);
      d2 += (d - bbd) * M0;
      d3 += (d - bd)  * M0 * num_t(bet1);
      d4 += (d - bbd) * M0;
      M0  = p.next(dd) - dd;
      if(d2 <= d1) {
        bet0 = 0;
        d1   = d2 = num_t(0);
      }
      if(d3 <= d4) {
        bet1 = 0;
        d3   = d4 = num_t(0);
      }
      bet0 ++;
      bet1 ++;
      M   = M0 * num_t(bet0 - bet1);
      if(! isfinite(M) || isnan(M) || t ++ <= range * 2) M = num_t(0);
      bbd = bd;
    }
    std::cout << d0 << ", " << M << std::endl;
  }
  return 0;
}

