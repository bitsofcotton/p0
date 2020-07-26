#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
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

template <typename T> T verbose(const T& d0, const T& b, const bool& flag) {
        T    res(0);
        auto d(d0);
  const auto bb(flag ? b * b : b);
  const auto bn(flag ? b : b * b);
  for(int i = 1; i < 80 && T(1) <= d; i ++) {
    const auto r(int((d - T(int(d / pow(bb, T(i)))) * pow(bb, T(i))) / pow(bb, T(i - 1))));
    res += (flag ? T(int(sqrt(T(r)))) : T(r * r)) * pow(bn, T(i - 1));
    d   -= T(r) * pow(bb, T(i - 1));
  }
  for(int i = 0; i < 80 && T(0) < d; i ++) {
    const auto r(int((d - T(int(d * pow(bb, T(i)))) / pow(bb, T(i))) * pow(bb, T(i + 1))));
    res += (flag ? T(int(sqrt(T(r)))) : T(r * r)) * pow(bn, - T(i + 1));
    d   -= T(r) * pow(bb, - T(i + 1));
  }
  return res;
}

int main(int argc, const char* argv[]) {
  int range(12);
  int base(7);
  int div(30);
  if(1 < argc)
    range = std::atoi(argv[1]);
  if(2 < argc)
    base  = std::atoi(argv[2]);
  if(3 < argc)
    div   = std::atoi(argv[3]);
  std::vector<P0B<num_t> > p;
  p.resize(div, P0B<num_t>(range));
  auto  q(p);
  std::string s;
  num_t d(0);
  auto  d0(d);
  auto  bd(d);
  auto  MM(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      d0 += (d - bd) * MM;
      MM = num_t(0);
      int cnt(0);
      for(int i = 0; i < p.size(); i ++) {
        const auto dd0(d + num_t(base) * num_t(i) / num_t(int(p.size())));
        const auto dd(dd0 < num_t(0) ? - verbose(- dd0, num_t(base), false)
                                     :   verbose(  dd0, num_t(base), false));
        const auto sd(sqrt(dd));
        try {
          const auto lM0(p[i].next(sd) / q[i].next(num_t(1) / sd));
          if(isfinite(lM0) && ! isnan(lM0)) {
            MM += (lM0 < num_t(0) ? - verbose(- lM0, num_t(base), true)
                                  :   verbose(  lM0, num_t(base), true)) - dd0;
            cnt ++;
          }
        } catch(...) {
          ;
        }
      }
      if(cnt)
        MM /= num_t(cnt);
    }
    std::cout << d0 << ", " << MM << std::endl << std::flush;
    bd = d;
  }
  return 0;
}

