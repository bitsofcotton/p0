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
  int   range(30);
  if(1 < argc) range = std::atoi(argv[1]);
  std::vector<P0B<num_t, true> > p;
  std::vector<P0C<num_t, true> > q;
  p.resize(std::atoi(argv[2]) * 2, P0B<num_t, true>(abs(range)));
  q.resize(std::atoi(argv[2]) * 2, P0C<num_t, true>(abs(range)));
  num_t d(0);
  auto  s0(d);
  auto  s1(d);
  auto  M0(d);
  std::vector<num_t> dd(std::atoi(argv[2]), num_t(0));
  auto  rr(dd);
  auto  M(dd);
  int   t(0);
  while(std::getline(std::cin, s, '\n')) {
    const auto bd0(d);
    std::stringstream ins(s);
    ins >> d;
    if(d != bd0) {
      if(bd0 != num_t(0) && M0 != num_t(0)) {
        s0 += (d - bd0) - M0;
        s1 += (d - bd0) * M0;
      }
      const auto bd(dd);
      for(int i = 0; i < dd.size(); i ++) {
        dd[i] += (d - bd0) * rr[i];
        rr[i] += num_t(arc4random_uniform(0x10000) + arc4random_uniform(0x10000) - 0x8000 * 2) / num_t(0x8000);
        if(dd[i] != num_t(0)) {
          M[i] = (range < 0 ? q[(t & 1) * dd.size() + i].next(dd[i])
                            : p[(t & 1) * dd.size() + i].next(dd[i])) - dd[i];
          if(num_t(0) < (dd[i] - bd[i]) * M[i])
            M0 += M[i];
        }
      }
      M0 /= num_t(dd.size());
    }
    std::cout << M0 << "," << s0 << ", " << s1 << std::endl << std::flush;
  }
  return 0;
}

