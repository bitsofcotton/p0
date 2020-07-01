#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include "assert.h"

#if !defined(_FLOAT_BITS_)
  #include <complex>
  #include <cmath>
  using namespace std;
  typedef long double num_t;
#else
  #include "ifloat.hh"
  template <typename T> using complex = Complex<T>;
# if _FLOAT_BITS_ == 32
    typedef uint32_t myuint;
    typedef int32_t  myint;
    typedef SimpleFloat<myuint, uint64_t, 32, int64_t> num_t;
# elif _FLOAT_BITS_ == 64
    typedef uint64_t myuint;
    typedef int64_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 64>, 64, int64_t> num_t;
# elif _FLOAT_BITS_ == 128
    typedef DUInt<uint64_t, 64> uint128_t;
    typedef Signed<uint128_t, 128> int128_t;
    typedef uint128_t myuint;
    typedef int128_t  myint;
    typedef SimpleFloat<myuint, DUInt<myuint, 128>, 128, int64_t> num_t;
# else
#   error cannot handle float
# endif
#endif

#include "simplelin.hh"
#include "p0.hh"

int main(int argc, const char* argv[]) {
  int range(30);
  int loop(200);
  int shift(16);
  if(1 < argc)
    range = std::atoi(argv[1]);
  if(2 < argc)
    loop  = std::atoi(argv[2]);
  if(3 < argc)
    shift = std::atoi(argv[3]);
  std::vector<P0B<num_t> > p;
  p.resize(loop ? loop : 1, P0B<num_t>(range));
  std::string s;
  num_t d(0);
  auto  d0(d);
  auto  bd(d);
  auto  MM(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      const auto delta(d - bd);
      d0 += delta * MM;
      if(loop) {
        auto dd(delta);
        MM  = num_t(1);
        for(int i = 0; i < p.size(); i ++) {
          MM += p[i].next(dd);
          dd *= delta / num_t(i + 2);
        }
#if !defined(_FLOAT_BITS_)
        MM = MM < num_t(0) ? - log(- MM) : log(MM);
#else
        MM = MM < num_t(0) ? - floor(log(- MM)) : floor(log(MM));
#endif
      } else
        MM = p[0].next(delta);
    }
#if !defined(_FLOAT_BITS_)
    std::cout << d0 * pow((long double)2, - (long double)shift) << ", " << MM << std::endl;
#else
    std::cout << myint((d0 >> (long long)shift).operator myuint()) << ", " << myint(MM.operator myuint()) << std::endl;
#endif
    bd = d;
  }
  return 0;
}

