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
  int denom(10000);
  int loop(8000);
  if(1 < argc)
    range = std::atoi(argv[1]);
  if(2 < argc)
    denom = std::atoi(argv[2]);
  if(3 < argc)
    loop  = std::atoi(argv[3]);
  std::vector<P0B<num_t> > p;
  p.resize(loop, P0B<num_t>(range));
  std::string s;
  num_t d(0);
  auto  d0(d);
  auto  bd(d);
  auto  MM(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    if(d != bd && bd != num_t(0)) {
      d0 += (d - bd) * MM;
      const auto  d0(log((d - bd) / num_t(denom) + num_t(1)));
            num_t dd(1);
      MM  = num_t(0);
      for(int i = 0; i < p.size(); i ++)
        MM += p[i].next(dd *= d0 / num_t(i + 1));
      MM *= num_t(denom);
    }
    std::cout << d0 << ", " << MM << std::endl << std::flush;
    bd = d;
  }
  return 0;
}

