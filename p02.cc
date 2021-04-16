#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <assert.h>

#include "ifloat.hh"
typedef myfloat num_t;
#include "simplelin.hh"
#include "p0.hh"

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int   range(30);
  if(1 < argc) range = std::atoi(argv[1]);
  std::vector<P0<num_t, true> > p;
  std::vector<P0C<num_t, true> > q;
  p.resize(std::atoi(argv[2]) * 2, P0<num_t, true>(abs(range)));
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
            M0 += M[i] * rr[i];
        }
      }
      t ++;
      M0 /= num_t(dd.size());
    }
    std::cout << M0 << "," << s0 << ", " << s1 << std::endl << std::flush;
  }
  return 0;
}

