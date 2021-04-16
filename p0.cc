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
  P0<num_t, true> p(abs(range));
  P0C<num_t, true> q(abs(range));
  num_t d(0);
  auto  s0(d);
  auto  s1(d);
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    const auto bd(d);
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      if(bd != num_t(0) && M != num_t(0)) {
        s0 += (d - bd) - M;
        s1 += (d - bd) * M;
      }
      M = (range < 0 ? q.next(d) : p.next(d)) - d;
    }
    std::cout << M << "," << s0 << ", " << s1 << std::endl << std::flush;
  }
  return 0;
}

