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
  int var(3);
  if(argc < 2)
    std::cerr << "p0 <var>?" << std::endl;
  if(1 < argc) var = std::atoi(argv[1]);
  std::cerr << "continue with p0 " << var << std::endl;
  P0<num_t,  true> p0(abs(var));
  P0W<num_t, true> q0(abs(var));
  auto  p1(p0);
  auto  q1(q0);
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
      const auto pn((var < 0 ? q0.next(d) : p0.next(d)) - d);
      // original function lower and higher frequency part, middle is ignored.
      M = (pn + num_t(1) / (var < 0 ? q1.next(num_t(1) / d)
                                    : p1.next(num_t(1) / d)) - d) / num_t(2);
      if(! isfinite(M) || isnan(M)) M = pn;
    }
    std::cout << M << "," << s0 << ", " << s1 << std::endl << std::flush;
  }
  return 0;
}

