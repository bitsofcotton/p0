#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>

#include "lieonn.hh"
typedef myfloat num_t;
#include "p0.hh"

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int   var(3);
  if(argc < 2)
    std::cerr << "p0 <var>?" << std::endl;
  if(1 < argc) var = std::atoi(argv[1]);
  std::cerr << "continue with p0 " << var <<  std::endl;
  P0<num_t, linearFeeder<num_t> > p(abs(var));
  P0<num_t, arctanFeeder<num_t> > q(abs(var));
  num_t d(0);
  auto  M(d);
  auto  s0(d);
  auto  s1(d);
  while(std::getline(std::cin, s, '\n')) {
    const auto bd(d);
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      s0 += (d - bd) - M;
      s1 += (d - bd) * M;
      // to make any sub-sequences to be the same meaning, no inverse condition,
      // this causes middle and high frequency parts to be ignored.
      M   = (var < 0 ? q.next(d) : p.next(d)) - d;
    }
    std::cout << M << "," << s0 << ", " << s1 << std::endl << std::flush;
  }
  return 0;
}

