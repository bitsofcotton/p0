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
  int  var(3);
  if(argc < 2)
    std::cerr << "p0 <var>?" << std::endl;
  if(1 < argc) var = std::atoi(argv[1]);
  std::cerr << "continue with p0 " << var <<  std::endl;
  P0<num_t, continuousFeeder<num_t, continuousFeeder<num_t, linearFeeder<num_t> > > > p(abs(var));
  P0<num_t, continuousFeeder<num_t, continuousFeeder<num_t, arctanFeeder<num_t> > > > q(abs(var));
  num_t d(0);
  auto  s0(d);
  auto  s1(d);
  auto  s2(d);
  auto  s3(d);
  auto  M(d);
  int   t(0);
  while(std::getline(std::cin, s, '\n')) {
    const auto bd(d);
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      if(M != num_t(0)) {
        s0 += (s3 = (d - bd) - (M - d));
        s1 += (s2 = (d - bd) * (M - d));
      }
      // original function lower and higher frequency part, middle is ignored.
      // to make any sub differential be the same meaning, no inverse condition.
      M = var < 0 ? q.next(d) : p.next(d);
      if(t ++ < abs(var)) M = num_t(0);
    }
    std::cout << M - d << "," << s0 << ", " << s1 << ", " << s2 << ", " << s3 << ", " << std::endl << std::flush;
  }
  return 0;
}

