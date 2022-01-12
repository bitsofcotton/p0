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
  auto var(3);
  auto rnd(1);
  if(argc < 3)
    std::cerr << argv[0] << " <len>? <rnd?>" << std::endl;
  if(1 < argc) var = std::atoi(argv[1]);
  if(2 < argc) rnd = std::atoi(argv[2]);
  std::cerr << "continue with " << argv[0] << " " << var << " " << rnd << std::endl;
  P0D<num_t, P0<num_t, idFeeder<num_t> > > p(abs(var), 1, rnd);
  P0D<num_t, P0<num_t, deltaFeeder<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > > > q(abs(var), 1, rnd);
  num_t d(0);
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    std::cout << D << ", " << (M = var < 0 ? q.next(d) : p.next(d)) << ", " << s.substr((int)ins.tellg() + 1) << std::endl << std::flush;
  }
  return 0;
}

