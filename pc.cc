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
  if(argc < 2)
    std::cerr << argv[0] << " <len>?" << std::endl;
  if(1 < argc) var = std::atoi(argv[1]);
  std::cerr << "continue with " << argv[0] << " " << var << std::endl;
  PC<num_t, P0D<num_t, P0<num_t, idFeeder<num_t> > > > p(P0D<num_t, P0<num_t, idFeeder<num_t> > >(abs(var)));
  PC<num_t, P0D<num_t, P0<num_t, deltaFeeder<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > > > > q(P0D<num_t, P0<num_t, deltaFeeder<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > > >(abs(var)));
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(M * d);
    std::cout << D << ", " << (M = var < 0 ? q.next(d) : p.next(d)) << ", " << s.substr((int)ins.tellg() + 1) << std::endl << std::flush;
  }
  return 0;
}

