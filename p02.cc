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
  P0D<num_t, P0<num_t, idFeeder<num_t> > > p1(abs(var));
  P0D<num_t, P0<num_t, idFeeder<num_t> > > p2(abs(var), 2);
  P0D<num_t, P0<num_t, deltaFeeder<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > > > q1(abs(var));
  P0D<num_t, P0<num_t, deltaFeeder<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > > > q2(abs(var), 2);
  num_t d(int(0));
  auto  M(d);
  auto  bM(d);
  auto  S(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * bM);
    S += d; bM = M;
    std::cout << D << ", " <<
      (M = var < 0 ? q2.next(S) - q1.next(S)
                   : p2.next(S) - p1.next(S)) << ", " <<
      s.substr((int)ins.tellg() + 1) << std::endl << std::flush;
  }
  return 0;
}

