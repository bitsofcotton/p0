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
  auto var(1);
  if(argc < 2)
    std::cerr << argv[0] << " <len>?" << std::endl;
  if(1 < argc) var  = std::atoi(argv[1]);
  std::cerr << "continue with " << argv[0] << " " << var << std::endl;
  shrinkMatrix<num_t, P0D<num_t, P0<num_t, idFeeder<num_t> > > > p(P0D<num_t, P0<num_t, idFeeder<num_t> > >(var < 0 ? - var : 3), var < 0 ? 1 : var);
  num_t d(0);
  auto  M(d);
  auto  S(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    const auto S0(S);
    std::cout << D << ", " << (M = p.next(S += d) - S0) << std::endl << std::flush;
  }
  return 0;
}

