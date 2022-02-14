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
typedef P0<num_t, idFeeder<num_t> > p0_0t;
typedef shrinkMatrix<num_t, p0_0t, true> p0_1t;
typedef northPole<num_t, p0_1t> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  auto var(1);
  auto le(12);
  if(argc < 2)
    std::cerr << argv[0] << " <len>? <Marctan>?" << std::endl;
  if(1 < argc) var = std::atoi(argv[1]);
  if(2 < argc) le  = std::atoi(argv[2]);
  std::cerr << "continue with " << argv[0] << " " << var << " " << le << std::endl;
  assert(0 <= le);
  p0_t  p(p0_1t(p0_0t(var < 0 ? - var : 3), var < 0 ? 1 : var), pow(num_t(int(2)), num_t(int(le))));
  int   t;
  num_t d(t ^= t);
  auto  M(d);
  auto  A(d);
  auto  S(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    const auto S0(S);
    A += d; t ++;
    std::cout << D << ", " << (M = p.next(S += d - A / num_t(t)) - S0 + A / num_t(t)) << std::endl << std::flush;
  }
  return 0;
}

