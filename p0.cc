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
typedef northPole<num_t, p0_0t> p0_1t;
typedef northPole<num_t, p0_1t> p0_2t;
typedef shrinkMatrix<num_t, p0_2t, true> p0_3t;
typedef northPole<num_t, p0_3t> p0_4t;
typedef northPole<num_t, p0_4t> p0_5t;
typedef compressIllegal<num_t, p0_5t> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  auto var(1);
  auto le(0);
  if(argc < 2)
    std::cerr << argv[0] << " <len>? <logexp>?" << std::endl;
  if(1 < argc) var = std::atoi(argv[1]);
  if(2 < argc) le  = std::atoi(argv[2]);
  std::cerr << "continue with " << argv[0] << " " << var << " " << le << std::endl;
  p0_t  p(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(var < 0 ? - var : 3))), var < 0 ? 1 : var))), le, le < 0 ? num_t(int(10000)) : num_t(int(20)));
  int   t;
  num_t d(t ^= t);
  auto  M(d);
  auto  A(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    A += d; t ++;
    std::cout << D << ", " << (M = p.next(d - A / num_t(t)) + A / num_t(t)) << std::endl << std::flush;
  }
  return 0;
}

