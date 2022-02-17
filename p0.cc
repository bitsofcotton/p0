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
typedef northPole<num_t, p0_1t> p0_2t;
typedef sumChain< num_t, p0_2t> p0_3t;
typedef northPole<num_t, p0_3t> p0_4t;
typedef shrinkMatrix<num_t, p0_4t, true> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  auto var(1);
  auto mabs2(0);
  if(argc < 2)
    std::cerr << argv[0] << " <len>? <Mabs in 2^x>?" << std::endl;
  if(1 < argc) var   = std::atoi(argv[1]);
  if(2 < argc) mabs2 = std::atoi(argv[2]);
  std::cerr << "continue with " << argv[0] << " " << var << " " << mabs2 << std::endl;
  p0_t  p(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(var < 0 ? - var : 3), max(1, var))))), max(1, abs(var)));
  P0<num_t, idFeeder<num_t> > q(3);
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    std::cout << D << ", " << (M = var ? p.next(d) : q.next(d)) << std::endl << std::flush;
  }
  return 0;
}

