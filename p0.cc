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
typedef shrinkMatrix<num_t, p0_2t> p0_3t;
typedef northPole<num_t, p0_3t> p0_4t;
typedef northPole<num_t, p0_4t> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  auto var(1);
  if(argc < 2)
    std::cerr << argv[0] << " <len>?" << std::endl;
  if(1 < argc) var  = std::atoi(argv[1]);
  std::cerr << "continue with " << argv[0] << " " << var << std::endl;
  p0_t p0(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(var < 0 ? - var : 3))), var < 0 ? 1 : var)));
  int   t;
  num_t zero(t ^= t);
  auto  d(zero);
  auto  M(d);
  std::vector<num_t> h;
  h.resize(abs(var) * 2 + 3, d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    if(d == zero) goto next;
    h[(t ++) % h.size()] = d;
    if(t < h.size()) goto next;
    {
      auto A(zero);
      auto S(zero);
      for(int i = 0; i < h.size(); i ++)
        A += h[i];
      A /= num_t(int(h.size()));
      auto p(p0);
      for(int j = 0; j < 2; j ++) {
        for(int i = 0; i < h.size(); i ++)
          (void)p.next(S += h[(t + i) % h.size()] - A);
        S = zero;
      }
      for(int i = 0; i < h.size() - 1; i ++)
        (void)p.next(S += h[(t + i) % h.size()] - A);
      const auto S0(S);
      auto pn(p.next(S += d - A));
      if(pn != zero) M = pn - S0 + A;
    }
   next:
    std::cout << D << ", " << M << std::endl << std::flush;
  }
  return 0;
}

