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
typedef shrinkMatrix<num_t, P0<num_t, idFeeder<num_t> > > p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  auto var(1);
  if(argc < 2)
    std::cerr << argv[0] << " <len>?" << std::endl;
  if(1 < argc) var  = std::atoi(argv[1]);
  std::cerr << "continue with " << argv[0] << " " << var << std::endl;
  northPole<num_t, northPole<num_t, p0_t> > p(northPole<num_t, p0_t>(p0_t(P0<num_t, idFeeder<num_t> >(var < 0 ? - var : 3), var < 0 ? 1 : var)));
  int   t;
  num_t d(t ^= t);
  auto  M(d);
  auto  S(d);
  auto  Snote(d);
  vector<num_t> h;
  h.resize(abs(var) * 2 + 3, num_t(int(0)));
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    const auto S0(S);
    Snote += d;
    if(d != num_t(int(0))) {
      h[(t ++) % h.size()] = d;
      if(! (t % h.size())) {
        S = num_t(int(0));
        for(int i = 0; i < h.size() - 1; i ++)
          (void)p.next(S += h[i] - Snote / num_t(int(t)));
      }
      auto pn(p.next(S += d - Snote / num_t(int(t))));
      std::cout << D << ", " << (pn == num_t(int(0)) ? M : M = pn - S0 + Snote / num_t(int(t))) << std::endl << std::flush;
    } else
      std::cout << D << ", " << M << std::endl << std::flush;
  }
  return 0;
}

