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
  northPole<num_t, shrinkMatrix<num_t, northPole<num_t, P0<num_t, idFeeder<num_t> > > > > p(shrinkMatrix<num_t, northPole<num_t, P0<num_t, idFeeder<num_t> > > >(northPole<num_t, P0<num_t, idFeeder<num_t> > >(P0<num_t, idFeeder<num_t> >(var < 0 ? - var : 3)), var < 0 ? 1 : var));
  auto  p0(p);
  int   t;
  num_t d(t ^= t);
  auto  M(d);
  auto  S(d);
  auto  A(d);
  vector<num_t> h;
  h.resize(abs(var) * 2 + 3, num_t(int(0)));
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    h[(t ++) % h.size()] = d;
    if(h.size() <= t) {
      A = S = num_t(int(0));
      p = p0;
      for(int i = 0; i < h.size(); i ++)
        A += h[i];
      A /= num_t(int(h.size()));
      for(int j = 0; j < 2; j ++)
        for(int i = 0; i < h.size(); i ++)
          (void)p.next(S += h[i] - A);
      for(int i = 0; i < h.size() - 1; i ++)
        (void)p.next(S += h[i] - A);
      t ^= t;
    }
    const auto S0(S);
    auto pn(p.next(S += d - A));
    std::cout << D << ", " << (pn == num_t(int(0)) ? M : M = pn - S0 + A) << std::endl << std::flush;
  }
  return 0;
}

