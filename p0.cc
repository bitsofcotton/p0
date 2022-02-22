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
typedef shrinkMatrix<num_t, p0_4t, true> p0_5t;
typedef sumChain< num_t, p0_5t, true> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  p0_t  p(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(3)) ))) ));
  auto  q(p);
  int   t;
  num_t d(t ^= t);
  auto  M( d);
  auto  MM(d);
  auto  Mx(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D( d * M);
    const auto DD(D * MM);
    if(Mx < abs(d)) Mx = abs(d) * num_t(int(2));
    M = max(- Mx, min(Mx, p.next(d)));
    if(4 < t) {
      MM = q.next(D);
      if(Mx != num_t(int(0))) MM /= Mx;
      MM = max(- Mx, min(Mx, MM));
    }
    if(t ++ <= 8)
      std::cout << num_t(int(0)) << ", " << num_t(int(0)) << std::endl << std::flush;
    else {
      std::cout << DD << ", " << M * MM << std::endl << std::flush;
      -- t;
    }
  }
  return 0;
}

