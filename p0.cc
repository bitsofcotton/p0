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
// N.B. on existing taylor series.
//      on the other hand, hypothesis 1~3-markov to be 3-markov prediction.
//         if we make hypothesis markov, >1-markov is needed because
//         if the status is null, only constant is accepted.
typedef P0<num_t, idFeeder<num_t> > p0_0t;
// N.B. on any R to R with sectional measurement.
//      on the other hand, hypothesis 1~3-markov to be 3-markov.
typedef shrinkMatrix<num_t, p0_0t, true> p0_1t;
typedef northPole<num_t, p0_1t> p0_2t;
typedef sumChain< num_t, p0_2t> p0_3t;
typedef northPole<num_t, p0_3t> p0_4t;
typedef sumChain< num_t, p0_4t, true> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  p0_t   p(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(3)) ))) );
  num_t  d(int(0));
  auto   M(d);
  auto   Mx(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    if(Mx < abs(d)) Mx = abs(d) * num_t(int(2));
          auto pn(p.next(d));
    if(isfinite(pn) && - Mx / sqrt(sqrt(SimpleMatrix<num_t>().epsilon)) < pn
                 && pn < Mx / sqrt(sqrt(SimpleMatrix<num_t>().epsilon)) )
      M = max(- Mx, min(Mx, move(pn)));
    else M = num_t(int(0));
    std::cout << D << ", " << (M /= Mx != num_t(int(0)) ? Mx : num_t(int(1))) << std::endl << std::flush;
  }
  return 0;
}

