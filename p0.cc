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
typedef northPole<num_t, p0_3t> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  p0_t   p(p0_3t(p0_2t(p0_1t(p0_0t(3)) )));
  myuint t;
  num_t  d(t ^= t);
  auto   M(d);
  auto   Mx(d);
  auto   S(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    if(Mx < abs(d)) Mx = abs(d) * num_t(int(2));
    const auto A((S += d) / num_t(++ t));
    // N.B. predict any R to R on 1~3 markov.
          auto pn(p.next(d - A) + A);
    if(isfinite(pn) &&
       isfinite(pn) && - Mx / sqrt(sqrt(SimpleMatrix<num_t>().epsilon)) &&
                    pn < Mx / sqrt(sqrt(SimpleMatrix<num_t>().epsilon)) )
      M = max(- Mx, min(Mx, move(pn)));
    else M = num_t(int(0));
    std::cout << D << ", " << M << std::endl << std::flush;
  }
  return 0;
}

