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
//      if we make hypothesis markov, >1-markov is needed because
//      if the status is null, only constant is accepted.
// N.B. if the sampling frequency is not enough, middle range of the original
//      function frequency (enough large bands) will effect prediction fail.
//      this is because we only observes highest and lowest frequency on
//      sampling points, so omitted part exists.
//      even if the parameter on P0 is large, situation unchange.
//      so we should use shrinkMatrix for them.
typedef P0<num_t, idFeeder<num_t> > p0_0t;
// N.B. on any of R to R with sectional measurement.
//      on the other hand, hypothesis 1~3-markov predict with 3-markov.
typedef sumChain< num_t, p0_0t> p0_1t;
typedef northPole<num_t, p0_1t> p0_2t;
typedef sumChain< num_t, p0_2t, true> p0_3t;
typedef sumChain< num_t, p0_3t> p0_4t;
typedef northPole<num_t, p0_4t> p0_5t;
typedef sumChain< num_t, p0_5t, true> p0_t;
// N.B. we apply them into probability.
//      if original function is lebesgue integrable and if the result is
//      continuous enough (without gulf), it's riemann integrable in probability.
//      on the other hand, for 0-markov's constant pred.
typedef shrinkMatrix<num_t, p0_t> p0_jt;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int   middle(3);
  int   step(1);
  if(argc <= 1) std::cerr << argv[0] << " <size> <step> : continue with ";
  if(1 < argc) middle = abs(std::atoi(argv[1]));
  if(2 < argc) step   = abs(std::atoi(argv[2]));
  std::cerr << argv[0] << " " << middle << " " << step << std::endl;
  p0_jt p(p0_t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(middle, step))))))), step);
  auto  q(p);
  num_t d(int(0));
  auto  M(d);
  auto  Mx(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    if(Mx < abs(d)) Mx = abs(d) * num_t(int(2));
          auto pn(p.next(d));
    if(isfinite(pn) && - Mx < pn && pn < Mx)
      M = max(- Mx, min(Mx, move(pn)));
    else M = num_t(int(0));
    std::cout << D << ", " << (M /= Mx != num_t(int(0)) ? Mx : num_t(int(1))) << std::endl << std::flush;
  }
  return 0;
}

