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
typedef shrinkMatrix<num_t, p0_0t> p0_1t;
typedef northPole<num_t, p0_1t> p0_2t;
typedef northPole<num_t, p0_2t> p0_3t;
// N.B. we apply them into probability.
//      if original function is lebesgue integrable and if the result is
//      continuous enough (without gulf), it's riemann integrable in probability.
//      on the other hand, for 0-markov's constant pred.
typedef shrinkMatrix<num_t, p0_3t> p0_t;
typedef avgOrigin<num_t, p0_t>    p0_avgo_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int   var(4);
  if(argc <= 1) std::cerr << argv[0] << " <size> : continue with ";
  if(1 < argc) var = abs(std::atoi(argv[1]));
  std::cerr << argv[0] << " " << var << std::endl;
  // N.B. this is not optimal but we use this:
  const auto step(int(exp(sqrt(log(num_t(var))))));
  // N.B. We need both p, q because of pnext step result s.t. pextend:
  p0_avgo_t  p(p0_t(p0_3t(p0_2t(p0_1t(p0_0t(var, step), step / 2))), (step / 2) + (step & 1)));
  p0_avgo_t  q(p0_t(p0_3t(p0_2t(p0_1t(p0_0t(var, (step - 1)), (step - 1) / 2))), ((step - 1) / 2) + ((step - 1) & 1)));
  num_t d(int(0));
  auto  M(d);
  auto  Mx(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    if(Mx < abs(d)) Mx = abs(d) * num_t(int(2));
          auto pn((p.next(d) + q.next(d)) / num_t(int(2)));
    if(isfinite(pn) && - Mx < pn && pn < Mx)
      M = max(- Mx, min(Mx, move(pn)));
    else M = num_t(int(0));
    std::cout << D << ", " << (M /= Mx != num_t(int(0)) ? Mx : num_t(int(1))) << std::endl << std::flush;
  }
  return 0;
}

