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
//      if the sampling frequency is not enough, middle range of the original
//      function frequency (enough large bands) will effect prediction fail.
//      this is because we only observes highest and lowest frequency on
//      sampling points, so omitted part exists.
//      even if the parameter on P0 is large, situation unchange.
//      so we should use shrinkMatrix for them.
typedef P0<num_t, idFeeder<num_t> > p0_0t;
// N.B. make information-rich not to associative/commutative.
typedef P0DFT<num_t, p0_0t, idFeeder<num_t> > p0_1t;
typedef P0DFT<num_t, p0_1t, idFeeder<num_t> > p0_2t;
// N.B. lg(7) < 3, zero divisor with var == 2 sectional.
typedef P0DFT<num_t, p0_2t, idFeeder<num_t> > p0_3t;
typedef P0DFT<num_t, p0_3t, idFeeder<num_t> > p0_4t;
// N.B. on any of R to R with sectional measurement.
typedef shrinkMatrix<num_t, p0_4t> p0_5t;
typedef northPole<num_t, p0_5t> p0_6t;
typedef northPole<num_t, p0_6t> p0_7t;
// N.B. we apply them into probability.
//    if original function is lebesgue integrable and if the result is
//    continuous enough (without gulf), it's riemann integrable in probability.
typedef shrinkMatrix<num_t, p0_7t>  p0_8t;
typedef sumChain<num_t, p0_8t>      p0_t;
typedef sumChain<num_t, p0_t, true> p0_at;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  const auto var(2);
  // N.B. up to 7-markov.
  const auto step(7);
  p0_t  p(p0_8t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, abs(var) * 2), abs(step)), abs(step)), abs(step)), abs(step)), abs(var)))), abs(var)) );
  p0_at q(p0_t(p0_8t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, abs(var) * 2), abs(step)), abs(step)), abs(step)), abs(step)), abs(var)))), abs(var)) ));
  num_t d(int(0));
  auto  M(d);
  auto  Mx(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    Mx = max(Mx, abs(d) * num_t(int(abs(var) * 2)));
    std::cout << D << ", " << (M = max(- Mx, min(Mx, var < 0 ? q.next(d) : p.next(d) ) ) ) << std::endl << std::flush;
  }
  return 0;
}

