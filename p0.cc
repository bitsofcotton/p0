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
// N.B. sectional measurement, also with probability.
typedef shrinkMatrix<num_t, p0_0t> p0_1t;
// N.B. make information-rich not to associative/commutative.
//      2 dimension semi-order causes (x, status) from input as sedenion.
typedef P0DFT<num_t, p0_1t, idFeeder<num_t> > p0_2t;
typedef P0DFT<num_t, p0_2t, idFeeder<num_t> > p0_3t;
typedef P0DFT<num_t, p0_3t, idFeeder<num_t> > p0_4t;
typedef P0DFT<num_t, p0_4t, idFeeder<num_t> > p0_5t;
// N.B. on any R to R into reasonable taylor.
typedef northPole<num_t, p0_5t> p0_6t;
typedef northPole<num_t, p0_6t> p0_7t;
// N.B. we make the prediction on (moved origin point) delta summation.
typedef sumChain<num_t, p0_7t>      p0_t;
typedef sumChain<num_t, p0_t, true> p0_at;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int var(- 2);
  if(argc <= 1) std::cerr << argv[0] << " <size> : continue with ";
  if(1 < argc) var = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << var << std::endl;
  // N.B. this is not optimal but we use this:
  const int step(max(num_t(3), exp(log(num_t(abs(var))) * log(num_t(abs(var))))));
  p0_t  p;
  p0_at q;
  if(var < 0)
    q = p0_at(p0_t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, abs(var)), abs(var)), step), step), step), step) ))) );
  else
    p = p0_t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, abs(var)), abs(var)), step), step), step), step) )));
  num_t d(int(0));
  auto  M(d);
  auto  Mx(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    Mx = max(Mx, abs(d) * num_t(int(abs(var) * 2)));
    std::cout << D << ", " << (M = max(- Mx, min(Mx, var < 0 ? q.next(d) : p.next(d) )) ) << std::endl << std::flush;
  }
  return 0;
}

