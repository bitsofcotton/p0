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
//      so we should use sectional measurement for them.
typedef P0<num_t, idFeeder<num_t> > p0_0t;
// N.B. sectional measurement, also expected value.
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
// N.B. we make the prediction on (delta) summation.
typedef sumChain<num_t, p0_7t>  p0_8t;
// N.B. we treat periodical part as non aligned complex arg part.
typedef logChain<num_t, p0_8t>  p0_9t;
typedef logChain<num_t, p0_9t>  p0_10t;
// N.B. we take average as origin of input.
typedef sumChain<num_t, p0_10t, true> p0_11t;
// N.B. if original sample lebesgue integrate is not enough continuous,
//      imitate original function by some of sample points,
//      but move origin point to average one, so a little better
//      original function estimation.
// N.B. frequency space *= 2 causes nyquist frequency ok.
typedef P0Expect<num_t, p0_11t> p0_t;

// N.B. plain complex form.
typedef northPole<num_t, p0_1t>  p0_s2t;
typedef northPole<num_t, p0_s2t> p0_s3t;
typedef sumChain<num_t, p0_s3t>  p0_s4t;
typedef logChain<num_t, p0_s4t>  p0_s5t;
typedef logChain<num_t, p0_s5t>  p0_s6t;
typedef sumChain<num_t, p0_s6t, true> p0_s7t;
typedef P0Expect<num_t, p0_s7t>  p0_st;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int var(2);
  if(argc <= 1) std::cerr << argv[0] << " <markov>? : continue with ";
  if(1 < argc) var = std::atoi(argv[1]);
  // N.B. this is not optimal but we use this:
  const int step(max(num_t(3), exp(log(num_t(abs(var))) * log(num_t(abs(var))))));
  p0_t  p;
  p0_st q;
  p0_0t r;
  num_t d(int(0));
  auto  M(d);
  auto  S(d);
  bool  need_init(true);
  while(std::getline(std::cin, s, '\n')) {
    if(need_init) {
      std::cerr << argv[0] << " " << var << std::endl;
      if(var < 0)
        q = p0_st(p0_s7t(p0_s6t(p0_s5t(p0_s4t(p0_s3t(p0_s2t(p0_1t(p0_0t(step, - var), - var) )) ) )) ), step);
      else if(0 < var)
        p = p0_t(p0_11t(p0_10t(p0_9t(p0_8t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, var), var), step), step), step), step) )) ) )) ), step);
      else
        r = p0_0t(3);
      need_init = false;
    }
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    std::cout << D << ", " << (M = num_t(var ? (var < 0 ? q.next(d) : p.next(d)) : r.next(d)) ) << ", " << (S += D) << std::endl << std::flush;;
  }
  return 0;
}

