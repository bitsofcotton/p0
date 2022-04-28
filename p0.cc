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
typedef sumChain<num_t, p0_10t, true> p0_t;
// N.B. if original sample lebesgue integrate is not enough continuous,
//      imitate original function by some of sample points,
//      but move origin point to average one, so a little better
//      original function estimation.
// N.B. frequency space *= 2 causes nyquist frequency ok.
// N.B. but this is equivalent from jammer on PRNG, and probe on some
//      measurable phenomenon.
//typedef P0Expect<num_t, p0_11t> p0_t;

// N.B. plain complex form.
typedef northPole<num_t, p0_1t>  p0_s2t;
typedef northPole<num_t, p0_s2t> p0_s3t;
typedef sumChain<num_t, p0_s3t>  p0_s4t;
typedef logChain<num_t, p0_s4t>  p0_s5t;
typedef logChain<num_t, p0_s5t>  p0_s6t;
typedef sumChain<num_t, p0_s6t, true> p0_st;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int var(2);
  if(argc <= 1) std::cerr << argv[0] << " <npred>? : continue with ";
  if(1 < argc) var = std::atoi(argv[1]);
  assert(0 <= var);
  // N.B. this is not optimal but we use this:
  const int step(max(num_t(3), exp(log(num_t(abs(var))) * log(num_t(abs(var))))));
  p0_t  p0(p0_10t(p0_9t(p0_8t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, var), var), var), var), var), var) )) ) )) );
  auto  p1(p0);
  p0_st q0(p0_s6t(p0_s5t(p0_s4t(p0_s3t(p0_s2t(p0_1t(p0_0t(step, var), var) )) ) )) );
  auto  q1(q0);
  p0_0t r(3);
  num_t d(int(0));
  auto  dS(d);
  auto  M(d);
  auto  S(d);
  SimpleMatrix<num_t> Mstore(step, 4);
  Mstore.O();
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    if(var) {
      // N.B. compete with dimension the original function might have.
      //      (x, f(x), status, const.) is eliminated in this method,
      //      but not stable in practical ones.
      const auto bdS(dS);
      dS += d;
      for(int i = 0; i < Mstore.rows() - 1; i ++)
        Mstore.row(i) = Mstore.row(i + 1);
      Mstore(Mstore.rows() - 1, 0) = p0.next(d);
      Mstore(Mstore.rows() - 1, 1) = q0.next(d);
      if(bdS != num_t(int(0)) && dS != num_t(int(0))) {
        const auto ddS(num_t(int(1)) / dS - num_t(int(1)) / bdS);
        Mstore(Mstore.rows() - 1, 2) = num_t(int(1)) / (p1.next(ddS) + num_t(int(1)) / dS);
        Mstore(Mstore.rows() - 1, 3) = num_t(int(1)) / (q1.next(ddS) + num_t(int(1)) / dS);
      } else
        Mstore(Mstore.rows() - 1, 2) = Mstore(Mstore.rows() - 1, 3) = num_t(int(0));
      const auto MM(Mstore.QR().col(Mstore.rows() - 1));
      M = MM[0] + MM[1] + MM[2] + MM[3];
    } else
      M = r.next(d);
    std::cout << D << ", " << M << ", " << (S += D) << std::endl << std::flush;
  }
  return 0;
}

