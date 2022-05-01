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
  if(argc <= 1) std::cerr << argv[0] << " <var>? : continue with ";
  if(1 < argc) var = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << var << std::endl;
  assert(0 <= var);
  num_t d(int(0));
  auto  M(d);
  auto  S(d);
  if(! var) {
    p0_0t p(3);
    while(std::getline(std::cin, s, '\n')) {
      std::stringstream ins(s);
      ins >> d;
      const auto D(d * M);
      std::cout << D << ", " << (M = p.next(d)) << ", " << (S += D) << std::endl << std::flush;
    }
    return 0;
  }
  // N.B. this is not optimal but we use this:
  const int step(max(num_t(3), exp(log(num_t(var)) * log(num_t(var)))));
  p0_t  p, pp;
  p0_st q, qq;
  auto  dS(d);
  bool  need_init(true);
  SimpleMatrix<num_t> Mstore(8, max(8, step));
  Mstore.O();
  while(std::getline(std::cin, s, '\n')) {
    if(need_init) {
      p  = p0_t(p0_10t(p0_9t(p0_8t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, var), var), var), var), var), var) )) ) )) );
      q  = p0_st(p0_s6t(p0_s5t(p0_s4t(p0_s3t(p0_s2t(p0_1t(p0_0t(step, var), var) )) ) )) );
      pp = p;
      qq = q;
      need_init = false;
    }
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    const auto bdS(dS);
    dS += d;
    if(bdS == num_t(int(0)) || dS == num_t(int(0))) {
      std::cout << D << ", " << M << ", " << (S += D) << ", " << 0 << std::endl << std::flush;
      continue;
    }
    const auto idS(num_t(int(1)) / dS);
    const auto dd(idS - num_t(int(1)) / bdS);
    if(! isfinite(dd) || isnan(dd)) {
      std::cout << D << ", " << M << ", " << (S += D) << ", " << 0 << std::endl << std::flush;
      continue;
    }
    // N.B. compete with dimension the original function might have.
    //      (x, f(x), status, const.) is eliminated twice in this method
    //      for anti-symmetric ones.
    const auto pd(pp.next(d));
    const auto qd(qq.next(d));
    const auto pdd(p.next(dd));
    const auto qdd(q.next(dd));
    const auto pp2(pdd + idS);
    const auto pq2(qdd + idS);
    const auto pp3(pd  +  dS);
    const auto pq3(qd  +  dS);
    for(int i = 0; i < Mstore.cols() - 1; i ++)
      Mstore.setCol(i, Mstore.col(i + 1));
    Mstore(0, Mstore.cols() - 1) =   pd;
    Mstore(1, Mstore.cols() - 1) =   qd;
    Mstore(2, Mstore.cols() - 1) = - pdd;
    Mstore(3, Mstore.cols() - 1) = - qdd;
    Mstore(4, Mstore.cols() - 1) = pp2 == num_t(int(0)) ? pp2 :    num_t(int(1)) / std::move(pp2) - dS;
    Mstore(5, Mstore.cols() - 1) = pq2 == num_t(int(0)) ? pq2 :    num_t(int(1)) / std::move(pq2) - dS;
    Mstore(6, Mstore.cols() - 1) = pp3 == num_t(int(0)) ? pp3 : - (num_t(int(1)) / std::move(pp3) - idS);
    Mstore(7, Mstore.cols() - 1) = pq3 == num_t(int(0)) ? pq3 : - (num_t(int(1)) / std::move(pq3) - idS);
    const auto lsvd(Mstore.SVD());
    const auto svd(lsvd * Mstore);
    std::vector<num_t> stat;
    stat.reserve(svd.rows());
    int rank;
    M = num_t(rank ^= rank);
    for(int i = 0; i < svd.rows(); i ++)
      stat.emplace_back(sqrt(svd.row(i).dot(svd.row(i))));
    auto sstat(stat);
    std::sort(sstat.begin(), sstat.end());
    for(int i = 0; i < svd.rows(); i ++)
      if(sstat[sstat.size() - 1] * sqrt(sqrt(SimpleMatrix<num_t>().epsilon)) <
         stat[i]) {
        auto sum(lsvd(i, 0));
        for(int j = 1; j < lsvd.cols(); j ++)
          sum += lsvd(i, j);
        M += svd(i, svd.cols() - 1) / stat[i] * sgn<num_t>(sum);
        rank ++;
      }
    if(! isfinite(M)) M = num_t(int(0));
    std::cout << D << ", " << (M /= num_t(rank ? rank : 1)) << ", " << (S += D) << ", " << rank << std::endl << std::flush;
  }
  return 0;
}

