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
// N.B. sectional measurement.
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
typedef sumChain<num_t, p0_7t>      p0_8t;
// N.B. we take probability on whole:
typedef shrinkMatrix<num_t, p0_8t> p0_t;
typedef sumChain<num_t, p0_t, true> p0_at;

// N.B. plain complex form.
typedef northPole<num_t, p0_1t> p0_s2t;
typedef northPole<num_t, p0_s2t> p0_s3t;
typedef sumChain<num_t, p0_s3t>      p0_s4t;
typedef shrinkMatrix<num_t, p0_s4t> p0_st;
typedef sumChain<num_t, p0_st, true> p0_ast;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int var(1);
  int look(1);
  if(argc <= 1) std::cerr << argv[0] << " <size>? <look>? : continue with ";
  if(1 < argc) var   = std::atoi(argv[1]);
  if(2 < argc) look  = std::atoi(argv[2]);
  std::cerr << argv[0] << " " << var << " " << look << std::endl;
  look += look + (look < 0 ? - 1 : 1);
  assert(look);
  // N.B. this is not optimal but we use this:
  const int  step(max(num_t(3), exp(log(num_t(abs(var * 2) + abs(look) - 1)) * log(num_t(abs(var * 2) + abs(look) - 1)))));
  // N.B. make periodical into complex plain arg periodical, not aligned.
  const auto recur(2);
  p0_t   p;
  p0_at  q;
  p0_st  pp;
  p0_ast qq;
  if(look < 0) {
    if(var < 0)
      qq = p0_ast(p0_st(p0_s4t(p0_s3t(p0_s2t(p0_1t(p0_0t(step, abs(var) + abs(look) - 1), abs(var)) )) ), abs(var)));
    else
      pp = p0_st(p0_s4t(p0_s3t(p0_s2t(p0_1t(p0_0t(step, abs(var) + abs(look) - 1), abs(var)) )) ), abs(var));
  } else {
    if(var < 0)
      q = p0_at(p0_t(p0_8t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, abs(var) + abs(look) - 1), abs(var)), step), step), step), step) )) ), abs(var)) );
    else
      p = p0_t(p0_8t(p0_7t(p0_6t(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(step, abs(var) + abs(look) - 1), abs(var)), step), step), step), step) )) ), abs(var));
  }
  num_t din(int(0));
  auto  Mx(din);
  std::vector<num_t> M;
  M.resize(abs(look), din);
  std::vector<num_t> d0;
  d0.resize(recur, din);
  auto  d1(d0);
  while(std::getline(std::cin, s, '\n')) {
    const auto bd(din);
    std::stringstream ins(s);
    ins >> din;
    num_t D(int(0));
    for(int i = 0; i < 2; i ++) {
      const auto  d(i ? din + din : din + bd);
            num_t dd(int(0));
      D += d * M[0];
      Mx = max(Mx, abs(d) * num_t(int(2)));
      if(recur) {
        d0[0] += d;
        for(int i = 1; i < d0.size(); i ++) {
          if(d1[i - 1] == num_t(int(0))) break;
          const auto work(d0[i - 1] / d1[i - 1] - num_t(int(1)));
          if(isfinite(work)) d0[i] += work;
          else break;
        }
        if(d1[d1.size() - 1] != num_t(int(0)) &&
           ! isfinite(dd = d0[d0.size() - 1] / d1[d1.size() - 1] - num_t(int(1))))
          dd = num_t(int(0));
      } else dd = d;
      for(int i = 0; i < M.size() - 1; i ++) M[i] = std::move(M[i + 1]);
      if(dd != num_t(int(0)))
        M[M.size() - 1] = (look < 0 ? (var < 0 ? qq.next(dd) : pp.next(dd)) : (var < 0 ? q.next(dd) : p.next(dd)));
      else M[M.size() - 1] = num_t(int(0));
      if(recur) for(int i = 0; i < d0.size(); i ++) M[M.size() - 1] *= d0[i];
      M[M.size() - 1] = max(- Mx, min(Mx, M[M.size() - 1]));
      d1 = d0;
    }
    std::cout << D << ", " << M[M.size() - 1] << ", " << M[M.size() - 2] << std::endl << std::flush;
  }
  return 0;
}

