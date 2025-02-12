#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include <stdint.h>
#include <sys/resource.h>

#include "lieonn.hh"
typedef myfloat num_t;

template <typename T> static inline T expscale(const T& x) {
  return sgn<T>(x) * (exp(abs(x)) - T(int(1)));
}

template <typename T> static inline T logscale(const T& x) {
  return sgn<T>(x) * log(abs(x) + T(int(1)));
}

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  if(1 < argc && argv[1][0] == 'r') {
    assert(3 < argc);
    const auto  len(std::atoi(argv[2]));
    const auto  step(std::atoi(argv[3]));
    const auto  d(diff<num_t>(len).row(len - 1));
    const auto& n(pnextcacher<num_t>(len, step));
    num_t sumn(int(0));
    auto  sumd(sumn);
    auto  n2n(sumn);
    auto  n2d(sumn);
    for(int i = 0; i < len; i ++) {
      std::cout << n[i] << ", " << d[i] << std::endl;
      sumn += n[i];
      sumd += d[i];
      n2n  += n[i] * n[i];
      n2d  += d[i] * d[i];
    }
    std::cout << sumn << ", " << sumd << ", " << sqrt(n2n) << ", " << sqrt(n2d) << std::endl;
    return 0;
  }
  int step(0);
  int length(0);
  if(argc < 2) std::cerr << argv[0] << " <step>? <length>? : continue with ";
  if(1 < argc) step   = std::atoi(argv[1]);
  if(2 < argc) length = std::atoi(argv[2]);
  assert(0 <= step);
  std::cerr << argv[0] << " " << step << " " << length << std::endl;
  PBond<num_t, P0maxRank<num_t> > p(max(abs(length), 2), P0maxRank<num_t>(step));
  SimpleVector<num_t> b;
  idFeeder<num_t> f(step ? step : 1);
  vector<num_t> heavy;
  std::string s;
  num_t zero(int(0));
  auto  d(zero);
  while(std::getline(std::cin, s, '\n')) {
    const auto& M(heavy.size() ? heavy[0] : (f.full ? f.res[0] : zero));
    std::stringstream ins(s);
    ins >> d;
#if defined(_CHAIN_)
    // std::cout << pseudoerfscale<num_t>((d = pseudoierfscale<num_t>(d)) - M) << ", ";
    std::cout << d - M << ", ";
#else
    std::cout << d * M << ", ";
#endif
    if(! length) {
      b.entity.emplace_back(d);
      if(step)
        std::cout << f.next(2 < b.size() ? P0maxRank<num_t>(step).next(b) : num_t(int(0)))[f.res.size() - 1] << std::endl << std::flush;
      else {
        if(heavy.size()) {
          for(int i = 1; i < heavy.size(); i ++)
            heavy[i - 1] = move(heavy[i]);
          heavy[heavy.size() - 1] = num_t(int(0));
        }
        if(1 < b.entity.size()) {
          heavy.emplace_back(num_t(int(0)));
          for(int i = 1; i <= heavy.size(); i ++)
            heavy[i - 1] += P0maxRank<num_t>(i).next(b);
        }
        std::cout << (heavy.size() ? heavy[0] /= num_t(int(heavy.size())) : num_t(int(0))) << std::endl << std::flush;
      }
    } else if(step)
      std::cout << f.next(0 < length ? p.next(d) : expscale<num_t>(p.next(logscale<num_t>(d))) )[f.res.size() - 1] << std::endl << std::flush;
    else assert(step);
  }
  return 0;
}

