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

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  if(1 < argc && argv[1][0] == 'r') {
    assert(2 < argc);
    const auto  len(std::atoi(argv[2]));
    const auto  d(diff<num_t>(len).row(len - 1));
    const auto& n(pnextcacher<num_t>(len, 1));
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
  int step(1);
  int length(0);
  if(argc < 2) std::cerr << argv[0] << " <step>? <length>? : continue with ";
  if(1 < argc) step   = std::atoi(argv[1]);
  if(2 < argc) length = std::atoi(argv[2]);
  std::cerr << argv[0] << " " << step << " " << length << std::endl;
  assert(0 < step);
  PBond<num_t, P0maxRank<num_t> > p(max(abs(length), 2), P0maxRank<num_t>(step));
  SimpleVector<num_t> b;
  idFeeder<num_t> f(step);
  std::string s;
  num_t zero(int(0));
  auto  d(zero);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    const auto& M(f.res[0]);
    ins >> d;
#if defined(_CHAIN_)
    // std::cout << pseudoerfscale<num_t>(((d = pseudoierfscale<num_t>(d / num_t(int(2)))) - M) * num_t(int(2))) << ", ";
    std::cout << d - M << ", ";
#else
    std::cout << d * M << ", ";
#endif
    if(! length) {
      b.entity.emplace_back(d);
      std::cout << f.next(2 < b.size() ? P0maxRank<num_t>().next(b) : num_t(int(0)))[step - 1] << std::endl << std::flush;
    } else
      std::cout << f.next(0 < length ? p.next(d) : expscale<num_t>(p.next(logscale<num_t>(d))) )[step - 1] << std::endl << std::flush;
  }
  return 0;
}

