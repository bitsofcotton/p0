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
#include <random>
#include <assert.h>
#include <stdint.h>
#include <sys/resource.h>

#include "lieonn.hh"
typedef myfloat num_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  int length(0);
  if(argc < 2) std::cerr << argv[0] << " <length>? : continue with ";
  if(1 < argc) length = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << length << std::endl;
  if(length < - 1) {
    const auto  len(abs(length));
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
  idFeeder<num_t> p(0);
  idFeeder<num_t> p3(length == 1 ? 3 : (length == 2 ? 4 : 0));
  std::string s;
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
#if defined(_CHAIN_)
    std::cout << d - M << ", ";
#else
    std::cout << d * M << ", ";
#endif
#if defined(_NONLIN_)
    std::cout << (M = expscale<num_t>(length ? (length < 3 ? (length == - 1 ? logscale<num_t>(d) : (length == 1 ? p0maxNext<num_t>(p3.next(logscale<num_t>(d))) : p3.next(logscale<num_t>(d))[0]) ) : deep<num_t, p0maxNext<num_t> >(p.next(logscale<num_t>(d)), length) ) ) : p0maxNext<num_t>(p.next(logscale<num_t>(d)) )) ) ) << std::endl << std::flush;
#else
    std::cout << (M = length ? (length < 3 ? (length == - 1 ? d : (length == 1 ? p0maxNext<num_t>(p3.next(d)) : p3.next(d)[0]) ) : deep<num_t, p0maxNext<num_t> >(p.next(d), length) ) : p0maxNext<num_t>(p.next(d)) ) << std::endl << std::flush;
#endif
  }
  return 0;
}

