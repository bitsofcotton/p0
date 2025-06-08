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
  idFeeder<num_t> p(max(0, - length));
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
    std::cout << (M = expscale<num_t>(length <= 0 ? p0maxNext<num_t>(p.next(logscale<num_t>(d))) : deep<num_t, p0maxNext<num_t> >(p.next(logscale<num_t>(d), length)) )) << std::endl << std::flush;
#else
    std::cout << (M = length <= 0 ? p0maxNext<num_t>(p.next(d)) : deep<num_t, p0maxNext<num_t> >(p.next(d), length) ) << std::endl << std::flush;
#endif
  }
  return 0;
}

