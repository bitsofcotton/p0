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
  int length(3);
  if(argc < 2) std::cerr << argv[0] << " <length>? : continue with ";
  if(1 < argc) length = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << length << std::endl;
  PBond<num_t, P0maxRank<num_t> > p(max(abs(length), 2));
  std::string s;
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    std::cout << d * M << ", ";
    std::cout << (length ? (0 < length ? M = p.next(d) : M = expscale<num_t>(p.next(logscale<num_t>(d)))) : M = d) << std::endl << std::flush;
  }
  return 0;
}

