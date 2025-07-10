#if !defined(_ONEBINARY_)
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

#define _COMPILE_PRED_
#include "lieonn.hh"
typedef myfloat num_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  int length(0);
  if(argc < 2) std::cerr << argv[0] << " <length>? : continue with ";
  if(1 < argc) length = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << length << std::endl;
  std::string s;
# if defined(_CHAIN_)
  const bool chain(true);
# else
  const bool chain(false);
# endif
#endif
  idFeeder<num_t> p(max(int(0), length));
  num_t d(int(0));
  num_t M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    std::cout << (chain ? d - M : d * M) << ", ";
#if defined(_NONLIN_)
    std::cout << (M = expscale<num_t>(p0maxNext<num_t>(p.next(logscale<num_t>(d))) )) << std::endl << std::flush;
#else
    std::cout << (M = p0maxNext<num_t>(p.next(d)) ) << std::endl << std::flush;
#endif
  }
#if !defined(_ONEBINARY_)
  return 0;
}
#endif

