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
  int length(0);
  if(argc < 2) std::cerr << argv[0] << " <length>? : continue with ";
  if(1 < argc) length = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << length << std::endl;
  idFeeder<num_t> p(0);
  std::string s;
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    std::cout << d * M << ", " << d - M << ", ";
    std::cout << (M = p0p<num_t, p0maxNext<num_t> >(p.next(d)) ) << std::endl << std::flush;
  }
  return 0;
}

