#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include "ifloat.hh"
#include "simplelin.hh"
#include "p0.hh"

const int range(20);
typedef SimpleFloat<uint64_t, DUInt<uint64_t, 64>, 64, int16_t> num_t;
typedef P0<num_t, Complex<num_t> > p0_t;

int main() {
  std::string s;
  p0_t  p0(range);
  num_t M(0);
  num_t bd(0);
  num_t d0(0);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    if(bd != d) {
      d0 += (d - bd) * M;
      try {
        M = p0.next(bd = d) - d;
      } catch(const char* e) {
        std::cerr << e << std::endl;
      }
    }
    std::cout << M << ", " << d0 << std::endl;
  }
  return 0;
}

