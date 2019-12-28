#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include "ifloat.hh"
#include "simplelin.hh"

#include <complex>
#include <cmath>
using namespace std;
// template <typename T> using complex = Complex<T>;
#include "p0.hh"
// typedef SimpleFloat<uint64_t, DUInt<uint64_t, 64>, 64, int16_t> num_t;
typedef double num_t;

int main(int argc, const char* argv[]) {
  std::string s;
  int range(20);
  if(1 < argc)
    range = std::atoi(argv[1]);
  P0<num_t>  p(range, 1);
  num_t M(0);
  num_t bd(0);
  num_t bbd(0);
  num_t d0(0);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    if(isfinite(M))
      d0 += (d - bd) * M;
    if(bd != d) {
      d0 += M * (d - bd);
      M   = p.next(d) - d;
    }
    bbd = bd;
    bd  = d;
    std::cout << M << ", " << d0 << ", " << std::endl << std::flush;
  }
  return 0;
}

