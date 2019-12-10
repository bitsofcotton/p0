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
using std::sqrt;
using std::atan2;
using std::log;
using std::exp;
using std::abs;
using std::sin;
using std::cos;
#include "p0.hh"
//typedef SimpleFloat<uint64_t, DUInt<uint64_t, 64>, 64, int16_t> num_t;
//typedef P0<num_t, Complex<num_t> > p0_t;
typedef double num_t;
typedef P0<num_t, std::complex<num_t> > p0_t;

int main(int argc, const char* argv[]) {
  std::string s;
  int range(20);
  if(1 < argc)
    range = std::atoi(argv[1]);
  assert(2 < range);
  int reset(80);
  if(2 < argc)
    reset = std::atoi(argv[2]);
  assert(2 < reset);
  p0_t  p(range, 2, 1);
  num_t M(0);
  num_t bd(0);
  num_t bbd(0);
  num_t d0(0);
  int   hit(0), nit(0);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    if(bd != d) {
      if((d - bbd) * M < num_t(0))
        hit ++;
      else if(num_t(0) < (d - bbd) * M)
        nit ++;
      d0 += M * (d - bbd);
      M   = p.next(d + bd) - d * num_t(2);
      bbd = bd;
      bd  = d;
    }
    std::cout << d0 << ", " << M << ", " << d << ", " << hit << ", " << nit << std::endl;
  }
  return 0;
}

