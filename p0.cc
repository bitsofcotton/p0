#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "assert.h"

#include <complex>
#include <cmath>
using namespace std;
typedef long double num_t;
/*
#include "ifloat.hh"
template <typename T> using complex = Complex<T>;
typedef SimpleFloat<uint64_t, DUInt<uint64_t, 64>, 64, int32_t> num_t;
*/
#include "simplelin.hh"
#include "p0.hh"

typedef P0C<num_t, P0B<num_t> > p_t;
p_t p;

int main(int argc, const char* argv[]) {
  std::string s;
  int range(20);
  int loop(14);
  if(1 < argc)
    range = std::atoi(argv[1]);
  if(2 < argc)
    loop  = std::atoi(argv[2]);
  p = p_t(range, loop);
  num_t d0(0);
  auto  MM(d0);
  auto  Md(d0);
  auto  bd(d0);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    if(bd != 0)
      Md = max(Md, abs(d - bd) * num_t(2));
    if(d != bd && Md != 0) {
      d0 += (d - bd) * MM;
      MM  = p.next(d / Md) * Md - d;
      if(! isfinite(MM) || isnan(MM))
        MM = num_t(0);
    }
    std::cout << d0 << ", " << MM << ", " << Md << std::endl;
    bd = d;
  }
  return 0;
}

