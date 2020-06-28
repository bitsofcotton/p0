#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "assert.h"

#include <complex>
#include <cmath>
using namespace std;
typedef long double num_t;
/*
#include "ifloat.hh"
template <typename T> using complex = Complex<T>;
typedef SimpleFloat<DUInt<uint64_t, 64>, DUInt<DUInt<uint64_t, 64>, 128>, 128, int32_t> num_t;
*/

#include "simplelin.hh"
#include "p0.hh"

int main(int argc, const char* argv[]) {
  std::string s;
  int range(12);
  int loop(80);
  if(1 < argc)
    range = std::atoi(argv[1]);
  if(2 < argc)
    loop  = std::atoi(argv[2]);
  std::vector<P0B<num_t> > p;
  p.resize(loop, P0B<num_t>(range));
  std::vector<num_t> rr;
  rr.resize(p.size(), num_t(0));
  auto  dd(rr);
  num_t d0(0);
  auto  d00(d0);
  auto  bd(d0);
  auto  MM(d0);
  int   t(0);
  num_t d;
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    if(d00 == num_t(0))
      d00 = d;
    if(d != bd) {
      d0 += (d - bd) * MM;
      MM  = num_t(0);
      for(int i = 0; i < p.size(); i ++) {
        auto MM0i(p[i].next(dd[i] += (d - bd) * rr[i]));
        MM += MM0i *= (rr[i] = i ? num_t(arc4random() & 0x7fff) / num_t(0x10000) + num_t(1) / num_t(2) : num_t(1));
      }
      if(t ++ < range * 2)
        MM = num_t(0);
      MM /= num_t(int(p.size()));
    }
    std::cout << d0 << ", " << MM << std::endl;
    bd = d;
  }
  return 0;
}

