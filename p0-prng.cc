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

int main(int argc, const char* argv[]) {
  std::string s;
  int range(40);
  int loop(8);
  if(1 < argc)
    range = std::atoi(argv[1]);
  if(2 < argc)
    loop  = std::atoi(argv[2]);
  P0C<num_t, P0B<num_t> > p(range, loop);
  num_t d0(0);
  auto  MM(d0);
  auto  MM0(d0);
  auto  Md(d0);
  auto  d1(d0);
  auto  d2(d0);
  auto  d3(d0);
  auto  bd(d0);
  auto  bbd(d0);
  auto  bd10(d0);
  auto  bd11(d0);
  auto  bd20(d0);
  auto  bd21(d0);
  int   t(0);
  auto  bet1(t);
  auto  bet2(t);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    if(bd != 0)
      Md = max(Md, abs(d - bd) * num_t(4));
    if(d != bd && Md != 0) {
      if(! isnan(MM) && MM != 0) {
        d0 += (d - bd) * MM * num_t(bet2 - bet1);
        d1 += (d - bd) * MM * num_t(bet1);
        d2 += (d - bd) * MM * num_t(bet2);
        d3 += (d - bbd) * MM;
      }
      MM  = p.next(d / Md) * Md - d;
      if(! isfinite(MM) || isnan(MM) || t ++ < range * 2)
        MM = num_t(0) / num_t(0);
      if((d3 - bd11) < (d1 - bd10)) {
        bet1 = 0;
        bd10 = d1;
        bd11 = d3;
      }
      if((d2 - bd20) < (d3 - bd21)) {
        bet2 = 0;
        bd20 = d2;
        bd21 = d3;
      }
      MM0 = MM * num_t((++ bet2) - (++ bet1));
      bbd = bd;
    }
    std::cout << d0 << ", " << MM0 << std::endl;
    bd = d;
  }
  return 0;
}

