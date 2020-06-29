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
  if(1 < argc)
    range = std::atoi(argv[1]);
  P0B<num_t> p(range);
  num_t d0(0);
  auto  MM(d0);
  auto  MM0(d0);
  auto  Md(d0);
  auto  d00(d0);
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
    if(d != bd) {
      d0 += (d - bd) * MM0;
      d1 += (d - bd)  * MM * num_t(bet1);
      d2 += (d - bd)  * MM * num_t(bet2);
      d3 += (d - bbd) * MM / num_t(2);
      if(d00 == num_t(0))
        d00 = d;
      MM = p.next(d - d00);
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

