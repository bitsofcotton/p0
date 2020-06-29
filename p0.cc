#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
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

template <typename T> const T& sgn(const T& x) {
  const static T zero(0);
  const static T one(1);
  const static T mone(- 1);
  return x == zero ? zero : (x < zero ? mone : one);
}

int main(int argc, const char* argv[]) {
  std::string s;
  int range(12);
  if(1 < argc)
    range = std::atoi(argv[1]);
  P0B<num_t> p(range);
  num_t d(0);
  auto  d0(d);
  auto  d1(d);
  auto  bd(d);
  auto  MM(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      d0 += d - MM;
      d1 += (d - bd) * (MM - bd);
      MM  = p.next(d);
    }
    std::cout << d0 << ", " << d1 << ", " << MM << std::endl;
    bd = d;
  }
  return 0;
}

