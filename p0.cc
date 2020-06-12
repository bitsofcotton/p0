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
typedef SimpleFloat<uint64_t, DUInt<uint64_t, 64>, 64, Signed<DUInt<DUInt<DUInt<DUInt<uint64_t, 64>, 128>, 256>, 512>, 1024> > num_t;
*/
#include "simplelin.hh"
#include "p0.hh"

int main(int argc, const char* argv[]) {
  std::string s;
  int range(4);
  if(1 < argc)
    range = std::atoi(argv[1]);
  P0B<num_t> p(range * 3);
  P0B<num_t> q(range * 4);
  P0B<num_t> r(range * 5);
  num_t d0(0);
  auto  bd(d0);
  auto  MM(d0);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      d0 += (d - bd) * MM;
      MM  = (p.next(d) + q.next(d) + r.next(d)) / num_t(3) - d;
    }
    std::cout << d0 << "," << MM << std::endl;
    bd = d;
  }
  return 0;
}

