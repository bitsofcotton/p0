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

const auto range(20);
int main(int argc, const char* argv[]) {
  std::string s;
  P0B<num_t, 6> p(range);
  num_t d0(0);
  auto  bd(d0);
  auto  MM(d0);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      d0 += (d - bd) * MM;
      MM  = p.next(d) - d;
    }
    std::cout << d0 << "," << MM << std::endl;
    bd = d;
  }
  return 0;
}

