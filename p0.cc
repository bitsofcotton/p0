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
typedef SimpleFloat<uint64_t, DUInt<uint64_t, 64>, 64, int32_t> num_t;
*/
#include "simplelin.hh"
#include "p0.hh"

typedef P0C<num_t, P0C<num_t, P0C<num_t, P0C<num_t, P0B<num_t> > > > > p_t;

int main(int argc, const char* argv[]) {
  std::string s;
  int range(30);
  if(1 < argc)
    range = std::atoi(argv[1]);
  p_t p(range);
  num_t d0(0);
  auto  bd(d0);
  auto  MM(d0);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    if(d != bd) {
      d0 += d - bd - MM;
      MM  = p.next(d) - d;
      bd  = d;
    }
    std::cout << d0 << ", " << MM << std::endl;
  }
  return 0;
}

