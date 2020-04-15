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
// typedef SimpleFloat<uint64_t, DUInt<uint64_t, 64>, 64, int32_t> num_t;
typedef SimpleFloat<DUInt<uint64_t, 64>, DUInt<DUInt<uint64_t, 64>, 128>, 128, int16_t> num_t;
*/
#include "simplelin.hh"
#include "p0.hh"

template <typename T> const T& sgn(const T& x) {
  const static T zero(0);
  const static T one(1);
  const static T mone(- 1);
  if(zero < x)
    return one;
  if(x < zero)
    return mone;
  return zero;
}

int main(int argc, const char* argv[]) {
  std::string s;
  int range(8);
  int look(1);
  if(1 < argc)
    range = std::atoi(argv[1]);
  if(2 < argc)
    look  = std::atoi(argv[2]);
  P0<num_t> p(range, look);
  SimpleVector<num_t> buf(range);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = num_t(0);
  auto  ibuf(buf);
  num_t d0(0);
  auto  MM(d0);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    const auto& bd(buf[buf.size() - 1]);
    if(d != bd) {
      d0 += (d - bd) * MM;
      for(int i = 0; i < buf.size() - 1; i ++) {
        buf[ i] = buf[ i + 1];
        ibuf[i] = ibuf[i + 1];
      }
      buf[ buf.size()  - 1] = d;
      ibuf[ibuf.size() - 1] = num_t(1) / d;
      MM  = (p.next(buf) + num_t(1) / p.next(ibuf)) / num_t(2) - p.lpf(buf.size()).row(buf.size() - 1).dot(buf);
    }
    std::cout << d0 << "," << MM << std::endl;
  }
  return 0;
}

