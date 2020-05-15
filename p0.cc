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
#include <stdlib.h>

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

template <typename T> T expscale(const T& x) {
  return sgn(x) * (exp(abs(x)) - T(1));
}

template <typename T> T logscale(const T& x) {
  return sgn(x) * log(abs(x) + T(1));
}

int main(int argc, const char* argv[]) {
  std::string s;
  int range(8);
  int recur(200);
  if(1 < argc)
    range = std::atoi(argv[1]);
  if(2 < argc)
    recur = std::atoi(argv[2]);
  P0<num_t> p(range);
  std::vector<SimpleVector<num_t> > buf;
  SimpleVector<num_t> wbuf(range);
  for(int i = 0; i < wbuf.size(); i ++)
    wbuf[i] = num_t(0);
  buf.resize(recur, wbuf);
  auto  bufe(buf);
  auto  bufl(buf);
  num_t d00(0);
  auto  MM0(d00);
  std::vector<num_t> d0;
  d0.resize(recur, num_t(0));
  auto  d1(d0);
  auto  d2(d0);
  auto  d3(d0);
  auto  bd10(d0);
  auto  bd11(d0);
  auto  bd20(d0);
  auto  bd21(d0);
  auto  MM(d0);
  auto  MMs(d0);
  auto  rr(d0);
  std::vector<int> bet1;
  bet1.resize(d0.size(), 0);
  auto  bet2(bet1);
  while(std::getline(std::cin, s, '\n')) {
    num_t d;
    std::stringstream ins(s);
    ins >> d;
    d = num_t(1) / d;
    if(! isfinite(d) || isnan(d)) {
      std::cout << d00 << ", " << MM0 << std::endl;
      continue;
    }
    const auto& bd(buf[0][buf[0].size() - 1]);
    const auto& bbd(buf[0][buf[0].size() - 2]);
    if(d != bd) {
      d00 += (d - bd) * MM0;
      MM0  = num_t(0);
      for(int j = 0; j < buf.size(); j ++) {
        d0[j] += (d - bd) * MM[j] * num_t(bet2[j] - bet1[j]);
        d1[j] += (d - bd) * MM[j] * num_t(bet1[j]);
        d2[j] += (d - bd) * MM[j] * num_t(bet2[j]);
        d3[j] += (d - bbd) * MM[j];
        for(int i = 0; i < buf[j].size() - 1; i ++) {
          buf[j][ i] = buf[j][ i + 1];
          bufe[j][i] = bufe[j][i + 1];
          bufl[j][i] = bufl[j][i + 1];
        }
        buf[j][ buf[j].size()  - 1] += (d - bd) * rr[j];
        bufe[j][bufe[j].size() - 1]  = expscale(buf[j][buf[j].size() - 1]);
        bufl[j][bufl[j].size() - 1]  = logscale(buf[j][buf[j].size() - 1]);
        rr[j] = j ? (num_t(arc4random() & 0xfffff) + num_t(0xfffff)) / num_t(0xfffff) : num_t(1);
        MM[j]  = ((p.next(buf[j]) + logscale(p.next(bufe[j])) +
                                    expscale(p.next(bufl[j]))) / num_t(3) -
                  buf[j][buf[j].size() - 1]) * rr[j];
        if(!isfinite(MM[j]) || isnan(MM[j]))
          MM[j] = num_t(0);
        if((d3[j] - bd11[j]) < (d1[j] - bd10[j])) {
          bet1[j] = 0;
          bd10[j] = d1[j];
          bd11[j] = d3[j];
        }
        if((d2[j] - bd20[j]) < (d3[j] - bd21[j])) {
          bet2[j] = 0;
          bd20[j] = d2[j];
          bd21[j] = d3[j];
        }
        bet1[j] ++;
        bet2[j] ++;
        MM0 += MM[j] * num_t(bet2[j] - bet1[j]);
      }
      if(!isfinite(MM0) || isnan(MM0))
        MM0 = num_t(0);
    }
    std::cout << d00 << ", " << MM0 << std::endl;
  }
  return 0;
}

