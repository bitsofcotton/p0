#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>

#include "lieonn.hh"
typedef myfloat num_t;
#include "p0.hh"

int main(int argc, const char* argv[]) {
  assert(1 < argc);
  const auto len(std::atoi(argv[1]));
  assert(0 < len);
  const auto step(argc <= 2 ? max(int(1), int(exp(sqrt(log(num_t(abs(len)))))) ) : std::atoi(argv[2]));
  P0<num_t, idFeeder<num_t> > p(len);
  const auto d(diff<num_t>(len).row(len - 1));
  const auto n(pnext<num_t>(len, step));
  const auto m(pnext<num_t>(len, step + 1));
  num_t sumn(int(0));
  auto  summ(sumn);
  auto  sumd(sumn);
  auto  n2n(sumn);
  auto  n2m(sumn);
  auto  n2d(sumn);
  for(int i = 0; i < len; i ++) {
    std::cout << n[i] << ", " << m[i] << ", " << d[i] << std::endl;
    sumn += n[i];
    summ += m[i];
    sumd += d[i];
    n2n  += n[i] * n[i];
    n2m  += m[i] * m[i];
    n2d  += d[i] * d[i];
  }
  std::cout << sumn << ", " << summ << ", " << sumd << ", " << sqrt(n2n) << ", " << sqrt(n2m) << ", " << sqrt(n2d) << std::endl;
  return 0;
}

