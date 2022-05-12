#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include <random>
#include <sys/resource.h>
#include "lieonn.hh"
#include "p0.hh"

/*
#if defined(_FLOAT_BITS_)
#define int int64_t
#endif
*/
typedef myfloat num_t;
/*
#if defined(_FLOAT_BITS_)
#undef int
#endif
*/
int main(int argc, const char* argv[]) {
/*
#if defined(_FLOAT_BITS_)
#define int int64_t
#endif
*/
  const auto method(std::atoi(argv[1]));
  const auto sum(std::atoi(argv[2]));
  const auto stat(std::atoi(argv[3]));
        auto line(std::atoi(argv[4]));
  std::random_device rd;
  std::mt19937_64 mt(rd());
  // cf. knuth_b for shuffle 128.
  std::shuffle_order_engine<std::linear_congruential_engine<unsigned int, 16807, 0, 2147483647>, 8> kb(rd());
  std::ranlux48 rl48(rd());
  std::string outbuf;
  P0recur<num_t, P0maxRank<num_t> > p(stat);
  int   t;
  num_t d(t ^= t);
  auto  M(d);
  auto  S(d);
  while(0 < line) {
    switch(method) {
    case 0:
      d += num_t(arc4random() & 0x7fffff) / (num_t(int(0x7fffff)) / num_t(int(2))) - num_t(int(1));
      break;
    case 1:
      d += num_t(int(mt()) & 0x7fffff) / (num_t(int(0x7fffff)) / num_t(int(2))) - num_t(int(1));
      break;
    case 2:
      d += num_t(int(kb()) & 0x7fffff) / (num_t(int(0x7fffff)) / num_t(int(2))) - num_t(int(1));
      break;
    case 3:
      d += num_t(int(rl48()) & 0x7fffff) / (num_t(int(0x7fffff)) / num_t(int(2))) - num_t(int(1));
      break;
    default:
      assert(0 && "Should not be reached.");
    }
    if(sum <= ++ t) {
      std::stringstream D(d * M);
      std::stringstream MM(M = p.next(d));
      std::stringstream SS(S += d * M);
      outbuf += D.str() + std::string(", ") + MM.str() + std::string(", ") + SS.str() + std::string("\n");
      d = num_t(t ^= t);
      -- line;
    }
  }
  std::cout << std::setprecision(30);
  std::cout << outbuf;
  return 0;
}

