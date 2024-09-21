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
#include <sys/resource.h>

#include "lieonn.hh"
typedef myfloat num_t;

template <typename T> inline T logscale(const T& x) {
  return sgn<T>(x) * log(T(int(1)) + abs(x));
}

template <typename T> inline T expscale(const T& x) {
  return sgn<T>(x) * (exp(abs(x)) - T(int(1)) );
}

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  int progression(2);
  if(argc < 2) std::cerr << argv[0] << " <raw progression num>? : continue with ";
  if(1 < argc) progression = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << progression << std::endl;
  Pprogression<num_t, P0maxRank<num_t> > p0(progression, int(exp(num_t(2 * 2))) );
  Pprogression<num_t, P0maxRank<num_t> > p1(progression, int(exp(num_t(2 * 2))) );
  Pprogression<num_t, P0maxRank<num_t> > p2(progression, int(exp(num_t(2 * 2))) );
  Pprogression<num_t, P0maxRank<num_t> > p3(progression, int(exp(num_t(2 * 2))) );
  Pprogression<num_t, P0maxRank<num_t> > p4(progression, int(exp(num_t(2 * 2))) );
  std::string s;
  num_t d(int(0));
  auto  M0(d);
  auto  M1(d);
  auto  M2(d);
  auto  M3(d);
  auto  M4(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    std::cout << d * M0 << ", " << d * M1 << ", ";
    std::cout << d * M2 << ", " << d * M3 << ", ";
    std::cout << d * M4 << ", ";
    std::cout << (M0 = expscale<num_t>(expscale<num_t>(p1.next(logscale<num_t>(logscale<num_t>(d)) / num_t(int(abs(progression))) ) )) ) << ", ";
    std::cout << (M1 = expscale<num_t>(p0.next(logscale<num_t>(d)) / num_t(int(abs(progression))) )) << ", ";
    std::cout << (M2 = p2.next(d) / num_t(int(abs(progression))) ) << ", ";
    std::cout << (M3 = logscale<num_t>(p1.next(expscale<num_t>(d)) / num_t(int(abs(progression))))) << ", ";
    std::cout << (M4 = logscale<num_t>(logscale<num_t>(p1.next(expscale<num_t>(expscale<num_t>(d))) / num_t(int(abs(progression))) )) ) << std::endl << std::flush;
  }
  return 0;
}

