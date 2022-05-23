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
#include "p0.hh"

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int status(20);
  if(1 < argc) status = std::atoi(argv[1]);
  assert(status);
  P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > > p0, p1, p2;
  P0recur<num_t, P0maxRank<num_t> > q0, q1, q2;
  if(0 < status) {
    p0 = P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > >(P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > >(P0recur<num_t, P0maxRank<num_t> >(abs(status)), abs(status)), abs(status));
    p1 = P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > >(P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > >(P0recur<num_t, P0maxRank<num_t> >(abs(status) * 2), abs(status) * 2), abs(status) * 2);
    p2 = P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > >(P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > >(P0recur<num_t, P0maxRank<num_t> >(abs(status) * 3), abs(status) * 3), abs(status) * 3);
  } else {
    q0 = P0recur<num_t, P0maxRank<num_t> >(abs(status));
    q1 = P0recur<num_t, P0maxRank<num_t> >(abs(status) * 2);
    q2 = P0recur<num_t, P0maxRank<num_t> >(abs(status) * 3);
  }
  num_t d(int(0));
  auto  M0(d);
  auto  M1(d);
  auto  M2(d);
  auto  MM(d);
  auto  S0(d);
  auto  S1(d);
  auto  S2(d);
  auto  SS(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D0(d * M0);
    const auto D1(d * M1);
    const auto D2(d * M2);
    const auto DD(d * MM);
    std::cout << D0 << ", " << D1 << ", " << D2 << ", " << DD;
    std::cout << ", " << (M0 = status < 0 ? q0.next(d) : p0.next(d));
    std::cout << ", " << (M1 = status < 0 ? q1.next(d) : p1.next(d));
    std::cout << ", " << (M2 = status < 0 ? q2.next(d) : p2.next(d));
    std::cout << ", " << (MM = (M0 + M1 + M2) / num_t(int(3)));
    std::cout << ", " << (S0 += D0) << ", " << (S1 += D1);
    std::cout << ", " << (S2 += D2) << ", " << (SS += DD);
    std::cout << std::endl << std::flush;
  }
  return 0;
}

