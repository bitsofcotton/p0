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
  assert(0 < status);
  P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > > p0, p1, p2;
  P0recur<num_t, P0maxRank<num_t> > q0, q1, q2;
  p0 = P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > >(P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > >(P0recur<num_t, P0maxRank<num_t> >(status), status), status);
  p1 = P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > >(P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > >(P0recur<num_t, P0maxRank<num_t> >(status * 2), status * 2), status * 2);
  p2 = P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > >(P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > >(P0recur<num_t, P0maxRank<num_t> >(status * 3), status * 3), status * 3);
  q0 = P0recur<num_t, P0maxRank<num_t> >(status);
  q1 = P0recur<num_t, P0maxRank<num_t> >(status * 2);
  q2 = P0recur<num_t, P0maxRank<num_t> >(status * 3);
  const num_t zero(int(0));
  auto  d(zero);
  auto  M0(d);
  auto  M1(d);
  auto  M2(d);
  auto  M3(d);
  auto  M4(d);
  auto  M5(d);
  auto  U0(d);
  auto  U1(d);
  auto  U2(d);
  auto  U3(d);
  auto  U4(d);
  auto  U5(d);
  auto  S0(d);
  auto  S1(d);
  auto  S2(d);
  auto  S3(d);
  auto  S4(d);
  auto  S5(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D0(d * M0);
    const auto D1(d * M1);
    const auto D2(d * M2);
    const auto D3(d * M3);
    const auto D4(d * M4);
    const auto D5(d * M5);
    const auto r0(p0.next(d));
    const auto r1(p1.next(d));
    const auto r2(p2.next(d));
    const auto r3(q0.next(d));
    const auto r4(q1.next(d));
    const auto r5(q2.next(d));
    U0 = max(U0, abs(r0));
    U1 = max(U1, abs(r1));
    U2 = max(U2, abs(r2));
    U3 = max(U3, abs(r3));
    U4 = max(U4, abs(r4));
    U5 = max(U5, abs(r5));
    std::cout << D0 << ", " << D1 << ", " << D2;
    std::cout << ", " << D3 << ", " << D4 << ", " << D5;
    std::cout << ", " << (M0 = U0 == zero ? U0 : r0 / U0);
    std::cout << ", " << (M1 = U1 == zero ? U1 : r1 / U1);
    std::cout << ", " << (M2 = U2 == zero ? U2 : r2 / U2);
    std::cout << ", " << (M3 = U3 == zero ? U3 : r3 / U3);
    std::cout << ", " << (M4 = U4 == zero ? U4 : r4 / U4);
    std::cout << ", " << (M5 = U5 == zero ? U5 : r5 / U5);
    std::cout << ", " << (S0 += D0) << ", " << (S1 += D1);
    std::cout << ", " << (S2 += D2) << ", " << (S3 += D3);
    std::cout << ", " << (S4 += D4) << ", " << (S5 += D5);
    std::cout << std::endl << std::flush;
  }
  return 0;
}

