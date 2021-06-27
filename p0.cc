#include <cstdio>
#include <cstring>
#include <iostream>
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
  std::cout << std::setprecision(30);
  std::string s;
  int  var(3);
  bool whole(false);
  if(argc < 2)
    std::cerr << "p0 <var>? <whole|partial>?" << std::endl;
  if(1 < argc) var   = std::atoi(argv[1]);
  if(2 < argc) whole = argv[2][0] == 'w';
  std::cerr << "continue with p0 " << var <<  " " << (const char*)(whole ? "whole" : "partial") << std::endl;
  P0<num_t, linearFeeder<num_t>, true>  pp0(abs(var));
  P0<num_t, linearFeeder<num_t>, false> qp0(abs(var));
  P0<num_t, arctanFeeder<num_t, true>, true>  pw0(abs(var));
  P0<num_t, arctanFeeder<num_t, true>, false> qw0(abs(var));
  auto  pp1(pp0);
  auto  qp1(qp0);
  auto  pw1(pw0);
  auto  qw1(qw0);
  num_t d(0);
  auto  s0(d);
  auto  s1(d);
  auto  s2(d);
  auto  s3(d);
  auto  M(d);
  auto  MM(d);
  auto  bd(d);
  int   t(0);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto bdd(d);
    const auto bMM(MM);
    if(MM < abs(d))
      MM = abs(d) * num_t(2);
    d -= bd;
    if(bMM != num_t(0))
      d /= bMM;
    if(d != num_t(0)) {
      if(M != num_t(0)) {
        s0 += (s3 = d - M);
        s1 += (s2 = d * M);
      }
      const auto pn((var < 0 ? (whole ? qw0.next(d) : qp0.next(d))
                             : (whole ? pw0.next(d) : pp0.next(d))));
      // original function lower and higher frequency part, middle is ignored.
      if(d == num_t(0)) M = pn;
      else if(whole) {
         const auto denom(var < 0 ? qw1.next(num_t(1) / d)
                                  : pw1.next(num_t(1) / d));
         M = denom == num_t(0) ? pn :
           (pn + num_t(1) / denom) / num_t(2);
      } else {
         const auto denom(var < 0 ? qp1.next(num_t(1) / d)
                                  : pp1.next(num_t(1) / d));
         M = denom == num_t(0) ? pn :
           (pn + num_t(1) / denom) / num_t(2);
      }
      if(! isfinite(M) || isnan(M)) M = pn;
      if(t ++ < abs(var)) M = num_t(0);
    }
    std::cout << (MM != num_t(0) ? M / MM : M) << "," << s0 << ", " << s1 << ", " << s2 << ", " << s3 << ", " << std::endl << std::flush;
    bd = bdd;
  }
  return 0;
}

