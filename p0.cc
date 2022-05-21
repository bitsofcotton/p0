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
  int status(60);
  if(1 < argc) status = std::atoi(argv[1]);
  assert(status);
  P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > > p;
  P0recur<num_t, P0maxRank<num_t> > q;
  if(0 < status)
    p = P0avg<num_t, P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > > >(P0measure<num_t, P0recur<num_t, P0maxRank<num_t> > >(P0recur<num_t, P0maxRank<num_t> >(abs(status)), abs(status)), abs(status));
  else
    q = P0recur<num_t, P0maxRank<num_t> >(abs(status));
  num_t d(int(0));
  auto  M(d);
  auto  S(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    std::cout << D << ", " << (M = status < 0 ? q.next(d) : p.next(d)) << ", " << (S += D) << std::endl << std::flush;
  }
  return 0;
}

