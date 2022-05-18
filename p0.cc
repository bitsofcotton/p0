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
  int status(6);
  if(1 < argc) status = std::atoi(argv[1]);
  if(status < 0)
    status = int(max(num_t(3), ceil(exp(log(num_t(- status)) * log(num_t(- status))))));
  assert(status);
  P0recur<num_t, P0maxRank<num_t> > p(status);
  num_t d(int(0));
  auto  M(d);
  auto  S(d);
  struct rusage t0;
  struct rusage t1;
  getrusage(RUSAGE_SELF, &t0);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    getrusage(RUSAGE_SELF, &t1);
    std::cout << D << ", " << (M = p.next(d)) << ", " << (S += D) << ", " << (t1.ru_utime.tv_sec - t0.ru_utime.tv_sec) << ", " << (t1.ru_utime.tv_usec - t0.ru_utime.tv_usec) << ", " << (t1.ru_stime.tv_sec - t0.ru_stime.tv_sec) << ", " << (t1.ru_stime.tv_usec - t0.ru_stime.tv_usec) << std::endl << std::flush;
    t0 = t1;
  }
  return 0;
}

