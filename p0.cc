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
  int status(3);
  if(argc < 2) std::cerr << argv[0] << " <status>? : continue with ";
  if(1 < argc) status = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << status << std::endl;
  assert(status);
  std::vector<P0maxRank<num_t> > p;
  const auto var0(max(int(1), int(exp(sqrt(log(num_t(abs(status))))))));
  p.reserve(abs(status) + abs(var0) * 2);
  for(int i = 0; i <= abs(var0) * 2; i ++)
    p.emplace_back(P0maxRank<num_t>(3, 1));
  for(int i = 1; i <= abs(status); i ++) {
    const auto var(max(int(1), int(exp(sqrt(log(num_t(i)))))));
    p.emplace_back(P0maxRank<num_t>(status < 0 ? - i : i, var));
  }
  int   t;
  num_t d(t ^= t);
  auto  M(d);
  auto  S(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    std::cout << D << ", " << (M = p[t ++].next(d)[0]) << ", " << (S += D) << std::endl << std::flush;
    for(int i = t; i < p.size(); i ++) p[i].next(d);
    if(p.size() <= t) {
      t ^= t;
      for(int i = 1; i <= abs(status); i ++) {
        const auto var(max(int(1), int(exp(sqrt(log(num_t(i)))))));
        p[i + abs(var0) * 2 - 1] = P0maxRank<num_t>(status < 0 ? - i : i, var);
      }
    }
  }
  return 0;
}

