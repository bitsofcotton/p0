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
  std::cout << std::setprecision(30);
  std::string s;
  auto var(3);
  auto step(1);
  if(argc < 2)
    std::cerr << argv[0] << " <len>? <step>?" << std::endl;
  if(1 < argc) var  = std::atoi(argv[1]);
  if(2 < argc) step = std::atoi(argv[2]);
  std::cerr << "continue with " << argv[0] << " " << var << " " << step << std::endl;
  P0D<num_t, P0<num_t, idFeeder<num_t> > > p(abs(var), abs(step));
  P0D<num_t, P0<num_t, deltaFeeder<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > > > q(abs(var), abs(step));
  int   t;
  num_t S(t ^= t);
  vector<num_t> M;
  M.resize(abs(step), num_t(int(0)));
  auto  d(M);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d[t % d.size()];
    auto D(d[t % d.size()] * M[t % M.size()]);
    if(0 < step)
      for(int i = 1; i < M.size(); i ++)
        D += d[(t + i) % d.size()] * M[t % M.size()];
    const auto S0(S);
    const auto t0((t ++) % d.size());
    std::cout << D << ", " << (M[t % M.size()] = step < 0 ? (var < 0 ? q.next(d[t0]) : p.next(d[t0])) : (var < 0 ? q.next(S += d[t0]) - S0 : p.next(S += d[t0]) - S0)) << ", " << s.substr((int)ins.tellg() + 1) << std::endl << std::flush;
  }
  return 0;
}

