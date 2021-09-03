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
  int   var(3);
  int   step(1);
  if(argc < 2)
    std::cerr << "p0 <var>? <step>?" << std::endl;
  if(1 < argc) var  = std::atoi(argv[1]);
  if(2 < argc) step = std::atoi(argv[2]);
  std::cerr << "continue with p0 " << var << " " << step << std::endl;
  P0<num_t, linearFeeder<num_t, idFeeder<num_t> > > p(abs(var), step);
  P0<num_t, arctanFeeder<num_t, idFeeder<num_t> > > q(abs(var), step);
  num_t d(0);
  auto  D(d);
  std::vector<num_t> M;
  M.resize(step, num_t(0));
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    D  = d * M[0];
    // to make any sub-sequences to be the same meaning, no inverse condition,
    // this causes middle and high frequency parts to be ignored.
    for(int i = 1; i < M.size(); i ++)
      M[i - 1] = M[i];
    M[M.size() - 1] = var < 0 ? q.next(d) : p.next(d);
    std::cout << D << ", " << M[M.size() - 1] << ", " << M[0] << ", " << d << std::endl << std::flush;
  }
  return 0;
}

