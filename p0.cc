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
  const auto var(10);
        int  step(1);
  if(argc < 2)
    std::cerr << "p0 <step>?" << std::endl;
  if(1 < argc) step = std::atoi(argv[1]);
  std::cerr << "continue with p0 " << step << std::endl;
  shrinkMatrix<num_t, P0<num_t, linearFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > > p(P0<num_t, linearFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > >(abs(var), abs(step)), abs(step));
  shrinkMatrix<num_t, P0<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > > q(P0<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > >(abs(var), abs(step)), abs(step));
  num_t d(0);
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    std::cout << D << ", " << (M = step < 0 ? q.next(d) : p.next(d)) << ", " << d << std::endl << std::flush;
  }
  return 0;
}

