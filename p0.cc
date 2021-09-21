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
  const auto var(25);
        int  step(1);
  if(argc < 2)
    std::cerr << "p0 <step>?" << std::endl;
  if(1 < argc) step = std::atoi(argv[1]);
  std::cerr << "continue with p0 " << step << std::endl;
  P0<num_t, linearFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > p(abs(var), abs(step));
  P0<num_t, arctanFeeder<num_t, sumFeeder<num_t, idFeeder<num_t> > > > q(abs(var), abs(step));
  num_t d(0);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    // to make any sub-sequences to be the same meaning, no inverse condition,
    // this causes middle and high frequency parts to be ignored.
    std::cout << (step < 0 ? q.next(d) : p.next(d)) << ", " << d << std::endl << std::flush;
  }
  return 0;
}

