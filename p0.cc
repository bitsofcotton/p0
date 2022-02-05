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
  auto var(1);
  auto step(1);
  if(argc < 2)
    std::cerr << argv[0] << " <len>?" << std::endl;
  if(1 < argc) var  = std::atoi(argv[1]);
  std::cerr << "continue with " << argv[0] << " " << var << std::endl;
  if(var < 0) var = - var;
  else {
    step = var;
    var  = 3;
  }
  std::vector<P0D<num_t, P0<num_t, idFeeder<num_t> > > > p;
  p.resize(step, P0D<num_t, P0<num_t, idFeeder<num_t> > >(var));
  int   t;
  num_t d(t ^= t);
  auto  S(d);
  auto  MM(d);
  std::vector<num_t> M;
  M.resize(step, d);
  auto  A(M);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * MM);
    S += d;
    for(int i = 0; i < A.size(); i ++)
      A[i] += S;
    const auto tt((t ++) % M.size());
    t    %= M.size();
    M[tt] = p[tt].next(A[tt]) - A[tt] * num_t(int(2)) + S * num_t(step);;
    A[tt] = num_t(int(0));
    MM    = M[0];
    for(int i = 1; i < M.size(); i ++) MM += M[i];
    std::cout << D << ", " << (MM /= num_t(step)) << std::endl << std::flush;
  }
  return 0;
}

