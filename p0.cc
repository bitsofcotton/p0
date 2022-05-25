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
  assert(0 < status);
  P0maxRank<num_t> p(status, max(int(1), int(exp(sqrt(log(num_t(status)))))));
  num_t d(int(0));
  std::vector<num_t> D;
  D.resize(3, d);
  auto  M(D);
  auto  S(D);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    for(int i = 0; i < D.size(); i ++) std::cout << (D[i]  = d * M[i]) << ", ";
    M = p.next(d);
    for(int i = 0; i < D.size(); i ++) std::cout <<  M[i] << ", ";
    for(int i = 0; i < S.size(); i ++) std::cout << (S[i] += D[i]) << ", ";
    std::cout << std::endl << std::flush;
  }
  return 0;
}

