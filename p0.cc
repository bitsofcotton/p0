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
  const auto var(max(int(1), int(exp(sqrt(log(num_t(status)))))));
  P0maxRank<num_t> p(status, var);;
  const num_t zero(int(0));
  auto  d(zero);
  std::vector<num_t> M;
  M.resize(4, d);
  auto  S(M);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    std::vector<num_t> D;
    D.reserve(M.size());
    for(int i = 0; i < M.size(); i ++) std::cout << (D[i] = d * M[i]) << ", ";
    const auto rp(p.next(d));
    assert(rp.size() == M.size());
    for(int i = 0; i < rp.size(); i ++) std::cout << (M[i] = rp[i]) << ", ";
    for(int i = 0; i < S.size(); i ++) std::cout << (S[i] += D[i]) << ", ";
    std::cout << std::endl << std::flush;
  }
  return 0;
}

