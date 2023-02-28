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
  P0maxRank<num_t> p(status ? abs(status) : 1);
  idFeeder<num_t> f(abs(status));
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    if(! status) M = num_t(int(1));
    else if(status < 0) {
      const auto& ff(f.next(d));
      if(f.full) {
        M = num_t(int(0));
        for(int i = 0; i < ff.size(); i ++) M += ff[i];
        M = - M;
      }
    } else M = p.next(d);
    std::cout << D << ", " << M << std::endl << std::flush;
  }
  return 0;
}

