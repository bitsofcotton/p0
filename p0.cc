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

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  int status(3);
  if(argc < 2) std::cerr << argv[0] << " <status>? : continue with ";
  if(1 < argc) status = std::atoi(argv[1]);
  // N.B. p0 is valid when input is continuous.
  //      the condition dimension up to 3 is from v2v tanglement
  //      with separated input, so original cont. input isn't affect.
  //      the condition dimension up to 7 is also from v2v but non commutative.
  if(0 < status) status = max(3, min(7, status));
  std::cerr << argv[0] << " " << status << std::endl;
  idFeeder<num_t> f(max(1, abs(status)));
  PBond<num_t, P0maxRank<num_t> > p(P0maxRank<num_t>(), max(1, abs(status)));
  std::string s;
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    std::cout << d * M << ", ";
    if(! status) M = num_t(int(1));
    else if(status < 0) {
      const auto& ff(f.next(d));
      if(f.full) {
        M = num_t(int(0));
        for(int i = 0; i < ff.size(); i ++) M += ff[i];
        M = - M;
      }
    } else M = p.next(d);
    std::cout << M << std::endl << std::flush;
  }
  return 0;
}

