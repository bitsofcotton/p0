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
  // XXX:
  // P0normalizeStat<num_t, P0midLin<num_t, P0alignStart<num_t> > > p(P0midLin<num_t, P0alignStart<num_t> >(P0alignStart<num_t>(status) ) );
  P0normalizeStat<num_t, P0midLin<num_t, P0alignStart<num_t> > > p((P0midLin<num_t, P0alignStart<num_t> >(P0alignStart<num_t>(status) ) ));
  num_t d(int(0));
  auto  M(d);
  auto  S(d);
  auto  Mx(M);
  while(std::getline(std::cin, s, '\n')) {
    const auto bM(M);
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M);
    Mx = max(Mx, abs(d) * num_t(int(2)));
    std::cout << D << ", " << (M = max(- Mx, min(Mx, p.next(d) )) ) << ", " << (S += D) << std::endl << std::flush;
  }
  return 0;
}

