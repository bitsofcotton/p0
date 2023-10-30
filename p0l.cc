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

template <typename T> inline T logscale(const T& x) {
  return sgn<T>(x) * log(T(int(1)) + abs(x));
}

template <typename T> inline T expscale(const T& x) {
  return sgn<T>(x) * (exp(abs(x)) - T(int(1)) );
}

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  int status(2);
  if(argc < 2) std::cerr << argv[0] << " <status>? : continue with ";
  if(1 < argc) status = std::atoi(argv[1]);
  assert(0 < status && status <= 3);
  // N.B. p0 is valid when input is continuous.
  //      the condition dimension up to 3 is from v2v tanglement
  //      with separated input, so original cont. input isn't affect.
  //      the condition dimension up to 7 is also from v2v but non commutative.
  std::cerr << argv[0] << " " << status << std::endl;
  idFeeder<num_t> in(status);
  idFeeder<num_t> out(status);
  PBond<num_t, P0maxRank<num_t> > p(P0maxRank<num_t>(status), exp(num_t(status * status)) );
  std::string s;
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    std::cout << d * M << ", ";
    const auto& inn(in.next(logscale<num_t>(logscale<num_t>(d)) ));
          auto  dd(inn[0]);
    for(int i = 1; i < inn.size(); i ++) dd += inn[i];
    const auto& outn(out.next(p.next(dd / num_t(int(inn.size())) )) );
    M = outn[0];
    for(int i = 1; i < outn.size(); i ++) M += outn[i];
    std::cout << (M = expscale<num_t>(expscale<num_t>(M / num_t(int(outn.size())) )) ) << std::endl << std::flush;
  }
  return 0;
}

