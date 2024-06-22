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
  int progression(1);
  if(argc < 2) std::cerr << argv[0] << " <progression num>? : continue with ";
  if(1 < argc) progression = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << progression << std::endl;
  if(1 < abs(progression)) {
    const auto& pn(pnTinySingle(abs(progression) - 1));
    progression = sgn<int>(progression);
    for(int i = 0; i < pn.size(); i ++) progression *= pn[i];
    std::cerr << "using raw progression: " << progression << std::endl;
  }
  // N.B. p0 is valid when input is continuous.
  //      the condition dimension up to 3 is from v2v tanglement
  //      with separated input, so original cont. input isn't affect.
  //      the condition dimension up to 7 is also from v2v but non commutative.
  // N.B. since P0maxRank isn't linear, Pprogression works.
  Pprogression<num_t, P0maxRank<num_t> > p(progression ? progression : 1, 3);
  idFeeder<num_t> in(1);
  idFeeder<num_t> out(1);
  PBond<num_t, P0maxRank<num_t> > q(int(exp(num_t(2 * 2))), P0maxRank<num_t>(2));
  std::string s;
  num_t d(int(0));
  auto  M(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    std::cout << d * M << ", ";
    if(! progression) {
      const auto& inn(in.next(logscale<num_t>(logscale<num_t>(d)) ));
            auto  dd(inn[0]);
      for(int i = 1; i < inn.size(); i ++) dd += inn[i];
      const auto& outn(out.next(q.next(dd / num_t(int(inn.size())) )) );
      M = outn[0];
      for(int i = 1; i < outn.size(); i ++) M += outn[i];
      M = expscale<num_t>(expscale<num_t>(M / num_t(int(outn.size())) ));
    } else if(progression == - 1) M = num_t(int(1));
    else M = p.next(d);
    std::cout << M << std::endl << std::flush;
  }
  return 0;
}

