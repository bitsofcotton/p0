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
  int status(- 3);
  if(argc < 2) std::cerr << argv[0] << " <status>? : continue with ";
  if(1 < argc) status = std::atoi(argv[1]);
  const int recur(status < 0 ? num_t(int(1)) : ceil(- log(SimpleMatrix<num_t>().epsilon()) / num_t(int(2)) ) );
  std::cerr << argv[0] << " " << status << " (" << recur << ")" << std::endl;
  assert(status && 0 < recur);
  std::vector<std::vector<P0maxRank<num_t> > > p;
  const auto var0(max(int(1), int(exp(sqrt(log(num_t(abs(status))))))));
  std::vector<P0maxRank<num_t> > p0;
  p0.reserve(abs(status));
  for(int i = 1; i <= abs(status); i ++)
    p0.emplace_back(P0maxRank<num_t>(i,
                      max(int(1), int(exp(sqrt(log(num_t(i)))))) ));
  p.resize(recur, p0);
  int   t;
  num_t d(t ^= t);
  auto  M0(d);
  auto  S(d);
  std::vector<num_t> M;
  M.resize(p.size(), d);
  auto  Mx(M);
  std::vector<idFeeder<num_t> > bb;
  bb.resize(p.size(), idFeeder<num_t>(var0 * 2));
  t -= var0 * 2;
  while(std::getline(std::cin, s, '\n')) {
    const auto bM(M);
    std::stringstream ins(s);
    ins >> d;
    const auto D(d * M0);
    M0 = num_t(int(1));
    if(t < 0)
      for(int i = 0; i < p.size(); i ++)
        for(int j = 0; j < p[i].size(); j ++)
          p[i][j].next(d);
    else
      for(int i = 0; i < p.size(); i ++) {
        Mx[i] = max(Mx[i], abs(d) * num_t(int(2)));
        M[i]  = p[i][t].next(d)[0];
        M[i] /= num_t(max(int(1), int(exp(sqrt(log(num_t(t + 1))))) ));
        if(abs(M[i] = max(- Mx[i], min(Mx[i], M[i]))) == Mx[i])
          M[i] = num_t(int(0));
        M0    *= (M[i] = sgn<num_t>(M[i]) * pow(abs(M[i]), num_t(int(1)) / num_t(int(p.size()))) );
        bb[i].next(d);
        d     *= bM[i];
      }
    std::cout << D << ", " << M0 << ", " << (S += D) << std::endl << std::flush;
    if(t ++ < 0) continue;
    for(int i = 0; i < p.size(); i ++)
      for(int j = t; j < p0.size(); j ++)
        p[i][j].next(bb[i].res[bb[i].res.size() - 1]);
    if(p0.size() <= t) {
      t ^= t;
      p.resize(0);
      p.resize(recur, p0);
      for(int i = 0; i < p.size(); i ++) {
        const auto& bbb(bb[i].res);
        for(int j = 0; j < p[i].size(); j ++)
          for(int k = 0; k < bbb.size(); k ++)
            p[i][j].next(bbb[k]);
      }
    }
  }
  return 0;
}

