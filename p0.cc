#if !defined(_ONEBINARY_)
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>
#include <stdint.h>
#include <sys/resource.h>

#include "lieonn.hh"
typedef myfloat num_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  int length(0);
  int step(0);
  if(argc < 2) std::cerr << argv[0] << " <length>? <step>? : continue with ";
  if(1 < argc) length = std::atoi(argv[1]);
  if(2 < argc) step   = std::atoi(argv[2]);
  assert(0 < step);
  std::cerr << argv[0] << " " << length << " " << step << std::endl;
  std::string s;
# if defined(_CHAIN_)
  const bool chain(true);
# else
  const bool chain(false);
# endif
#endif
  idFeeder<SimpleVector<num_t> > p(length * step);
  SimpleVector<num_t> d;
  idFeeder<SimpleVector<num_t> > MM(step);
  while(std::getline(std::cin, s, '\n')) {
    SimpleVector<num_t> M(MM.res[0]);
    int cnt(1);
    for(int i = 0; i < s.size(); i ++) if(s[i] == ',') cnt ++;
    d.resize(cnt);
    int i, j;
    for(i = 0, j = 0; i < s.size(); i ++) {
      std::stringstream ins(s.substr(i, s.size() - i));
      ins >> d[j ++];
      for( ; s[i] != ',' && i < s.size(); i ++) ;
    }
    if(M.size() < d.size()) M.entity.resize(d.size(), num_t(int(0)) );
    for(int i = 0; i < d.size(); i ++)
      std::cout << (chain ? d[i] - M[i] : d[i] * M[i]) << ", ";
    p.next(d);
    SimpleVector<num_t> bM(M);
    if(p.full) {
      const SimpleVector<SimpleVector<num_t> > pp(
        skipX<SimpleVector<num_t> >(p.res, step));
      for(int i = 0; i < d.size(); i ++) {
        idFeeder<num_t> buf(pp.size());
        for(int j = 0; j < pp.size(); j ++) buf.next(pp[j][i]);
        assert(buf.full);
        M[i] = p0maxNext<num_t>(buf.res);
      }
      MM.next(M);
    }
    for(int i = 0; i < M.size() - 1; i ++)
      std::cout << (chain ? bM[i] : M[i]) << ", ";
    std::cout << (chain ? bM[bM.size() - 1] : M[M.size() - 1]) << std::endl << std::flush;
    M = MM.res[0];
  }
#if !defined(_ONEBINARY_)
  return 0;
}
#endif

