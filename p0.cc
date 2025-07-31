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

#define _COMPILE_PRED_
#include "lieonn.hh"
typedef myfloat num_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  int length(0);
  if(argc < 2) std::cerr << argv[0] << " <length>? : continue with ";
  if(1 < argc) length = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << length << std::endl;
  std::string s;
# if defined(_CHAIN_)
  const bool chain(true);
# else
  const bool chain(false);
# endif
#endif
  idFeeder<std::vector<num_t> > p(length);
  std::vector<num_t> d;
  std::vector<num_t> M;
  while(std::getline(std::cin, s, '\n')) {
    int cnt(1);
    for(int i = 0; i < s.size(); i ++) if(s[i] == ',') cnt ++;
    d.resize(cnt);
    int i, j;
    for(i = 0, j = 0; i < s.size(); i ++) {
      std::stringstream ins(s);
      ins >> d[j ++];
      for( ; s[i] != ',' && i < s.size(); i ++) ;
    }
    if(M.size() < d.size()) M.resize(d.size(), num_t(int(0)) );
    for(int i = 0; i < d.size(); i ++)
      std::cout << (chain ? d[i] - M[i] : d[i] * M[i]) << ", ";
    p.next(d);
    if(p.full)
      for(int i = 0; i < d.size(); i ++) {
        idFeeder<num_t> buf(p.res.size());
        for(int j = 0; j < p.res.size(); j ++) buf.next(p.res[j][i]);
        assert(buf.full);
        M[i] = p0maxNext<num_t>(buf.res);
      }
    for(int i = 0; i < M.size() - 1; i ++)
      std::cout << M[i] << ", ";
    std::cout << M[M.size() - 1] << std::endl << std::flush;
  }
#if !defined(_ONEBINARY_)
  return 0;
}
#endif

