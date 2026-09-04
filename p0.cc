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
#include <limits>
#include <assert.h>
#include <stdint.h>
#include <sys/resource.h>

#include "lieonn.hh"
typedef myfloat num_t;
lieonn_t lieonn;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  int length(1);
  if(1 < argc) length = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << length << std::endl;
  lieonnStaticInit();
  std::string s;
# if defined(_CHAIN_)
  const bool chain(true);
# else
  const bool chain(false);
# endif
#endif
  idFeeder<SimpleVector<num_t> > p(length);
  SimpleVector<num_t> M;
  while(std::getline(std::cin, s, '\n')) {
    SimpleVector<num_t> d(s2sv<num_t>(s));
    if(M.size() < d.size()) M.entity.resize(d.size(), num_t(int(0)) );
    for(int i = 0; i < d.size(); i ++)
      std::cout << (chain ? d[i] - M[i] : d[i] * M[i]) << ", ";
    p.next(d);
    SimpleVector<num_t> bM(M);
    if(p.full) {
      const SimpleVector<SimpleVector<num_t> > pp(p.res);
      for(int i = 0; i < d.size(); i ++) {
        idFeeder<num_t> buf(pp.size());
        for(int j = 0; j < pp.size(); j ++) buf.next(pp[j][i]);
        assert(buf.full);
        M[i] = p0maxNext<num_t>(buf.res);
      }
    }
    for(int i = 0; i < M.size() - 1; i ++)
      std::cout << (chain ? bM[i] : M[i]) << ", ";
    std::cout << (chain ? bM[bM.size() - 1] : M[M.size() - 1]) << std::endl << std::flush;
  }
#if !defined(_ONEBINARY_)
  lieonnStaticDestroy();
  return 0;
}
#endif

