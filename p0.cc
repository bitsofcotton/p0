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

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(_MIMALLOC_)
#define MIMALLOC_OVERRIDE_H
#define MIMALLOC_NEW_DELETE_H
#include <mimalloc.h>
#endif

#if defined(_FLOAT_BITS_)
#define int int64_t
#endif
#include "lieonn.hh"
typedef myfloat num_t;
lieonn_t lieonn;

#if defined(_FLOAT_BITS_)
#undef int
#endif
int main(int argc, const char* argv[]) {
#if defined(_FLOAT_BITS_)
#define int int64_t
#endif
  std::cout << std::setprecision(30);
  int length(1);
  if(1 < argc) length = std::atoi(argv[1]);
  std::cerr << argv[0] << " " << length << std::endl;
  lieonnStaticInit();
  std::string s;
# if defined(_P_RAW_)
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
    if(chain) {
      for(int i = 0; i < M.size() - 1; i ++)
        std::cout << M[i] << ", ";
      std::cout << M[M.size() - 1] << endl;
    } else std::cout << flush;
    if(p.full) M = p0maxNext<num_t>(p.res);
    if(! chain) {
      for(int i = 0; i < M.size() - 1; i ++)
        std::cout << M[i] << ", ";
      std::cout << M[M.size() - 1] << endl;
    }
  }
#if !defined(_ONEBINARY_)
  lieonnStaticDestroy();
  return 0;
}
#endif

