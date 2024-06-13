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

template <typename T> T sqrtscale(const T& x) {
  return sgn<T>(x) * sqrt(abs(x));
}

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  assert(1 < argc);
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  vector<vector<SimpleMatrix<num_t> > > in;
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    in.emplace_back(move(work));
    assert(in[0].size() == in[in.size() - 1].size() &&
           in[0][0].rows() == in[in.size() - 1][0].rows() &&
           in[0][0].cols() == in[in.size() - 1][0].cols());
  }
  assert(~ (in.size() & 1));
  vector<SimpleMatrix<num_t> > res;
  res.resize(in[0].size(), SimpleMatrix<num_t>(in[0][0].rows(), in[0][0].cols()).O());
  for(int i = 0; i < in[0].size(); i ++) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static,1)
#endif
    for(int j = 0; j < in[0][i].rows(); j ++) {
      for(int k = 0; k < in[0][i].cols(); k ++) {
        PBond<num_t, P0maxRank<num_t> > p(P0maxRank<num_t>(), in.size() / 2);
        for(int m = 0; m < in.size() / 2; m ++)
          res[i](j, k) = sqrtscale<num_t>(p.next((in[m * 2][i](j, k) - num_t(int(1)) / num_t(int(2))) * (in[m * 2 + 1][i](j, k) - num_t(int(1)) / num_t(int(2))))) + num_t(int(1)) / num_t(int(2));
      }
    }
  }
  if(! savep2or3<num_t>("pout.ppm", normalize<num_t>(res)) )
    cerr << "failed to save." << endl;
  return 0;
}

