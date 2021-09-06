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

#include "lieonn.hh"
typedef myfloat num_t;
#include "p0.hh"

int main(int argc, const char* argv[]) {
  assert(1 < argc);
  const auto len(std::atoi(argv[1]));
  P0<num_t, idFeeder<num_t> > p(len);
  for(int i = 0; i < len; i ++)
    std::cout << p.p[i] << ", " << diff<num_t>(len)(len - 1, i) << std::endl;
  return 0;
}

