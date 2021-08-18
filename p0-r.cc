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
  for(int i = 0; i < 80; i ++)
    std::cout << nextP0<num_t>(80)[i] << ", " << diff<num_t>(80)(79, i) << std::endl;
  return 0;
}

