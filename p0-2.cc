#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <assert.h>

#include "ifloat.hh"
typedef myfloat num_t;
#include "simplelin.hh"
#include "p0.hh"

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  int range(20);
  if(1 < argc)
    range = std::atoi(argv[1]);
/*
  const auto pn(nextP0<num_t, true>(range));
*/
  SimpleVector<num_t> pn(range);
  for(int i = 0; i < pn.size(); i ++)
    pn[i] =
      cos((num_t(1) / num_t(8) + num_t(i) / num_t(pn.size())) * num_t(2) * num_t(4) * atan2(num_t(1), num_t(1)) * num_t(std::atoi(argv[2])));
//  const auto qn(diff<num_t>(- pn.size()) * diff<num_t>(pn.size()) * pn);
  const auto qn(diff<num_t>(pn.size()) * pn);
  for(int i = 0; i < pn.size(); i ++)
    std::cout << qn[i] << ", " << pn[i] << std::endl;
  std::cout << std::endl;
  return 0;
}

