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
typedef P0<num_t, idFeeder<num_t> > p0_0t;
typedef shrinkMatrix<num_t, p0_0t, true> p0_1t;
typedef northPole<num_t, p0_1t> p0_2t;
typedef sumChain< num_t, p0_2t> p0_3t;
typedef northPole<num_t, p0_3t> p0_4t;
typedef shrinkMatrix<num_t, p0_4t, true> p0_5t;
typedef sumChain< num_t, p0_5t, true> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  p0_t  p0(p0_5t(p0_4t(p0_3t(p0_2t(p0_1t(p0_0t(3)) ))) ));
  std::vector<p0_t> p;
  auto  q(p);
  num_t d(int(0));
  auto  Mp(d);
  auto  Mq(d);
  auto  Mx(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto Dp(d  * Mp);
    const auto Dq(Dp * Mq);
    if(Mx < abs(d) * num_t(int(2))) Mx = abs(d) * num_t(int(2));
    if(Mx == num_t(int(0))) {
      std::cout << Dq << ", " << Mp * Mq << std::endl << std::endl;
      continue;
    }
    p.emplace_back(p0);
    q.emplace_back(p0);
    Mp = Mq = num_t(int(0));
    // N.B. vanish first period.
    for(int i = 0; i < p.size(); i ++) Mp += max(- Mx, min(Mx, p[i].next(d )));
    // N.B. vanish double period causes no period on the stream.
    //      this results the right hand will no return in probability 1 ideally.
    for(int i = 0; i < q.size(); i ++) Mq += max(- Mx, min(Mx, q[i].next(Dp)));
    Mp /= num_t(int(p.size())) * Mx;
    Mq /= num_t(int(q.size())) * Mx;
    std::cout << Dq << ", " << Mp * Mq << std::endl << std::flush;
  }
  return 0;
}

