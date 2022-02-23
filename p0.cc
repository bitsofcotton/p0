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
// N.B. on existing taylor series.
typedef P0<num_t, idFeeder<num_t> > p0_0t;
// N.B. on sectional measurement.
typedef shrinkMatrix<num_t, p0_0t, true> p0_1t;
// N.B. on any of L2(R).
typedef northPole<num_t, p0_1t> p0_2t;
// N.B. on probability.
typedef shrinkMatrix<num_t, p0_2t, true> p0_t;

int main(int argc, const char* argv[]) {
  std::cout << std::setprecision(30);
  std::string s;
  p0_t   p0(p0_2t(p0_1t(p0_0t(3)) ));
  std::vector<p0_t> p;
  auto   q(p);
  num_t  d(int(0));
  auto   Mp(d);
  auto   Mq(d);
  auto   Mx(d);
  while(std::getline(std::cin, s, '\n')) {
    std::stringstream ins(s);
    ins >> d;
    const auto Dp(d  * Mp);
    const auto Dq(Dp * Mq);
    if(Mx < abs(d)) Mx = abs(d) * num_t(int(2));
    if(Mx == num_t(int(0))) {
      std::cout << Dq << ", " << Mp * Mq << std::endl << std::endl;
      continue;
    }
    p.emplace_back(p0);
    q.emplace_back(p0);
    int  cntp;
    auto cntq(cntp ^= cntp);
    Mp = Mq = num_t(cntp);
    // N.B. vanish periods.
    for(int i = 0; i < p.size(); i ++) {
      auto pi(p[i].next(d));
      if(isfinite(pi) && - Mx / sqrt(sqrt(SimpleMatrix<num_t>().epsilon)) < pi
                   && pi < Mx / sqrt(sqrt(SimpleMatrix<num_t>().epsilon)) 
                   && pi != num_t(int(0)) ) {
        Mp += max(- Mx, min(Mx, move(pi))); cntp ++; }
    }
    // N.B. vanish double period causes no period on the stream.
    //      this results the right hand will no return in probability 1 ideally.
    for(int i = 0; i < q.size(); i ++) {
      auto qi(q[i].next(Dp));
      if(isfinite(qi) && - Mx / sqrt(sqrt(SimpleMatrix<num_t>().epsilon)) < qi
                   && qi < Mx / sqrt(sqrt(SimpleMatrix<num_t>().epsilon))
                   && qi != num_t(int(0)) ) {
        Mq += max(- Mx, min(Mx, move(qi))); cntq ++; }
    }
    Mp /= num_t(cntp ? cntp : 1) * Mx;
    Mq /= num_t(cntq ? cntq : 1) * Mx;
    std::cout << Dq << ", " << Mp * Mq << std::endl << std::flush;
  }
  return 0;
}

