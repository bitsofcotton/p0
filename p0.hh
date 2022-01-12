/*
 BSD 3-Clause License

Copyright (c) 2019-2021, bitsofcotton (kazunobu watatsu)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#if !defined(_P0_)

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::flush;
using std::endl;

template <typename T> SimpleVector<T> pnext(const int& size, const int& step = 1) {
  SimpleVector<T> p;
  if(size <= 1) {
    p.resize(1);
    p[0] = T(1);
  } else if(size <= 2) {
    p.resize(2);
    p[0] = T(1) / T(2);
    p[1] = T(1) / T(2) + T(1);
    p   /= T(2);
  } else {
    const auto file(std::string("./.cache/lieonn/next-") +
      std::to_string(size) + std::string("-") + std::to_string(step) +
#if defined(_FLOAT_BITS_)
      std::string("-") + std::to_string(_FLOAT_BITS_)
#else
      std::string("-ld")
#endif
    );
    ifstream cache(file.c_str());
    if(cache.is_open()) {
      cache >> p;
      cache.close();
    } else {
      p = taylor(size, T(size + step - 1));
      cerr << "." << flush;
      const auto pp(pnext<T>(size - 1, step));
      for(int j = 0; j < pp.size(); j ++)
        p[j - pp.size() + p.size()] += pp[j] * T(size - 1);
      p /= T(size);
      ofstream ocache(file.c_str());
      ocache << p << endl;
      ocache.close();
    }
  }
  assert(p.size() == size);
  return p;
}

template <typename T, typename feeder> class P0 {
public:
  typedef SimpleVector<T> Vec;
  inline P0() { ; }
  inline P0(const int& size, const int& step = 1) {
    f = feeder(size);
    p = pnext<T>(size, step);
  }
  inline ~P0() { ; };
  inline T next(const T& in) {
    const auto& ff(f.next(in));
    static const T zero(int(0));
    return f.full ? p.dot(ff) : zero;
  }
  Vec p;
  feeder f;
};

template <typename T, typename P> class P0Dsgn {
public:
  typedef SimpleVector<T> Vec;
  inline P0Dsgn() { ; }
  inline P0Dsgn(const int& size, const int& step = 1, const int& recur = 1) {
    p.resize(recur * recur, P(size, step));
    q.resize(p.size(), P(size, step));
    r.resize(p.size(), P(size, step));
    brnd.resize(p.size(), b = T(int(0)));
    vbrnd.resize(recur, b);
  }
  inline ~P0Dsgn() { ; };
  inline T next(const T& in) {
    assert(p.size() == q.size() && q.size() == r.size() &&
           p.size() == brnd.size());
    auto rnd(brnd);
    auto vrnd(vbrnd);
    for(int i = 0; i < rnd.size(); i ++)
      rnd[i]  = T(int(arc4random_uniform(0x80000001))) /
                                   T(int(0x80000000 + 1));
    for(int i = 0; i < vrnd.size(); i ++)
      vrnd[i] = T(int(arc4random_uniform(0x80000001))) /
                                   T(int(0x80000000 + 1));
    T res(int(0));
    for(int i = 0; i < vrnd.size(); i ++) {
      T pp(int(0));
      T qq(int(0));
      T rr(int(0));
      for(int j = 0; j < vrnd.size(); j ++) {
        const auto din(in * rnd[i * vrnd.size() + j] * vrnd[i]);
        const auto ddelta(din - b * brnd[i * vrnd.size() + j] * vbrnd[i]);
        pp += sgn<T>(p[i * vrnd.size() + j].next( ddelta));
        qq += abs(q[i * vrnd.size() + j].next(abs(ddelta)));
        rr += abs(r[i * vrnd.size() + j].next(abs(din)));
      }
      res += sgn<T>(- sgn<T>(pp) * abs(qq) / T(int(vrnd.size())) * T(int(2)) + in * vrnd[i]) * abs(rr) / T(int(vrnd.size())) * T(int(2));
    }
    b = in;
    brnd  = rnd;
    vbrnd = vrnd;
    return res /= T(int(vrnd.size())) / T(int(2));
  }
  vector<P> p;
  vector<P> q;
  vector<P> r;
  vector<T> brnd;
  vector<T> vbrnd;
  T b;
};

#define _P0_
#endif

