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

template <typename T, typename P> class P0D {
public:
  typedef SimpleVector<T> Vec;
  inline P0D() { ; }
  inline P0D(const int& size, const int& step = 1) {
    p = P(size, step);
    q = P(size, step);
    r[0] = r[1] = T(int(1));
    M[0] = M[1] = s[0] = s[1] = bin = T(t ^= t);
  }
  inline ~P0D() { ; };
  inline pair<T, T> next(const T& in) {
    t ++;
    s[t & 1] += in - bin;
    const auto MM((t & 1 ? p : q).next(s[t & 1]));
    if(M[t & 1] != T(int(0))) r[t & 1] *= M[t & 1] / (M[t & 1] + MM);
    M[t & 1] += MM;
    bin = in;
    return make_pair(r[t & 1] * MM, r[t & 1]);
  }
  int t;
  T s[2];
  T r[2];
  T M[2];
  T bin;
  P p;
  P q;
};

template <typename T, typename P> class P3 {
public:
  typedef SimpleVector<T> Vec;
  inline P3() { ; }
  inline P3(const int& size, const int& step = 1) {
    pd = P0D<T, P>(size, step);
    p  = P(size, step);
    Md = M = s = T(t ^= t);
    -- t;
  }
  inline ~P3() { ; };
  // we have only per 3 step glitch.
  inline T next(const T& in) {
    static const T zero(int(0));
    const auto bMd(Md);
    Md = pd.next(in).first;
    if(! ((++ t) % 3)) return zero;
    s += bMd * M * in;
    if(t % 3 == 2) {
      M = p.next(s);
      s = T(int(0));
      return Md * M;
    }
    return zero;
  }
  P0D<T, P> pd;
  P p;
  T M;
  T Md;
  T s;
  int t;
};

#define _P0_
#endif

