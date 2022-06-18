/*
 BSD 3-Clause License

Copyright (c) 2019-2022, bitsofcotton (kazunobu watatsu)
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
    p[0] = T(1) / T(2) + T(step < 0 ? 1 : 0);
    p[1] = T(1) / T(2) + T(step < 0 ? 0 : 1);
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
      p = taylor(size, T(step < 0 ? step : size + step - 1));
      cerr << "." << flush;
      if(abs(step) * 2 < size) {
        const auto pp(pnext<T>(size - 1, step));
        for(int j = 0; j < pp.size(); j ++)
          p[step < 0 ? j : j - pp.size() + p.size()] += pp[j] * T(size - 1);
        p /= T(size);
        ofstream ocache(file.c_str());
        ocache << p << endl;
        ocache.close();
      }
    }
  }
  assert(p.size() == size);
  return p;
}

template <typename T> SimpleVector<T> minsq(const int& size) {
  assert(1 < size);
  const T xsum(size * (size - 1) / 2);
  const T xdot(size * (size - 1) * (2 * size - 1) / 6);
  const auto denom(xdot * T(size) - xsum * xsum);
  SimpleVector<T> s(size);
  for(int i = 0; i < s.size(); i ++)
    s[i] = (T(i) * T(size) - xsum) / denom;
  return s;
}

template <typename T> const SimpleVector<T>& pnextcache(const int& size, const int& step) {
  assert(0 < size && 0 <= step);
  static vector<vector<SimpleVector<T> > > cp;
  if(cp.size() <= size) cp.resize(size + 1, vector<SimpleVector<T> >());
  if(cp[size].size() <= step) cp[size].resize(step + 1, SimpleVector<T>());
  if(cp[size][step].size()) return cp[size][step];
  return cp[size][step] = (pnext<T>(size, step) + pnext<T>(size, step + 1)) / T(int(2));
}

template <typename T> const SimpleVector<T>& mscache(const int& size) {
  assert(0 < size);
  static vector<SimpleVector<T> > ms;
  if(ms.size() <= size) ms.resize(size + 1, SimpleVector<T>());
  if(ms[size].size()) return ms[size];
  return ms[size] = minsq<T>(size);
}

template <typename T> class Pnull {
public:
  inline Pnull() { ; }
  inline ~Pnull() { ; };
  inline T next(const T& in) { return T(int(0)); }
};

template <typename T, typename feeder> class P0 {
public:
  inline P0() { ; }
  inline P0(const int& size, const int& step = 1) {
    f = feeder(size);
    this->step = step;
  }
  inline ~P0() { ; };
  inline T next(const T& in) {
    const auto& ff(f.next(in));
/*
    if(f.full) {
      auto avg(ff);
      T avgf(avg[0]);
      for(int i = 1; i < avg.size(); i ++) avgf += avg[i];
      for(int i = 0; i < avg.size(); i ++) avg[i] -= avgf / T(int(avg.size()));
      std::cerr << sqrt(avg.dot(avg) / ff.dot(ff)) << std::endl;
    }
*/
    return f.full ? pnextcache<T>(ff.size(), step).dot(ff) : T(int(0));
  }
  int step;
  feeder f;
};

template <typename T, typename P> class P0inv {
public:
  inline P0inv() { ; }
  inline P0inv(P&& p) { this->p = p; }
  inline ~P0inv() { ; }
  inline T next(const T& in) {
    static const T zero(int(0));
    static const T one(int(1));
    if(in == zero) return zero;
    const auto pn(p.next(one / in));
    if(pn == zero) return in;
    return one / pn;
  }
  P p;
};

template <typename T> const SimpleMatrix<complex<T> >& dftcache(const int& size) {
  assert(size != 0);
  static vector<SimpleMatrix<complex<T> > > cdft;
  static vector<SimpleMatrix<complex<T> > > cidft;
  if(0 < size) {
    if(cdft.size() <= size) cdft.resize(size + 1, SimpleMatrix<complex<T> >());
    if(cdft[size].rows() && cdft[size].cols()) return cdft[size];
    return cdft[size] = dft<T>(size);
  }
  if(cidft.size() <= abs(size)) cidft.resize(abs(size) + 1, SimpleMatrix<complex<T> >());
  if(cidft[abs(size)].rows() && cidft[abs(size)].cols()) return cidft[abs(size)];
  return cidft[abs(size)] = dft<T>(size);
}

template <typename T, typename P, typename feeder> class P0DFT {
public:
  inline P0DFT() { ; }
  inline P0DFT(P&& p, const int& size) {
    f = feeder(size);
    (this->p).resize(size, p);
    q = this->p;
  }
  inline ~P0DFT() { ; };
  inline T next(const T& in) {
    const auto& fn(f.next(in));
    if(! f.full) return T(int(0));
    auto ff(dftcache<T>(fn.size()) * fn.template cast<complex<T> >());
    assert(ff.size() == p.size() && p.size() == q.size());
    for(int i = 0; i < ff.size(); i ++)
      if(! (ff[i].real() == T(int(0)) && ff[i].imag() == T(int(0)) ) )
        ff[i] = abs(p[i].next(abs(ff[i]))) * exp(complex<T>(T(int(0)), q[i].next(arg(ff[i]))));
    return dftcache<T>(- fn.size()).row(fn.size() - 1).dot(ff).real();
  }
  vector<P> p;
  vector<P> q;
  feeder f;
};

template <typename T, typename P> class northPole {
public:
  inline northPole() { ; }
  inline northPole(P&& p) { this->p = p; }
  inline ~northPole() { ; }
  inline T next(const T& in) {
    static const T zero(int(0));
    static const T one(int(1));
    static const T M(atan(one / sqrt(SimpleMatrix<T>().epsilon())));
    if(! isfinite(in) || in == zero) return in;
    auto work(atan(in));
    // assert(- M < work && work < M);
    work = atan(one / work);
    // assert(- M < work && work < M);
    work = p.next(work);
    if(! isfinite(work) || work == zero) return in;
    work = tan(max(- M, min(M, one / tan(max(- M, min(M, work))))));
    if(isfinite(work)) return move(work);
    return in;
  }
  P p;
};

template <typename T, typename P, bool avg = false> class sumChain {
public:
  inline sumChain() { S = T(t ^= t); }
  inline sumChain(P&& p) { this->p = p; S = T(t ^= t); }
  inline ~sumChain() { ; }
  inline T next(const T& in) {
    S += in;
    if(! avg) return p.next(S) - S;
    const auto A(S / T(++ t));
    return p.next(in - A) + A;
  }
  T S;
  P p;
  myuint t;
};

template <typename T, typename P> class logChain {
public:
  inline logChain() { ; }
  inline logChain(P&& p) { this->p = p; S = T(int(0)); }
  inline ~logChain() { ; }
  inline T next(const T& in) {
    static const T zero(int(0));
    static const T one(int(1));
    const auto bS(S);
    S += in;
    if(bS == zero) return zero;
    const auto dd(S / bS - one);
    if(! isfinite(dd)) return zero;
    return p.next(dd) * S;
  }
  T S;
  P p;
};

template <typename T> class P0maxRank {
public:
  inline P0maxRank() { ; }
  inline P0maxRank(const int& status, const int& var) {
    assert(0 < status && 0 < var);
    p = p0_st(p0_s6t(p0_s5t(p0_s4t(p0_s3t(p0_s2t(p0_1t(p0_0t(status, var), var) )) ) )) );
    q = p0_it(p0_i6t(p0_i5t(p0_i4t(p0_i3t(p0_i2t(p0_i1t(p0_i0t(p0_0t(status, var) ), var) )) ) )) );
    r = p0_at(p0_a6t(p0_a5t(p0_a4t(p0_a3t(p0_a2t(p0_a1t() )) ) )) );
  }
  inline ~P0maxRank() { ; }
  inline T next(const T& in) {
    return (p.next(in) + q.next(in) + r.next(in)) / T(int(3));
  }
  // N.B. on existing taylor series.
  //      if the sampling frequency is not enough, middle range of the original
  //      function frequency (enough large bands) will effect prediction fail.
  //      this is because we only observes highest and lowest frequency on
  //      sampling points, so omitted part exists.
  //      even if the parameter on P0 is large, situation unchange.
  //      so we should use sectional measurement for them.
  typedef P0<T, idFeeder<T> > p0_0t;
  // N.B. sectional measurement, also expected value.
  typedef shrinkMatrix<T, p0_0t> p0_1t;
/*
  // N.B. make information-rich not to associative/commutative.
  //      2 dimension semi-order causes (x, status) from input as sedenion.
  typedef P0DFT<T, p0_1t, idFeeder<T> > p0_2t;
  typedef P0DFT<T, p0_2t, idFeeder<T> > p0_3t;
  typedef P0DFT<T, p0_3t, idFeeder<T> > p0_4t;
  typedef P0DFT<T, p0_4t, idFeeder<T> > p0_5t;
  // N.B. on any R to R into reasonable taylor.
  typedef northPole<T, p0_5t> p0_6t;
  typedef northPole<T, p0_6t> p0_7t;
  // N.B. we treat periodical part as non aligned complex arg part.
  typedef logChain<T, p0_7t>  p0_8t;
  typedef logChain<T, p0_8t>  p0_9t;
  // N.B. we make the prediction on (delta) summation.
  typedef sumChain<T, p0_9t>  p0_10t;
  // N.B. we take average as origin of input.
  typedef sumChain<T, p0_10t, true> p0_t;
  // N.B. if original sample lebesgue integrate is not enough continuous,
  //      imitate original function by some of sample points,
  //      but move origin point to average one, so a little better
  //      original function estimation.
  // N.B. frequency space *= 2 causes nyquist frequency ok.
  // N.B. but this is equivalent from jammer on PRNG, and probe on some
  //      measurable phenomenon.
  //typedef P0Expect<T, p0_11t> p0_t;
  // N.B. this needs huge memory to run.
*/
  // N.B. plain complex form.
  typedef northPole<T, p0_1t>  p0_s2t;
  typedef northPole<T, p0_s2t> p0_s3t;
  typedef logChain<T, p0_s3t>  p0_s4t;
  typedef logChain<T, p0_s4t>  p0_s5t;
  typedef sumChain<T, p0_s5t>  p0_s6t;
  typedef sumChain<T, p0_s6t, true> p0_st;

  typedef P0inv<T, p0_0t> p0_i0t;
  typedef shrinkMatrix<T, p0_i0t> p0_i1t;
  typedef northPole<T, p0_i1t> p0_i2t;
  typedef northPole<T, p0_i2t> p0_i3t;
  typedef logChain<T, p0_i3t> p0_i4t;
  typedef logChain<T, p0_i4t> p0_i5t;
  typedef sumChain<T, p0_i5t>  p0_i6t;
  typedef sumChain<T, p0_i6t, true> p0_it;

  typedef sumChain<T, Pnull<T>, true> p0_a1t;
  typedef northPole<T, p0_a1t>    p0_a2t;
  typedef northPole<T, p0_a2t>    p0_a3t;
  typedef logChain<T, p0_a3t>     p0_a4t;
  typedef logChain<T, p0_a4t>     p0_a5t;
  typedef sumChain<T, p0_a5t>     p0_a6t;
  typedef sumChain<T, p0_a6t, true> p0_at;
  p0_st p;
  p0_it q;
  p0_at r;
};

template <typename T, typename P> class P0normalizeStat {
public:
  inline P0normalizeStat() { ; }
  inline P0normalizeStat(P&& p) { this->p = p; }
  inline ~P0normalizeStat() { ; }
  inline T next(T in) {
    static const auto L(pow(SimpleMatrix<T>().epsilon(), - T(int(1)) / T(int(8))));
    static const int recur(max(T(int(1)), exp(sqrt(- log(SimpleMatrix<T>().epsilon()) ))) );
    for(int i = 0; i <= recur; i ++) in  = logscale(in * L);
    auto res(p.next(in));
    for(int i = 0; i <= recur; i ++) res = expscale(res) / L;
    return res;
  }
  inline T expscale(const T& in) {
    return sgn<T>(in) * (exp(abs(in)) - T(int(1)));
  }
  inline T logscale(const T& in) {
    return sgn<T>(in) * log(abs(in) + T(int(1)));
  }
  P p;
};

template <typename T, typename P> class P0midLin {
public:
  inline P0midLin() { ; }
  inline P0midLin(P&& p) { this->p = p; b = T(int(0)); }
  inline ~P0midLin() { ; }
  inline T next(const T& in) {
    p.next(b + in);
    return p.next(in + in) - (b = in);
  }
  T b;
  P p;
};

template <typename T> class P0alignStart {
public:
  inline P0alignStart() { ; }
  inline P0alignStart(const int& status) {
    assert(0 < status);
    p0.reserve(status);
    for(int i = 1; i <= status; i ++)
      p0.emplace_back(P0maxRank<T>(i,
        max(int(1), int(exp(sqrt(log(T(i)))))) ));
    p = p0;
    M = T(t ^= t);
    b = idFeeder<T>(abs(t -= max(int(1), int(exp(sqrt(log(T(status)))))) ));
  }
  inline ~P0alignStart() { ; }
  inline const T& next(const T& in) {
    if(t < 0)
      for(int i = 0; i < p.size(); i ++) p[i].next(in);
    else {
      M = p[t].next(in) / T(max(int(1), int(exp(sqrt(log(T(t + 1)))))));
      for(int i = 0; i < p.size(); i ++) if(t != i) p[i].next(in);
    }
    b.next(in);
    if(t ++ < 0) return M;
    if(p0.size() <= t) {
      t ^= t;
      p  = p0;
      const auto& bb(b.res);
      for(int i = 0; i < p.size(); i ++)
        for(int j = 0; j < bb.size(); j ++)
          p[i].next(bb[j]);
    }
    return M;
  }
  T M;
  int t;
  idFeeder<T> b;
  vector<P0maxRank<T> > p0;
  vector<P0maxRank<T> > p;
};

template <typename T, typename P, typename Q> class PthenQ {
public:
  inline PthenQ() { ; }
  inline PthenQ(P&& p, Q&& q) { this->p = p; this->q = q; M = T(int(0)); }
  inline ~PthenQ() { ; }
  inline T next(const T& in) {
    auto res(q.next(in * M));
    M = p.next(in);
    return res;
  }
  P p;
  Q q;
  T M;
};

#define _P0_
#endif

