/*
 BSD 3-Clause License

Copyright (c) 2019-2020, bitsofcotton
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

template <typename T> class P0 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleVector<complex<T> > VecU;
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<complex<T> > MatU;
  inline P0();
  inline ~P0();
  const Mat&  nextP(const int& size);
private:
  const MatU& seed(const int& size, const bool& idft);
  const Mat&  diff(const int& size);
  inline Vec  taylor(const int& size, const T& step);
  const T&    Pi() const;
  const complex<T>& J() const;
};

template <typename T> inline P0<T>::P0() {
  ;
}

template <typename T> inline P0<T>::~P0() {
  ;
}

template <typename T> const T& P0<T>::Pi() const {
  const static auto pi(atan2(T(1), T(1)) * T(4));
  return pi;
}

template <typename T> const complex<T>& P0<T>::J() const {
  const static auto i(complex<T>(T(0), T(1)));
  return i;
}

template <typename T> const typename P0<T>::MatU& P0<T>::seed(const int& size, const bool& f_idft) {
  assert(0 < size);
  static vector<MatU> dft;
  static vector<MatU> idft;
  if(dft.size() <= size)
    dft.resize(size + 1, MatU());
  if(idft.size() <= size)
    idft.resize(size + 1, MatU());
  if((!f_idft) &&  dft[size].rows() == size &&  dft[size].cols() == size)
    return dft[size];
  if(  f_idft  && idft[size].rows() == size && idft[size].cols() == size)
    return idft[size];
  auto& edft( dft[size]);
  auto& eidft(idft[size]);
  edft.resize(size, size);
  eidft.resize(size, size);
  for(int i = 0; i < edft.rows(); i ++)
    for(int j = 0; j < edft.cols(); j ++) {
      const auto theta(- T(2) * Pi() * T(i) * T(j) / T(edft.rows()));
      edft(i, j)  = complex<T>(cos(  theta), sin(  theta));
      eidft(i, j) = complex<T>(cos(- theta), sin(- theta)) / complex<T>(T(size));
    }
  if(f_idft)
    return eidft;
  return edft;
}

template <typename T> const typename P0<T>::Mat& P0<T>::diff(const int& size) {
  assert(0 < size);
  static vector<Mat> D;
  if(D.size() <= size)
    D.resize(size + 1, Mat());
  if(D[size].rows() == size && D[size].cols() == size)
    return D[size];
  auto& d(D[size]);
  d.resize(size, size);
  auto  DD(seed(size, false));
  DD.row(0) *= complex<T>(T(0));
  T nd(0);
  T ni(0);
  for(int i = 1; i < DD.rows(); i ++) {
    const auto phase(- J() * T(2) * Pi() * T(i) / T(DD.rows()));
    const auto phase2(complex<T>(T(1)) / phase);
    DD.row(i) *= phase;
    nd += abs(phase)  * abs(phase);
    ni += abs(phase2) * abs(phase2);
  }
  return d = (seed(size, true) * DD).template real<T>() * sqrt(sqrt(T(DD.rows() - 1) / (nd * ni)));
}

template <typename T> inline typename P0<T>::Vec P0<T>::taylor(const int& size, const T& step) {
  const auto  spt(min(size - 1, std::max(0, int(std::ceil(step)))));
  const auto  residue(step - T(spt));
        Vec   tayl(size);
  const auto& D(diff(size));
        auto  dt(D * residue);
  for(int i = 0; i < tayl.size(); i ++)
    tayl[i] = T(0);
  tayl[spt] = T(1);
  for(int i = 2; 0 <= i; i ++) {
    const auto bt(tayl);
    tayl += dt.row(spt);
    if(bt == tayl)
      break;
    dt = (D * dt) * (residue / T(i));
  }
  return tayl;
}

template <typename T> const typename P0<T>::Mat& P0<T>::nextP(const int& size) {
  assert(1 < size);
  static vector<Mat> P;
  if(P.size() <= size)
    P.resize(size + 1, Mat());
  if(P[size].rows() == size && P[size].cols() == size)
    return P[size];
  auto& p(P[size]);
  p.resize(size, size);
  Mat   extends(p.cols() * 2 - 1, p.cols());
  Mat   revextends(extends.rows(), extends.cols());
  for(int i = 0; i < extends.rows(); i ++)
    for(int j = 0; j < extends.cols(); j ++)
      extends(i, j) = i / 2 == j && i % 2 == 0 ? T(1) : T(0);
  for(int i = 0; i < p.cols() - 1; i ++)
    extends.row(i * 2 + 1) = taylor(p.cols(), T(i) + T(1) / T(2));
  for(int i = 0; i < extends.rows(); i ++)
    for(int j = 0; j < extends.cols(); j ++)
      revextends(i, j) = extends(extends.rows() - 1 - i,
                                 extends.cols() - 1 - j);
  const auto reverse(revextends.transpose() * taylor(p.cols() * 2 - 1, - T(2)));
  p.row(0) = extends.transpose() * taylor(p.cols() * 2 - 1, T(p.cols() * 2));
  for(int i = 0; i < reverse.size(); i ++)
    p(0, i) += reverse[reverse.size() - i - 1];
  p.row(0) /= T(2);
  for(int i = 1; i < p.rows(); i ++)
    for(int j = 0; j < p.cols(); j ++)
      p(i, j) = p(0, (p.cols() - i + j) % p.cols());
  return p;
}


template <typename T> class P0B {
public:
  typedef SimpleVector<T> Vec;
  inline P0B();
  inline P0B(const int& size);
  inline ~P0B();
  inline T next(const T& in);
private:
  P0<T> p;
  Vec   buf;
  int   t;
};

template <typename T> inline P0B<T>::P0B() {
  t = 0;
}

template <typename T> inline P0B<T>::P0B(const int& size) {
  buf.resize(size);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  t = 0;
}

template <typename T> inline P0B<T>::~P0B() {
  ;
}

template <typename T> inline T P0B<T>::next(const T& in) {
  buf[(t ++) % buf.size()] = in;
  return p.nextP(buf.size()).row(t % buf.size()).dot(buf);
}


template <typename T, typename U> class P0C {
public:
  typedef SimpleVector<T> Vec;
  inline P0C();
  inline P0C(const int& size, const int& loop);
  inline ~P0C();
  T next(const T& in, const int& idx = 0);
private:
  std::vector<U> p;
};

template <typename T, typename U> inline P0C<T,U>::P0C() {
  ;
}

template <typename T, typename U> inline P0C<T,U>::P0C(const int& size, const int& loop) {
  assert(1 < size && 1 < loop);
  p.resize(pow(2, loop) - 1, U(size));
}

template <typename T, typename U> inline P0C<T,U>::~P0C() {
  ;
}

template <typename T, typename U> T P0C<T,U>::next(const T& in, const int& idx) {
  const static T halfPi(atan2(T(1), T(1)) * T(4) / T(2));
  const auto inpi(in * halfPi);
  // 0 -> 1, 2 -> ... -> 2^n - 1, ..., 2^(n + 1) - 2
  const auto stt(((p.size() + 2) / 2 - 2 + 2) / 2 - 1);
  if(idx < stt) {
    const auto M(atan2(next(cos(inpi), 2 * idx + 1), next(sin(inpi), 2 * idx + 2)) - inpi);
    return atan2(cos(M), sin(M)) / halfPi + in;
  }
  const auto M(atan2(p[2 * idx - stt].next(cos(inpi)), p[2 * idx + 1 - stt].next(sin(inpi))) - inpi);
  return atan2(cos(M), sin(M)) / halfPi + in;
}

#define _P0_
#endif

