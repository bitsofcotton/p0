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

template <typename T, bool denoise = false> const SimpleVector<T>& nextP0(const int& size) {
  static vector<SimpleVector<T> > P;
  if(P.size() <= size)
    P.resize(size + 1, SimpleVector<T>());
  auto& p(P[size]);
  if(p.size() != size) {
    if(size <= 1) {
      p    = SimpleVector<T>(1);
      p[0] = T(1);
    } else {
      p = taylor<T>(size, T(size));
      if(denoise) {
        std::cerr << "." << std::flush;
        const auto& pp(nextP0<T, denoise>(size - 1));
        for(int i = 0; i < pp.size(); i ++)
          p[i - pp.size() + p.size()] += pp[i] * T(size - 1);
        p /= T(size);
      }
    }
  }
  return p;
}

template <typename T, bool denoise = false> class P0 {
public:
  typedef SimpleVector<T> Vec;
  inline P0();
  inline P0(const int& size);
  inline ~P0();
  inline T next(const T& in);
private:
  Vec buf;
};

template <typename T, bool denoise> inline P0<T, denoise>::P0() {
  ;
}

template <typename T, bool denoise> inline P0<T, denoise>::P0(const int& size) {
  assert(0 < size);
  buf.resize(size);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
}

template <typename T, bool denoise> inline P0<T, denoise>::~P0() {
  ;
}

template <typename T, bool denoise> inline T P0<T, denoise>::next(const T& in) {
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i]  = buf[i + 1];
  buf[buf.size() - 1] = atan(in);
  return tan(nextP0<T, denoise>(buf.size()).dot(buf));
}


template <typename T, bool denoise = false> class P0W {
public:
  typedef SimpleVector<T> Vec;
  inline P0W();
  inline P0W(const int& size);
  inline ~P0W();
  inline T next(const T& in);
private:
  Vec buf;
};

template <typename T, bool denoise> inline P0W<T, denoise>::P0W() {
  ;
}

template <typename T, bool denoise> inline P0W<T, denoise>::P0W(const int& size) {
  assert(0 < size);
  buf.resize(size);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
}

template <typename T, bool denoise> inline P0W<T, denoise>::~P0W() {
  ;
}

template <typename T, bool denoise> inline T P0W<T, denoise>::next(const T& in) {
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i]  = buf[i + 1];
  buf[buf.size() - 1] = in;
  auto avg(buf[0]);
  for(int i = 1; i < buf.size(); i ++)
    avg += buf[i];
  avg /= T(buf.size());
  auto buf2(buf);
  for(int i = 0; i < buf2.size(); i ++)
    buf2[i] -= avg;
  for(int i = 1; i < buf2.size(); i ++)
    buf2[i] += buf2[i - 1];
  for(int i = 0; i < buf2.size(); i ++)
    buf2[i] = atan(buf2[i]);
  return tan(nextP0<T, denoise>(buf2.size()).dot(buf2)) - buf2[buf2.size() - 1] + avg;
}

#define _P0_
#endif

