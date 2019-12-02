#if !defined(_P0_)

template <typename T, typename U> class P0 {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleVector<U> VecU;
  typedef SimpleMatrix<T> Mat;
  typedef SimpleMatrix<U> MatU;
  inline P0(const int& range, const int& shrink = 4);
  inline ~P0();
  inline T next(const T& in);
private:
  Vec pred;
  const MatU& seed(const int& size, const bool& idft);
  const Mat&  diff(const int& size, const bool& integrate);
  const Vec&  nextTaylor(const int& size);
  const Vec&  nextDeepTaylor(const int& size);
  const T& Pi() const;
  const U& J()  const;
  Vec buf;
};

template <typename T, typename U> P0<T,U>::P0(const int& range, const int& shrink) {
  assert(1 < range && 0 < shrink);
  buf.resize(range);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  const auto& pred0(nextDeepTaylor(range * shrink));
  pred.resize(range);
  for(int i = 0; i < pred.size(); i ++)
    pred[i] = T(0);
  for(int i = 0; i < pred0.size(); i ++)
    pred[i / shrink] += pred0[i];
  pred /= T(shrink);
}

template <typename T, typename U> P0<T,U>::~P0() {
  ;
}

template <typename T, typename U> inline T P0<T,U>::next(const T& in) {
  for(int i = 1; i < buf.size(); i ++)
    buf[i - 1] = buf[i];
  buf[buf.size() - 1] = in;
  return pred.dot(buf);
}

template <typename T, typename U> const T& P0<T,U>::Pi() const {
  const static auto pi(atan2(T(1), T(1)) * T(4));
  return pi;
}

template <typename T, typename U> const U& P0<T,U>::J() const {
  const static auto i(sqrt(U(T(- 1))));
  return i;
}

template <typename T, typename U> const typename P0<T,U>::MatU& P0<T,U>::seed(const int& size, const bool& f_idft) {
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
  dft[ size].resize(size, size);
  idft[size].resize(size, size);
  auto& edft( dft[size]);
  auto& eidft(idft[size]);
  for(int i = 0; i < edft.rows(); i ++)
    for(int j = 0; j < edft.cols(); j ++) {
      const auto theta(- T(2) * Pi() * T(i) * T(j) / T(edft.rows()));
      edft(i, j)  = U(cos(  theta), sin(  theta));
      eidft(i, j) = U(cos(- theta), sin(- theta)) / U(T(size));
    }
  if(f_idft)
    return eidft;
  return edft;
}

template <typename T, typename U> const typename P0<T,U>::Mat& P0<T,U>::diff(const int& size, const bool& integrate) {
  assert(0 < size);
  static vector<Mat> D;
  static vector<Mat> I;
  if(D.size() <= size)
    D.resize(size + 1, Mat());
  if(I.size() <= size)
    I.resize(size + 1, Mat());
  if((!integrate) && D[size].rows() == size && D[size].cols() == size)
    return D[size];
  if(  integrate  && I[size].rows() == size && I[size].cols() == size)
    return I[size];
  D[size].resize(size, size);
  I[size].resize(size, size);
  auto& d(D[size]);
  auto& i(I[size]);
  auto  DD(seed(size, false));
  DD.row(0) *= U(T(0));
  auto  II(DD);
  T nd(0);
  T ni(0);
  for(int i = 1; i < DD.rows(); i ++) {
    const auto phase(- J() * T(2) * Pi() * T(i) / T(DD.rows()));
    const auto phase2(U(T(1)) / phase);
    DD.row(i) *= phase;
    II.row(i) /= phase;
    nd += abs(phase)  * abs(phase);
    ni += abs(phase2) * abs(phase2);
  }
  d = (seed(size, true) * DD).template real<T>() * sqrt(sqrt(T(DD.rows() - 1) / (nd * ni)));
  i = (seed(size, true) * II).template real<T>() * sqrt(sqrt(T(II.rows() - 1) / (nd * ni)));
  if(integrate)
    return i;
  return d;
}

template <typename T, typename U> const typename P0<T,U>::Vec& P0<T,U>::nextTaylor(const int& size) {
  assert(0 < size);
  static vector<Vec> ntayl;
  if(ntayl.size() <= size)
    ntayl.resize(size + 1, Vec());
  if(ntayl[size].size() == size)
    return ntayl[size];
  ntayl[size].resize(size);
        auto& v(ntayl[size]);
  const auto& D(diff(size, false));
        auto  dd(D.row(D.rows() - 1));
  for(int i = 0; i < v.size() - 1; i ++)
    v[i] = T(0);
  v[v.size() - 1] = T(1);
  for(int i = 2; 0 <= i; i ++) {
    const auto bv(v);
    v += dd;
    if(bv == v)
      break;
    dd = (D * dd) / T(i);
  }
  return v;
}

template <typename T, typename U> const typename P0<T,U>::Vec& P0<T,U>::nextDeepTaylor(const int& size) {
  assert(0 < size);
  static vector<Vec> dtayl;
  if(dtayl.size() <= size)
    dtayl.resize(size + 1, Vec());
  if(dtayl[size].size() == size)
    return dtayl[size];
  auto& p(dtayl[size]);
  p = nextTaylor(size);
  p[p.size() - 1] += T(1);
  for(int i = 2; i < p.size(); i ++) {
    const auto& q(nextTaylor(i));
    for(int j = 0; j < q.size(); j ++)
      p[p.size() + j - q.size()] += q[j];
  }
  return p /= T(p.size());
}

#define _P0_
#endif

