// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random> 

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

template<typename Real>
std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> gramSchmidt(
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> & A 
){
  size_t m = A.rows();
  size_t n = A.cols();
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Q =
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Zero(m, n);
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> R =
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
  for(size_t k=0; k<n; k++) {
    Q.col(k) = A.col(k);
    if(k > 0) {
      for(size_t i=0; i<k-1; i++){
        R(i,k) = Q.col(i).dot(Q.col(k));
        Q.col(k) = Q.col(k) - R.coeff(i,k) * Q.col(i);
      }
    }
    R(k,k) = sqrt(Q.col(k).dot(Q.col(k)));
  }
  return {Q,R};
} 

template <typename Real>
Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> tsolveAndMultiply(
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>& C) {
  Eigen::TriangularView<
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Upper>
    T = A.template triangularView<Eigen::Upper>();
  return T.template solve<Eigen::OnTheRight>(C);  // C %*% solve(A)
}

template <typename Real>
std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> QRdecomp(
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>&
      M) {  // for nrows >= ncols
  Eigen::HouseholderQR<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> qr =
    M.householderQr();
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> R_ =
    qr.matrixQR().template triangularView<Eigen::Upper>();
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Q_ = qr.householderQ();
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> R =
    R_.block(0, 0, M.cols(), M.cols());
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Q =
    Q_.block(0, 0, M.rows(), M.cols());
  return {Q, R};
}

template<typename Real>
Eigen::Matrix<Real, Eigen::Dynamic, 1> solve(
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> & A, 
    const Eigen::Matrix<Real, Eigen::Dynamic, 1> & b
){
  return A.colPivHouseholderQr().solve(b);
} 

template<typename Real>
int rankMatrix(const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> & M){
  Eigen::ColPivHouseholderQR<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> qr(M);
  return qr.rank();
}

std::default_random_engine generator; 
std::normal_distribution<double> gaussian(0.0,1.0);

template<typename Real>
Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> gmatrix(size_t nrows, size_t ncols){
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> G(nrows,ncols);
  for(size_t i=0; i<nrows; i++){
    for(size_t j=0; j<ncols; j++){
      G(i,j) = (Real)(gaussian(generator));
    }
  }
  return G;
}

template<typename Real>
int sgn(Real x){
  return (0 < x) - (x < 0);
}

std::uniform_real_distribution<double> runif(0.0,1.0);

template<typename Real>
Real approx(Real x, unsigned n){
  return round(x * pow(10, n)) / pow(10, n);
}

template<typename Real>
std::vector<Real> fidSample(
    Eigen::Matrix<Real, Eigen::Dynamic, 1> VT2, 
    Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum, 
    Real L, Real U){
  Real ZZ, wt; // outputs
  size_t p = VTsum.size(); // = VT2.size()
  std::vector<size_t> high, low, zero, zeronot;
  size_t lhigh=0, llow=0, lzero=0, lzeronot;
  for(size_t i=0; i<p; i++){
    if(VT2.coeff(i)==0){
      lzero += 1;
      zero.push_back(i);
    }else{
      zeronot.push_back(i);
      if(VT2.coeff(i)>0){
        lhigh += 1;
        high.push_back(i);
      }else{
        llow += 1;
        low.push_back(i);
      }
    }
  }
  lzeronot = p-lzero;
  Real MAX, MIN;
  Real infty = std::numeric_limits<Real>::infinity();
  if((lhigh>0 && llow>0) || lzero>0){
    Real temp;
    std::vector<int> UU(p);
    std::vector<int> LL(p);
    std::vector<int> SS(p);
    for(size_t i=0; i<p; i++){
      UU[i] = sgn<Real>(U-VTsum.coeff(i));
      LL[i] = sgn<Real>(L-VTsum.coeff(i));
      SS[i] = sgn<Real>(VT2.coeff(i));
    }
    if(lzero==p){
      MAX = infty;
      MIN = -infty;
      bool anyUUpos = false, anyLLneg = false;
      size_t i = 0;
      while(i<p && !anyUUpos){
        if(UU[i]>0){
          anyUUpos = true;
        }
        i += 1;
      }
      i = 0;
      while(i<p && !anyLLneg){
        if(LL[i]<0){
          anyLLneg = true;
        }
        i += 1;
      }
      temp = (Real)(anyUUpos && anyLLneg);
    }else{
      size_t i = 0;
      bool c1 = true, c2 = true, d1 = true, d2 = true;
      while(i < lzero && c1){
        if(UU[zero[i]]==-1 || LL[zero[i]]==-1){
          c1 = false;
        }
        i += 1;
      }
      i = 0;
      while(i < lzero && c2){
        if(UU[zero[i]]==1 || LL[zero[i]]==1){
          c2 = false;
        }
        i += 1;
      }
      i = 0;
      while(i < p && d1){
        if(SS[i] == -1){
          d1 = false;
        }
        i += 1;
      }
      i = 0;
      while(i < p && d2){
        if(SS[i] == 1){
          d2 = false;
        }
        i += 1;
      }
      if((d1 && c1) || (d2 && c2)){
        MAX = infty;
        MIN = infty;
        for(size_t i=0; i<lzeronot; i++){
          size_t zni = zeronot[i];
          MIN = std::min(MIN, std::min((U-VTsum.coeff(zni))/VT2.coeff(zni), 
                                       (L-VTsum.coeff(zni))/VT2.coeff(zni)));
        }
        temp = 1-(atan(MIN)/PI+0.5);
      }else if((d2 && c1) || (d1 && c2)){
        MIN = -infty;
        MAX = -infty;
        for(size_t i=0; i<lzeronot; i++){
          size_t zni = zeronot[i];
          MAX = std::max(MAX, 
                         std::max((U-VTsum.coeff(zni))/VT2.coeff(zni), 
                                  (L-VTsum.coeff(zni))/VT2.coeff(zni)));
        }
        temp = atan(MAX)/PI+0.5; 				
      }else{
        Real Hmax = -infty;
        Real Hmin = infty;
        for(size_t i=0; i<lhigh; i++){
          size_t hi = high[i];
          Real xu = (U-VTsum.coeff(hi))/VT2.coeff(hi);
          Real xl = (L-VTsum.coeff(hi))/VT2.coeff(hi);
          Hmax = std::max(Hmax, std::max(xu,xl));
          Hmin = std::min(Hmin, std::min(xu,xl));
        }
        Real Lmax = -infty;
        Real Lmin = infty;
        for(size_t i=0; i<llow; i++){
          size_t li = low[i];
          Real xu = (U-VTsum.coeff(li))/VT2.coeff(li);
          Real xl = (L-VTsum.coeff(li))/VT2.coeff(li);
          Lmax = std::max(Lmax, std::max(xu,xl));
          Lmin = std::min(Lmin, std::min(xu,xl));
        }
        Real bpos, tpos, bneg, tneg;
        if(approx<Real>(Lmin-Hmax,12)>=0){
          bpos = -infty;
          tpos = Hmax;
          bneg = Lmin;
          tneg = infty;
        }else if(approx<Real>(Hmin-Lmax,12)>=0){
          bpos = Hmin;
          tpos = infty;
          bneg = -infty;
          tneg = Lmax;
        }else{
          bpos = -infty;
          tpos = infty;
          bneg = -infty;
          tneg = infty;
        }
        Real Pprob, Nprob;
        if(tpos==infty){
          Pprob = 1-(atan(bpos)/PI+0.5);
        }else{
          Pprob = atan(tpos)/PI+0.5;
        }
        if(tneg==infty){
          Nprob = 1-(atan(bneg)/PI+0.5);
        }else{
          Nprob = atan(tneg)/PI+0.5;
        }
        temp = Pprob+Nprob;
        Pprob = Pprob/temp;
        Nprob = 1-Pprob;
        if(runif(generator) <= Nprob){
          MIN = bneg;
          MAX = tneg;
        }else{
          MIN = bpos;
          MAX = tpos;
        }
      }
    }
    Real y = atan(MAX)/PI+0.5;
    Real x = atan(MIN)/PI+0.5;
    Real u = x+(y-x)*runif(generator);
    ZZ = tan(PI*(u-0.5));
    Real ZZ2 = ZZ*ZZ;
    wt = exp(-ZZ2/2)*(1+ZZ2)*temp;
  }else{
    MAX = -infty;
    MIN = infty;
    for(size_t i=0; i<p; i++){
      Real xu = (U-VTsum.coeff(i))/VT2.coeff(i);
      Real xl = (L-VTsum.coeff(i))/VT2.coeff(i);
      MAX = std::max(MAX, std::max(xu,xl));
      MIN = std::min(MIN, std::min(xu,xl));
    }
    Real y = atan(MAX)/PI+0.5; 
    Real x = atan(MIN)/PI+0.5;
    Real u = x + (y-x)*runif(generator);
    ZZ = tan(PI*(u-0.5));
    Real ZZ2 = ZZ*ZZ;
    wt = exp(-ZZ2/2)*(1+ZZ2)*(y-x);
  }
  
  std::vector<Real> out = {ZZ, wt};
  return out;
}

Eigen::VectorXi cppunique(Eigen::VectorXi v){
  int size = v.size();
  for (int i = 0; i < size; i++) {
    for (int j = i + 1; j < size;) {
      if (v.coeff(j) == v.coeff(i)) {
        for (int k = j; k+1 < size; k++) {
          v(k) = v.coeff(k + 1);
        }
        size--;
      } else {
        j++;
      }
    }
  }
  Eigen::VectorXi out = v.head(size);
  return out;
} 

std::vector<std::vector<int>> cartesianProduct(
    const std::vector<std::vector<int>>& v){
  std::vector<std::vector<int>> s = {{}};
  for(auto& u : v){
    std::vector<std::vector<int>> r;
    for(auto& x : s){
      for(auto y : u){
        r.push_back(x);
        r.back().push_back(y);
      }
    }
    s.swap(r);
  }
  return s;
}

std::vector<std::vector<int>> combinations(std::vector<int> C, int n){
  std::vector<std::vector<int>> sets;
  for(size_t i=0; i<C.size(); i++){
    sets.push_back({C[i],C[i]+n});
  }
  return cartesianProduct(sets);
}

Eigen::MatrixXi vv2matrix(std::vector<std::vector<int>> U, 
                          size_t nrow, size_t ncol){
  Eigen::MatrixXi out(nrow,ncol);
  for(size_t i=0; i<nrow; i++){
    for(size_t j=0; j<ncol; j++){
      out(i,j) = U[j][i];
    }
  }
  return out;
}

size_t spow(size_t base, size_t exp){
  size_t result = 1;
  while(exp){
    if(exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}

Eigen::VectorXi Vsort(Eigen::VectorXi V){ 
  std::sort(V.data(),V.data()+V.size());
  return V;
}

template<typename Real>
std::vector<Real> Vcumsum(const Eigen::Matrix<Real, Eigen::Dynamic, 1> & vec){
  std::vector<Real> out(vec.size());
  out[0] = vec.coeff(0);
  for(int i=1; i<vec.size(); i++) {
    out[i] = out[i-1] + vec.coeff(i);
  }
  return out;
}

template <class T>
std::vector<T> zero2n(T n){
  std::vector<T> out(n);
  for(T i=0; i<n; i++){
    out[i] = i;
  }
  return out;
}

std::vector<size_t> sample_int(size_t n){
  std::vector<size_t> elems = zero2n(n);
  std::shuffle(elems.begin(), elems.end(), generator);
  return elems;
}

template<typename Real>
Real rchisq(int df){
  Real out = 0;
  for(int i=0; i<df; i++){
    Real toadd = (Real)(gaussian(generator));
    out += toadd*toadd;
  }
  return out;
}

Eigen::MatrixXd umatrix(size_t nrows, size_t ncols){
  Eigen::MatrixXd U(nrows,ncols);
  for(size_t i=0; i<nrows; i++){
    for(size_t j=0; j<ncols; j++){
      U(i,j) = runif(generator);
    }
  }
  return U;
}

template<typename Real>
Eigen::MatrixXd pickCoordinates(
    size_t Dim, size_t N, size_t fe, 
    std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> VT, 
    Eigen::MatrixXd U){
  Eigen::MatrixXd VTend(Dim,N);
  for(size_t i=0; i<N; i++){
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VTi = VT[i]; 
    for(size_t j=0; j<Dim; j++){
      if(U.coeff(j,i) < 0.5){
        if(j < fe){
          VTend(j,i) = VTi.row(j).minCoeff();
        }else{
          Real x = VTi.row(j).minCoeff();
          if(x < 0){
            x = 0;
          }
          VTend(j,i) = (double)x;
        }
      }else{
        if(j < fe){
          VTend(j,i) = (double)VTi.row(j).maxCoeff();
        }else{
          Real x = VTi.row(j).maxCoeff();
          if(x < 0){
            x = 0;
          }
          VTend(j,i) = (double)x;
        }
      }
    }
  }
  return VTend;
}

typedef struct GFI {
  Eigen::MatrixXd vertices;
  Eigen::VectorXd weights;
  Rcpp::NumericVector ess;
} GFI;

template<typename Real>
GFI gfilmm_(
    Eigen::Matrix<Real, Eigen::Dynamic, 1> L, 
    Eigen::Matrix<Real, Eigen::Dynamic, 1> U, 
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> FE, 
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> RE, 
    Eigen::MatrixXi RE2, 
    Rcpp::IntegerVector E,
    size_t N, double thresh){
  Eigen::Matrix<Real, Eigen::Dynamic, 1> WTnorm(N); // output:weights
  const size_t n = L.size();
  const size_t fe = FE.cols(); // si FE=NULL, passer une matrice n x 0
  const size_t re = RE2.cols();
  const size_t Dim = fe+re;
  const size_t Dimm1 = Dim-1;
  Eigen::MatrixXd VERTEX(Dim, N); // output
  const Rcpp::IntegerVector Esum = Rcpp::cumsum(E);
  
  //-------- SET-UP ALGORITHM OBJECTS ------------------------------------------
  std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> Z(re);
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> weight = 
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Ones(n,N);
  Rcpp::NumericVector ESS(n, (double)N);
  std::vector<int> C; // initial constraints 
  std::vector<int> K; // complement of C
  std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> VT(N); // vertices
  
  //-------- SAMPLE ALL Z's / SET-UP WEIGHTS -----------------------------------
  std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> A(N);
  for(size_t k=0; k<N; k++){
    A[k].resize(n,re);
  }
  for(size_t j=0; j<re; j++){  
    Z[j] = gmatrix<Real>(E(j),N); 
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> 
      M = RE.block(0, Esum(j)-E(j), n, E(j)) * Z[j];
    for(size_t k=0; k<N; k++){
      for(size_t i=0; i<n; i++){
        A[k](i,j) = M.coeff(i,k);
      }
    } 
  }
  
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> AA(n, Dim);
  AA << FE, A[0]; 
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> AT(0, Dim);
  int r = 0;
  for(size_t i=0; i<n; i++){
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Atemp(AT.rows()+1,Dim);
    Atemp << AT, AA.row(i);
    if(rankMatrix<Real>(Atemp) > r){
      AT = Atemp;
      r += 1;
      C.push_back((int)i);
    }else{
      K.push_back((int)i);
    }
  }
  

  Eigen::VectorXi K_start = 
    Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(C.data(), Dim);
  for(size_t i=0; i<n-Dim; i++){
    for(size_t j=0; j<N; j++){
      Z[re-1](K[i],j) = 0.0;
    }
  }
  
  //-------- FIND INITIAL VERTICES ---------------------------------------------
  std::vector<std::vector<int>> USE = combinations(C,n);
  size_t twoPowerDim = spow(2, Dim);
  Eigen::MatrixXi tUSE = vv2matrix(USE, Dim, twoPowerDim);
  std::vector<Eigen::MatrixXi> CC(N, tUSE); // constraints
  Eigen::Matrix<Real, Eigen::Dynamic, 1> b(2*n);
  b << U, -L;
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> FEFE(2*n,fe);
  FEFE << FE, -FE;

  for(size_t k=0; k<N; k++){ 
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> V(Dim, twoPowerDim);
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> AkAk(2*n,re);
    AkAk << A[k], -A[k];
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> AA(2*n,Dim);
    AA << FEFE, AkAk;
    for(size_t i=0; i<twoPowerDim; i++){
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> AAuse(Dim,Dim);
      Eigen::Matrix<Real, Eigen::Dynamic, 1> buse(Dim);
      for(size_t j=0; j<Dim; j++){
        buse(j) = b.coeff(USE[i][j]);
        for(size_t l=0; l<Dim; l++){
          AAuse(j,l) = AA.coeff(USE[i][j],l);
        }
      }  
      V.col(i) = solve<Real>(AAuse, buse);
    }
    VT[k] = V; 
  }
  std::vector<size_t> VC(N, twoPowerDim); // number of vertices
  
  //-------- MAIN ALGORITHM ----------------------------------------------------
  //double break_point = 10;
  double lengthK = (double)(n-Dim);
  size_t K_n = (size_t)(ceil(lengthK/10.0));
  std::vector<Eigen::VectorXi> K_temp(K_n);
  Eigen::VectorXi KV = Eigen::Map<Eigen::VectorXi, 
                                  Eigen::Unaligned>(K.data(), n-Dim);
  for(size_t i=0; i<K_n-1; i++){
    K_temp[i] = KV.segment(i*10, 10);
  }
  K_temp[K_n-1] = KV.tail(n-Dim-(K_n-1)*10);
  Eigen::VectorXi K1(0);
  for(size_t k_n=0; k_n<K_n; k_n++){
    K1.conservativeResize(K1.size()+K_temp[k_n].size());
    K1.tail(K_temp[k_n].size()) = K_temp[k_n];
    for(int ki=0; ki<K_temp[k_n].size(); ki++){
      int k = K_temp[k_n].coeff(ki);
      if(k_n>0){
        for(size_t i=0; i<re; i++){
          if(E(i)>Z[i].rows()){ 
            int nrowsZi = Z[i].rows();
            Z[i].conservativeResize(E(i), Eigen::NoChange);
            Z[i].bottomRows(E(i)-nrowsZi) = gmatrix<Real>(E(i)-nrowsZi, N);
          }
        }
      }

      for(size_t i=0; i<N; i++){
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VTi = VT[i]; 
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1 = VTi.topRows(Dimm1);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> VT2 = VTi.row(Dimm1);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1t(re-1);
        for(size_t j=0; j<re-1; j++){
          Z1t(j) = Z[j](RE2.coeff(k,j),i);
        }
        Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1(Dimm1); 
        Z1 << FE.row(k).transpose(), Z1t;
        Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum = VT1.transpose() * Z1;
        
        std::vector<Real> sample = 
          fidSample<Real>(VT2, VTsum, L.coeff(k), U.coeff(k));

        Real ZZ = sample[0];
        Real wt = sample[1];
        Z[re-1](k,i) = ZZ;
        weight(k,i) = wt;
        VTsum += ZZ * VT2;
        
        // fidVertex
        Eigen::MatrixXi CCi = CC[i];
        Real Lk = L.coeff(k);
        Real Uk = U.coeff(k);
        size_t p = VTsum.size();
        std::vector<size_t> whichl, checkl, whichu, checku, both;
        for(size_t h=0; h<p; h++){
          bool b = false;
          if(VTsum.coeff(h) >= Lk){
            whichl.push_back(h);
            b = true;
          }else{
            checkl.push_back(h);
          }
          if(VTsum.coeff(h) <= Uk){
            whichu.push_back(h);
            if(b){
              both.push_back(h);
            }
          }else{
            checku.push_back(h);
          }
        }
        Eigen::MatrixXi CCtemp(Dim, 0);
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VTtemp(Dim, 0);
        int vert = 0;
        size_t lcheckl = checkl.size();
        size_t lwhichl = whichl.size();
        if(lcheckl != 0){
          Eigen::MatrixXi CA(Dim, lcheckl); 
          Eigen::MatrixXi CB(Dim, lwhichl);
          for(size_t h=0; h < Dim; h++){
            for(size_t j=0; j < lcheckl; j++){
              CA(h,j) = CCi.coeff(h, checkl[j]);
            }
            for(size_t j=0; j < lwhichl; j++){
              CB(h,j) = CCi.coeff(h, whichl[j]);
            }
          }
          Eigen::MatrixXi INT = Eigen::MatrixXi::Zero(2*n,lcheckl);
          for(size_t ll=0; ll<lcheckl; ll++){
            for(size_t h=0; h<Dim; h++){
              INT(CA.coeff(h,ll),ll) = 1;
            }
          }
          Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum_cl(lcheckl);
          Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum_wl(lwhichl);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1_cl(Dim, lcheckl);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1_wl(Dim, lwhichl);
          for(size_t h=0; h<lcheckl; h++){
            VTsum_cl(h) = VTsum.coeff(checkl[h]);
            for(size_t j=0; j<Dim; j++){
              VT1_cl(j,h) = VTi.coeff(j,checkl[h]);
            }
          }
          for(size_t h=0; h<lwhichl; h++){
            VTsum_wl(h) = VTsum.coeff(whichl[h]);
            for(size_t j=0; j<Dim; j++){
              VT1_wl(j,h) = VTi.coeff(j,whichl[h]);
            }
          }
          for(size_t ii=0; ii<p-lcheckl; ii++){
            Eigen::MatrixXi INT2(Dim, lcheckl);
            for(size_t h=0; h<Dim; h++){
              for(size_t j=0; j<lcheckl; j++){
                INT2(h,j) = INT(CB.coeff(h,ii),j);
              }
            }
            for(size_t j=0; j<lcheckl; j++){
              int colSum = 0;
              for(size_t h=0; h<Dim; h++){
                colSum += INT2.coeff(h,j);
              }
              if(colSum == (int)Dimm1){
                vert += 1;
                std::vector<int> inter(Dim);
                size_t m = 0;
                for(size_t h=0; h<Dim; h++){
                  if(INT2(h,j)==1){
                    inter[m] = CB.coeff(h,ii);
                    m += 1;
                  }
                }
                inter[Dimm1] = k+(int)n;
                CCtemp.conservativeResize(Eigen::NoChange, vert); 
                for(size_t h=0; h<Dim; h++){
                  CCtemp(h,vert-1) = inter[h]; 
                }
                Real lambda = (Lk-VTsum_wl(ii))/(VTsum_cl(j)-VTsum_wl(ii));
                VTtemp.conservativeResize(Eigen::NoChange, vert);
                for(size_t h=0; h<Dim; h++){
                  VTtemp(h,vert-1) = 
                    lambda*VT1_cl.coeff(h,j) + (1-lambda)*VT1_wl.coeff(h,ii);
                }
              }
            }
          }
        }
        size_t lchecku = checku.size();
        size_t lwhichu = whichu.size();
        if(lchecku != 0){
          Eigen::MatrixXi CA(Dim, lchecku); 
          Eigen::MatrixXi CB(Dim, lwhichu);
          for(size_t h=0; h < Dim; h++){
            for(size_t j=0; j < lchecku; j++){
              CA(h,j) = CCi.coeff(h, checku[j]);
            }
            for(size_t j=0; j < lwhichu; j++){
              CB(h,j) = CCi.coeff(h, whichu[j]);
            }
          }
          Eigen::MatrixXi INT = Eigen::MatrixXi::Zero(2*n,lchecku);
          for(size_t ll=0; ll<lchecku; ll++){
            for(size_t h=0; h<Dim; h++){
              INT(CA.coeff(h,ll),ll) = 1;
            }
          }
          Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum_cu(lchecku);
          Eigen::Matrix<Real, Eigen::Dynamic, 1> VTsum_wu(lwhichu);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1_cu(Dim, lchecku);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VT1_wu(Dim, lwhichu);
          for(size_t h=0; h<lchecku; h++){
            VTsum_cu(h) = VTsum.coeff(checku[h]);
            for(size_t j=0; j<Dim; j++){
              VT1_cu(j,h) = VTi.coeff(j,checku[h]);
            }
          }
          for(size_t h=0; h<lwhichu; h++){
            VTsum_wu(h) = VTsum.coeff(whichu[h]);
            for(size_t j=0; j<Dim; j++){
              VT1_wu(j,h) = VTi.coeff(j,whichu[h]);
            }
          }
          for(size_t ii=0; ii<p-lchecku; ii++){
            Eigen::MatrixXi INT2(Dim, lchecku);
            for(size_t h=0; h<Dim; h++){
              for(size_t j=0; j<lchecku; j++){
                INT2(h,j) = INT(CB.coeff(h,ii),j);
              }
            }
            for(size_t j=0; j<lchecku; j++){
              int colSum = 0;
              for(size_t h=0; h<Dim; h++){
                colSum += INT2.coeff(h,j);
              }
              if(colSum == (int)Dimm1){
                vert += 1;
                std::vector<int> inter(Dim);
                size_t m = 0;
                for(size_t h=0; h<Dim; h++){
                  if(INT2.coeff(h,j)==1){
                    inter[m] = CB.coeff(h,ii);
                    m += 1;
                  }
                }
                inter[Dimm1] = k;
                CCtemp.conservativeResize(Eigen::NoChange, vert); 
                for(size_t h=0; h<Dim; h++){
                  CCtemp(h,vert-1) = inter[h]; 
                }
                Real lambda = 
                  (Uk-VTsum_wu.coeff(ii))/(VTsum_cu.coeff(j)-VTsum_wu.coeff(ii));
                VTtemp.conservativeResize(Eigen::NoChange, vert);
                for(size_t h=0; h<Dim; h++){
                  VTtemp(h,vert-1) = 
                    lambda*VT1_cu.coeff(h,j) + (1-lambda)*VT1_wu.coeff(h,ii);
                }
              }
            }
          }
        }
        size_t lboth = both.size();
        if(lboth>0){
          for(size_t j=0; j<lboth; j++){
            vert += 1;
            CCtemp.conservativeResize(Eigen::NoChange, vert); 
            VTtemp.conservativeResize(Eigen::NoChange, vert); 
            for(size_t h=0; h<Dim; h++){
              CCtemp(h,vert-1) = CCi.coeff(h,both[j]);
              VTtemp(h,vert-1) = VTi.coeff(h,both[j]);
            }
          }
        }
        VC[i] = vert;
        CC[i] = CCtemp;
        VT[i] = VTtemp;
      }

      Eigen::Matrix<Real, Eigen::Dynamic, 1> WT(N);
      for(size_t j=0; j<N; j++){
        WT(j) = weight.col(j).prod();
      }
      Real WTsum = WT.sum();
      WTnorm = WT/WTsum;
      ESS(k) = (double)(1/(WTnorm.dot(WTnorm)));

      if(ESS(k)<thresh && k < K[n-Dim-1]){
        std::vector<size_t> N_sons(N,0);
        std::vector<Real> dist = Vcumsum<Real>(WTnorm);
        double aux = runif(generator);
        std::vector<double> u(N);
        for(size_t i=0; i<N; i++){
          u[i] = (aux + (double)i) / (double)N;
        }
        size_t j = 0; 
        for(size_t i=0; i<N; i++){
          while(u[i]>dist[j]){ 
            j = j+1; 
          }
          N_sons[j] = N_sons[j]+1; 
        }
        std::vector<int> zero2k_ = zero2n<int>(k+1);
        Eigen::VectorXi zero2k = 
          Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(zero2k_.data(), k+1);
        Eigen::VectorXi JJ0(k+1+Dim);
        JJ0 << zero2k, K_start;
        Eigen::VectorXi JJ = cppunique(JJ0);
        int lenJJ = JJ.size();
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> REJJ(lenJJ,Esum(re-1)); 
        for(int jj=0; jj<lenJJ; jj++){
          REJJ.row(jj) = RE.row(JJ.coeff(jj));
        }
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> FEJJ(lenJJ,fe); 
        if(fe>0){
          for(int jj=0; jj<lenJJ; jj++){
            FEJJ.row(jj) = FE.row(JJ.coeff(jj));
          }
        }
        std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> ZZ(re);
        for(size_t ii=0; ii<re; ii++){
          ZZ[ii].resize(E(ii), 0);
        }
        std::vector<size_t> VCVC(N,0);
        std::vector<Eigen::MatrixXi> CCCC(N);
        std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> VTVT(N);
        for(size_t i=0; i<N; i++){
          size_t Nsons_i = N_sons[i];
          if(Nsons_i){
            std::vector<size_t> VCtemp(Nsons_i, VC[i]);
            std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> Ztemp(re);
            size_t copy = Nsons_i-1;
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VTi = VT[i];
            std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> 
              VTtemp(Nsons_i, VTi);
            for(size_t ii=0; ii<re; ii++){
              Ztemp[ii].resize(E(ii), Nsons_i);
              for(size_t j=0; j<Nsons_i; j++){
                Ztemp[ii].col(j) = Z[ii].col(i);
              }
            }			
            if(copy){
              std::vector<size_t> ord = sample_int(re);
              for(size_t j=0; j<re; j++){
                size_t kk = ord[j];
                for(size_t ii=0; ii<copy; ii++){
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> XX(lenJJ, 0);
                  for(size_t jj=0; jj<re; jj++){
                    if(jj != kk){
                      XX.conservativeResize(Eigen::NoChange, XX.cols()+1);
                      Eigen::Matrix<Real, Eigen::Dynamic, 1> newcol(lenJJ);
                      newcol = 
                        REJJ.block(0, Esum(jj)-E(jj), lenJJ, E(jj)) *
                        Ztemp[jj].col(ii);
                      XX.rightCols(1) = newcol; 
                    }
                  }
                  Eigen::VectorXi v(lenJJ);
                  for(int jj=0; jj<lenJJ; jj++){
                    v(jj) = RE2(JJ.coeff(jj),kk);
                  }
                  Eigen::VectorXi vv = cppunique(v);
                  Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1_(vv.size());
                  for(int jj=0; jj<vv.size(); jj++){
                    Z1_(jj) = Ztemp[kk].coeff(vv(jj),ii);
                  }
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> CO2_ = 
                    REJJ.block(0, Esum(kk)-E(kk), lenJJ, E(kk));
                  std::vector<int> pcolsums;
                  for(int jj=0; jj<E(kk); jj++){
                    Real colsum = CO2_.col(jj).sum();
                    if(colsum > 0){
                      pcolsums.push_back(jj);
                    }
                  }
                  size_t ncolsCO2 = pcolsums.size();
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> CO2(lenJJ, ncolsCO2);
                  for(size_t jj=0; jj<ncolsCO2; jj++){
                    CO2.col(jj) = CO2_.col(pcolsums[jj]);
                  }
                  std::vector<int> Z00;
                  for(int jj=0; jj<Z1_.size(); jj++){
                    if(Z1_.coeff(jj) != 0){
                      Z00.push_back(jj);
                    }
                  }
                  size_t lenZ1 = Z00.size();
                  Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1(lenZ1);
                  for(size_t jj=0; jj<lenZ1; jj++){
                    Z1(jj) = Z1_.coeff(Z00[jj]);
                  }
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> XXX(lenJJ, Dimm1);
                  if(fe>0){
                    XXX << FEJJ, XX;
                  }else{
                    XXX = XX;
                  }
                  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> MAT(lenJJ,((int)Dim)-1+ncolsCO2);
                  MAT << -XXX, CO2;
                  Eigen::FullPivLU<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> lu(MAT);
                  int rk = lu.rank();
                  if(rk < MAT.cols()){
                    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> NUL = lu.kernel();
                    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> n1 = NUL.topRows(NUL.rows()-ncolsCO2);
                    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> n2 = NUL.bottomRows(ncolsCO2);
                    std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> QR = QRdecomp<Real>(n2);
                //    std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> QR = gramSchmidt<Real>(n2);
                    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> O2 = QR[0];
                    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> R = QR[1];
                    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> O1 = tsolveAndMultiply<Real>(R, n1);
                    //Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> O1 = 
                    //  n1 * R.colPivHouseholderQr().inverse();
                    int rankO2 = O2.cols();
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> a(rankO2);
                    a = O2.transpose() * Z1;
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> O2a(lenZ1); 
                    O2a = O2 * a;
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> tau_(lenZ1);
                    tau_ = Z1 - O2a;
                    Real b = sqrt(tau_.dot(tau_));
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> tau = tau_/b;
                    Real bb = sqrt(rchisq<Real>(lenZ1-rankO2));
                    Real bbb = b/bb;
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> aa(rankO2);
                    for(int jj=0; jj<rankO2; jj++){
                      aa(jj) = (Real)(gaussian(generator));
                    }
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> O2aa(lenZ1);
                    O2aa = O2*aa;
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> bbtau(lenZ1);
                    bbtau = bb*tau;
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> MM3(lenZ1);
                    MM3 = O2aa + bbtau;
                    for(size_t jj=0; jj<lenZ1; jj++){
                      Ztemp[kk](Z00[jj],ii) = MM3(jj); 
                    }
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> M2(Dimm1);
                    M2 = O1 * (bbb*aa - a);
                    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VTtemp_ii = VTtemp[ii];
                    for(size_t jj=0; jj<VC[i]; jj++){
                      for(size_t vert=0; vert<Dim; vert++){
                        if(vert < fe+kk){
                          VTtemp[ii](vert,jj) = VTtemp_ii.coeff(vert,jj) - 
                            VTtemp_ii.coeff(fe+kk,jj) * M2.coeff(vert);
                        }else if(vert > fe+kk){
                          VTtemp[ii](vert,jj) = VTtemp_ii.coeff(vert,jj) - 
                            VTtemp_ii.coeff(fe+kk,jj) * M2.coeff(vert-1);
                        }
                      }
                      VTtemp[ii](fe+kk,jj) = bbb*VTtemp_ii.coeff(fe+kk,jj);
                    }
                  }else{
                    Real b = sqrt(Z1.dot(Z1));
                    Eigen::Matrix<Real, Eigen::Dynamic, 1> tau = Z1/b;
                    Real bb = sqrt(rchisq<Real>(lenZ1));
                    for(size_t jj=0; jj<lenZ1; jj++){
                      Ztemp[kk](Z00[jj],ii) = bb*tau.coeff(jj);
                    }
                    VTtemp[ii].row(fe+kk) *= b/bb;
                  }
                }
              }
            }
            for(size_t ii=0; ii<re; ii++){
              ZZ[ii].conservativeResize(Eigen::NoChange, 
                                        ZZ[ii].cols()+(int)Nsons_i);
              ZZ[ii].rightCols(Nsons_i) = Ztemp[ii];
            }
            size_t d = 0;
            for(size_t ii=0; ii<i; ii++){
              d += N_sons[ii];
            }
            size_t dd = d + Nsons_i;
            for(size_t ii=d; ii<dd; ii++){
              VCVC[ii] = VCtemp[ii-d];
            }
            for(size_t kk=0; kk<Nsons_i; kk++){
              VTVT[kk+d] = VTtemp[kk];
              CCCC[kk+d] = CC[i];
            }
          }
        }
        Z = ZZ;
        VT = VTVT;
        VC = VCVC;
        CC = CCCC;
        weight = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Ones(n,N);
      }
    }

    //------------determine signs ----------------------------------------------
    Eigen::MatrixXi signs = Eigen::MatrixXi::Zero(re,N);
    for(size_t i=0; i<N; i++){
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VTi = VT[i];
      for(size_t j=0; j<re; j++){
        Eigen::Matrix<Real, Eigen::Dynamic, 1> row = VTi.row(fe+j);
        bool allneg = true;
        size_t jj = 0;
        while(jj<VC[i] && allneg){
          allneg = row.coeff(jj) < 0;
          jj += 1;
        }
        if(allneg){
          signs(j,i) = -1;
        }
      }
    }

    //--------FINAL RESAMPLE ---------------------------------------------------
    std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> ZZ(re); // resized below
    std::vector<Eigen::VectorXi> nn(re);
    std::vector<int> lengths_nn(re);
    std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> VTVT(N);
    Eigen::VectorXi n1_(K1.size()+(int)Dim);
    n1_ << K1,K_start;
    Eigen::VectorXi n1 = Vsort(n1_);
    for(size_t ii=0; ii<re; ii++){
      Eigen::VectorXi vec(n1.size());
      for(int jj=0; jj<n1.size(); jj++){
        vec(jj) = RE2.coeff(n1(jj),ii);
      }
      nn[ii] = cppunique(vec);
      lengths_nn[ii] = nn[ii].size();
      ZZ[ii].resize(lengths_nn[ii], 0);
    }

    for(size_t i=0; i<N; i++){
      Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> VTtemp = VT[i];
      std::vector<Eigen::Matrix<Real, Eigen::Dynamic, 1>> Ztemp(re);
      for(size_t ii=0; ii<re; ii++){
        Ztemp[ii].resize(lengths_nn[ii]);
        for(int iii=0; iii<lengths_nn[ii]; iii++){
          Ztemp[ii](iii) = Z[ii].coeff(nn[ii].coeff(iii),i);
        }
      }
      std::vector<size_t> ord = sample_int(re);
      for(size_t j=0; j<re; j++){
        size_t kk = ord[j];
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> XX(n,0);
        for(size_t jj=0; jj<re; jj++){
          if(jj != kk){
            XX.conservativeResize(Eigen::NoChange, XX.cols()+1);
            Eigen::Matrix<Real, Eigen::Dynamic, 1> newcol(n);
            newcol = RE.block(0, Esum(jj)-E(jj), n, 
                              lengths_nn[jj]) * Ztemp[jj];
            XX.rightCols(1) = newcol;
          }
        }
        Eigen::Matrix<Real, Eigen::Dynamic, 1> Z1 = Ztemp[kk];
        int lenZ1 = Z1.size();
        int ncolsCO2 = lengths_nn[kk];
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> CO2 = RE.block(0, Esum(kk)-E(kk), n, ncolsCO2);
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> XXX(n, Dimm1);
        if(fe>0){
          XXX << FE, XX;
        }else{
          XXX = XX;
        }
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> MAT(n, Dimm1+ncolsCO2);
        MAT << -XXX, CO2;
        Eigen::FullPivLU<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> lu(MAT);
        int rk = lu.rank();
        if(rk < MAT.cols()){
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> NUL = lu.kernel();
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> n1 = NUL.topRows(NUL.rows()-ncolsCO2);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> n2 = NUL.bottomRows(ncolsCO2);
          std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> QR = QRdecomp<Real>(n2);
          //std::vector<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> QR = gramSchmidt<Real>(n2);
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> O2 = QR[0];
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> R = QR[1];
          Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> O1 = tsolveAndMultiply<Real>(R, n1);
          //Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> O1 = 
          //  n1 * R.colPivHouseholderQr().inverse();
          int rankO2 = O2.cols();
          Eigen::Matrix<Real, Eigen::Dynamic, 1> a(rankO2);
          a = O2.transpose() * Z1;
          Eigen::Matrix<Real, Eigen::Dynamic, 1> O2a(lenZ1); 
          O2a = O2 * a;
          Eigen::Matrix<Real, Eigen::Dynamic, 1> tau_(lenZ1);
          tau_ = Z1 - O2a;
          Real b = sqrt(tau_.dot(tau_));
          Eigen::Matrix<Real, Eigen::Dynamic, 1> tau = tau_/b;
          Real bb = sqrt(rchisq<Real>(lenZ1-rankO2));
          Real bbb = b/bb;
          Eigen::Matrix<Real, Eigen::Dynamic, 1> aa(rankO2);
          for(int jj=0; jj<rankO2; jj++){
            aa(jj) = (Real)(gaussian(generator));
          }
          Eigen::Matrix<Real, Eigen::Dynamic, 1> O2aa(lenZ1);
          O2aa = O2*aa;
          Eigen::Matrix<Real, Eigen::Dynamic, 1> bbtau(lenZ1);
          bbtau = bb*tau;
          Ztemp[kk] = O2aa + bbtau;
          Eigen::Matrix<Real, Eigen::Dynamic, 1> M2(Dimm1);
          M2 = O1 * (bbb*aa - a);
          for(size_t jj=0; jj<VC[i]; jj++){
            for(size_t vert=0; vert<Dim; vert++){
              if(vert < fe+kk){
                VTtemp(vert,jj) -= VTtemp.coeff(fe+kk,jj) * M2.coeff(vert);
              }else if(vert > fe+kk){
                VTtemp(vert,jj) -= VTtemp.coeff(fe+kk,jj) * M2.coeff(vert-1);
              }
            }
            VTtemp(fe+kk,jj) *= bbb;
          }
        }else{
          Real b = sqrt(Z1.dot(Z1));
          Eigen::Matrix<Real, Eigen::Dynamic, 1> tau = Z1/b;
          Real bb = sqrt(rchisq<Real>(lenZ1));
          Ztemp[kk] = bb*tau;
          VTtemp.row(fe+kk) *= b/bb;
        }
      }
      for(size_t ii=0; ii<re; ii++){
        ZZ[ii].conservativeResize(Eigen::NoChange, ZZ[ii].cols()+1);
        ZZ[ii].rightCols(1) = Ztemp[ii];
      }
      VTVT[i] = VTtemp;
    }

    Z = ZZ; 
    VT = VTVT;
    
    //---------flip negatives --------------------------------------------------			
    for(size_t i=0; i<N; i++){
      for(size_t j=0; j<re; j++){
        if(signs(j,i) == -1){ 
          VT[i].row(fe+j) *= -1;
          Z[j].col(i) *= -1;
        }
      }
    }

    if(k_n == K_n-1){ //if finished pick coordinates
      Eigen::MatrixXd unif = umatrix(Dim, N); 
      VERTEX = pickCoordinates<Real>(Dim, N, fe, VT, unif);
    }
    
  }

  GFI out;
  Eigen::VectorXd ow = Eigen::VectorXd::Zero(N);
  out.weights = ow;
  out.vertices = VERTEX; //.resize(Dim, N);
  out.ess = ESS;
  
  for(auto j=0; j<N; j++){
    out.weights(j) = (double)WTnorm.coeff(j);
  //  for(auto i=0; i<Dim; i++){
  //    out.vertices(i,j) = (double)VERTEX.coeff(i,j);
  //  }
  }
  
  return out;
}

Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> d2l(
    Eigen::MatrixXd M
){
  Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> out(M.rows(), M.cols());
  for(auto i=0; i<M.rows(); i++){
    for(auto j=0; j<M.cols(); j++){
      out(i,j) = (long double)M.coeff(i,j);
    }
  }
  return out;
}

Eigen::Matrix<long double, Eigen::Dynamic, 1> d2lVector(
    Eigen::VectorXd V
){
  Eigen::Matrix<long double, Eigen::Dynamic, 1> out(V.rows());
  for(auto i=0; i<V.size(); i++){
    out(i) = (long double)V.coeff(i);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List gfilmm_double(
  Eigen::VectorXd L, Eigen::VectorXd U, 
  Eigen::MatrixXd FE, Eigen::MatrixXd RE, 
  Eigen::MatrixXi RE2, 
  Rcpp::IntegerVector E,
  size_t N, double thresh)
{
  
  GFI gfi = gfilmm_<double>(L, U, FE, RE, RE2, E, N, thresh);
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("VERTEX") = gfi.vertices,
    Rcpp::Named("WEIGHT") = gfi.weights
  );
  out.attr("ESS") = gfi.ess;
  return out;
  
}

// [[Rcpp::export]]
Rcpp::List gfilmm_long(
    Eigen::VectorXd L, Eigen::VectorXd U, 
    Eigen::MatrixXd FE, Eigen::MatrixXd RE, 
    Eigen::MatrixXi RE2, 
    Rcpp::IntegerVector E,
    size_t N, double thresh)
{
  
  GFI gfi = gfilmm_<long double>(d2lVector(L), d2lVector(U), d2l(FE), d2l(RE), RE2, E, N, thresh);
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("VERTEX") = gfi.vertices,
    Rcpp::Named("WEIGHT") = gfi.weights
  );
  out.attr("ESS") = gfi.ess;
  return out;
  
}
