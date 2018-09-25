// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

Eigen::VectorXd solve(const Eigen::MatrixXd & A, const Eigen::VectorXd & b){
  return A.colPivHouseholderQr().solve(b);
} 

Eigen::MatrixXd tsolveAndMultiply(const Eigen::MatrixXd & A, const Eigen::MatrixXd & C){
  Eigen::MatrixXd M = A.triangularView<Eigen::Upper>().solve<Eigen::OnTheRight>(C);
  return M; // C %*% solve(A)
} 

Rcpp::List nullSpace(const Eigen::MatrixXd & M){ 
  Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
  Eigen::MatrixXd nspace = lu.kernel(); // not orthonormal
  int r = lu.rank();
  return Rcpp::List::create(Rcpp::Named("kernel") = nspace,
                            Rcpp::Named("rank") = r);
}

Rcpp::List QRdecomp(const Eigen::MatrixXd & M){ // for nrows >= ncols
  Eigen::HouseholderQR<Eigen::MatrixXd> qr = M.householderQr();
  Eigen::MatrixXd R_ = qr.matrixQR().triangularView<Eigen::Upper>();
  Eigen::MatrixXd Q_ = qr.householderQ();
  Eigen::MatrixXd R = R_.block(0,0,M.cols(),M.cols());
  Eigen::MatrixXd Q = Q_.block(0,0,M.rows(),M.cols());
  return Rcpp::List::create(Rcpp::Named("Q") = Q,
                            Rcpp::Named("R") = R);
}

int rankMatrix(const Eigen::MatrixXd & M){
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(M);
  return qr.rank();
}

std::default_random_engine generator; 
std::normal_distribution<double> gaussian(0.0,1.0);

Eigen::MatrixXd gmatrix(size_t nrows, size_t ncols){
  Eigen::MatrixXd G(nrows,ncols);
  for(size_t i=0; i<nrows; i++){
    for(size_t j=0; j<ncols; j++){
      G(i,j) = gaussian(generator);
    }
  }
  return G;
}



int sgn(double x){
  return (0 < x) - (x < 0);
}


std::uniform_real_distribution<double> runif(0.0,1.0);

double approx(double x, unsigned n){
  return round(x * pow(10, n)) / pow(10, n);
}

// [[Rcpp::export]]
std::vector<double> fidSample(Eigen::VectorXd VT2, Eigen::VectorXd VTsum, 
                     double L, double U){
  double ZZ, wt; // outputs
  size_t p = VTsum.size(); // = VT2.size()
  std::vector<size_t> high, low, zero, zeronot;
  size_t lhigh=0, llow=0, lzero=0, lzeronot;
  for(size_t i=0; i<p; i++){
    if(VT2(i)==0){
      lzero += 1;
      zero.push_back(i);
    }else{
      zeronot.push_back(i);
      if(VT2(i)>0){
        lhigh += 1;
        high.push_back(i);
      }else{
        llow += 1;
        low.push_back(i);
      }
    }
  }
  lzeronot = p-lzero;
  double MAX, MIN;
  double infty = std::numeric_limits<double>::infinity();
  if((lhigh>0 && llow>0) || lzero>0){
    double temp;
    std::vector<int> UU(p);
    std::vector<int> LL(p);
    std::vector<int> SS(p);
    for(size_t i=0; i<p; i++){
      UU[i] = sgn(U-VTsum(i));
      LL[i] = sgn(L-VTsum(i));
      SS[i] = sgn(VT2(i));
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
      temp = (double)(anyUUpos && anyLLneg);
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
          MIN = std::min(MIN, std::min((U-VTsum(zni))/VT2(zni),(L-VTsum(zni))/VT2(zni)));
        }
        temp = 1-(atan(MIN)/PI+0.5);
      }else if((d2 && c1) || (d1 && c2)){
        MIN = -infty;
        MAX = -infty;
        for(size_t i=0; i<lzeronot; i++){
          size_t zni = zeronot[i];
          MAX = std::max(MAX, std::max((U-VTsum(zni))/VT2(zni),(L-VTsum(zni))/VT2(zni)));
        }
        temp = atan(MAX)/PI+0.5; 				
      }else{
        double Hmax = -infty;
        double Hmin = infty;
        for(size_t i=0; i<lhigh; i++){
          size_t hi = high[i];
          double xu = (U-VTsum(hi))/VT2(hi);
          double xl = (L-VTsum(hi))/VT2(hi);
          Hmax = std::max(Hmax, std::max(xu,xl));
          Hmin = std::min(Hmin, std::min(xu,xl));
        }
        double Lmax = -infty;
        double Lmin = infty;
        for(size_t i=0; i<llow; i++){
          size_t li = low[i];
          double xu = (U-VTsum(li))/VT2(li);
          double xl = (L-VTsum(li))/VT2(li);
          Lmax = std::max(Lmax, std::max(xu,xl));
          Lmin = std::min(Lmin, std::min(xu,xl));
        }
        double bpos, tpos, bneg, tneg;
        if(approx(Lmin-Hmax,12)>=0){
          bpos = -infty;
          tpos = Hmax;
          bneg = Lmin;
          tneg = infty;
        }else if(approx(Hmin-Lmax,12)>=0){
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
        double Pprob, Nprob;
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
    double y = atan(MAX)/PI+0.5;
    double x = atan(MIN)/PI+0.5;
    double u = x+(y-x)*runif(generator);
    ZZ = tan(PI*(u-0.5));
    double ZZ2 = ZZ*ZZ;
    wt = exp(-ZZ2/2)*(1+ZZ2)*temp;
  }else{
    MAX = -infty;
    MIN = infty;
    for(size_t i=0; i<p; i++){
      double xu = (U-VTsum(i))/VT2(i);
      double xl = (L-VTsum(i))/VT2(i);
      MAX = std::max(MAX, std::max(xu,xl));
      MIN = std::min(MIN, std::min(xu,xl));
    }
    double y = atan(MAX)/PI+0.5; 
    double x = atan(MIN)/PI+0.5;
    double u = x + (y-x)*runif(generator);
    ZZ = tan(PI*(u-0.5));
    double ZZ2 = ZZ*ZZ;
    wt = exp(-ZZ2/2)*(1+ZZ2)*(y-x);
  }
  
  std::vector<double> out = {ZZ, wt};
  return out;
}

Eigen::VectorXi cppunique(Eigen::VectorXi v){
  int size = v.size();
  for (int i = 0; i < size; i++) {
    for (int j = i + 1; j < size;) {
      if (v(j) == v(i)) {
        for (int k = j; k+1 < size; k++) {
          v(k) = v(k + 1);
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

std::vector<std::vector<int>> cartesianProduct(const std::vector<std::vector<int>>& v){
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

Eigen::MatrixXi vv2matrix(std::vector<std::vector<int>> U, size_t nrow, size_t ncol){
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

Eigen::VectorXi Vsort(Eigen::VectorXi V){ // dans Eigen ?
  std::sort(V.data(),V.data()+V.size());
  return V;
}

std::vector<double> Vcumsum(const Eigen::VectorXd & vec){
  std::vector<double> out(vec.size());
  out[0] = vec(0);
  for(int i=1; i<vec.size(); i++) {
    out[i] = out[i-1] + vec(i);
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

double rchisq(int df){
  double out = 0;
  for(int i=0; i<df; i++){
    double toadd = gaussian(generator);
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

Eigen::MatrixXd pickCoordinates(size_t Dim, size_t N, size_t fe, std::vector<Eigen::MatrixXd> VT, Eigen::MatrixXd U){
  Eigen::MatrixXd VTend(Dim,N);
  for(size_t i=0; i<N; i++){
    Eigen::MatrixXd VTi = VT[i]; 
    for(size_t j=0; j<Dim; j++){
      if(U(j,i) < 0.5){
        if(j < fe){
          VTend(j,i) = VTi.row(j).minCoeff();
        }else{
          double x = VTi.row(j).minCoeff();
          if(x < 0){
            x = 0;
          }
          VTend(j,i) = x;
        }
      }else{
        if(j < fe){
          VTend(j,i) = VTi.row(j).maxCoeff();
        }else{
          double x = VTi.row(j).maxCoeff();
          if(x < 0){
            x = 0;
          }
          VTend(j,i) = x;
        }
      }
    }
  }
  return VTend;
}

// [[Rcpp::export]]
Rcpp::List gfilmm_(Eigen::VectorXd L, Eigen::VectorXd U, 
                  Eigen::MatrixXd FE, Eigen::MatrixXd RE, 
                  Eigen::MatrixXi RE2, Rcpp::IntegerVector E,
                  size_t N, size_t thresh){
  Eigen::VectorXd WTnorm(N); // output:weights
  const size_t n = L.size();
  const size_t fe = FE.cols(); // si FE=NULL, passer une matrice n x 0
  const size_t re = RE2.cols();
  const size_t Dim = fe+re;
  const size_t Dimm1 = Dim-1;
  Eigen::MatrixXd VERTEX(Dim, N); // output
  const Rcpp::IntegerVector Esum = Rcpp::cumsum(E);
  
  //-------- SET-UP ALGORITHM OBJECTS ------------------------------------------
  std::vector<Eigen::MatrixXd> Z(re);
  Eigen::MatrixXd weight = Eigen::MatrixXd::Ones(E(re-1),N);
  Rcpp::NumericVector ESS(n, (double)N);
  std::vector<int> C; // initial constraints 
  std::vector<int> K; // complement of C
  std::vector<Eigen::MatrixXd> VT(N); // vertices

  //-------- SAMPLE ALL Z's / SET-UP WEIGHTS -----------------------------------
  std::vector<Eigen::MatrixXd> A(N);
  for(size_t k=0; k<N; k++){
    A[k].resize(n,re);
  }
  for(size_t j=0; j<re; j++){  
    Z[j] = gmatrix(E(j),N); 
    Eigen::MatrixXd M = RE.block(0, Esum(j)-E(j), n, E(j)) * Z[j];
    for(size_t k=0; k<N; k++){
      for(size_t i=0; i<n; i++){
        A[k](i,j) = M(i,k);
      }
    } 
  }
  
  Eigen::MatrixXd AA(n, Dim);
  AA << FE,A[1];
  Eigen::MatrixXd AT(0, Dim);
  int r = 0;
  for(size_t i=0; i<n; i++){
    Eigen::MatrixXd Atemp(AT.rows()+1,Dim);
    Atemp << AT, AA.row(i);
    if(rankMatrix(Atemp) > r){
      AT = Atemp;
      r = rankMatrix(AT);
      C.push_back((int)i);
    }else{
      K.push_back((int)i);
    }
  }

  Eigen::VectorXi K_start = 
    Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(C.data(), Dim);
  for(size_t i=0; i<n-Dim; i++){
    for(size_t j=0; j<N; j++){
      Z[re-1](K[i],j) = 0;
    }
  }
  
  //-------- FIND INITIAL VERTICES ---------------------------------------------
  std::vector<std::vector<int>> USE = combinations(C,n);
  size_t twoPowerDim = spow(2, Dim);
  Eigen::MatrixXi tUSE = vv2matrix(USE, Dim, twoPowerDim);
  std::vector<Eigen::MatrixXi> CC(N, tUSE); // constraints
  Eigen::VectorXd b(2*n);
  b << U, -L;
  Eigen::MatrixXd FEFE(2*n,fe);
  FEFE << FE,-FE;
  for(size_t k=0; k<N; k++){
    Eigen::MatrixXd V(Dim, twoPowerDim);
    Eigen::MatrixXd AkAk(2*n,re);
    AkAk << A[k], -A[k];
    Eigen::MatrixXd AA(2*n,Dim);
    AA << FEFE,AkAk;
    for(size_t i=0; i<twoPowerDim; i++){
      Eigen::MatrixXd AAuse(Dim,Dim);
      Eigen::VectorXd buse(Dim);
      for(size_t j=0; j<Dim; j++){
        buse(j) = b(USE[i][j]);
        for(size_t l=0; l<Dim; l++){
          AAuse(j,l) = AA(USE[i][j],l);
        }
      }  
      V.col(i) = solve(AAuse, buse);
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
      int k = K_temp[k_n](ki);
      if(k_n>0){
        for(size_t i=0; i<re; i++){
          if(E[i]>Z[i].rows()){ 
            int nrowsZi = Z[i].rows();
            Z[i].conservativeResize(E(i), Eigen::NoChange);
            Z[i].bottomRows(E(i)-nrowsZi) = gmatrix(E(i)-nrowsZi, N);
          }
        }
      }
      for(size_t i=0; i<N; i++){
        Eigen::MatrixXd VTi = VT[i]; 
        Eigen::MatrixXd VT1 = VTi.topRows(Dimm1);
        Eigen::VectorXd VT2 = VTi.row(Dimm1);
        Eigen::VectorXd Z1t(re-1);
        for(size_t j=0; j<re-1; j++){
          Z1t(j) = Z[j](RE2(k,j),i);
        }
        Eigen::VectorXd Z1(Dimm1);
        Z1 << FE.row(k), Z1t;
        Eigen::VectorXd VTsum = VT1.transpose() * Z1;
        
        std::vector<double> sample = fidSample(VT2, VTsum, L(k), U(k));
        double ZZ = sample[0];
        double wt = sample[1];
        Z[re-1](k,i) = ZZ;
        weight(k,i) = wt;
        VTsum += ZZ * VT2;

        // Rcpp::List vertex = 
        //   fidVertex(VTi, CC[i], VTsum, L(k), U(k), Dim, (int)n, k);
        // VC[i] = vertex["vert"];
        // CC[i] = vertex["CCtemp"];
        // VT[i] = vertex["VTtemp"];
        
        // fidVertex
        Eigen::MatrixXi CCi = CC[i];
        double Lk = L(k);
        double Uk = U(k);
        size_t p = VTsum.size();
        std::vector<size_t> whichl, checkl, whichu, checku, both;
        for(size_t h=0; h<p; h++){
          bool b = false;
          if(VTsum(h) >= Lk){
            whichl.push_back(h);
            b = true;
          }else{
            checkl.push_back(h);
          }
          if(VTsum(h) <= Uk){
            whichu.push_back(h);
            if(b){
              both.push_back(h);
            }
          }else{
            checku.push_back(h);
          }
        }
        Eigen::MatrixXi CCtemp(Dim, 0);
        Eigen::MatrixXd VTtemp(Dim, 0);
        int vert = 0;
        size_t lcheckl = checkl.size();
        size_t lwhichl = whichl.size();
        if(lcheckl < p){
          Eigen::MatrixXi CA(Dim, lcheckl); 
          Eigen::MatrixXi CB(Dim, lwhichl);
          for(size_t h=0; h < Dim; h++){
            for(size_t j=0; j < lcheckl; j++){
              CA(h,j) = CCi(h, checkl[j]);
            }
            for(size_t j=0; j < lwhichl; j++){
              CB(h,j) = CCi(h, whichl[j]);
            }
          }
          Eigen::MatrixXi INT = Eigen::MatrixXi::Zero(2*n,lcheckl);
          for(size_t ll=0; ll<lcheckl; ll++){
            for(size_t h=0; h<Dim; h++){
              INT(CA(h,ll),ll) = 1;
            }
          }
          Eigen::VectorXd VTsum_cl(lcheckl);
          Eigen::VectorXd VTsum_wl(lwhichl);
          Eigen::MatrixXd VT1_cl(Dim, lcheckl);
          Eigen::MatrixXd VT1_wl(Dim, lwhichl);
          for(size_t h=0; h<lcheckl; h++){
            VTsum_cl(h) = VTsum(checkl[h]);
            for(size_t j=0; j<Dim; j++){
              VT1_cl(j,h) = VTi(j,checkl[h]);
            }
          }
          for(size_t h=0; h<lwhichl; h++){
            VTsum_wl(h) = VTsum(whichl[h]);
            for(size_t j=0; j<Dim; j++){
              VT1_wl(j,h) = VTi(j,whichl[h]);
            }
          }
          for(size_t ii=0; ii<p-lcheckl; ii++){
            Eigen::MatrixXi INT2(Dim, lcheckl);
            for(size_t h=0; h<Dim; h++){
              for(size_t j=0; j<lcheckl; j++){
                INT2(h,j) = INT(CB(h,ii),j);
              }
            }
            for(size_t j=0; j<lcheckl; j++){
              int colSum = 0;
              for(size_t h=0; h<Dim; h++){
                colSum += INT2(h,j);
              }
              if(colSum == (int)Dimm1){
                vert += 1;
                std::vector<int> inter(Dim);
                size_t m = 0;
                for(size_t h=0; h<Dim; h++){
                  if(INT2(h,j)==1){
                    inter[m] = CB(h,ii);
                    m += 1;
                  }
                }
                inter[Dimm1] = k+(int)n;
                CCtemp.conservativeResize(Eigen::NoChange, vert); 
                for(size_t h=0; h<Dim; h++){
                  CCtemp(h,vert-1) = inter[h]; 
                }
                double lambda = (Lk-VTsum_wl(ii))/(VTsum_cl(j)-VTsum_wl(ii));
                VTtemp.conservativeResize(Eigen::NoChange, vert);
                for(size_t h=0; h<Dim; h++){
                  VTtemp(h,vert-1) = 
                    lambda*VT1_cl(h,j) + (1-lambda)*VT1_wl(h,ii);
                }
              }
            }
          }
        }
        size_t lchecku = checku.size();
        size_t lwhichu = whichu.size();
        if(lchecku < p){
          Eigen::MatrixXi CA(Dim, lchecku); 
          Eigen::MatrixXi CB(Dim, lwhichu);
          for(size_t h=0; h < Dim; h++){
            for(size_t j=0; j < lchecku; j++){
              CA(h,j) = CCi(h, checku[j]);
            }
            for(size_t j=0; j < lwhichu; j++){
              CB(h,j) = CCi(h, whichu[j]);
            }
          }
          Eigen::MatrixXi INT = Eigen::MatrixXi::Zero(2*n,lchecku);
          for(size_t ll=0; ll<lchecku; ll++){
            for(size_t h=0; h<Dim; h++){
              INT(CA(h,ll),ll) = 1;
            }
          }
          Eigen::VectorXd VTsum_cu(lchecku);
          Eigen::VectorXd VTsum_wu(lwhichu);
          Eigen::MatrixXd VT1_cu(Dim, lchecku);
          Eigen::MatrixXd VT1_wu(Dim, lwhichu);
          for(size_t h=0; h<lchecku; h++){
            VTsum_cu(h) = VTsum(checku[h]);
            for(size_t j=0; j<Dim; j++){
              VT1_cu(j,h) = VTi(j,checku[h]);
            }
          }
          for(size_t h=0; h<lwhichu; h++){
            VTsum_wu(h) = VTsum(whichu[h]);
            for(size_t j=0; j<Dim; j++){
              VT1_wu(j,h) = VTi(j,whichu[h]);
            }
          }
          for(size_t ii=0; ii<p-lchecku; ii++){
            Eigen::MatrixXi INT2(Dim, lchecku);
            for(size_t h=0; h<Dim; h++){
              for(size_t j=0; j<lchecku; j++){
                INT2(h,j) = INT(CB(h,ii),j);
              }
            }
            for(size_t j=0; j<lchecku; j++){
              int colSum = 0;
              for(size_t h=0; h<Dim; h++){
                colSum += INT2(h,j);
              }
              if(colSum == (int)Dimm1){
                vert += 1;
                std::vector<int> inter(Dim);
                size_t m = 0;
                for(size_t h=0; h<Dim; h++){
                  if(INT2(h,j)==1){
                    inter[m] = CB(h,ii);
                    m += 1;
                  }
                }
                inter[Dimm1] = k;
                CCtemp.conservativeResize(Eigen::NoChange, vert); 
                for(size_t h=0; h<Dim; h++){
                  CCtemp(h,vert-1) = inter[h]; 
                }
                double lambda = (Uk-VTsum_wu(ii))/(VTsum_cu(j)-VTsum_wu(ii));
                VTtemp.conservativeResize(Eigen::NoChange, vert);
                for(size_t h=0; h<Dim; h++){
                  VTtemp(h,vert-1) = 
                    lambda*VT1_cu(h,j) + (1-lambda)*VT1_wu(h,ii);
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
              CCtemp(h,vert-1) = CCi(h,both[j]);
              VTtemp(h,vert-1) = VTi(h,both[j]);
            }
          }
        }
        VC[i] = vert;
        CC[i] = CCtemp;
        VT[i] = VTtemp;
      }
      
      Eigen::VectorXd WT(N);
      for(size_t j=0; j<N; j++){
        WT(j) = weight.col(j).prod();
      }
      double WTsum = WT.sum();
      WTnorm = WT/WTsum;
      ESS(k) = 1/(WTnorm.dot(WTnorm));
      if(ESS(k)<thresh && k < K[n-Dim-1]){
        std::vector<size_t> N_sons(N,0);
        std::vector<double> dist = Vcumsum(WTnorm);
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
        Eigen::MatrixXd REJJ(lenJJ,Esum(re-1)); 
        for(int jj=0; jj<lenJJ; jj++){
          REJJ.row(jj) = RE.row(JJ(jj));
        }
        Eigen::MatrixXd FEJJ(lenJJ,fe); 
        if(fe>0){
          for(int jj=0; jj<lenJJ; jj++){
            FEJJ.row(jj) = FE.row(JJ(jj));
          }
        }
        std::vector<Eigen::MatrixXd> ZZ(re);
        for(size_t ii=0; ii<re; ii++){
          ZZ[ii].resize(E(ii), 0);
        }
        std::vector<size_t> VCVC(N,0);
        std::vector<Eigen::MatrixXi> CCCC(N);
        std::vector<Eigen::MatrixXd> VTVT(N);
        for(size_t i=0; i<N; i++){
          size_t Nsons_i = N_sons[i];
          if(Nsons_i){
            std::vector<size_t> VCtemp(Nsons_i, VC[i]);
            std::vector<Eigen::MatrixXd> Ztemp(re);
            size_t copy = Nsons_i-1;
            Eigen::MatrixXd VTi = VT[i];
            std::vector<Eigen::MatrixXd> VTtemp(Nsons_i, VTi);
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
                  Eigen::MatrixXd XX(lenJJ, 0);
                  for(size_t jj=0; jj<re; jj++){
                    if(jj != kk){
                      XX.conservativeResize(Eigen::NoChange, XX.cols()+1);
                      Eigen::VectorXd newcol(lenJJ);
                      newcol = 
                        REJJ.block(0, Esum(jj)-E(jj), lenJJ, E(jj)) *
                        Ztemp[jj].col(ii);
                      XX.rightCols(1) = newcol; 
                    }
                  }
                  Eigen::VectorXi v(lenJJ);
                  for(int jj=0; jj<lenJJ; jj++){
                    v(jj) = RE2(JJ(jj),kk);
                  }
                  Eigen::VectorXi vv = cppunique(v);
                  Eigen::VectorXd Z1_(vv.size());
                  for(int jj=0; jj<vv.size(); jj++){
                    Z1_(jj) = Ztemp[kk](vv(jj),ii);
                  }
                  Eigen::MatrixXd CO2_ = 
                    REJJ.block(0, Esum(kk)-E(kk), lenJJ, E(kk));
                  std::vector<double> colsumsCO2(E(kk));
                  std::vector<int> pcolsums;
                  for(int jj=0; jj<E(kk); jj++){
                    double colsum = CO2_.col(jj).sum();
                    if(colsum > 0){
                      pcolsums.push_back(jj);
                    }
                  }
                  size_t ncolsCO2 = pcolsums.size();
                  Eigen::MatrixXd CO2(lenJJ, ncolsCO2);
                  for(size_t jj=0; jj<ncolsCO2; jj++){
                    CO2.col(jj) = CO2_.col(pcolsums[jj]);
                  }
                  std::vector<int> Z00;
                  for(int jj=0; jj<Z1_.size(); jj++){
                    if(Z1_(jj) != 0){
                      Z00.push_back(jj);
                    }
                  }
                  size_t lenZ1 = Z00.size();
                  Eigen::VectorXd Z1(lenZ1);
                  for(size_t jj=0; jj<lenZ1; jj++){
                    Z1(jj) = Z1_(Z00[jj]);
                  }
                  Eigen::MatrixXd XXX(lenJJ, Dimm1);
                  if(fe>0){
                    XXX << FEJJ, XX;
                  }else{
                    XXX = XX;
                  }
                  Eigen::MatrixXd MAT(lenJJ,((int)Dim)-1+ncolsCO2);
                  MAT << -XXX, CO2;
                  Rcpp::List kern = nullSpace(MAT);
                  int rk = kern["rank"];
                  if(rk < MAT.cols()){
                    Eigen::MatrixXd NUL = kern["kernel"];
                    Eigen::MatrixXd n1 = NUL.topRows(NUL.rows()-ncolsCO2);
                    Eigen::MatrixXd n2 = NUL.bottomRows(ncolsCO2);
                    Rcpp::List QR = QRdecomp(n2);
                    Eigen::MatrixXd O2 = QR["Q"];
                    Eigen::MatrixXd R = QR["R"];
                    Eigen::MatrixXd O1 = tsolveAndMultiply(R, n1);
                    //Rcpp::Rcout << "nrow O1 = Dimm1? " << O1.rows() << std::endl;
                    int rankO2 = O2.cols();
                    Eigen::VectorXd a(rankO2);
                    a = O2.transpose() * Z1;
                    //Rcpp::Rcout << "ncol O2 " << O2.cols() << std::endl;
                    //Rcpp::Rcout << "size Z1 " << Z1.size() << std::endl;
                    //Rcpp::Rcout << "nrow O2 = size Z1 ? " << O2.rows() << std::endl;
                    Eigen::VectorXd O2a(lenZ1); 
                    O2a = O2 * a;
                    Eigen::VectorXd tau_(lenZ1);
                    tau_ = Z1 - O2a;
                    double b = sqrt(tau_.dot(tau_));
                    Eigen::VectorXd tau = tau_/b;
                    double bb = sqrt(rchisq(lenZ1-rankO2));
                    double bbb = b/bb;
                    Eigen::VectorXd aa(rankO2);
                    for(int jj=0; jj<rankO2; jj++){
                      aa(jj) = gaussian(generator);
                    }
                    Eigen::VectorXd O2aa(lenZ1);
                    O2aa = O2*aa;
                    Eigen::VectorXd bbtau(lenZ1);
                    bbtau = bb*tau;
                    Eigen::VectorXd MM3(lenZ1);
                    MM3 = O2aa + bbtau;
                    //   Eigen::VectorXd MM3 = O2*aa + bb*tau;
                    //   //Rcpp::Rcout << "nrow O2 " << O2.rows() << std::endl;
                    //   //Rcpp::Rcout << "size tau " << tau.size() << std::endl;
                    //   ////Rcpp::Rcout << "size Z00 " << Z00.size() << std::endl;
                    //   //Rcpp::Rcout << "nrow Ztemp " << Ztemp[kk].rows() << std::endl;
                    //   //Rcpp::Rcout << "size MM3 " << MM3.size() << std::endl;
                    for(size_t jj=0; jj<lenZ1; jj++){
                      Ztemp[kk](Z00[jj],ii) = MM3(jj); // pkoi size Z00 ?
                    }
                    Eigen::VectorXd M2(Dimm1);
                    M2 = O1 * (bbb*aa - a);
                    Eigen::MatrixXd VTtemp_ii = VTtemp[ii];
                    for(size_t jj=0; jj<VC[i]; jj++){
                      for(size_t vert=0; vert<Dim; vert++){
                        if(vert < fe+kk){
                          VTtemp[ii](vert,jj) = VTtemp_ii(vert,jj) - 
                            VTtemp_ii(fe+kk,jj) * M2(vert);
                        }else if(vert > fe+kk){
                          VTtemp[ii](vert,jj) = VTtemp_ii(vert,jj) - 
                            VTtemp_ii(fe+kk,jj) * M2(vert-1);
                        }
                      }
                      VTtemp[ii](fe+kk,jj) = bbb*VTtemp_ii(fe+kk,jj);
                    }
                  }else{
                    double b = sqrt(Z1.dot(Z1));
                    Eigen::VectorXd tau = Z1/b;
                    double bb = sqrt(rchisq(lenZ1));
                    for(size_t jj=0; jj<lenZ1; jj++){
                      Ztemp[kk](Z00[jj],ii) = bb*tau(jj);
                    }
                    VTtemp[ii].block(fe+kk,0,1,VC[i]) *= b/bb;
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
        weight = Eigen::MatrixXd::Ones(E(re-1),N);
      }
    }
    
    //------------determine signs ----------------------------------------------
    Eigen::MatrixXi signs = Eigen::MatrixXi::Zero(re,N);
    for(size_t i=0; i<N; i++){
      Eigen::MatrixXd VTi = VT[i];
      for(size_t j=0; j<re ; j++){
        Eigen::VectorXd row = VTi.row(fe+j);
        bool allneg = true;
        size_t jj = 0;
        while(jj<VC[i] && allneg){
          allneg = row(jj) < 0;
          jj += 1;
        }
        if(allneg){
          signs(j,i) = -1;
        }
      }
    }
    
    //--------FINAL RESAMPLE ---------------------------------------------------
    std::vector<Eigen::MatrixXd> ZZ(re);
    std::vector<Eigen::VectorXi> nn(re);
    std::vector<int> lengths_nn(re);
    std::vector<Eigen::MatrixXd> VTVT(N);
    Eigen::VectorXi n1_(K1.size()+(int)Dim);
    n1_ << K1,K_start;
    Eigen::VectorXi n1 = Vsort(n1_);
    for(size_t ii=0; ii<re; ii++){
      Eigen::VectorXi vec(n1.size());
      for(int jj=0; jj<n1.size(); jj++){
        vec(jj) = RE2(n1(jj),ii);
      }
      nn[ii] = cppunique(vec);
      lengths_nn[ii] = nn[ii].size();
      
      ZZ[ii].resize(lengths_nn[ii], 0);
    }
    for(size_t i=0; i<N; i++){
      Eigen::MatrixXd VTtemp = VT[i];
      std::vector<Eigen::VectorXd> Ztemp(re);
      for(size_t ii=0; ii<re; ii++){
        Ztemp[ii].resize(lengths_nn[ii]);
      }
      for(size_t ii=0; ii<re; ii++){
        // Ztemp[ii].resize(lengths_nn[ii]);
        for(int iii=0; iii<lengths_nn[ii]; iii++){
          Ztemp[ii](iii) = Z[ii](nn[ii](iii),i);
        }
      }
      std::vector<size_t> ord = sample_int(re);
      for(size_t j=0; j<re; j++){
        size_t kk = ord[j];
        Eigen::MatrixXd XX(n,0);
        for(size_t jj=0; jj<re; jj++){
          if(jj != kk){
            XX.conservativeResize(Eigen::NoChange, XX.cols()+1);
            Eigen::VectorXd newcol(n);
            newcol = RE.block(0, Esum(jj)-E(jj), n, 
                              lengths_nn[jj]) * Ztemp[jj];
            XX.rightCols(1) = newcol;
          }
        }
        Eigen::VectorXd Z1 = Ztemp[kk];
        int lenZ1 = Z1.size();
        int ncolsCO2 = lengths_nn[kk];
        Eigen::MatrixXd CO2 = RE.block(0, Esum(kk)-E(kk), n, ncolsCO2);
        Eigen::MatrixXd XXX(n, Dimm1);
        if(fe>0){
          XXX << FE,XX;
        }else{
          XXX = XX;
        }
        Eigen::MatrixXd MAT(n, Dimm1+ncolsCO2);
        MAT << -XXX,CO2;
        Rcpp::List kern = nullSpace(MAT);
        int rk = kern["rank"];
        if(rk < MAT.cols()){
          Eigen::MatrixXd NUL = kern["kernel"];
          Eigen::MatrixXd n1 = NUL.topRows(NUL.rows()-ncolsCO2);
          Eigen::MatrixXd n2 = NUL.bottomRows(ncolsCO2);
          Rcpp::List QR = QRdecomp(n2);
          Eigen::MatrixXd O2 = QR["Q"];
          Eigen::MatrixXd R = QR["R"];
          Eigen::MatrixXd O1 = tsolveAndMultiply(R, n1);
          int rankO2 = O2.cols();
          Eigen::VectorXd a(rankO2);
          a = O2.transpose() * Z1;
          Eigen::VectorXd O2a(lenZ1); 
          O2a = O2 * a;
          Eigen::VectorXd tau_(lenZ1);
          tau_ = Z1 - O2a;
          double b = sqrt(tau_.dot(tau_));
          Eigen::VectorXd tau = tau_/b;
          double bb = sqrt(rchisq(lenZ1-rankO2));
          double bbb = b/bb;
          Eigen::VectorXd aa(rankO2);
          for(int jj=0; jj<rankO2; jj++){
            aa(jj) = gaussian(generator);
          }
          Eigen::VectorXd O2aa(lenZ1);
          O2aa = O2*aa;
          Eigen::VectorXd bbtau(lenZ1);
          bbtau = bb*tau;
          //Eigen::VectorXd MM3(O2.rows());
          //Rcpp::Rcout << "length Ztemp[kk] " << Ztemp[kk].size() << std::endl;
          //Rcpp::Rcout << "= nrow O2 ? " << O2.rows() << std::endl;
          Ztemp[kk] = O2aa + bbtau;
          Eigen::VectorXd M2(Dimm1);
          M2 = O1 * (bbb*aa - a);
          for(size_t jj=0; jj<VC[i]; jj++){
            for(size_t vert=0; vert<Dim; vert++){
              if(vert < fe+kk){
                VTtemp(vert,jj) -= VTtemp(fe+kk,jj) * M2(vert);
              }else if(vert > fe+kk){
                VTtemp(vert,jj) -= VTtemp(fe+kk,jj) * M2(vert-1);
              }
            }
            VTtemp(fe+kk,jj) = bbb*VTtemp(fe+kk,jj);
          }
        }else{
          double b = sqrt(Z1.dot(Z1));
          Eigen::VectorXd tau = Z1/b;
          double bb = sqrt(rchisq(lenZ1));
          Ztemp[kk] = bb*tau;
          VTtemp.row(fe+kk) *= b/bb;
          //VTtemp.block(fe+kk,0,1,VC[i]) *= b/bb; // c'est simplement row fe+kk
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
      VERTEX = pickCoordinates(Dim, N, fe, VT, unif);
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("VERTEX") = VERTEX,
                            Rcpp::Named("WEIGHT") = WTnorm, 
                            Rcpp::Named("ESS") = ESS);
}

// TODO: essaye fidVertex et fidSample dans gfimm - galère avec les 0-based... mais non
// virer weights[i<re-1]
// virer signs=1
// intégrer fidVertex et fidSample au code, et nullSpace et QR
// stocker N_sons[i] dans une variable, lengths_nn[ii], etc
// remplacer les O2.cols(), etc