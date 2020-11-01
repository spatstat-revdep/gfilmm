// fidVertex
Eigen::MatrixXi CCi = CC[i];
double Lk = L.coeff(k);
double Uk = U.coeff(k);
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
Eigen::MatrixXd VTtemp(Dim, 0);
int vert = 0;
size_t lcheckl = checkl.size();
size_t lwhichl = whichl.size();
if(lcheckl < p){
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
  Eigen::VectorXd VTsum_cl(lcheckl);
  Eigen::VectorXd VTsum_wl(lwhichl);
  Eigen::MatrixXd VT1_cl(Dim, lcheckl);
  Eigen::MatrixXd VT1_wl(Dim, lwhichl);
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
        double lambda = (Lk-VTsum_wl(ii))/(VTsum_cl(j)-VTsum_wl(ii));
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
if(lchecku < p){
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
  Eigen::VectorXd VTsum_cu(lchecku);
  Eigen::VectorXd VTsum_wu(lwhichu);
  Eigen::MatrixXd VT1_cu(Dim, lchecku);
  Eigen::MatrixXd VT1_wu(Dim, lwhichu);
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
        double lambda =
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
