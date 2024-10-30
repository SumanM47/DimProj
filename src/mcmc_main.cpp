// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

List Slicesto3D(arma::cube Y, arma::mat S, double bma, int M, int N,
                arma::vec mu_init, arma::vec rho_init, arma::mat A_init, arma::vec sig2_init,
                int niters, int nburn, int nthin, int adapt,
                arma::vec mu0, arma::vec rho0, arma::mat A0, double musig2, double sigsig2,
                double sig2mu, double sig2rho, double sig2A){

  // Bookkeeping and Processing inputs
  int I = Y.n_rows;
  int J = Y.n_cols;
  int n = Y.n_slices;
  int K = rho0.n_elem;

  arma::vec mu = mu_init;
  arma::mat A = A_init;
  arma::vec rho = rho_init;
  arma::vec sig2 = sig2_init;

  // Creating and reshaping the first row of the covariance matrices
  arma::vec d2 = arma::pow(S.col(0) - S(0,0),2) + arma::pow(S.col(1) - S(0,1),2);
  arma::mat D2 = arma::reshape(d2,N,M);
  arma::mat D21 = arma::join_rows(D2,arma::fliplr(D2.cols(1,M-2)));
  arma::mat D2ext = arma::join_cols(D21,arma::flipud(D21.rows(1,N-2)));

  int Next = D2ext.n_rows;
  int Mext = D2ext.n_cols;
  int nbig = Next*Mext;

  arma::mat bigtemp = arma::zeros(Next,Mext);
  bigtemp.submat(0,0,N-1,M-1) = arma::ones(N,M);
  arma::uvec uind = arma::find(arma::vectorise(bigtemp));

  //std::cout << "Checkpoint 1" << std::endl;

  // Creating variables and intermediary terms
  arma::mat Cs = arma::zeros(I,K);
  double c1 = sqrt(arma::datum::pi)*bma;
  double temp1=0.0;
  arma::vec temp2=arma::regspace(0,I-1);
  arma::cube res = arma::zeros(I,J,n);
  arma::cube bigres = arma::zeros(I,J,nbig);
  arma::mat quad = arma::zeros(I,J);
  arma::mat bigquad = quad;
  arma::cube U = arma::zeros(n,I,K);
  arma::mat ez = arma::zeros(Next,Mext);
  ez(0,0) = 1.0;
  arma::cx_mat tempfftbigU2 = arma::fft2(ez);
  arma::vec diaginv = arma::ones(K);

  arma::cube SigCube(Next,Mext,K);

  arma::mat Cmat = arma::zeros(I,I), Cmatchol = arma::zeros(I,I), iCmatchol = arma::zeros(I,I);


  //std::cout << "Checkpoint 1.1" << std::endl;

  // Computing the scale terms
  for(int kind=0;kind<K;kind++){
    temp1 = rho(kind);
    SigCube.slice(kind) = arma::exp(-D2ext/(2*temp1*temp1));
    Cs.col(kind) = c1*((temp2+1)%arma::normcdf((temp2+1)*bma/temp1) + (temp2-1)%arma::normcdf((temp2-1)*bma/temp1) - 2*temp2%arma::normcdf(temp2*bma/temp1)) + (temp1/sqrt(2.0))*(arma::exp(-0.5*arma::pow((temp2+1)*bma/temp1,2)) + arma::exp(-0.5*arma::pow((temp2-1)*bma/temp1,2)) - 2*arma::exp(-0.5*arma::pow((temp2)*bma/temp1,2)));
  }


  // Generating initial latent process and computing related terms
  arma::mat tempnorm = arma::zeros(Next,Mext);
  arma::cx_mat tempSigeig(tempnorm,tempnorm);
  arma::cx_mat tempfftbigU(tempnorm,tempnorm);
  arma::cx_mat tempbigUquadhalf(tempnorm,tempnorm);

  arma::cube bigU = arma::zeros(nbig,I,K);
  arma::cube bigUstar = bigU;
  arma::mat bigUquad = arma::zeros(I,K);
  arma::cube bigUgradpart = arma::zeros(nbig,I,K);
  arma::vec logdets = arma::zeros(K);
  arma::mat tempbigU = arma::zeros(Next,Mext);

  //std::cout << "Checkpoint 1.2" << std::endl;

  for(int kind=0; kind<K; kind++){
    tempSigeig = arma::fft2(SigCube.slice(kind));
    logdets(kind) = arma::accu(arma::log(arma::real(tempSigeig)));
    for(int iii=0; iii<I; iii++){
    tempnorm = arma::randn(Next,Mext);
    tempbigU = arma::real(arma::ifft2(arma::sqrt(tempSigeig)%arma::fft2(tempnorm)));
    U.slice(kind).col(iii) = arma::vectorise(tempbigU.submat(0,0,N-1,M-1));
    tempfftbigU = arma::fft2(tempbigU);
    tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
    bigUquad(iii,kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;
    tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU2;
    diaginv(kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;
    bigUgradpart.slice(kind).col(iii) = arma::vectorise(arma::real(arma::ifft2((1/tempSigeig)%tempfftbigU)));
    bigU.slice(kind).col(iii) = arma::vectorise(tempbigU);
    }
    bigUstar.slice(kind) = bigU.slice(kind)*arma::chol(arma::toeplitz(Cs.col(kind)));
  }

  //std::cout << "Checkpoint 1.3" << std::endl;

  // Computing intermediary terms
  arma::vec tempres = arma::zeros(n);
  arma::vec tempres2 = tempres, tempres3 = tempres2;
  arma::mat bigUstari = arma::zeros(nbig,K), bigUstark = arma::zeros(nbig,I), Ustari = arma::zeros(n,K);

  for(int iind=0;iind<I;iind++){
    bigUstari = bigUstar.col_as_mat(iind);
    Ustari = bigUstari.rows(uind);
  for(int jind=0;jind<J;jind++){
    tempres = mu(jind)*bma + Ustari*A.row(jind).t();
      tempres2 = Y(arma::span(iind),arma::span(jind),arma::span::all);
      tempres3 = tempres2 - tempres;
      res.tube(iind,jind) = tempres3;
      quad(iind,jind) = arma::sum(tempres3%tempres3);
    }
  }

  bigres.slices(uind) = res;
  bigquad = arma::sum(arma::pow(bigres,2),2);

  //std::cout << "Checkpoint 2" << std::endl;

  // Initiating intermediary terms required for the MCMC

  // Likelihoods
  double curlik_mu, canlik_mu, curlik_A, canlik_A, curlik_rho, canlik_rho, curlik_sig2, canlik_sig2, curlik_bigU, canlik_bigU;

  // Gradients and candidates
  arma::vec curgrad_mu = arma::zeros(J), cangrad_mu =  arma::zeros(J), can_mu = mu;
  // For Individual Updates of A
  //arma::vec curgrad_aj = arma::zeros(K), cangrad_aj =  arma::zeros(K), tempA = arma::zeros(K);
  // For Block Update of A
  arma::mat curgrad_a = arma::zeros(J,K), cangrad_a = arma::zeros(J,K);
  arma::vec tempA = arma::zeros(K);

  arma::vec curgrad_sig2 = arma::zeros(J), cangrad_sig2 =  arma::zeros(J), can_sig2 = sig2;
  arma::mat can_A = A;
  double can_rho = 0;
  arma::vec canCs = arma::zeros(I);
  //arma::mat curgrad_sig2 = arma::zeros(I,J), cangrad_sig2 =  arma::zeros(I,J), can_sig2 = sig2;

  //std::cout << "Checkpoint 2.1" << std::endl;

  //arma::mat curgrad_bigU = arma::zeros(Next,Mext), cangrad_bigU = arma::zeros(Next,Mext);
  //arma::vec curgrad_bigU = arma::zeros(nbig), cangrad_bigU= arma::zeros(nbig); //Individual Update
  arma::cube curgrad_bigU = arma::zeros(nbig,I,K), cangrad_bigU = arma::zeros(nbig,I,K); //Joint Update
  //arma::mat subres_row = arma::zeros(J,n), subres_col= arma::zeros(I,n), subres_slice= arma::zeros(I,J);
  arma::mat subres_row = arma::zeros(J,nbig), subres_col= arma::zeros(I,nbig), subres_slice= arma::zeros(I,J);

  // Initiating and calculating the inverse Hessian matrices
  // diagonal Hessian matrices turned to vectors
  // block diagonal matrices made into cube to use each slice independently
  //arma::vec nHessinv_mu = arma::ones(J), nHessinv_U= arma::ones(K);
  arma::vec nHessinv_mu = arma::ones(J), nHessinv_sig2 = arma::ones(J);
  //arma::mat nHessinv_A = arma::ones(K,J);
  arma::cube nHessinv_A = arma::ones(K,K,J), nHessinv_A_chol = arma::ones(K,K,J), nHessinv_A_cholinv = arma::ones(K,K,J);
  //arma::mat nHessinv_sig2 = arma::ones(I,J);
  //arma::cube nHessinv_U = arma::ones(nbig,I,K);
  //x//arma::vec nHessinv_U = arma::ones(K);
  arma::mat nHessinv_U = arma::ones(I*K,I*K), nHess_U= arma::ones(I*K,I*K);
  //arma::cube arma::mat nHessinv_A = arma::ones(K,K,J), nHessinv_U = arma::ones(Next,Mext,K);

  //std::cout << "Checkpoint 2.2" << std::endl;

  double qprop, qcur, la, canlogdets;
  //arma::vec canbigUquad = arma::zeros(I); // Individual U Update
  arma::mat canbigUquad = arma::zeros(I,K); // Joint U update
  //arma::mat canquad = quad, canSigCube = SigCube.slice(0), canbigUgradpart = bigUgradpart.slice(0), canU = U;
  //arma::mat canSigCube = SigCube.slice(0), canbigUgradpart = bigUgradpart.slice(0), canbigquad=bigquad; // Individual Update
  arma::mat canSigCube = SigCube.slice(0), canbigquad=bigquad; // Joint U update
  arma::cube canbigUgradpart = bigUgradpart; // Joint U update
  //arma::cube canres = res, can_bigU = bigU;
  arma::cube canbigres = bigres, can_bigU = bigU, canbigUstar = bigUstar;

  //std::cout << "Checkpoint 2.3" << std::endl;

  //for(int jind=0;jind<J;jind++){
    //nHessinv_mu(jind) = 1/((1/sig2mu) + n*bma*bma*arma::sum(1/(sig2.col(jind))));
    nHessinv_mu = 1/((1/sig2mu) + nbig*bma*bma*(I/sig2));
  //}

  //x//nHessinv_sig2 = 1/(1/sigsig2 + 0.5*quad/sig2);
  nHessinv_sig2 = 1/((1/sigsig2) + 0.5*(arma::sum(bigquad,0).t()/sig2));

  //std::cout << "Checkpoint 2.4" << std::endl;

  //x//arma::vec dUtU = arma::sum(arma::pow(U,2),0).t();
  //x//arma::vec dUtU = arma::zeros(K);
  //x//arma::mat dUtU = arma::zeros(K,K);
  //arma::mat dUtU = arma::zeros(K,I);
  arma::vec tempHessA = arma::zeros(K);
  //x//for(int kind=0;kind<K;kind++){
    //x//dUtU(kind) = arma::accu(arma::pow(bigUstar.slice(kind),2));
  //x//}
  //for(int iind=0; iind < I; iind++){
    //dUtU.col(iind) = arma::trans(arma::sum(arma::pow(bigUstar.col_as_mat(iind),2),0));
  //}
  //tempHessA = arma::zeros(K,1);
  //for(int iind=0;iind<I;iind++){
    //tempHessA += dUtU.col(iind);
  //}
  //for(int jind=0;jind<J;jind++){
      //nHessinv_A.col(jind) = 1/((1/sig2A) + (tempHessA/sig2(jind)));
  //}

  arma::mat UtU = arma::zeros(K,K);
  for(int iind=0;iind<I;iind++){
    bigUstari = bigUstar.col_as_mat(iind);
    UtU += bigUstari.t()*bigUstari;
  }

  for(int jind=0;jind<J;jind++){
    arma::mat tempmat = UtU/sig2(jind) + arma::eye(K,K)/sig2A;
    nHessinv_A.slice(jind) = arma::inv(tempmat);
    nHessinv_A_chol.slice(jind) = arma::chol(nHessinv_A.slice(jind),"lower");
    nHessinv_A_cholinv.slice(jind) = arma::inv(nHessinv_A_chol.slice(jind));
  }

  //std::cout << "Checkpoint 2.5" << std::endl;

  arma::mat dCmatchol = arma::zeros(I,K);

  for(int kind=0;kind<K;kind++){
    arma::mat Cm = arma::toeplitz(Cs.col(kind));
    bool icholsuc = false;
    while(icholsuc == false){
      icholsuc = arma::chol(iCmatchol,Cm);
      if(icholsuc == false){
        Cm += arma::eye(I,I)*1e-6;
        Cs(0,kind) += 1e-6;
      }
    }
    dCmatchol.col(kind) = diagvec(iCmatchol);
    for(int iind=0;iind<I;iind++){
    //double tempscalebigU = 0;
      //double tempnHessU = diaginv(kind) + dCmatchol(iind,kind)*dCmatchol(iind,kind)*arma::sum(arma::pow(A.col(kind),2)/sig2);
    //tempnHessU = diaginv(kind)*arma::ones(nbig,I);
    //for(int iind=0; iind<I; iind++){
      //tempnHessU.col(iind) += arma::sum(arma::pow(A.col(kind),2)/sig2);
    //}
      for(int k2ind=0; k2ind<kind; k2ind++){
        for(int i2ind=0; i2ind<=iind; i2ind++){
          nHess_U(kind*I+iind,k2ind*I+i2ind) = dCmatchol(iind,kind)*dCmatchol(i2ind,k2ind)*arma::sum(A.col(kind)%A.col(k2ind)/sig2);
          nHess_U(k2ind*I+i2ind,kind*I+iind) = dCmatchol(iind,kind)*dCmatchol(i2ind,k2ind)*arma::sum(A.col(kind)%A.col(k2ind)/sig2);
        }
      }

      for(int i2ind=0; i2ind<=iind; i2ind++){
        nHess_U(kind*I+iind,kind*I+i2ind) = diaginv(kind) + dCmatchol(iind,kind)*dCmatchol(i2ind,kind)*arma::sum(arma::pow(A.col(kind),2)/sig2);
        nHess_U(kind*I+i2ind,kind*I+iind) = diaginv(kind) + dCmatchol(iind,kind)*dCmatchol(i2ind,kind)*arma::sum(arma::pow(A.col(kind),2)/sig2);
      }
      //nHessinv_U(iind,kind) = 1/tempnHessU;
    //nHessinv_U.slice(kind).submat(0,0,N-1,M-1) = (1/(diaginv(kind)/cs(kind) + tempscalebigU))*arma::ones(N,M);
    //nHessinv_U.slice(kind) = (1/(diaginv(kind)/cs(kind) + tempscalebigU))*arma::ones(Next,Mext);
    }
  }

  nHessinv_U = inv(nHess_U);
  arma::mat nHessinv_U_ch = arma::eye(I*K,I*K);
  bool ucholsuc = false;
  while(ucholsuc==false){
    ucholsuc = arma::chol(nHessinv_U_ch,nHessinv_U);
    if(ucholsuc==false){
      nHessinv_U += arma::eye(I*K,I*K)*1e-6;
    }
  }
  arma::mat nHessinv_U_chinv = inv(nHessinv_U_ch);

  //nHessinv_U = arma::ones(K);

  //std::cout << "done with Hessians" << std::endl;
  // Storage
  arma::mat keep_mu = arma::zeros(niters,J);
  arma::cube keep_A = arma::zeros(niters,J,K);
  arma::mat keep_rho = arma::zeros(niters,K);
  arma::mat keep_sig2 = arma::zeros(niters,J);
  //arma::mat acc_U = arma::zeros(I,K);
  double acc_U = 0.0;

  double acc_mu = 0.0;
  //arma::vec acc_A = arma::zeros(J);
  double acc_Aall = 0;
  arma::vec acc_rho = arma::zeros(K);
  double acc_sig2 = 0.0;

  // proposal variances
  double h_mu0 = 2.7225, h_mu = h_mu0;
  //arma::vec h_A0 = (5.6644/(K))*arma::ones(J), h_A = h_A0;
  double h_A0all = 5.6644/(J*K), h_Aall = h_A0all;
  arma::vec h_rho0 = (5.6644)*arma::ones(K)/(0.5*rho*nbig*I + 1/sig2rho), h_rho = h_rho0;
  double h_sig20 = (5.6644/(J)), h_sig2 = h_sig20;
  //arma::mat h_bigU0 = (2.7225/pow(nbig,0.33))*arma::ones(I,K), h_bigU = h_bigU0;
  double h_bigU0 = (2.7225/pow(nbig*I*K,0.33)), h_bigU = h_bigU0;

  double adp = nburn/2;

  //double h_mu0 = 0.5, h_mu = h_mu0;
  //arma::vec h_A0 = 0.5*arma::ones(J), h_A = h_A0;
  //arma::vec h_rho0 = 0.025*arma::ones(K), h_rho = h_rho0;
  //double h_sig20 = 0.5, h_sig2 = h_sig20;
  //arma::mat h_bigU0 = 0.5*arma::ones(I,K), h_bigU = h_bigU0;

  //std::cout << "Checkpoint 3" << std::endl;

  //std::cout << "Here we go" << std::endl;

  // GO!
  for(int i=0;i<niters;i++){
    for(int ii=0;ii<nthin;ii++){

      // Update mu

      //curlik_mu = -0.5*arma::sum(arma::pow(mu-mu0,2))/sig2mu - 0.5*arma::accu(quad/sig2);
      curlik_mu = -0.5*arma::sum(arma::pow(mu-mu0,2))/sig2mu - 0.5*arma::sum(arma::sum(bigquad,0)/sig2.t());
      //subres_slice = arma::sum(res,2);
      subres_slice = arma::sum(bigres,2);
      curgrad_mu = -(mu - mu0)/sig2mu + bma*(arma::sum(subres_slice,0).t()/sig2);

      canbigres = bigres;
      canbigquad = bigquad;
      can_mu = mu + 0.5*h_mu*(nHessinv_mu%curgrad_mu) + sqrt(h_mu)*(arma::sqrt(nHessinv_mu)%arma::randn(J));

      for(int jind=0;jind<J;jind++){
        //tempres = can_mu(jind)*bma + U*A.row(jind).t();
        //for(int iind=0;iind<I;iind++){
          //tempres2 = Y(arma::span(iind),arma::span(jind),arma::span::all);
          //tempres3 = tempres2 - tempres;
          //canres.tube(iind,jind) = tempres3;
          //canquad(iind,jind) = arma::sum(tempres3%tempres3);
        //}
        canbigres(arma::span::all,arma::span(jind),arma::span::all) += (mu(jind) - can_mu(jind))*bma;
      }
      canbigquad = arma::sum(arma::pow(canbigres,2),2);

      //canlik_mu = -0.5*arma::sum(arma::pow(can_mu-mu0,2))/sig2mu - 0.5*arma::accu(canquad/sig2);
      canlik_mu = -0.5*arma::sum(arma::pow(can_mu-mu0,2))/sig2mu - 0.5*arma::sum(arma::sum(canbigquad,0)/sig2.t());
      //subres_slice = arma::sum(canres,2);
      subres_slice = arma::sum(canbigres,2);
      cangrad_mu = -(can_mu - mu0)/sig2mu + bma*(arma::sum(subres_slice,0).t()/sig2);

      qprop = -0.5*arma::sum(arma::pow(can_mu - mu - 0.5*h_mu*(nHessinv_mu%curgrad_mu),2)/nHessinv_mu)/h_mu;
      qcur = -0.5*arma::sum(arma::pow(mu - can_mu - 0.5*h_mu*(nHessinv_mu%cangrad_mu),2)/nHessinv_mu)/h_mu;

      la = canlik_mu + qcur - curlik_mu - qprop;
      if(log(arma::randu()) < la){
        //for(int jind=0;jind<J;jind++){
          //bigres(arma::span::all,arma::span(jind),arma::span::all) += mu(jind) - can_mu(jind);
        //}
        mu = can_mu;
        //res = canres;
        //bigres.slices(uind) = res;
        //quad = canquad;
        bigres=canbigres;
        bigquad=canbigquad;
        acc_mu += 1.0/nthin;
      }

      //std::cout << "Checkpoint 4" << std::endl;

      // Update A

      // Block Update
      curlik_A = 0.0;
      canlik_A = 0.0;
      qprop=0.0;
      qcur=0.0;
      for(int jind=0;jind<J;jind++){
        curlik_A += -0.5*arma::sum(arma::pow(A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(bigquad.col(jind))/sig2(jind);
        for(int kind=0;kind<K;kind++){
          double tempdA = 0.0;
          for(int iind=0;iind<I;iind++){
            arma::vec tempU = bigUstar.slice(kind).col(iind);
            arma::vec nbigtemp = bigres.tube(iind,jind);
            tempdA += arma::sum(tempU%nbigtemp);
          }
          curgrad_a(jind,kind) = -(A(jind,kind) - A0(jind,kind))/sig2A + tempdA/sig2(jind);
        }
      }

      can_A = A;
      canbigres = bigres;
      canbigquad = bigquad;

      for(int jind=0;jind<J;jind++){
        //tempA = A.row(jind).t() + 0.5*h_Aall*(nHessinv_A.col(jind)%curgrad_a.row(jind).t()) + sqrt(h_Aall)*(arma::sqrt(nHessinv_A.col(jind))%arma::randn(K));
        tempA = A.row(jind).t() + 0.5*h_Aall*(nHessinv_A.slice(jind)*curgrad_a.row(jind).t()) + sqrt(h_Aall)*(nHessinv_A_chol.slice(jind)*arma::randn(K));
        can_A.row(jind) = tempA.t();
      }

      for(int jind=0;jind<J;jind++){
        arma::mat bigUa = arma::zeros(nbig,I);
        for(int kind=0;kind<K;kind++){
          bigUa += bigUstar.slice(kind)*(A(jind,kind) - can_A(jind,kind));
        }

        for(int iind=0;iind<I;iind++){
          canbigres.tube(iind,jind) += bigUa.col(iind);
        }
      }
        canbigquad = arma::sum(arma::pow(canbigres,2),2);

      for(int jind=0;jind<J;jind++){
        for(int kind=0;kind<K;kind++){
          double tempdA = 0.0;
          for(int iind=0;iind<I;iind++){
            arma::vec tempU = bigUstar.slice(kind).col(iind);
            arma::vec nbigtemp = canbigres.tube(iind,jind);
            tempdA += arma::sum(tempU%nbigtemp);
          }
          cangrad_a(jind,kind) = -(can_A(jind,kind) - A0(jind,kind))/sig2A + tempdA/sig2(jind);
        }
        }

        for(int jind=0;jind<J;jind++){
          canlik_A += -0.5*arma::sum(arma::pow(can_A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(canbigquad.col(jind))/sig2(jind);

          //qprop += -0.5*arma::sum(arma::pow((can_A.row(jind).t() - A.row(jind).t() - 0.5*h_Aall*(nHessinv_A.col(jind)%curgrad_a.row(jind).t())),2)/nHessinv_A.col(jind))/h_Aall;
          //qcur += -0.5*arma::sum(arma::pow((A.row(jind).t() - can_A.row(jind).t() - 0.5*h_Aall*(nHessinv_A.col(jind)%cangrad_a.row(jind).t())),2)/nHessinv_A.col(jind))/h_Aall;

          qprop += -0.5*arma::sum(arma::pow(nHessinv_A_cholinv.slice(jind)*(can_A.row(jind).t() - A.row(jind).t() - 0.5*h_Aall*(nHessinv_A.slice(jind)*curgrad_a.row(jind).t())),2))/h_Aall;
          qcur += -0.5*arma::sum(arma::pow(nHessinv_A_cholinv.slice(jind)*(A.row(jind).t() - can_A.row(jind).t() - 0.5*h_Aall*(nHessinv_A.slice(jind)*cangrad_a.row(jind).t())),2))/h_Aall;


      }

      la = canlik_A + qcur - curlik_A - qprop;
      if(log(arma::randu()) < la){
        A = can_A;
        bigres = canbigres;
        bigquad = canbigquad;
        acc_Aall += 1.0/nthin;
      }

      // Individual Update
      //for(int jind=0;jind<J;jind++){
      //curlik_A = -0.5*arma::sum(arma::pow(A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(bigquad.col(jind))/sig2(jind);

      //for(int kind=0;kind<K;kind++){
          //double tempdA = 0.0;
          //for(int iind=0;iind<I;iind++){
            //arma::vec tempU = bigUstar.slice(kind).col(iind);
            //arma::vec nbigtemp = bigres.tube(iind,jind);
            //tempdA += arma::sum(tempU%nbigtemp);
          //}
          //curgrad_aj(kind) = -(A(jind,kind) - A0(jind,kind))/sig2A + tempdA/sig2(jind);
        //}

        //std::cout<<"Checkpoint 4.1"<< std::endl;

        //can_A = A;
        //canbigres = bigres;
        //canbigquad = bigquad;


        //tempA = A.row(jind).t() + 0.5*h_A(jind)*(nHessinv_A.col(jind)%curgrad_aj) + sqrt(h_A(jind))*(arma::sqrt(nHessinv_A.col(jind))%arma::randn(K));
        //can_A.row(jind) = tempA.t();


        //arma::mat bigUa = arma::zeros(nbig,I);
        //for(int kind=0;kind<K;kind++){
          //bigUa += bigUstar.slice(kind)*(A(jind,kind) - can_A(jind,kind));
        //}
        //for(int iind=0;iind<I;iind++){
          //canbigres.tube(iind,jind) += bigUa.col(iind);
        //}
        //canbigquad = arma::sum(arma::pow(canbigres,2),2);

        //for(int kind=0;kind<K;kind++){
          //double tempdA = 0.0;
          //for(int iind=0;iind<I;iind++){
            //arma::vec tempU = bigUstar.slice(kind).col(iind);
            //arma::vec nbigtemp = canbigres.tube(iind,jind);
            //tempdA += arma::sum(tempU%nbigtemp);
          //}
          //cangrad_aj(kind) = -(can_A(jind,kind) - A0(jind,kind))/sig2A + tempdA/sig2(jind);
        //}

        //std::cout<<"Checkpoint 4.2"<< std::endl;

        //qprop = -0.5*arma::sum(arma::pow((can_A.row(jind).t() - A.row(jind).t() - 0.5*h_A(jind)*(nHessinv_A.col(jind)%curgrad_aj)),2)/nHessinv_A.col(jind))/h_A(jind);
        //qcur = -0.5*arma::sum(arma::pow((A.row(jind).t() - can_A.row(jind).t() - 0.5*h_A(jind)*(nHessinv_A.col(jind)%cangrad_aj)),2)/nHessinv_A.col(jind))/h_A(jind);


      //canlik_A = -0.5*arma::sum(arma::pow(can_A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(canbigquad.col(jind))/sig2(jind);

      //la = canlik_A + qcur - curlik_A - qprop;
      //if(log(arma::randu()) < la){
        //A = can_A;
        //bigres = canbigres;
        //bigquad = canbigquad;
        //acc_A(jind) += 1.0/nthin;
      //}
      //}

      //std::cout << "size of A:" << arma::size(A) << std::endl;
      //std::cout << "Checkpoint 5" << std::endl;

      // Update rho
      // Likelihood now includes the residual part since Ustar has C_k in it. Need to update accordingly

      arma::mat can_dCmatchol = dCmatchol;

      for(int kind=0;kind<K;kind++){
        curlik_rho = -0.5*pow(log(rho(kind)/rho0(kind)),2)/sig2rho  - 0.5*I*logdets(kind) - 0.5*arma::accu(bigUquad.col(kind)) -0.5*arma::sum(arma::sum(bigquad,0)/sig2.t());


        can_rho = exp(log(rho(kind)) + sqrt(h_rho(kind))*arma::randn());

        canSigCube = arma::exp(-D2ext/(2*can_rho*can_rho));
        canCs = c1*((temp2+1)%arma::normcdf((temp2+1)*bma/can_rho) + (temp2-1)%arma::normcdf((temp2-1)*bma/can_rho) - 2*temp2%arma::normcdf(temp2*bma/can_rho)) + (can_rho/sqrt(2.0))*(arma::exp(-0.5*arma::pow((temp2+1)*bma/can_rho,2)) + arma::exp(-0.5*arma::pow((temp2-1)*bma/can_rho,2)) - 2*arma::exp(-0.5*arma::pow((temp2)*bma/can_rho,2)));


        tempSigeig = arma::fft2(canSigCube);
        canlogdets = arma::accu(arma::log(arma::real(tempSigeig)));


        canbigUquad.col(kind) = bigUquad.col(kind); // Joint U update
        canbigUgradpart.slice(kind) = bigUgradpart.slice(kind); // Joint U update
        //canbigUquad = bigUquad.col(kind); // Individual U update
        //canbigUgradpart = bigUgradpart.slice(kind); // Individual U update
        canbigUstar = bigUstar;
        canbigres = bigres;
        canbigquad = bigquad;

        for(int iii=0; iii<I; iii++){
        tempfftbigU = arma::fft2(arma::reshape(bigU.slice(kind).col(iii),Next,Mext));
        tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
        canbigUquad(iii,kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig; //Joint U Update
        canbigUgradpart.slice(kind).col(iii) = arma::vectorise(arma::real(arma::ifft2((1/tempSigeig)%tempfftbigU))); // Joint U Update
        //canbigUquad(iii) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig; //Individual U Update
        //canbigUgradpart.col(iii) = arma::vectorise(arma::real(arma::ifft2((1/tempSigeig)%tempfftbigU))); // Individual U Update
        }

        bool cholsuc = false;
        Cmat = arma::toeplitz(canCs);
        while(cholsuc == false){
          cholsuc = arma::chol(Cmatchol,Cmat);
          if(cholsuc == false){
            Cmat += arma::eye(I,I)*1e-6;
            canCs(0) += 1e-6;
          }
        }

        can_dCmatchol.col(kind) = diagvec(Cmatchol);

        canbigUstar.slice(kind) = bigU.slice(kind)*Cmatchol;
        for(int iind=0; iind<I; iind++){
          for(int jind=0; jind<J; jind++){
            canbigres.tube(iind,jind) += (bigUstar.col_as_mat(iind) - canbigUstar.col_as_mat(iind))*(A.row(jind).t());
          }
        }
        canbigquad = arma::sum(arma::pow(canbigres,2),2);

        canlik_rho = -0.5*pow(log(can_rho/rho0(kind)),2)/sig2rho - 0.5*I*canlogdets - 0.5*arma::accu(canbigUquad.col(kind)) -0.5*arma::sum(arma::sum(canbigquad,0)/sig2.t());

        la = canlik_rho - curlik_rho;
        if(log(arma::randu()) < la){
          rho(kind) = can_rho;
          Cs.col(kind) = canCs;
          SigCube.slice(kind) = canSigCube;
          logdets(kind) = canlogdets;
          bigUquad.col(kind) = canbigUquad.col(kind); // Joint U Update
          bigUgradpart.slice(kind) = canbigUgradpart.slice(kind); // Joint U Update
          //bigUquad.col(kind) = canbigUquad; // Individual U Update
          //bigUgradpart.slice(kind) = canbigUgradpart; // Individual U Update
          bigres=canbigres;
          bigquad=canbigquad;
          dCmatchol=can_dCmatchol;
          acc_rho(kind) += 1.0/nthin;
        }

      }

      //std::cout << "Checkpoint 6" << std::endl;

      // Update sig2

      //x//curlik_sig2 = - arma::accu(arma::pow(arma::log(sig2)-musig2,2)/(2*sigsig2)) -0.5*n*arma::accu(arma::log(sig2)) -0.5*arma::accu(quad/sig2);
      curlik_sig2 = - arma::sum(arma::pow(arma::log(sig2)-musig2,2)/(2*sigsig2)) -0.5*nbig*I*arma::sum(arma::log(sig2)) -0.5*arma::sum(arma::sum(bigquad,0).t()/sig2);
      //x//curgrad_sig2 = -0.5*n + 0.5*(quad/sig2) - (arma::log(sig2) - musig2)/sigsig2;
      curgrad_sig2 = -0.5*nbig*I + 0.5*arma::sum(bigquad,0).t()/sig2 - (arma::log(sig2) - musig2)/sigsig2;

      can_sig2 = arma::exp(arma::log(sig2) + 0.5*h_sig2*(nHessinv_sig2%curgrad_sig2) + sqrt(h_sig2)*(arma::sqrt(nHessinv_sig2)%arma::randn(J)));


      //x//canlik_sig2 = - arma::accu(arma::pow(arma::log(can_sig2)-musig2,2)/(2*sigsig2)) -0.5*n*arma::accu(arma::log(can_sig2)) -0.5*arma::accu(quad/can_sig2);
      canlik_sig2 = - arma::sum(arma::pow(arma::log(can_sig2)-musig2,2)/(2*sigsig2)) -0.5*nbig*I*arma::sum(arma::log(can_sig2)) -0.5*arma::sum(arma::sum(bigquad,0).t()/can_sig2);
      //x//cangrad_sig2 = -0.5*n + 0.5*(quad/can_sig2) - (arma::log(can_sig2) - musig2)/sigsig2;
      cangrad_sig2 = -0.5*nbig*I + 0.5*arma::sum(bigquad,0).t()/can_sig2 - (arma::log(can_sig2) - musig2)/sigsig2;

      qprop = -0.5*arma::sum(arma::pow(arma::log(can_sig2) - arma::log(sig2) - 0.5*h_sig2*(nHessinv_sig2%curgrad_sig2),2)/nHessinv_sig2)/h_sig2;
      qcur = -0.5*arma::sum(arma::pow(arma::log(sig2) - arma::log(can_sig2) - 0.5*h_sig2*(nHessinv_sig2%cangrad_sig2),2)/nHessinv_sig2)/h_sig2;

      la = canlik_sig2 + qcur - curlik_sig2 - qprop;
      if(log(arma::randu()) < la){
        sig2 = can_sig2;
        acc_sig2 += 1.0/nthin;
      }

      //std::cout << "Checkpoint 7" << std::endl;

      //std::cout << "Checkpoint 8" << std::endl;

      // Update bigU and U

      //Joint update

      curlik_bigU = -0.5*arma::accu(bigUquad) - 0.5*arma::sum(arma::sum(bigquad,0)/sig2.t());


      for(int kind=0;kind<K;kind++){
        for(int iind=0;iind<I;iind++){
          arma::mat tempgradpartU = arma::zeros(1,nbig);
          subres_row = arma::sum(bigres,0);
          tempgradpartU = ((A.col(kind)/sig2).t())*subres_row;

          curgrad_bigU.slice(kind).col(iind) = -bigUgradpart.slice(kind).col(iind) + dCmatchol(iind,kind)*tempgradpartU.t();
        }
      }

      can_bigU = bigU;
      canbigres=bigres;
      canbigUquad=bigUquad;
      canbigUstar = bigUstar;
      canbigUgradpart = bigUgradpart;
      canbigquad = bigquad;

      arma::cube rn = arma::randn(nbig,I,K);

      for(int kind=0;kind<K;kind++){
        for(int iind=0;iind<I;iind++){
          //can_bigU.slice(kind).col(iind) = bigU.slice(kind).col(iind) + 0.5*h_bigU*(nHessinv_U(iind,kind)*curgrad_bigU.slice(kind).col(iind)) + sqrt(h_bigU)*(sqrt(nHessinv_U(iind,kind))*arma::randn(nbig));
          //can_bigU.slice(kind).col(iind) = bigU.slice(kind).col(iind) + 0.5*h_bigU*(nHessinv_u(kind*I+iind,)*curgrad_bigU.slice(kind).col(iind)) + sqrt(h_bigU)*(nHessinv_U_ch(iind,kind))

          arma::vec tempumid = arma::zeros(nbig), tempumid2 = arma::zeros(nbig);
          for(int k2ind=0; k2ind<K; k2ind++){
            for(int i2ind=0; i2ind<I; i2ind++){
              tempumid += nHessinv_U(kind*I+iind,k2ind*I+i2ind)*curgrad_bigU.slice(k2ind).col(i2ind);
              tempumid2 += nHessinv_U_ch(kind*I+iind,k2ind*I+i2ind)*rn.slice(k2ind).col(i2ind);
            }
          }
          can_bigU.slice(kind).col(iind) = bigU.slice(kind).col(iind) + 0.5*h_bigU*tempumid + sqrt(h_bigU)*tempumid2;

          tempSigeig = arma::fft2(SigCube.slice(kind));
          tempfftbigU = arma::fft2(arma::reshape(can_bigU.slice(kind).col(iind),Next,Mext));
          tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
          canbigUquad(iind,kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;
          canbigUgradpart.slice(kind).col(iind) = arma::vectorise(arma::real(arma::ifft2((1/tempSigeig)%tempfftbigU)));

          canbigUstar.slice(kind) = can_bigU.slice(kind)*arma::chol(arma::toeplitz(Cs.col(kind)));

        //std::cout << "Checkpoint 8.2" << std::endl;
          for(int jind=0;jind<J;jind++){
            canbigres.tube(iind,jind) += (bigUstar.col_as_mat(iind) - canbigUstar.col_as_mat(iind))*(A.row(jind).t());
          }
        }
      }

      canbigquad = arma::sum(arma::pow(canbigres,2),2);



      canlik_bigU = -0.5*arma::accu(canbigUquad) - 0.5*arma::sum(arma::sum(canbigquad,0)/sig2.t());

      for(int kind=0;kind<K;kind++){
        for(int iind=0;iind<I;iind++){
          arma::mat tempgradpartU = arma::zeros(1,nbig);
          subres_row = arma::sum(canbigres,0);
          tempgradpartU = ((A.col(kind)/sig2).t())*subres_row;
          cangrad_bigU.slice(kind).col(iind)= -canbigUgradpart.slice(kind).col(iind) + dCmatchol(iind,kind)*tempgradpartU.t();
        }
      }

      qcur=0;
      qprop=0;


      arma::cube resc = arma::zeros(nbig,I,K), resp = arma::zeros(nbig,I,K);

      for(int kind=0;kind<K;kind++){
        for(int iind=0;iind<I;iind++){
          arma::vec rescur = arma::zeros(nbig), resprop = arma::zeros(nbig);
          //qprop += -0.5*arma::accu(arma::pow(can_bigU.slice(kind).col(iind) - bigU.slice(kind).col(iind) - 0.5*h_bigU*(nHessinv_U(iind,kind)*curgrad_bigU.slice(kind).col(iind)),2)/nHessinv_U(iind,kind))/h_bigU;
          //qcur += -0.5*arma::accu(arma::pow(bigU.slice(kind).col(iind) - can_bigU.slice(kind).col(iind) - 0.5*h_bigU*(nHessinv_U(iind,kind)*cangrad_bigU.slice(kind).col(iind)),2)/nHessinv_U(iind,kind))/h_bigU;
          for(int k2ind=0; k2ind<K; k2ind++){
            for(int i2ind=0; i2ind<I; i2ind++){
              rescur += nHessinv_U(kind*I+iind,k2ind*I+i2ind)*curgrad_bigU.slice(k2ind).col(i2ind);
              resprop += nHessinv_U_ch(kind*I+iind,k2ind*I+i2ind)*cangrad_bigU.slice(k2ind).col(i2ind);
            }
          }
          resp.slice(kind).col(iind) = can_bigU.slice(kind).col(iind) - bigU.slice(kind).col(iind) -0.5*h_bigU*rescur;
          resc.slice(kind).col(iind) = bigU.slice(kind).col(iind) - can_bigU.slice(kind).col(iind) -0.5*h_bigU*resprop;
        }
      }

      for(int kind=0; kind<K;kind++){
        for(int iind=0;iind<I;iind++){
          arma::vec rescur2 = arma::zeros(nbig), resprop2 = arma::zeros(nbig);
          for(int k2ind=0;k2ind<K;k2ind++){
            for(int i2ind=0; i2ind<I;i2ind++){
              rescur2 += nHessinv_U_chinv(kind*I+iind,k2ind*I+i2ind)*resc.slice(k2ind).col(i2ind);
              resprop2 += nHessinv_U_chinv(kind*I+iind,k2ind*I+i2ind)*resp.slice(k2ind).col(i2ind);
            }
          }

          qprop += arma::accu(arma::pow(resprop2,2))/(2*h_bigU);
          qcur += arma::accu(arma::pow(rescur2,2))/(2*h_bigU);
        }
      }

      la = canlik_bigU - curlik_bigU +qcur - qprop;
      if(log(arma::randu()) < la){
        bigU = can_bigU;
        bigres=canbigres;
        bigquad = canbigquad;
        bigUquad = canbigUquad;
        bigUgradpart = canbigUgradpart;
        bigUstar = canbigUstar;
        acc_U += 1.0/nthin;
      }




      //Individual Update
      //for(int kind=0;kind<K;kind++){
        //for(int iind=0; iind<I;iind++){
        //curlik_bigU = -0.5*arma::accu(bigUquad(iind,kind)) - 0.5*arma::sum(arma::sum(bigquad,0)/sig2.t());

        //arma::mat tempgradpartU = arma::zeros(1,nbig);
        //subres_row = arma::sum(bigres,0);
        //tempgradpartU = ((A.col(kind)/sig2).t())*subres_row;

        //curgrad_bigU= -bigUgradpart.slice(kind).col(iind) + tempgradpartU.t();

        //can_bigU = bigU;
        //canbigres=bigres;
        //canbigUquad=bigUquad.col(kind);
        //canbigUstar = bigUstar;
        //canbigUgradpart = bigUgradpart.slice(kind);
        //canbigquad = bigquad;
        //can_bigU.slice(kind).col(iind) = bigU.slice(kind).col(iind) + 0.5*h_bigU(iind,kind)*(nHessinv_U(kind)*curgrad_bigU) + sqrt(h_bigU(iind,kind))*(sqrt(nHessinv_U(kind))*arma::randn(nbig,1));

        //tempSigeig = arma::fft2(SigCube.slice(kind));
        //tempfftbigU = arma::fft2(arma::reshape(can_bigU.slice(kind).col(iind),Next,Mext));
        //tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
        //canbigUquad(iind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;
        //canbigUgradpart.col(iind) = arma::vectorise(arma::real(arma::ifft2((1/tempSigeig)%tempfftbigU)));

        //canbigUstar.slice(kind) = can_bigU.slice(kind)*arma::chol(arma::toeplitz(Cs.col(kind)));

        //for(int jind=0;jind<J;jind++){
          //canbigres.tube(iind,jind) += (bigUstar.col_as_mat(iind) - canbigUstar.col_as_mat(iind))*(A.row(jind).t());
        //}

        //canbigquad = arma::sum(arma::pow(canbigres,2),2);

        //canlik_bigU = -0.5*canbigUquad(iind) - 0.5*arma::sum(arma::sum(canbigquad,0)/sig2.t());

        //tempgradpartU = arma::zeros(1,nbig);
        //subres_row = arma::sum(canbigres,0);
        //tempgradpartU = ((A.col(kind)/sig2).t())*subres_row;
        //cangrad_bigU= -canbigUgradpart.col(iind) + tempgradpartU.t();

        //qprop = -0.5*arma::accu(arma::pow(can_bigU.slice(kind).col(iind) - bigU.slice(kind).col(iind) - 0.5*h_bigU(iind,kind)*(nHessinv_U(kind)*curgrad_bigU),2)/nHessinv_U(kind))/h_bigU(iind,kind);
        //qcur = -0.5*arma::accu(arma::pow(bigU.slice(kind).col(iind) - can_bigU.slice(kind).col(iind) - 0.5*h_bigU(iind,kind)*(nHessinv_U(kind)*cangrad_bigU),2)/nHessinv_U(kind))/h_bigU(iind,kind);

        //la = canlik_bigU - curlik_bigU +qcur - qprop;
        //if(log(arma::randu()) < la){
          //bigU.slice(kind).col(iind) = can_bigU.slice(kind).col(iind);
          //bigres=canbigres;
          //bigquad = canbigquad;
          //bigUquad(iind,kind) = canbigUquad(iind);
          //bigUgradpart.slice(kind).col(iind) = canbigUgradpart.col(iind);
          //bigUstar = canbigUstar;
          //acc_U(iind,kind) += 1.0/nthin;
        //}

        //}
      //}


      //std::cout << "Checkpoint 9" << std::endl;

      // done with all updates, go again
    }
    // store thinned samples only
    keep_mu.row(i) = mu.t();
    keep_A.row(i) = A;
    keep_rho.row(i) = rho.t();
    keep_sig2.row(i) = sig2.t();

    if(adapt!=0){
      if(i < adp){
      //double itr = (i+1)*nthin;
      //h_mu = h_mu0*exp(((acc_mu/itr) - 0.574)/sqrt(itr));
      //h_A = h_A0%arma::exp(((acc_A/itr) - 0.574)/sqrt(itr));
      //h_rho = h_rho0%arma::exp(((acc_rho/itr) - 0.23)/sqrt(itr));
      //h_sig2 = h_sig20*exp(((acc_sig2/itr) - 0.574)/sqrt(itr));
      //h_bigU = h_bigU0%arma::exp(((acc_U/itr) - 0.574)/sqrt(itr));

      if(acc_mu > 0.5){
        h_mu = h_mu*1.2;
      }
      if(acc_mu < 0.3){
        h_mu = h_mu*0.8;
      }

      // Block Updates
      if(acc_Aall > 0.5){
        h_Aall = h_Aall*1.2;
      }
      if(acc_Aall < 0.3){
        h_Aall = h_Aall*0.8;
      }

      // Individual Updates
      //for(int jind=0;jind<J;jind++){
        //if(acc_A(jind) > 0.5){
          //h_A(jind) = h_A(jind)*1.2;
        //}
        //if(acc_A(jind) < 0.3){
          //h_A(jind) = h_A(jind)*0.8;
        //}
      //}
      for(int kind=0;kind<K;kind++){
        if(acc_rho(kind) > 0.5){
          h_rho(kind) = h_rho(kind)*1.2;
        }
        if(acc_rho(kind) < 0.3){
          h_rho(kind) = h_rho(kind)*0.8;
        }
      }

      if(acc_sig2 > 0.5){
        h_sig2 = h_sig2*1.2;
      }
      if(acc_sig2 < 0.3){
        h_sig2 = h_sig2*0.8;
      }

      //Joint Update
      if(acc_U > 0.5){
        h_bigU = h_bigU*1.2;
      }
      if(acc_U < 0.3){
        h_bigU = h_bigU*0.8;
      }

      //Individual Update
      //for(int iind=0;iind<I;iind++){
        //for(int kind=0;kind<K;kind++){
          //if(acc_U(iind,kind) > 0.5){
            //h_bigU(iind,kind) = h_bigU(iind,kind)*1.2;
          //}
          //if(acc_U(iind,kind) < 0.3){
            //h_bigU(iind,kind) = h_bigU(iind,kind)*0.8;
          //}
        //}
      //}


      acc_mu = 0.0;
      //acc_A = arma::zeros(J);
      acc_Aall = 0.0;
      acc_rho = arma::zeros(K);
      acc_sig2 = 0.0;
      //acc_U = arma::zeros(I,K);
      acc_U = 0.0;

      }
    }

  }


  double eitr = nthin*(niters + nburn - adp);

  //std::cout << "Checkpoint 10" << std::endl;

  return List::create(Named("mu") = keep_mu,
                      Named("A") = keep_A,
                      Named("rho") = keep_rho,
                      Named("sigma2") = keep_sig2,
                      Named("acc_U") = acc_U/eitr,
                      Named("acc_mu") = acc_mu/eitr,
                      Named("acc_A") = acc_Aall/eitr,
                      Named("acc_rho") = acc_rho/eitr,
                      Named("acc_sig2") = acc_sig2/eitr,
                      Named("niters") = niters,
                      Named("nthin") = nthin,
                      Named("nburn") = nburn);


}
