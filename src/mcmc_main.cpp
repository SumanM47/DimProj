// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

List Slicesto3D(arma::cube Y, arma::mat S, double bma, int M, int N,
                arma::vec mu_init, arma::vec rho_init, arma::mat A_init, arma::mat sig2_init,
                int niters, int nburn, int nthin,
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
  arma::mat sig2 = sig2_init;

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
    bigUstari = bigUstar.col(iind);
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
  arma::vec curgrad_aj = arma::zeros(K), cangrad_aj =  arma::zeros(K), tempA = arma::zeros(K);
  arma::mat can_A = A;
  double can_rho = 0;
  arma::vec canCs = arma::zeros(I);
  arma::mat curgrad_sig2 = arma::zeros(I,J), cangrad_sig2 =  arma::zeros(I,J), can_sig2 = sig2;

  //std::cout << "Checkpoint 2.1" << std::endl;

  arma::mat curgrad_bigU = arma::zeros(Next,Mext), cangrad_bigU = arma::zeros(Next,Mext);
  //arma::mat subres_row = arma::zeros(J,n), subres_col= arma::zeros(I,n), subres_slice= arma::zeros(I,J);
  arma::mat subres_row = arma::zeros(J,nbig), subres_col= arma::zeros(I,nbig), subres_slice= arma::zeros(I,J);

  // Initiating and calculating the inverse Hessian matrices
  // diagonal Hessian matrices turned to vectors
  // block diagonal matrices made into cube to use each slice independently
  //arma::vec nHessinv_mu = arma::ones(J), nHessinv_U= arma::ones(K);
  arma::vec nHessinv_mu = arma::ones(J);
  arma::mat nHessinv_A = arma::ones(K,J), nHessinv_sig2 = arma::ones(I,J);
  //arma::mat nHessinv_sig2 = arma::ones(I,J);
  arma::cube nHessinv_U = arma::ones(nbig,I,K);
  //arma::cube arma::mat nHessinv_A = arma::ones(K,K,J), nHessinv_U = arma::ones(Next,Mext,K);

  //std::cout << "Checkpoint 2.2" << std::endl;

  double qprop, qcur, la, canlogdets;
  arma::vec canbigUquad = arma::zeros(I);
  //arma::mat canquad = quad, canSigCube = SigCube.slice(0), canbigUgradpart = bigUgradpart.slice(0), canU = U;
  arma::mat canSigCube = SigCube.slice(0), canbigUgradpart = bigUgradpart.slice(0), canbigquad=bigquad;
  //arma::cube canres = res, can_bigU = bigU;
  arma::cube canbigres = bigres, can_bigU = bigU, canbigUstar = bigUstar;

  //std::cout << "Checkpoint 2.3" << std::endl;

  for(int jind=0;jind<J;jind++){
    //nHessinv_mu(jind) = 1/((1/sig2mu) + n*bma*bma*arma::sum(1/(sig2.col(jind))));
    nHessinv_mu(jind) = 1/((1/sig2mu) + nbig*bma*bma*arma::sum(1/(sig2.col(jind))));
  }

  //nHessinv_sig2 = 1/(1/sigsig2 + 0.5*quad/sig2);
  nHessinv_sig2 = 1/(1/sigsig2 + 0.5*bigquad/sig2);

  //std::cout << "Checkpoint 2.4" << std::endl;

  //arma::vec dUtU = arma::sum(arma::pow(U,2),0).t();
  //arma::vec dUtU = arma::zeros(K);
  //arma::mat dUtU = arma::zeros(K,K);
  arma::mat dUtU = arma::zeros(K,I);
  arma::vec tempHessA = arma::zeros(K);
  //for(int kind=0;kind<K;kind++){
    //dUtU(kind) = arma::accu(arma::pow(bigUstar.slice(kind),2));
  //}
  for(int iind=0; iind < I; iind++){
    dUtU.col(iind) = arma::trans(arma::sum(arma::pow(bigUstar.col_as_mat(iind),2),0));
  }
  for(int jind=0;jind<J;jind++){
    tempHessA = arma::zeros(K,1);
      for(int iind=0;iind<I;iind++){
        tempHessA += dUtU.col(iind)/sig2(iind,jind);
      }
      nHessinv_A.col(jind) = 1/(1/sig2A + tempHessA);
    }

  //std::cout << "Checkpoint 2.5" << std::endl;

  for(int kind=0;kind<K;kind++){
    //double tempscalebigU = 0;
    arma::mat tempnHessU = arma::zeros(nbig,I);
    tempnHessU = diaginv(kind)*arma::ones(nbig,I);
    for(int iind=0; iind<I; iind++){
      tempnHessU.col(iind) += arma::sum(arma::pow(A.col(kind),2)/sig2.row(iind).t());
    }
    nHessinv_U.slice(kind) = 1/tempnHessU;
    //nHessinv_U.slice(kind).submat(0,0,N-1,M-1) = (1/(diaginv(kind)/cs(kind) + tempscalebigU))*arma::ones(N,M);
    //nHessinv_U.slice(kind) = (1/(diaginv(kind)/cs(kind) + tempscalebigU))*arma::ones(Next,Mext);
  }
  //nHessinv_U = arma::ones(K);

  //std::cout << "done with Hessians" << std::endl;
  // Storage
  arma::mat keep_mu = arma::zeros(niters,J);
  arma::cube keep_A = arma::zeros(niters,J,K);
  arma::mat keep_rho = arma::zeros(niters,K);
  arma::cube keep_sig2 = arma::zeros(niters,I,J);
  arma::mat acc_U = arma::zeros(I,K);

  // proposal variances
  double h_mu = 2.7225/pow(J,0.33);
  arma::vec h_A = (5.6644/(K))*arma::ones(J);
  arma::vec h_rho = 0.4*(5.6644)*arma::ones(K);
  double h_sig2 = (5.6644/(I*J));
  arma::vec h_bigU = (2.7225/(K*pow(nbig*I,0.33)))*arma::ones(K);

  //std::cout << "Checkpoint 3" << std::endl;

  // GO!
  for(int i=0;i<niters;i++){
    for(int ii=0;ii<nthin;ii++){

      // Update mu

      //curlik_mu = -0.5*arma::sum(arma::pow(mu-mu0,2))/sig2mu - 0.5*arma::accu(quad/sig2);
      curlik_mu = -0.5*arma::sum(arma::pow(mu-mu0,2))/sig2mu - 0.5*arma::accu(bigquad/sig2);
      //subres_slice = arma::sum(res,2);
      subres_slice = arma::sum(bigres,2);
      curgrad_mu = -(mu - mu0)/sig2mu + bma*arma::sum(subres_slice/sig2,0).t();

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
      canlik_mu = -0.5*arma::sum(arma::pow(can_mu-mu0,2))/sig2mu - 0.5*arma::accu(canbigquad/sig2);
      //subres_slice = arma::sum(canres,2);
      subres_slice = arma::sum(canbigres,2);
      cangrad_mu = -(can_mu - mu0)/sig2mu + bma*arma::sum(subres_slice/sig2,0).t();

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
      }

      //std::cout << "Checkpoint 4" << std::endl;

      // Update A

      for(int jind=0;jind<J;jind++){
      //curlik_A = -0.5*arma::sum(arma::pow(A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(quad.col(jind)/sig2.col(jind));
      curlik_A = -0.5*arma::sum(arma::pow(A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(bigquad.col(jind)/sig2.col(jind));

      //qprop = 0;
      //qcur = 0;

        //subres_col = res(arma::span::all,arma::span(jind),arma::span::all);
        //curgrad_aj = -(A.row(jind).t() - A0.row(jind).t())/sig2A + arma::sum((U.t()*subres_col.t())*arma::diagmat(1/sig2.col(jind)),1);
        for(int kind=0;kind<K;kind++){
          double tempdA = 0.0;
          for(int iind=0;iind<I;iind++){
            arma::vec tempU = bigUstar.slice(kind).col(iind);
            arma::vec nbigtemp = bigres.tube(iind,jind);
            tempdA += arma::sum(tempU%nbigtemp)/sig2(iind,jind);
          }
          curgrad_aj(kind) = -(A(jind,kind) - A0(jind,kind))/sig2A + tempdA;
        }

        //std::cout<<"Checkpoint 4.1"<< std::endl;

        can_A = A;
        //canres = res;
        //canquad = quad;
        canbigres = bigres;
        canbigquad = bigquad;

        tempA = A.row(jind).t() + 0.5*h_A(jind)*(nHessinv_A.col(jind)%curgrad_aj) + sqrt(h_A(jind))*(arma::sqrt(nHessinv_A.col(jind))%arma::randn(K));
        can_A.row(jind) = tempA.t();

        //tempres = mu(jind)*bma + U*(can_A.row(jind).t());
        //tempres = mu(jind)*bma + U*(tempA.t());
        for(int iind=0;iind<I;iind++){
          //tempres2 = Y(arma::span(iind),arma::span(jind),arma::span::all);
          //tempres3 = tempres2 - tempres;
          //canres.tube(iind,jind) = tempres3;
          //canquad(iind,jind) = arma::sum(tempres3%tempres3);
        }

        arma::mat bigUa = arma::zeros(nbig,I);
        for(int kind=0;kind<K;kind++){
          bigUa += bigUstar.slice(kind)*(A(jind,kind) - can_A(jind,kind));
        }
        //arma::vec bigUav = arma::vectorise(bigUa);
        for(int iind=0;iind<I;iind++){
          canbigres.tube(iind,jind) += bigUa.col(iind);
        }
        canbigquad = arma::sum(arma::pow(canbigres,2),2);

        //subres_col = canres(arma::span::all,arma::span(jind),arma::span::all);
        //cangrad_aj = -(can_A.row(jind).t() - A0.row(jind).t())/sig2A + arma::sum((U.t()*subres_col.t())*arma::diagmat(1/sig2.col(jind)),1);
        for(int kind=0;kind<K;kind++){
          //arma::vec tempU = arma::vectorise(bigU.slice(kind));
          double tempdA = 0.0;
          for(int iind=0;iind<I;iind++){
            arma::vec tempU = bigUstar.slice(kind).col(iind);
            arma::vec nbigtemp = canbigres.tube(iind,jind);
            tempdA += arma::sum(tempU%nbigtemp)/sig2(iind,jind);
          }
          cangrad_aj(kind) = -(can_A(jind,kind) - A0(jind,kind))/sig2A + tempdA;
        }

        //std::cout<<"Checkpoint 4.2"<< std::endl;

        qprop = -0.5*arma::sum(arma::pow((can_A.row(jind).t() - A.row(jind).t() - 0.5*h_A(jind)*(nHessinv_A.col(jind)%curgrad_aj)),2)/nHessinv_A.col(jind))/h_A(jind);
        qcur = -0.5*arma::sum(arma::pow((A.row(jind).t() - can_A.row(jind).t() - 0.5*h_A(jind)*(nHessinv_A.col(jind)%cangrad_aj)),2)/nHessinv_A.col(jind))/h_A(jind);


      //canlik_A = -0.5*arma::sum(arma::pow(can_A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(canquad.col(jind)/sig2.col(jind));
      canlik_A = -0.5*arma::sum(arma::pow(can_A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(canbigquad.col(jind)/sig2.col(jind));

      la = canlik_A + qcur - curlik_A - qprop;
      if(log(arma::randu()) < la){
        //arma::mat bigUa = arma::zeros(Next,Mext);
        //for(int kind=0;kind<K;kind++){
          //bigUa += bigU.slice(kind)*(A(jind,kind) - can_A(jind,kind));
        //}
        //for(int iind=0;iind<I;iind++){
          //bigres.tube(iind,jind) += arma::vectorise(bigUa);
        //}
        A = can_A;
        //res = canres;
        //quad = canquad;
        bigres = canbigres;
        bigquad = canbigquad;
      }
      //bigres.slices(uind) = res;
      }

      //std::cout << "size of A:" << arma::size(A) << std::endl;
      //std::cout << "Checkpoint 5" << std::endl;

      // Update rho
      // Likelihood now includes the residual part since Ustar has C_k in it. Need to update accordingly

      for(int kind=0;kind<K;kind++){
        curlik_rho = -0.5*pow(log(rho(kind)/rho0(kind)),2)/sig2rho  - 0.5*I*logdets(kind) - 0.5*arma::accu(bigUquad.col(kind)) -0.5*arma::accu(bigquad/sig2);


        can_rho = exp(log(rho(kind)) + sqrt(h_rho(kind))*arma::randn());

        canSigCube = arma::exp(-D2ext/(2*can_rho*can_rho));
        canCs = c1*((temp2+1)%arma::normcdf((temp2+1)*bma/can_rho) + (temp2-1)%arma::normcdf((temp2-1)*bma/can_rho) - 2*temp2%arma::normcdf(temp2*bma/can_rho)) + (can_rho/sqrt(2.0))*(arma::exp(-0.5*arma::pow((temp2+1)*bma/can_rho,2)) + arma::exp(-0.5*arma::pow((temp2-1)*bma/can_rho,2)) - 2*arma::exp(-0.5*arma::pow((temp2)*bma/can_rho,2)));


        tempSigeig = arma::fft2(canSigCube);
        canlogdets = arma::accu(arma::log(arma::real(tempSigeig)));
        for(int iii=0; iii<I; iii++){
        tempfftbigU = arma::fft2(arma::reshape(bigU.slice(kind).col(iii),Next,Mext));
        tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
        canbigUquad(iii) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;
        canbigUgradpart.col(iii) = arma::vectorise(arma::real(arma::ifft2((1/tempSigeig)%tempfftbigU)));
        }
        canbigUstar.slice(kind) = bigU.slice(kind)*arma::chol(arma::toeplitz(canCs));
        for(int iind=0; iind<I; iind++){
          for(int jind=0; jind<J; jind++){
            canbigres.tube(iind,jind) += (bigUstar.col_as_mat(iind) - canbigUstar.col_as_mat(iind))*(A.row(jind).t());
          }
        }
        canbigquad = arma::sum(arma::pow(canbigres,2),2);

        canlik_rho = -0.5*pow(log(can_rho/rho0(kind)),2)/sig2rho - 0.5*I*canlogdets - 0.5*arma::accu(canbigUquad) -0.5*arma::accu(canbigquad/sig2);

        la = canlik_rho - curlik_rho;
        if(log(arma::randu()) < la){
          rho(kind) = can_rho;
          Cs.col(kind) = canCs;
          SigCube.slice(kind) = canSigCube;
          logdets(kind) = canlogdets;
          bigUquad.col(kind) = canbigUquad;
          bigUgradpart.slice(kind) = canbigUgradpart;
          bigres=canbigres;
          bigquad=canbigquad;
        }

      }

      //std::cout << "Checkpoint 6" << std::endl;

      // Update sig2

      //curlik_sig2 = - arma::accu(arma::pow(arma::log(sig2)-musig2,2)/(2*sigsig2)) -0.5*n*arma::accu(arma::log(sig2)) -0.5*arma::accu(quad/sig2);
      curlik_sig2 = - arma::accu(arma::pow(arma::log(sig2)-musig2,2)/(2*sigsig2)) -0.5*nbig*arma::accu(arma::log(sig2)) -0.5*arma::accu(bigquad/sig2);
      //curgrad_sig2 = -0.5*n + 0.5*(quad/sig2) - (arma::log(sig2) - musig2)/sigsig2;
      curgrad_sig2 = -0.5*nbig + 0.5*(bigquad/sig2) - (arma::log(sig2) - musig2)/sigsig2;

      can_sig2 = arma::exp(arma::log(sig2) + 0.5*h_sig2*(nHessinv_sig2%curgrad_sig2) + sqrt(h_sig2)*(arma::sqrt(nHessinv_sig2)%arma::randn(I,J)));


      //canlik_sig2 = - arma::accu(arma::pow(arma::log(can_sig2)-musig2,2)/(2*sigsig2)) -0.5*n*arma::accu(arma::log(can_sig2)) -0.5*arma::accu(quad/can_sig2);
      canlik_sig2 = - arma::accu(arma::pow(arma::log(can_sig2)-musig2,2)/(2*sigsig2)) -0.5*nbig*arma::accu(arma::log(can_sig2)) -0.5*arma::accu(bigquad/can_sig2);
      //cangrad_sig2 = -0.5*n + 0.5*(quad/can_sig2) - (arma::log(can_sig2) - musig2)/sigsig2;
      cangrad_sig2 = -0.5*nbig + 0.5*(bigquad/can_sig2) - (arma::log(can_sig2) - musig2)/sigsig2;

      qprop = -0.5*arma::accu(arma::pow(arma::log(can_sig2) - arma::log(sig2) - 0.5*h_sig2*(nHessinv_sig2%curgrad_sig2),2)/nHessinv_sig2)/h_sig2;
      qcur = -0.5*arma::accu(arma::pow(arma::log(sig2) - arma::log(can_sig2) - 0.5*h_sig2*(nHessinv_sig2%cangrad_sig2),2)/nHessinv_sig2)/h_sig2;

      la = canlik_sig2 + qcur - curlik_sig2 - qprop;
      if(log(arma::randu()) < la){
        sig2 = can_sig2;
      }

      //std::cout << "Checkpoint 7" << std::endl;

      //std::cout << "Checkpoint 8" << std::endl;

      // Update bigU and U

      for(int kind=0;kind<K;kind++){
        for(int iind=0; iind<I;iind++){
        //curlik_bigU = -0.5*bigUquad(kind)/cs(kind) - 0.5*arma::accu(quad/sig2);
        curlik_bigU = -0.5*arma::accu(bigUquad(iind,kind)) - 0.5*arma::accu(bigquad.row(iind)/sig2.row(iind));

        arma::mat tempgradpartU = arma::zeros(1,n);
        //arma::mat tempgradpartbigUnew = arma::zeros(Next,Mext);
          //subres_row = res(arma::span(iind),arma::span::all,arma::span::all);
          subres_row = bigres(arma::span(iind),arma::span::all,arma::span::all);
          tempgradpartU = ((A.col(kind)/sig2.row(iind).t()).t())*subres_row;

        //tempgradpartbigUnew.submat(0,0,N-1,M-1) = arma::reshape(tempgradpartU.t(),N,M);
        //tempgradpartbigUnew = arma::reshape(tempgradpartU.t(),Next,Mext);
        curgrad_bigU= -bigUgradpart.slice(kind).col(iind) + tempgradpartU.t();

        //canU = U;
        can_bigU = bigU;
        canbigres=bigres;
        canbigquad=bigquad;
        canbigUstar = bigUstar;
        //canres=res;
        //canquad=quad;

        //std::cout << "Checkpoint 8.1" << std::endl;

        //can_bigU.slice(kind) = bigU.slice(kind) + 0.5*h_bigU(kind)*nHessinv_U(kind)*curgrad_bigU + sqrt(h_bigU(kind)*nHessinv_U(kind))*arma::randn(Next,Mext);
        can_bigU.slice(kind).col(iind) = bigU.slice(kind).col(iind) + 0.5*h_bigU(kind)*(nHessinv_U.slice(kind).col(iind)%curgrad_bigU) + sqrt(h_bigU(kind))*(arma::sqrt(nHessinv_U.slice(kind).col(iind))%arma::randn(nbig,1));

        tempSigeig = arma::fft2(SigCube.slice(kind));
        tempfftbigU = arma::fft2(arma::reshape(can_bigU.slice(kind).col(iind),Next,Mext));
        tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
        canbigUquad(iind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;
        canbigUgradpart.col(iind) = arma::vectorise(arma::real(arma::ifft2((1/tempSigeig)%tempfftbigU)));

        canbigUstar.slice(kind) = can_bigU.slice(kind)*arma::chol(arma::toeplitz(Cs.col(kind)));

        //canU.col(kind) = arma::vectorise(can_bigU.slice(kind).submat(0,0,N-1,M-1));

        //std::cout << "Checkpoint 8.2" << std::endl;
        for(int jind=0;jind<J;jind++){
          canbigres.tube(iind,jind) += (bigUstar.col_as_mat(iind) - canbigUstar.col_as_mat(iind))*(A.row(jind).t());
        }

        //for(int jind=0;jind<J;jind++){
          //tempres = mu(jind)*bma + canU*(A.row(jind).t());
          //arma::mat bigUa = arma::zeros(Next,Mext);
          //for(int kkind=0;kkind<K;kkind++){
            //bigUa += (bigU.slice(kkind) - can_bigU.slice(kkind))*A(jind,kkind);
          //}
          //for(int iind=0;iind<I;iind++){
            //tempres2 = Y(arma::span(iind),arma::span(jind),arma::span::all);
            //tempres3 = tempres2 - tempres;
            //canres.tube(iind,jind) = tempres3;
            //canquad(iind,jind) = arma::sum(tempres3%tempres3);
            //canbigres.tube(iind,jind) += arma::vectorise(bigUa);
          //}
        //}

        canbigquad = arma::sum(arma::pow(canbigres,2),2);

        //canbigres.slices(uind) = canres;

        //canlik_bigU = -0.5*canbigUquad/cs(kind) - 0.5*arma::accu(canquad/sig2);
        canlik_bigU = -0.5*canbigUquad(iind) - 0.5*arma::accu(canbigquad.row(iind)/sig2.row(iind));

        //std::cout << "Checkpoint 8.3" << std::endl;

        //tempgradpartU = arma::zeros(1,n);
        tempgradpartU = arma::zeros(1,nbig);
        //tempgradpartbigUnew = arma::zeros(Next,Mext);
        //for(int iind=0;iind<I;iind++){
          //subres_row = canres(arma::span(iind),arma::span::all,arma::span::all);
          subres_row = canbigres(arma::span(iind),arma::span::all,arma::span::all);
          tempgradpartU = ((A.col(kind)/sig2.row(iind).t()).t())*subres_row;
        //}
        //tempgradpartbigUnew.submat(0,0,N-1,M-1) = arma::reshape(tempgradpartU.t(),N,M);
        //tempgradpartbigUnew = arma::reshape(tempgradpartU.t(),Next,Mext);
        cangrad_bigU= -canbigUgradpart.col(iind) + tempgradpartU.t();

        //qprop = -0.5*arma::accu(arma::pow(can_bigU.slice(kind) - bigU.slice(kind) - 0.5*h_bigU(kind)*(nHessinv_U(kind)*curgrad_bigU),2)/nHessinv_U(kind))/h_bigU(kind);
        //qcur = -0.5*arma::accu(arma::pow(bigU.slice(kind) - can_bigU.slice(kind) - 0.5*h_bigU(kind)*(nHessinv_U(kind)*cangrad_bigU),2)/nHessinv_U(kind))/h_bigU(kind);

        qprop = -0.5*arma::accu(arma::pow(can_bigU.slice(kind).col(iind) - bigU.slice(kind).col(iind) - 0.5*h_bigU(kind)*(nHessinv_U.slice(kind).col(iind)%curgrad_bigU),2)/nHessinv_U.slice(kind).col(iind))/h_bigU(kind);
        qcur = -0.5*arma::accu(arma::pow(bigU.slice(kind).col(iind) - can_bigU.slice(kind).col(iind) - 0.5*h_bigU(kind)*(nHessinv_U.slice(kind).col(iind)%cangrad_bigU),2)/nHessinv_U.slice(kind).col(iind))/h_bigU(kind);



        //std::cout << "Checkpoint 8.4" << std::endl;

        la = canlik_bigU - curlik_bigU +qcur - qprop;
        if(log(arma::randu()) < la){
          bigU.slice(kind).col(iind) = can_bigU.slice(kind).col(iind);
          //U.col(kind) = canU.col(kind);
          //res = canres;
          bigres=canbigres;
          //quad = canquad;
          bigquad = canbigquad;
          bigUquad(iind,kind) = canbigUquad(iind);
          bigUgradpart.slice(kind).col(iind) = canbigUgradpart.col(iind);
          acc_U(iind,kind) += 1.0/(niters*nthin);
        }
        //bigres.slices(uind) = res;

        }
      }


      //std::cout << "Checkpoint 9" << std::endl;

      // done with all updates, go again
    }
    // store thinned samples only
    keep_mu.row(i) = mu.t();
    keep_A.row(i) = A;
    keep_rho.row(i) = rho.t();
    keep_sig2.row(i) = sig2;
  }


  //std::cout << "Checkpoint 10" << std::endl;

  return List::create(Named("mu") = keep_mu,
                      Named("A") = keep_A,
                      Named("rho") = keep_rho,
                      Named("sigma2") = keep_sig2,
                      Named("acc_U") = acc_U,
                      Named("niters") = niters,
                      Named("nthin") = nthin,
                      Named("nburn") = nburn);


}
