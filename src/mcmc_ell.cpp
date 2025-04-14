// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

List Slicesto3D_ell(arma::cube Y, arma::mat S, double bma, int M, int N,
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
  arma::mat lam = arma::zeros(Next,Mext);
  arma::mat tempfftiid1 = arma::zeros(Next,Mext);
  arma::mat tempfftiid2 = arma::zeros(Next,Mext);
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
    // Cs.col(kind) = c1*((temp2+1)%arma::normcdf((temp2+1)*bma/temp1) + (temp2-1)%arma::normcdf((temp2-1)*bma/temp1) - 2*temp2%arma::normcdf(temp2*bma/temp1)) + (temp1/sqrt(2.0))*(arma::exp(-0.5*arma::pow((temp2+1)*bma/temp1,2)) + arma::exp(-0.5*arma::pow((temp2-1)*bma/temp1,2)) - 2*arma::exp(-0.5*arma::pow((temp2)*bma/temp1,2)));
    Cs.col(kind) = c1*((temp2+1)%arma::normcdf((temp2+1)*bma/temp1) + (temp2-1)%arma::normcdf((temp2-1)*bma/temp1) - 2*temp2%arma::normcdf(temp2*bma/temp1)) + (temp1/(2.0))*(arma::exp(-0.5*arma::pow((temp2+1)*bma/temp1,2)) + arma::exp(-0.5*arma::pow((temp2-1)*bma/temp1,2)) - 2*arma::exp(-0.5*arma::pow((temp2)*bma/temp1,2)));
  }


  // Generating initial latent process and computing related terms
  arma::mat tempnorm = arma::zeros(Next,Mext);
  arma::cx_mat tempSigeig(tempnorm,tempnorm);
  arma::cx_mat tempfftbigU(tempnorm,tempnorm);
  arma::cx_mat tempbigUquadhalf(tempnorm,tempnorm);
  arma::cx_mat tempfftiid(tempnorm,tempnorm);

  arma::cube bigU = arma::zeros(nbig,I,K);
  arma::cube bigUstar = bigU;
  arma::mat bigUquad = arma::zeros(I,K);
  arma::vec logdets = arma::zeros(K);
  arma::mat tempbigU = arma::zeros(Next,Mext);

  //std::cout << "Checkpoint 1.2" << std::endl;

  arma::mat dCmatchol = arma::zeros(I,K);

  for(int kind=0; kind<K; kind++){
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
    tempSigeig = arma::fft2(SigCube.slice(kind));
    lam = arma::real(tempSigeig);
    lam.clean(0.0);
    lam.replace(0,1e-12);
    logdets(kind) = arma::accu(arma::log(lam));
    lam = arma::sqrt(lam);
    for(int iii=0; iii<I; iii++){
      tempfftiid1 = lam%arma::randn(Next,Mext);
      tempfftiid2 = lam%arma::randn(Next,Mext);
      tempfftiid = arma::cx_mat(tempfftiid1,tempfftiid2);
      //tempfftiid = arma::fft2(arma::randn(Next,Mext));
      tempbigU = lam%arma::real(tempfftiid)*2*sqrt(arma::datum::pi);
      //tempbigU = arma::real(arma::ifft2(tempfftiid))/nbig;

      U.slice(kind).col(iii) = arma::vectorise(tempbigU.submat(0,0,N-1,M-1));
      tempfftbigU = arma::fft2(tempbigU);
      tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
      bigUquad(iii,kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;
      tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU2;
      diaginv(kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;
      bigU.slice(kind).col(iii) = arma::vectorise(tempbigU);
    }
    bigUstar.slice(kind) = bigU.slice(kind)*iCmatchol;
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

  // Candidates

  arma::vec can_mu = arma::zeros(J), can_sig2 = arma::zeros(J);
  // For Individual Updates of A
  // For Block Update of A
  arma::vec tempA = arma::zeros(K);

  arma::mat can_A = A;
  double can_rho = 0;
  arma::vec canCs = arma::zeros(I);

  //std::cout << "Checkpoint 2.1" << std::endl;

  arma::mat subres_row = arma::zeros(J,nbig), subres_col= arma::zeros(I,nbig), subres_slice= arma::zeros(I,J);


  //std::cout << "Checkpoint 2.2" << std::endl;

  double la, canlogdets;
  arma::mat canbigUquad = arma::zeros(I,K); // Joint U update
  arma::mat canSigCube = SigCube.slice(0), canbigquad=bigquad; // Joint U update
  arma::cube canbigres = bigres, can_bigU = bigU, canbigUstar = bigUstar;

  //std::cout << "Checkpoint 2.3" << std::endl;


  arma::vec nu_mu = arma::zeros(J);
  arma::mat nu_A = arma::zeros(J,K);
  arma::vec nu_sig2 = arma::zeros(J);
  arma::cube nu_bigU = arma::zeros(nbig,I,K);

  double theta_mu=0, theta_mu_min=0,theta_mu_max=0;
  arma::vec theta_A=arma::zeros(J), theta_A_min=arma::zeros(J),theta_A_max=arma::zeros(J);
  double theta_sig2=0, theta_sig2_min=0,theta_sig2_max=0;
  arma::mat theta_bigU=arma::zeros(I,K), theta_bigU_min=arma::zeros(I,K),theta_bigU_max=arma::zeros(I,K);

  // Storage
  arma::mat keep_mu = arma::zeros(niters,J);
  arma::cube keep_A = arma::zeros(niters,J,K);
  arma::mat keep_rho = arma::zeros(niters,K);
  arma::mat keep_sig2 = arma::zeros(niters,J);
  arma::vec keep_logL = arma::zeros(niters);
  double lu = 0.0;

  arma::vec acc_rho = arma::zeros(K);

  arma::vec h_rho0 = (5.6644)*arma::ones(K)/(0.5*rho*nbig*I + 1/sig2rho), h_rho = h_rho0;

  int adp = ceil(nburn/2.0);

  int maxCount = 100;

  double logL = 0.0;

  //std::cout << "Checkpoint 3" << std::endl;

  //std::cout << "Here we go" << std::endl;

  // GO!
  for(int i=0;i<niters;i++){
    for(int ii=0;ii<nthin;ii++){

      // Set thetas
      theta_mu = 2*arma::datum::pi*arma::randu(), theta_mu_min = theta_mu - 2*arma::datum::pi, theta_mu_max = theta_mu;
      theta_A = 2*arma::datum::pi*arma::randu(J), theta_A_min = theta_A - 2*arma::datum::pi, theta_A_max = theta_A;
      theta_sig2 = 2*arma::datum::pi*arma::randu(), theta_sig2_min = theta_sig2 - 2*arma::datum::pi, theta_sig2_max = theta_sig2;
      theta_bigU = 2*arma::datum::pi*arma::randu(I,K), theta_bigU_min = theta_bigU - 2*arma::datum::pi, theta_bigU_max = theta_bigU;

      // Update mu

      nu_mu = sqrt(sig2mu)*arma::randn(J);
      curlik_mu = -0.5*arma::sum(arma::pow(mu-mu0,2))/sig2mu - 0.5*arma::sum(arma::sum(bigquad,0)/sig2.t());

      can_mu = mu0 + (mu-mu0)*cos(theta_mu) + nu_mu*sin(theta_mu);

      canbigres = bigres;
      canbigquad = bigquad;
      for(int jind=0;jind<J;jind++){
        canbigres(arma::span::all,arma::span(jind),arma::span::all) += (mu(jind) - can_mu(jind))*bma;
      }
      canbigquad = arma::sum(arma::pow(canbigres,2),2);
      canlik_mu = -0.5*arma::sum(arma::pow(can_mu-mu0,2))/sig2mu - 0.5*arma::sum(arma::sum(canbigquad,0)/sig2.t());

      la = canlik_mu - curlik_mu;
      lu = log(arma::randu());

      int count_mu = 0;

      while(lu >= la){
        if(count_mu >= maxCount){break;}
        if(theta_mu < 0){theta_mu_min = theta_mu;}
        if(theta_mu >= 0){theta_mu_max = theta_mu;}
        theta_mu = theta_mu_min + (theta_mu_max - theta_mu_min)*arma::randu();

        can_mu = mu0 + (mu-mu0)*cos(theta_mu) + nu_mu*sin(theta_mu);

        canbigres = bigres;
        canbigquad = bigquad;
        for(int jind=0;jind<J;jind++){
          canbigres(arma::span::all,arma::span(jind),arma::span::all) += (mu(jind) - can_mu(jind))*bma;
        }
        canbigquad = arma::sum(arma::pow(canbigres,2),2);
        canlik_mu = -0.5*arma::sum(arma::pow(can_mu-mu0,2))/sig2mu - 0.5*arma::sum(arma::sum(canbigquad,0)/sig2.t());

        la = canlik_mu - curlik_mu;

        count_mu += 1;

      }
      if(count_mu < maxCount){
        mu = can_mu;
        bigres = canbigres;
        bigquad = canbigquad;
      }


      //std::cout << "Checkpoint 4" << std::endl;

      // Update A

      // curlik_A = 0.0;
      // canlik_A = 0.0;
      //
      // nu_A = sqrt(sig2A)*arma::randn(J,K);
      //
      // for(int jind=0;jind<J;jind++){
      //   curlik_A += -0.5*arma::sum(arma::pow(A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(bigquad.col(jind))/sig2(jind);
      // }
      //
      // can_A = A0 + A*cos(theta_A) + nu_A*sin(theta_A);
      //
      // canbigres = bigres;
      // canbigquad = bigquad;
      //
      // for(int jind=0;jind<J;jind++){
      //   arma::mat bigUa = arma::zeros(nbig,I);
      //   for(int kind=0;kind<K;kind++){
      //     bigUa += bigUstar.slice(kind)*(A(jind,kind) - can_A(jind,kind));
      //   }
      //
      //   for(int iind=0;iind<I;iind++){
      //     canbigres.tube(iind,jind) += bigUa.col(iind);
      //   }
      // }
      // canbigquad = arma::sum(arma::pow(canbigres,2),2);
      //
      // for(int jind=0;jind<J;jind++){
      //   canlik_A += -0.5*arma::sum(arma::pow(can_A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(canbigquad.col(jind))/sig2(jind);
      // }
      //
      // la = canlik_A - curlik_A;
      // lu = log(arma::randu());
      //
      // int count_A = 0;
      //
      // while(lu >= la){
      //   if(count_A >= maxCount){break;}
      //   if(theta_A < 0){theta_A_min = theta_A;}
      //   if(theta_A >= 0){theta_A_max = theta_A;}
      //
      //   theta_A = theta_A_min + (theta_A_max - theta_A_min)*arma::randu();
      //
      //   can_A = A0 + A*cos(theta_A) + nu_A*sin(theta_A);
      //
      //   canbigres = bigres;
      //   canbigquad = bigquad;
      //
      //   for(int jind=0;jind<J;jind++){
      //     arma::mat bigUa = arma::zeros(nbig,I);
      //     for(int kind=0;kind<K;kind++){
      //       bigUa += bigUstar.slice(kind)*(A(jind,kind) - can_A(jind,kind));
      //     }
      //
      //     for(int iind=0;iind<I;iind++){
      //       canbigres.tube(iind,jind) += bigUa.col(iind);
      //     }
      //   }
      //   canbigquad = arma::sum(arma::pow(canbigres,2),2);
      //
      //   for(int jind=0;jind<J;jind++){
      //     canlik_A += -0.5*arma::sum(arma::pow(can_A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(canbigquad.col(jind))/sig2(jind);
      //   }
      //
      //   la = canlik_A - curlik_A;
      //
      //   count_A += 1;
      //
      // }
      //
      // if(count_A < maxCount){
      //   A = can_A;
      //   bigres = canbigres;
      //   bigquad = canbigquad;
      // }

      for(int jind=0;jind<J;jind++){
        curlik_A = 0.0;
        canlik_A = 0.0;

        nu_A = sqrt(sig2A)*arma::randn(K);

        curlik_A = -0.5*arma::sum(arma::pow(A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(bigquad.col(jind))/sig2(jind);

        can_A = A;
        can_A.row(jind) = A0.row(jind) + (A.row(jind)-A0.row(jind))*cos(theta_A(jind)) + nu_A.t()*sin(theta_A(jind));

        canbigres = bigres;
        canbigquad = bigquad;

        arma::mat bigUa = arma::zeros(nbig,I);
        for(int kind=0;kind<K;kind++){
          bigUa += bigUstar.slice(kind)*(A(jind,kind) - can_A(jind,kind));
        }

        for(int iind=0;iind<I;iind++){
          canbigres.tube(iind,jind) += bigUa.col(iind);
        }

        canbigquad = arma::sum(arma::pow(canbigres,2),2);

        canlik_A = -0.5*arma::sum(arma::pow(can_A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(canbigquad.col(jind))/sig2(jind);

        la = canlik_A - curlik_A;
        lu = log(arma::randu());

        int count_A = 0;

        while(lu >= la){
          if(count_A >= maxCount){break;}
          if(theta_A(jind) < 0){theta_A_min(jind) = theta_A(jind);}
          if(theta_A(jind) >= 0){theta_A_max(jind) = theta_A(jind);}

          theta_A(jind) = theta_A_min(jind) + (theta_A_max(jind) - theta_A_min(jind))*arma::randu();

          can_A = A;
          can_A.row(jind) = A0.row(jind) + (A.row(jind)-A0.row(jind))*cos(theta_A(jind)) + nu_A.t()*sin(theta_A(jind));

          canbigres = bigres;
          canbigquad = bigquad;

          arma::mat bigUa = arma::zeros(nbig,I);
          for(int kind=0;kind<K;kind++){
            bigUa += bigUstar.slice(kind)*(A(jind,kind) - can_A(jind,kind));
          }

          for(int iind=0;iind<I;iind++){
            canbigres.tube(iind,jind) += bigUa.col(iind);
          }

          canbigquad = arma::sum(arma::pow(canbigres,2),2);

          canlik_A = -0.5*arma::sum(arma::pow(can_A.row(jind).t() - A0.row(jind).t(),2))/sig2A - 0.5*arma::sum(canbigquad.col(jind))/sig2(jind);

          la = canlik_A - curlik_A;
          count_A += 1;

        }

        if(count_A < maxCount){
          A = can_A;
          bigres = canbigres;
          bigquad = canbigquad;
        }

      }



      //std::cout << "Checkpoint 5" << std::endl;

      // Update rho
      // Likelihood now includes the residual part since Ustar has C_k in it. Need to update accordingly

      arma::mat can_dCmatchol = dCmatchol;

      for(int kind=0;kind<K;kind++){
        curlik_rho = -0.5*pow(log(rho(kind)/rho0(kind)),2)/sig2rho  - 0.5*I*logdets(kind) - 0.5*arma::accu(bigUquad.col(kind)) -0.5*arma::sum(arma::sum(bigquad,0)/sig2.t());


        can_rho = exp(log(rho(kind)) + sqrt(h_rho(kind))*arma::randn());

        canSigCube = arma::exp(-D2ext/(2*can_rho*can_rho));
        // canCs = c1*((temp2+1)%arma::normcdf((temp2+1)*bma/can_rho) + (temp2-1)%arma::normcdf((temp2-1)*bma/can_rho) - 2*temp2%arma::normcdf(temp2*bma/can_rho)) + (can_rho/sqrt(2.0))*(arma::exp(-0.5*arma::pow((temp2+1)*bma/can_rho,2)) + arma::exp(-0.5*arma::pow((temp2-1)*bma/can_rho,2)) - 2*arma::exp(-0.5*arma::pow((temp2)*bma/can_rho,2)));
        canCs = c1*((temp2+1)%arma::normcdf((temp2+1)*bma/can_rho) + (temp2-1)%arma::normcdf((temp2-1)*bma/can_rho) - 2*temp2%arma::normcdf(temp2*bma/can_rho)) + (can_rho/(2.0))*(arma::exp(-0.5*arma::pow((temp2+1)*bma/can_rho,2)) + arma::exp(-0.5*arma::pow((temp2-1)*bma/can_rho,2)) - 2*arma::exp(-0.5*arma::pow((temp2)*bma/can_rho,2)));


        tempSigeig = arma::fft2(canSigCube);
        canlogdets = arma::accu(arma::log(arma::real(tempSigeig)));


        canbigUquad.col(kind) = bigUquad.col(kind); // Joint U update
        canbigUstar = bigUstar;
        canbigres = bigres;
        canbigquad = bigquad;

        for(int iii=0; iii<I; iii++){
          tempfftbigU = arma::fft2(arma::reshape(bigU.slice(kind).col(iii),Next,Mext));
          tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
          canbigUquad(iii,kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig; //Joint U Update
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
        //la = canlik_rho - curlik_rho;
        if(log(arma::randu()) < la){
          rho(kind) = can_rho;
          Cs.col(kind) = canCs;
          SigCube.slice(kind) = canSigCube;
          logdets(kind) = canlogdets;
          bigUquad.col(kind) = canbigUquad.col(kind); // Joint U Update
          bigres=canbigres;
          bigquad=canbigquad;
          dCmatchol=can_dCmatchol;
          acc_rho(kind) += 1.0/nthin;
        }

      }

      //std::cout << "Checkpoint 6" << std::endl;

      // Update sig2

      nu_sig2 = sqrt(sigsig2)*arma::randn(J);

      curlik_sig2 = - arma::sum(arma::pow(arma::log(sig2)-musig2,2)/(2*sigsig2)) -0.5*nbig*I*arma::sum(arma::log(sig2)) -0.5*arma::sum(arma::sum(bigquad,0).t()/sig2);

      can_sig2 = arma::exp(musig2 + (arma::log(sig2)-musig2)*cos(theta_sig2) + nu_sig2*sin(theta_sig2));

      canlik_sig2 = - arma::sum(arma::pow(arma::log(can_sig2)-musig2,2)/(2*sigsig2)) -0.5*nbig*I*arma::sum(arma::log(can_sig2)) -0.5*arma::sum(arma::sum(bigquad,0).t()/can_sig2);

      la = canlik_sig2 - curlik_sig2;
      lu = log(arma::randu());

      int count_sig2 = 0;

      while(lu >= la){
        if(count_sig2 > maxCount){break;}
        if(theta_sig2 < 0){theta_sig2_min = theta_sig2;}
        if(theta_sig2 >= 0){theta_sig2_max = theta_sig2;}

        theta_sig2 = theta_sig2_min + (theta_sig2_max - theta_sig2_min)*arma::randu();

        can_sig2 = arma::exp(musig2 + (arma::log(sig2)-musig2)*cos(theta_sig2) + nu_sig2*sin(theta_sig2));

        canlik_sig2 = - arma::sum(arma::pow(arma::log(can_sig2)-musig2,2)/(2*sigsig2)) -0.5*nbig*I*arma::sum(arma::log(can_sig2)) -0.5*arma::sum(arma::sum(bigquad,0).t()/can_sig2);

        la = canlik_sig2 - curlik_sig2;
        count_sig2 += 1;
      }

      if(count_sig2 < maxCount){
        sig2= can_sig2;
      }

      //std::cout << "Checkpoint 7" << std::endl;

      //std::cout << "Checkpoint 8" << std::endl;

      // Update bigU and U

      //Joint update

      for(int kind=0;kind<K;kind++){
        tempSigeig = arma::fft2(SigCube.slice(kind));
        lam = arma::real(tempSigeig);
        lam.clean(0.0);
        lam.replace(0,1e-12);
        lam = arma::sqrt(lam);
        Cmatchol = arma::chol(arma::toeplitz(Cs.col(kind)));
        for(int iind=0;iind<I;iind++){
          tempfftiid1 = lam%arma::randn(Next,Mext);
          tempfftiid2 = lam%arma::randn(Next,Mext);
          tempfftiid = arma::cx_mat(tempfftiid1,tempfftiid2);
          //tempfftiid = arma::fft2(arma::randn(Next,Mext));
          nu_bigU.slice(kind).col(iind) = arma::vectorise(lam%arma::real(tempfftiid))*2*sqrt(arma::datum::pi);
          //nu_bigU.slice(kind).col(iind) = arma::vectorise(arma::real(arma::ifft2(tempfftiid)))/nbig;

          curlik_bigU = -0.5*arma::accu(bigUquad) - 0.5*arma::sum(arma::sum(bigquad,0)/sig2.t());

          can_bigU = bigU;
          can_bigU.slice(kind).col(iind) = bigU.slice(kind).col(iind)*cos(theta_bigU(iind,kind)) + nu_bigU.slice(kind).col(iind)*sin(theta_bigU(iind,kind));

          canbigres=bigres;
          canbigUquad=bigUquad;
          canbigUstar = bigUstar;
          canbigquad = bigquad;


          tempfftbigU = arma::fft2(arma::reshape(can_bigU.slice(kind).col(iind),Next,Mext));
          tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
          canbigUquad(iind,kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;

          canbigUstar.slice(kind) = can_bigU.slice(kind)*Cmatchol;

          for(int i2ind=0;i2ind<I;i2ind++){
            for(int j2ind=0;j2ind<J;j2ind++){
              canbigres.tube(i2ind,j2ind) += (bigUstar.col_as_mat(i2ind) - canbigUstar.col_as_mat(i2ind))*(A.row(j2ind).t());
            }
          }

          canbigquad = arma::sum(arma::pow(canbigres,2),2);

          canlik_bigU = -0.5*arma::accu(canbigUquad) - 0.5*arma::sum(arma::sum(canbigquad,0)/sig2.t());

          la = canlik_bigU - curlik_bigU;
          lu = log(arma::randu());

          int count_bigU = 0;

          while(lu >= la){
            if(count_bigU >= maxCount){break;}
            if(theta_bigU(iind,kind) < 0){theta_bigU_min(iind,kind) = theta_bigU(iind,kind);}
            if(theta_bigU(iind,kind) > 0){theta_bigU_max(iind,kind) = theta_bigU(iind,kind);}

            theta_bigU(iind,kind) = theta_bigU_min(iind,kind) + (theta_bigU_max(iind,kind) - theta_bigU_min(iind,kind))*arma::randu();

            can_bigU = bigU;
            can_bigU.slice(kind).col(iind) = bigU.slice(kind).col(iind)*cos(theta_bigU(iind,kind)) + nu_bigU.slice(kind).col(iind)*sin(theta_bigU(iind,kind));

            canbigres=bigres;
            canbigUquad=bigUquad;
            canbigUstar = bigUstar;
            canbigquad = bigquad;


            tempfftbigU = arma::fft2(arma::reshape(can_bigU.slice(kind).col(iind),Next,Mext));
            tempbigUquadhalf = (1/arma::sqrt(tempSigeig))%tempfftbigU;
            canbigUquad(iind,kind) = real(arma::cdot(tempbigUquadhalf,tempbigUquadhalf))/nbig;

            canbigUstar.slice(kind) = can_bigU.slice(kind)*Cmatchol;

            for(int i2ind=0;i2ind<I;i2ind++){
              for(int j2ind=0;j2ind<J;j2ind++){
                canbigres.tube(i2ind,j2ind) += (bigUstar.col_as_mat(i2ind) - canbigUstar.col_as_mat(i2ind))*(A.row(j2ind).t());
              }
            }

            canbigquad = arma::sum(arma::pow(canbigres,2),2);

            canlik_bigU = -0.5*arma::accu(canbigUquad) - 0.5*arma::sum(arma::sum(canbigquad,0)/sig2.t());

            la = canlik_bigU - curlik_bigU;
            count_bigU += 1;
          }
          if(count_bigU < maxCount){
            bigU = can_bigU;
            bigres = canbigres;
            bigquad = canbigquad;
            bigUstar= canbigUstar;
            bigUquad = canbigUquad;
          }

        }
      }

      logL = -0.5*nbig*I*arma::sum(arma::log(sig2)) - 0.5*arma::sum(arma::sum(bigquad,0)/sig2.t()) - 0.5*I*arma::accu(logdets) - 0.5*arma::accu(bigUquad) - 0.5*arma::sum(arma::pow(mu-mu0,2))/sig2mu -0.5*arma::accu(arma::pow(A - A0,2))/sig2A - 0.5*arma::sum(pow(log(rho/rho0),2))/sig2rho - arma::sum(arma::pow(arma::log(sig2)-musig2,2)/(2*sigsig2));


      //std::cout << "Checkpoint 9" << std::endl;

      // done with all updates, go again
    }
    // store thinned samples only
    keep_mu.row(i) = mu.t();
    keep_A.row(i) = A;
    keep_rho.row(i) = rho.t();
    keep_sig2.row(i) = sig2.t();
    keep_logL(i) = logL;

    if(adapt!=0){
      if(i < adp){

        for(int kind=0;kind<K;kind++){
          if(acc_rho(kind) > 0.5){
            h_rho(kind) = h_rho(kind)*1.2;
          }
          if(acc_rho(kind) < 0.3){
            h_rho(kind) = h_rho(kind)*0.8;
          }
        }

        acc_rho = arma::zeros(K);
      }
    }

  }


  double eitr = nthin*(niters + nburn - adp);

  //std::cout << "Checkpoint 10" << std::endl;

  return List::create(Named("mu") = keep_mu,
                      Named("A") = keep_A,
                      Named("rho") = keep_rho,
                      Named("sigma2") = keep_sig2,
                      Named("logPosterior") = keep_logL,
                      Named("acc_rho") = acc_rho/eitr,
                      Named("niters") = niters,
                      Named("nthin") = nthin,
                      Named("nburn") = nburn);


}
