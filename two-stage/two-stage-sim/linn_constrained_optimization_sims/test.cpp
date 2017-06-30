  /* @file test.cpp
   @brief Code to test the constrained opt simulation.  
 */
   
  // some of these might not be necessary 
  #include <Rcpp.h>
  #include <math.h>
  #include <ctime>
  #include <boost/numeric/ublas/vector.hpp>
  #include <boost/numeric/ublas/vector_proxy.hpp>
  #include <boost/numeric/ublas/matrix.hpp>
  #include <gsl/gsl_randist.h>
  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_blas.h>
   
  #include "genModel.hpp"
  #include "genDist.hpp"
  #include "estPars.hpp"
  #include "estFnls.hpp"
  #include "gslFnls.hpp"
  #include "learn.hpp"
   
   
   RcppExport SEXP genEx (SEXP n, SEXP m, SEXP mu, SEXP beta10, 
    SEXP beta11, SEXP gamma0, SEXP gamma1, 
    SEXP betaY20, SEXP betaY21, SEXP betaZ20, SEXP
    betaZ21, SEXP sqrtSigma, SEXP J, SEXP N, SEXP kappa, SEXP alpha,
    SEXP ss){

    using namespace std;
    namespace ublas = boost::numeric::ublas;
    // _ss is the seed to use
    Rcpp::NumericVector Rss (ss);
    unsigned long int _ss = static_cast<unsigned long int>(Rss (0)); 
    // _n is the training set size
    Rcpp::NumericVector Rn (n);
    int _n = static_cast<int>(Rn (0));
    // _m is the test set size for evaluation of policies
    Rcpp::NumericVector Rm (m);
    int _m = static_cast<int>(Rm (0));
    Rcpp::NumericVector Rkappa (kappa);
    double _kappa = Rkappa (0);
    Rcpp::NumericVector Ralpha (alpha);
    double _alpha = Ralpha (0);
    Rcpp::NumericVector RN (N);
    int _N = static_cast<int>(RN (0));
    Rcpp::NumericVector RJ (J);
    int _J = static_cast<int>(RJ (0));
    Rcpp::NumericVector Rmu (mu);
    Rcpp::NumericVector Rbeta10 (beta10);
    Rcpp::NumericVector Rbeta11 (beta11);
    Rcpp::NumericVector Rgamma0 (gamma0);
    Rcpp::NumericVector Rgamma1 (gamma1);
    Rcpp::NumericVector RbetaY20 (betaY20);
    Rcpp::NumericVector RbetaY21 (betaY21);
    Rcpp::NumericVector RbetaZ20 (betaZ20);
    Rcpp::NumericVector RbetaZ21 (betaZ21);
    Rcpp::NumericMatrix RsqrtSigma (sqrtSigma);
    
    ublas::vector<double> _mu (1);
    _mu (0) = Rmu (0);
    ublas::vector<double> _beta10 (2);
    ublas::vector<double> _beta11 (2);
    ublas::vector<double> _gamma0 (2);
    ublas::vector<double> _gamma1 (2);
    ublas::vector<double> _betaY20 (2);
    ublas::vector<double> _betaY21 (2);
    ublas::vector<double> _betaZ20 (2);
    ublas::vector<double> _betaZ21 (2);
    for (int i = 0; i < 2; ++i){
      _beta10 (i) = Rbeta10 (i);
      _beta11 (i) = Rbeta11 (i);
      _gamma0 (i) = Rgamma0 (i);
      _gamma1 (i) = Rgamma1 (i);
      _betaY20 (i) = RbetaY20 (i);
      _betaY21 (i) = RbetaY21 (i);
      _betaZ20 (i) = RbetaZ20 (i);
      _betaZ21 (i) = RbetaZ21 (i);
    }
    ublas::matrix<double> _sqrtSigma (2, 2);
    for (int i = 0; i < 2; ++i){
      for (int j = 0; j < 2; ++j){
        _sqrtSigma (i, j) = RsqrtSigma (i, j);
      }
    }
    
    ublas::vector<double> _X1ts (_m);
    ublas::vector<double> _xits (_m);
    ublas::vector<double> _X2ts (_m);
    ublas::vector<int> _A1ts (_m);
    ublas::vector<int> _A2ts (_m);
    ublas::matrix<double> _H1ts (_m, 2);
    ublas::matrix<double> _H2ts (_m, 2);
    ublas::vector<double> _epsYts (_m);
    ublas::vector<double> _epsZts (_m);
    ublas::vector<double> _Yts (_m);
    ublas::vector<double> _Zts (_m);
    _X1ts.assign (ublas::zero_vector<double> (_m));
    _xits.assign (ublas::zero_vector<double> (_m));
    _X2ts.assign (ublas::zero_vector<double> (_m));
    _A1ts.assign (ublas::zero_vector<int> (_m));
    _A2ts.assign (ublas::zero_vector<int> (_m));
    _H1ts.assign (ublas::zero_matrix<double> (_m, 2));
    _H2ts.assign (ublas::zero_matrix<double> (_m, 2));
    _epsYts.assign (ublas::zero_vector<double> (_m));
    _epsZts.assign (ublas::zero_vector<double> (_m));
    _Yts.assign (ublas::zero_vector<double> (_m));
    _Zts.assign (ublas::zero_vector<double> (_m));

    ublas::vector<double> _X1 (_n);
    ublas::vector<double> _xi (_n);
    ublas::vector<double> _X2 (_n);
    ublas::vector<int> _A1 (_n);
    ublas::vector<int> _A2 (_n);
    ublas::matrix<double> _H1 (_n, 2);
    ublas::matrix<double> _H2 (_n, 2);
    ublas::vector<double> _epsY (_n);
    ublas::vector<double> _epsZ (_n);
    ublas::vector<double> _Y (_n);
    ublas::vector<double> _Z (_n);
    _X1.assign (ublas::zero_vector<double> (_n));
    _xi.assign (ublas::zero_vector<double> (_n));
    _X2.assign (ublas::zero_vector<double> (_n));
    _A1.assign (ublas::zero_vector<int> (_n));
    _A2.assign (ublas::zero_vector<int> (_n));
    _H1.assign (ublas::zero_matrix<double> (_n, 2));
    _H2.assign (ublas::zero_matrix<double> (_n, 2));
    _epsY.assign (ublas::zero_vector<double> (_n));
    _epsZ.assign (ublas::zero_vector<double> (_n));
    _Y.assign (ublas::zero_vector<double> (_n));
    _Z.assign (ublas::zero_vector<double> (_n));
    
  // declare a random number generator of type mt19937 and allocate
  //  memory to it 
    gsl_rng *_r=gsl_rng_alloc(gsl_rng_mt19937);  
    
  // set the seed 
    gsl_rng_set(_r, (unsigned long int)_ss);
    
  // Generate a test set for evaluation of policies
  // function located in genModel.hpp
    genModel::generate (_r, _m, _mu, _beta10, _beta11, _gamma0, 
      _gamma1, _betaY20, _betaY21, _betaZ20, _betaZ21, _sqrtSigma, 
      _X1ts, _A1ts, _xits, _X2ts, _A2ts, _H1ts, _H2ts, _epsYts, 
      _epsZts, _Yts, _Zts); 

  // estimate E(Y) and E(Z) under randomization
    double ranY = 0.0;
    double ranZ = 0.0;
    for (int i = 0; i < _m; ++i){
      ranY += _Yts (i);
      ranZ += _Zts (i);
    }  
    ranY /= static_cast<double>(_m);
    ranZ /= static_cast<double>(_m);

  // Generate a training set
  // function located in genModel.hpp
    genModel::generate (_r, _n, _mu, _beta10, _beta11, _gamma0, 
      _gamma1, _betaY20, _betaY21, _betaZ20, _betaZ21,
      _sqrtSigma, _X1, _A1, _xi, _X2, _A2, _H1, _H2, _epsY, _epsZ, 
      _Y, _Z);   
    
    ublas::vector<double> _etaOpt1 (2);
    ublas::vector<double> _etaOpt2 (2);
    _etaOpt1.assign (ublas::zero_vector<double> (2));
    _etaOpt2.assign (ublas::zero_vector<double> (2));
    double _meanY = 0.0;
    double _meanZ = 0.0;
    double _uncstY = 0.0;
    double _uncstZ = 0.0;
    double _mxYmnZ = 0.0;
    double _feasible = 0.0;
    int _lq2 = 200;
    // function located in learn.hpp
    cqlearn::etaLearn (_r, _X1, _A1, _X2, _A2, _H1, _H2, _Y, _Z, 
      _kappa, _alpha, _N, _J, _lq2, _meanY, _meanZ, _feasible, 
      _uncstY, _uncstZ, _mxYmnZ, _etaOpt1, _etaOpt2); 

    ublas::vector<double> _Yeval (_m);
    ublas::vector<double> _Zeval (_m);
    _Yeval.assign (ublas::zero_vector<double> (_m));
    _Zeval.assign (ublas::zero_vector<double> (_m));
    double meanYopt = 0.0;
    double meanZopt = 0.0;

    if (_feasible > 0.0){
      // Estimate E^{\etaOpt}(Y) and E^{\etaOpt}(Z)
      genEval::genOpt (_mu, _beta10, _beta11, _gamma0, 
      _gamma1, _betaY20, _betaY21, _betaZ20, _betaZ21,
      _sqrtSigma, _X1ts, _xits, _epsYts, _epsZts, _etaOpt1, 
      _etaOpt2, _Yeval, _Zeval);

      for (int i = 0; i < _m; ++i){
        meanYopt += _Yeval (i);
        meanZopt += _Zeval (i);
      }
      meanYopt /= static_cast<double>(_m);
      meanZopt /= static_cast<double>(_m);
    }  
    
  // free the memory allocated to r 
    gsl_rng_free(_r);
    
    Rcpp::NumericVector Rfeasible (1);
    Rfeasible = _feasible;
    
    Rcpp::NumericVector RetaHat1 (2);
    Rcpp::NumericVector RetaHat2 (2);
    for (int i = 0; i < 2; ++i){
      RetaHat1 (i) = _etaOpt1 (i);
      RetaHat2 (i) = _etaOpt2 (i);
    }
    Rcpp::NumericVector RmeanY (1);
    RmeanY (0) = _meanY;
    Rcpp::NumericVector RmeanYopt (1);
    RmeanYopt (0) = meanYopt;
    Rcpp::NumericVector RmeanZ (1);
    RmeanZ (0) = _meanZ;
    Rcpp::NumericVector RmeanZopt (1);
    RmeanZopt (0) = meanZopt;
    Rcpp::NumericVector RuncstY (1);
    RuncstY (0) = _uncstY;
    Rcpp::NumericVector RuncstZ (1);
    RuncstZ (0) = _uncstZ;
    Rcpp::NumericVector RmxYmnZ (1);
    RmxYmnZ (0) = _mxYmnZ;
    
    return Rcpp::List::create (Rcpp::Named ("etaHat1", RetaHat1),
      Rcpp::Named ("etaHat2", RetaHat2),
      Rcpp::Named ("feasible", Rfeasible),
      Rcpp::Named ("meanYopt", RmeanYopt),
      Rcpp::Named ("meanZopt", RmeanZopt),
      Rcpp::Named ("meanY", RmeanY),
      Rcpp::Named ("meanZ", RmeanZ),
      Rcpp::Named ("uncstY", RuncstY),
      Rcpp::Named ("uncstZ", RuncstZ),
      Rcpp::Named ("mxYmnZ", RmxYmnZ));
    
  }
  