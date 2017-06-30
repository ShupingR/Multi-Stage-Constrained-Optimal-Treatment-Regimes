/* @file estPars.hpp
 */
#ifndef _ESTPARS_HPP
#define _ESTPARS_HPP

#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "ols.hpp"
#include "genDist.hpp"

/* @fn estS2
   @brief Function to estimate parameters needed from the second stage
   regression of a response R on H2, A2.
   @param s2design Design matrix where each row is (1, X2, A2, A2*X2);
   @param R a response variable to be maximized
   @return betaHat20, betaHat21 parameter vectors
   @return residR vector of residuals from the regression
   @return sdR standard deviation of the residuals
*/

namespace stage2 {
  namespace ublas = boost::numeric::ublas;
  bool estS2 (const ublas::matrix<double>& s2design, 
		 const ublas::vector<double>& R, 
		 ublas::vector<double>& betaHat20,
		 ublas::vector<double>& betaHat21,
		 ublas::vector<double>& residsR,
		 double& sdR){

    ublas::vector<double> bHat2 (2*betaHat20.size()); 
    ols::lm (s2design, R, bHat2, residsR, sdR);

    betaHat20 (0) = bHat2 (0);
    betaHat20 (1) = bHat2 (1);
    betaHat21 (0) = bHat2 (2);
    betaHat21 (1) = bHat2 (3);

    return true;
  }
} // end namespace


/* @fn mvEst 
   @brief Function to implement mean/variance modeling for X2 given
          first-stage quantities in constrained Q-learning
   @param X2 vector whose conditional mean and variance functions we
             want to estimate
   @param s1design design matrix with each row (1, X1, A1, A1*X1)
   @param thetaHat estimated mean function parameters
   @param gammaHat estimated variance function parameters
*/
namespace meanVar{
  namespace ublas = boost::numeric::ublas;
  bool mvEst (const ublas::matrix<double>& s1design,
	      const ublas::vector<double>& X2, 
	      ublas::vector<double>& thetaHat,
	      ublas::vector<double>& gammaHat){

	      int _d1 = X2.size ();

	      // MEAN MODELING 
	      /* Fit the mean function and return estimated parameters
		 in vector _thetaHat */  
	      ublas::vector<double> _muResids (_d1);   
	      double _muSD = 0.0;   
	      ols::lm (s1design, X2, thetaHat, _muResids, _muSD);
	      //	      _muFitted = prod (_s1opt, _thetaHat); 

	      // VARIANCE MODELING
	      /* Create new response which is the log of the squared
		 residuals */  
	      ublas::vector<double> _log2Res (_d1);
	      for (int i = 0; i < _d1; ++i){
		_log2Res (i) = log (_muResids (i)*_muResids (i));
	      }

	      /* Regress log of squared resids on _s1design */
	      ublas::vector<double> _sResids (_d1);   
	      double _sSD = 0.0;   
	      ols::lm (s1design, _log2Res, gammaHat, _sResids, _sSD);

	      // STANDARDIZED RESIDUALS
	      /* Create vector of standardized residuals with
		 incorrect intercept term */
	      ublas::vector<double> _s1xgam = prod (s1design, gammaHat);
	      ublas::vector<double> _wrongInt (_d1);
	      for (int i = 0; i < _d1; ++i){
		_wrongInt (i) = _muResids (i) / exp (_s1xgam (i) / 2.0);
	      }
	      
	      /* Compute the sd of the incorrect standardized residuals,
		 correct the intercept */ 
	      double _sd2 = 0.0;
	      sdev::sdCmpt (_wrongInt, _sd2);
	      gammaHat[0] += 2*log (_sd2);
	      
	      // ublas::vector<double> _s1xgamOpt = prod (_s1opt,
	      // _gammaHat); 
	      // for (int i = 0; i < _d1; ++i){
	      //    _sFitted (i) = exp (_s1xgamOpt (i) / 2.0);
	      // }

	      // genNorm::ranNorm (_J, _mysamp);
	      
	      return true;
  }
} // end namespace	      

namespace meanVarEval {
  namespace ublas = boost::numeric::ublas;
  bool mvEval (const double& X1,
	       const int& A1, 
	       const ublas::vector<double>& thetaHat,
	       const ublas::vector<double>& gammaHat,
	       double& muEval, 
	       double& sigmaEval){

    ublas::vector<double> s1design (4);
    s1design (0) = 1.0;
    s1design (1) = X1; 
    s1design (2) = A1*1.0;
    s1design (3) = A1*X1;

    muEval = inner_prod (s1design, thetaHat);
    sigmaEval = sqrt (exp (inner_prod (s1design, gammaHat)));

    return true;
  }
} // end namespace

namespace mdls{
  namespace ublas = boost::numeric::ublas;
  bool learn (const gsl_rng*& _r, 
	      const ublas::vector<double>& X1,
	      const ublas::vector<int>& A1, 
	      const ublas::vector<double>& X2,
	      const ublas::vector<int>& A2, 
	      const ublas::vector<double>& Y, 
	      ublas::vector<double>& bHatY20,
	      ublas::vector<double>& bHatY21,
	      ublas::vector<double>& residsR,
	      double& sdR,
	      ublas::vector<double>& _thetaHat,
	      ublas::vector<double>& _gammaHat,
	      const int& _J,
	      ublas::vector<double>& _normSamp,
	      ublas::vector<double>& _unif){

    int _n = static_cast<int>(X1.size ());
    int _N = static_cast<int>(_unif.size ());
    // Create the second-stage design matrix
    ublas::matrix<double> _s2design (_n, 4);
    for (int i = 0; i < _n; ++i){
      _s2design (i, 0) = 1.0;
      _s2design (i, 1) = X2 (i);
      _s2design (i, 2) = A2 (i)*1.0;
      _s2design (i, 3) = A2 (i)*X2 (i);
    }

    // Regress Y on H2, A2
    stage2::estS2 (_s2design, Y, bHatY20, bHatY21, residsR, sdR);

    // Create the first-stage design matrix
    ublas::matrix<double> _s1design (_n, 4);
    for (int i = 0; i < _n; ++i){
      _s1design (i, 0) = 1.0;
      _s1design (i, 1) = X1 (i);
      _s1design (i, 2) = A1 (i)*1.0;
      _s1design (i, 3) = A1 (i)*X1 (i);
    }

    // Estimate parameters of the mean and variance functions of X2
    // given X1 and A1
    meanVar::mvEst (_s1design, X2, _thetaHat, _gammaHat);
    
    // Generate iid N(0, 1) samples to estimate the inner integral with
    // respect to the standardized residuals obtained from mean/variance
    // modeling of X2 given X1, A1. For now we assume these standardized
    // residuals are normal, which is correctly specified for our
    // current generative model
    // function located in genDist.hpp
    genNorm::ranNorm (_J, _r, _normSamp);

    // Generate Unif[0,1] samples
    // function located in genDist.hpp
    genUnif::ranUnif (_N, _r, _unif);
    return true;
  }
} // end namespace
	     
#endif
