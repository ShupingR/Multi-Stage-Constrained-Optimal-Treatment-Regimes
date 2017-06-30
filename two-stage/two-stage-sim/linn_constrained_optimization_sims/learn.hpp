 /* @file test.cpp 
	@brief Code to test the constrained opt simulation.  
  */

 // some of these might not be necessary 

 #include <Rcpp.h>
 #include <math.h>
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include "genDist.hpp"
 #include "estPars.hpp"
 #include "estFnls.hpp"
 #include "gslFnls.hpp"

	namespace cqlearn{
		using namespace std;
		namespace ublas = boost::numeric::ublas;
		bool etaLearn (const gsl_rng* _r,
			const ublas::vector<double>& _X1, 
			const ublas::vector<int>& _A1, 	 
			const ublas::vector<double>& _X2, 
			const ublas::vector<int>& _A2, 
			const ublas::matrix<double>& _H1, 
			const ublas::matrix<double>& _H2,
			const ublas::vector<double>& _Y, 
			const ublas::vector<double> _Z, 
			const double& _kappa, 
			const double& _alpha,
			const int& _N, 
			const int& _J, 
      const int& _lq2,
			double& meanY, 
			double& meanZ, 
			double& feasible, 
      double& uncstY,
      double& uncstZ,
      double& mxYmnZ,
			ublas::vector<double>& etaOpt1, 
			ublas::vector<double>& etaOpt2){   
			int _n = _X1.size ();

      // Estimate E(Y) under randomization
      double ranY = 0.0;
      for (int i = 0; i < _n; ++i){
        ranY += _Y (i);
      }
      ranY /= static_cast<double>(_n);

  // Estimate parameters 
			ublas::vector<double> _bHatY20 (2);
			ublas::vector<double> _bHatY21 (2);
			ublas::vector<double> _residsY (_n);
			double _sdY = 0.0;
			ublas::vector<double> _thetaHatY (4);
			ublas::vector<double> _gammaHatY (4);
			ublas::vector<double> _normSampY (_J);
			_normSampY.assign (ublas::zero_vector<double> (_J));
			ublas::vector<double> _unifY (_N);
			_unifY.assign (ublas::zero_vector<double> (_N));

  // function located in estPars.hpp
      mdls::learn (_r, _X1, _A1, _X2, _A2, _Y, _bHatY20, 
        _bHatY21, _residsY, _sdY, _thetaHatY, _gammaHatY, _J, 
        _normSampY, _unifY); 

      ublas::vector<double> eHat1 (2);
      ublas::vector<double> eHat2 (2);
      eHat1.assign (ublas::zero_vector<double> (2));
      eHat2.assign (ublas::zero_vector<double> (2));

  // Now we have all parameters needed to max E^{\eta}(Y) wrt
  // \eta 
  // using GSL, which is equivalent to min -E^{\eta}(Y).  
      int _lenpar = 3 + _n + 2 + 2 + 1 + 4 + 4 + _J + _lq2 + 2 + _N;  
      double _max = 1.0;
      ublas::vector<double> _myparsY (_lenpar);
      _myparsY.assign (ublas::zero_vector<double> (_lenpar));
  // function located in gslFnls.hpp
      prep::preOpt (_Y, _J, _X1, _bHatY20, _bHatY21, _sdY,
        _thetaHatY, _gammaHatY, _normSampY, _unifY, _N, _max, _lq2, 
        _myparsY);  

  // function located in genDist.hpp
      genNorm::ranNorm (2, _r, eHat1);
      genNorm::ranNorm (2, _r, eHat2);

	// function located in gslFnls.hpp
   		gslmain::main (_myparsY, eHat1, eHat2);

      // Estimate mean of Y under the est optimal regime
      meanY = 0.0;
      fnl::fnlEval (_H1, _Y, _bHatY20, _bHatY21, _sdY, 
        _thetaHatY, _gammaHatY, eHat1, eHat2, _normSampY, 
        _unifY, _N, _lq2, meanY); 

      uncstY = meanY; 
      double diffY = meanY - ranY;

    // Check mean of Z under estimated optimal regime
    // Estimate parameters 
   		ublas::vector<double> _bHatZ20 (2);
   		ublas::vector<double> _bHatZ21 (2);
   		ublas::vector<double> _residsZ (_X1.size ());
   		double _sdZ = 0.0;
   		ublas::vector<double> _thetaHatZ (4);
   		ublas::vector<double> _gammaHatZ (4);
   		ublas::vector<double> _normSampZ (_J);
   		_normSampZ.assign (ublas::zero_vector<double> (_J));
   		ublas::vector<double> _unifZ (_N);
   		_unifZ.assign (ublas::zero_vector<double> (_N));

   		mdls::learn (_r, _X1, _A1, _X2, _A2, _Z, _bHatZ20, _bHatZ21,
   			_residsZ, _sdZ, _thetaHatZ, _gammaHatZ, _J,
   			_normSampZ, _unifZ); 

   		meanZ = 0.0;
   		fnl::fnlEval (_H1, _Z, _bHatZ20, _bHatZ21, _sdZ, 
            _thetaHatZ, _gammaHatZ, eHat1, eHat2, _normSampZ, 
            _unifZ, _N, _lq2, meanZ); 

      mxYmnZ = meanZ; 
   		feasible = 0.0;  

   		if (meanZ < _kappa)
   		{
	// if estimated E^{\eta}(Z) < kappa, the constraint is
	// met for 
	// the unconstrained optimal \eta and we are done
 
   			etaOpt1 (0) = eHat1 (0);
   			etaOpt1 (1) = eHat1 (1);
   			etaOpt2 (0) = eHat2 (0);
   			etaOpt2 (1) = eHat2 (1);
   			feasible = 1.0; 
   			return true;
   		}
   		else
   		{
   			_max = -1.0;
   			ublas::vector<double> _myparsZ (_lenpar);
   			_myparsZ.assign (ublas::zero_vector<double> (_lenpar));
   			prep::preOpt (_Z, _J, _X1, _bHatZ20, _bHatZ21, _sdZ,
   				_thetaHatZ, _gammaHatZ, _normSampZ, _unifZ, _N, _max,
   				_lq2, _myparsZ);  

   			genNorm::ranNorm (2, _r, eHat1);
   			genNorm::ranNorm (2, _r, eHat2);
   			gslmain::main (_myparsZ, eHat1, eHat2);

	// check min E^{\eta}(Z) to see if any feasible regime exists
   			meanZ = 0.0;
   			fnl::fnlEval (_H1, _Z, _bHatZ20, _bHatZ21, _sdZ, _thetaHatZ, 
   				_gammaHatZ, eHat1, eHat2, _normSampZ, _unifZ, _N, _lq2,
   				meanZ);  

        ublas::vector<double> etaOptZ1 (2);
        ublas::vector<double> etaOptZ2 (2);
        for (int i = 0; i < 2; ++i){
          etaOptZ1 (i) = eHat1 (i);
          etaOptZ2 (i) = eHat2 (i);
        }

        double _minAdZ = 0.0;
        for (int i = 0; i < 24; ++i){
          genNorm::ranNorm (2, _r, eHat1);
          genNorm::ranNorm (2, _r, eHat2);
          gslmain::main (_myparsZ, eHat1, eHat2);
          fnl::fnlEval (_H1, _Z, _bHatZ20, _bHatZ21, _sdZ, _thetaHatZ, 
            _gammaHatZ, eHat1, eHat2, _normSampZ, _unifZ, _N, _lq2,
            _minAdZ);  
          if (_minAdZ < meanZ){
            for (int j = 0; j < 2; ++j){
              etaOptZ1 (j) = eHat1 (j);
              etaOptZ2 (j) = eHat2 (j);
            }
            meanZ = _minAdZ;
          }
          _minAdZ = 0.0;
        }
        uncstZ = meanZ;

   			if (meanZ >= _kappa)
   			{
   				feasible = -1.0;

   				return true;
   			}
   			else
   			{
          feasible = 2.0;
          double logterm = log (_kappa - meanZ);
   				double lambda = diffY/logterm;
          ublas::vector<double> lamvec (6);
          if (lambda < 0.0){
            lambda *= -1.0;
          }
          lamvec (0) = lambda;
   				int lenparB = 3 + _n + 8 + 2 + 16 + 2*_J + 2*_lq2 + 4 + 2*_N; 
   				ublas::vector<double> _myparsB (lenparB);
   				_myparsB.assign (ublas::zero_vector<double> (lenparB));
   				prepB::preOptB (_Y, _Z, _J, _X1, _bHatY20, _bHatY21,
   					_bHatZ20, _bHatZ21, _sdY, _sdZ, _thetaHatY,
   					_gammaHatY, _thetaHatZ,
   					_gammaHatZ, _normSampY, _normSampZ, _unifY, 
   					_unifZ, _N, lambda, _kappa, _lq2, _myparsB); 
          double _minEst = 0.0;
   				gslmainB::mainB (_myparsB, eHat1, eHat2, _minEst); 
   				for (int i = 0; i < 2; ++i){
   					etaOpt1 (i) = eHat1 (i);
   					etaOpt2 (i) = eHat2 (i);
   				}

          double _minAd = 0.0;
          for (int i = 0; i < 8; ++i){
//            genNorm::ranNorm (2, _r, eHat1);
//            genNorm::ranNorm (2, _r, eHat2);
            gslmainB::mainB (_myparsB, eHat1, eHat2, _minAd);
            if (_minAd < _minEst){
              for (int j = 0; j < 2; ++j){
                etaOpt1 (j) = eHat1 (j);
                etaOpt2 (j) = eHat2 (j);
              }
              _minEst = _minAd;
            }
            _minAd = 0.0;
          }

     // Estimate mean of Y under this regime
          ublas::vector<double> mnYB (6);
          meanY = 0.0;
          fnl::fnlEval (_H1, _Y, _bHatY20, _bHatY21, _sdY, 
            _thetaHatY, _gammaHatY, etaOpt1, etaOpt2, _normSampY, 
            _unifY, _N, _lq2, meanY); 
          mnYB (0) = meanY;

          ublas::vector<double> _mins (6);
   				_mins (0) = _minEst;
   				// double _minopt = _minEst;

          ublas::vector<double> etaTemp1 (2);
          ublas::vector<double> etaTemp2 (2);
          for (int i = 0; i < 2; ++i){
            etaTemp1 (i) = etaOpt1 (i);
            etaTemp2 (i) = etaOpt2 (i);
          }

	    // Loop over the seq of decreasing (to zero) lambdas
   				for (int i = 0; i < 5; ++i){
   					lambda *= _alpha;
            lamvec (i + 1) = lambda;
   					prepB::preOptB (_Y, _Z, _J, _X1, _bHatY20, _bHatY21,
   						_bHatZ20, _bHatZ21, _sdY, _sdZ, _thetaHatY,
   						_gammaHatY, _thetaHatZ,
   						_gammaHatZ, _normSampY, _normSampZ, 
   						_unifY, _unifZ, _N,
   						lambda, _kappa, _lq2, _myparsB);
   					gslmainB::mainB (_myparsB, etaTemp1, etaTemp2, _minEst); 
            _minAd = 0.0;
            for (int j = 0; j < 5; ++j){
              genNorm::ranNorm (2, _r, eHat1);
              genNorm::ranNorm (2, _r, eHat2);
              for (int k = 0; k < 2; ++k){
                eHat1 (k) *= 0.25;
                eHat2 (k) *= 0.25;
                eHat1 (k) += etaTemp1 (k);
                eHat2 (k) += etaTemp2 (k);
              }
              gslmainB::mainB (_myparsB, eHat1, eHat2, _minAd);
              if (_minAd < _minEst){
                for (int j = 0; j < 2; ++j){
                  etaTemp1 (j) = eHat1 (j);
                  etaTemp2 (j) = eHat2 (j);
                }
                _minEst = _minAd;
              }
              _minAd = 0.0;
            }

      // Estimate mean of Y under this regime
//            meanY = 0.0;
//            fnl::fnlEval (_H1, _Y, _bHatY20, _bHatY21, _sdY, 
//              _thetaHatY, _gammaHatY, etaTemp1, etaTemp2, 
//              _normSampY, 
//              _unifY, _N, _lq2, meanY); 
//            mnYB (i + 1) = meanY; 

//  					_mins (i + 1) = _minEst;
//   					if (_minEst < _mins (i))
//   					{
//   						_minopt = _minEst;
//   					}
          }

          for (int j = 0; j < 2; ++j){
            etaOpt1 (j) = etaTemp1 (j);
            etaOpt2 (j) = etaTemp2 (j);
          }
            return true;

          }
        }
      }
	} // end namespace
	