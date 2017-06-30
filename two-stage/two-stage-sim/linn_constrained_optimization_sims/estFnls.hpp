  /* @file estFnls.hpp
 */
  #ifndef _ESTFNLS_HPP
  #define _ESTFNLS_HPP
  
  #include <math.h>
  #include <boost/numeric/ublas/vector.hpp>
  #include <boost/numeric/ublas/vector_proxy.hpp>
  #include <boost/numeric/ublas/matrix.hpp>
  #include <boost/math/distributions/normal.hpp>
  using boost::math::normal;
  #include <boost/range/algorithm/sort.hpp>
  #include <iostream>


  /* @fn distEval
   @brief estimates P^{\eta}_{R}(r) for response R and fixed \eta, r
  */

   namespace dist {
    namespace ublas = boost::numeric::ublas;
    bool distEval (const ublas::matrix<double>& H1,
     const ublas::vector<double>& bHat20,
     const ublas::vector<double>& bHat21,
     const double& sdR, 
     const ublas::vector<double>& thetaHat,
     const ublas::vector<double>& gammaHat,
     const double& eta10,
     const double& eta11,
     const double& eta20,
     const double& eta21,
     const ublas::vector<double>& normSamp,
     const double& q, 
     double& cdfEval){

      ublas::vector<double> eta1 (2);
      ublas::vector<double> eta2 (2);
      eta1 (0) = eta10;
      eta1 (1) = eta11;
      eta2 (0) = eta20;
      eta2 (1) = eta21;

      int _n = H1.size1 ();
      int _J = normSamp.size ();
      ublas::matrix<double> hatX2 (_n, _J);
      ublas::matrix<double> prbs (_n, _J);
      double mu = 0.0;
      double sig = 0.0;
      int a1 = 0; 
      int a2 = 0;
      double X1 = 0.0;
      double inside = 0.0;
      ublas::vector<double> h1te1 = prod (H1, eta1); 
      normal d;

      cdfEval = 0.0;

      for (int i = 0; i < _n; ++i){
        X1 = static_cast<double>(H1 (i, 1));
        a1 = static_cast<int>((h1te1 (i) >= 0) - (h1te1 (i) < 0));
        meanVarEval::mvEval (X1, a1, thetaHat, gammaHat, mu, sig); 

        for (int j = 0; j < _J; ++j){
         hatX2 (i, j) = mu + sig*normSamp (j); 
         a2 = static_cast<int>((eta2 (0) + hatX2 (i, j)*eta2 (1) >= 0) - (eta2 (0) + hatX2 (i, j)*eta2 (1) < 0)); 
         inside = q - bHat20 (0) - hatX2 (i, j)*bHat20 (1) 
         - a2*(bHat21 (0) + hatX2 (i, j)*bHat21 (1));
         inside /= sdR; 
         prbs (i, j) = cdf (d, inside);
         cdfEval += prbs (i, j);
       }
     }

     cdfEval /= static_cast<double>(_n*_J); 
     return true;
   }
  } // end namespace

  /* @fn estQnts
   @brief estimates quantiles of a vector
  */
   namespace qnts {
    namespace ublas = boost::numeric::ublas;
    bool estQnts (const ublas::vector<double>& R,
      const ublas::vector<double>& prbs,
      ublas::vector<double>& qts){

      int _n = R.size ();
      int _lq = prbs.size ();
    /* Get quantiles of R*/
 //     int _lq1 = _lq + 1;
 //     int _lq2 = _lq + 2;
      ublas::vector<double> newR (_n);
      for (int i = 0; i < _n; ++i){
        newR (i) = R (i);
      }
      ublas::vector<double> _sortR = boost::sort (newR);
      for (int i = 0; i < _lq; ++i){
        qts (i) = _sortR (static_cast<int>(floor (_n*prbs (i))));
      }
    // std::cout << _qntR << "This is _qntR \n";

      return true;
    }
  } // end namespace


  /* @fn distEval
   @brief estimates P^{\eta}_{R}(r) for response R and fixed \eta, r
  */

   namespace fnl {
    namespace ublas = boost::numeric::ublas;
    bool fnlEval (const ublas::matrix<double>& H1,
      const ublas::vector<double>& R,
      const ublas::vector<double>& bHat20,
      const ublas::vector<double>& bHat21,
      const double& sdR, 
      const ublas::vector<double>& thetaHat,
      const ublas::vector<double>& gammaHat,
      const ublas::vector<double>& eta1,
      const ublas::vector<double>& eta2,
      const ublas::vector<double>& normSamp,
      const ublas::vector<double>& unif,
      const int& _N, 
      const int& _lq2,
      double& meanR){

 //     int _n = H1.size1 ();
 //     int _J = normSamp.size ();
      ublas::vector<double> _prSeq (2);
      _prSeq (0) = 0.01;
      _prSeq (1) = 0.99;
      ublas::vector<double> _init (2);
      _init.assign (ublas::zero_vector<double> (2));
    // function located in estFnls.hpp
      qnts::estQnts (R, _prSeq, _init);

      double _dist = _init (1) - _init (0);
      double _lower = _init (0) - _dist/4.0;
      double _upper = _init (1) + _dist/4.0;

    // create equally spaced grid of 200 points from _lower to
    // _upper 
      ublas::vector<double> _qts (_lq2);
      _qts (0) = _lower;
      double _steps = (_upper - _lower)/static_cast<double>(_lq2 - 1);
      for (int i = 1; i < _lq2; ++i){
        _qts (i) = _qts (i - 1) + _steps;
      }

    /* Evaluate the distribution function of R under (eta1, eta2)
       regime at each grid point. */  
      ublas::vector<double> _prbsR (_lq2);
      double cdfVal = 0.0;
      double eta10 = eta1 (0);
      double eta11 = eta1 (1);
      double eta20 = eta2 (0);
      double eta21 = eta2 (1);

      for (int i = 0; i < _lq2; ++i){
        dist::distEval (H1, bHat20, bHat21, sdR, thetaHat, gammaHat,
          eta10, eta11, eta20, eta21, normSamp, _qts (i),
          cdfVal);   
        _prbsR (i) = cdfVal;
      }

  //      std::cout << _prbsR << "This is _prbsR \n";
    /* Get interval probabilities and the quantile midpoints of those
       intervals */  
      int _lq1 = _lq2 - 1;
      int _lq3 = _lq2 + 1;
      ublas::vector<double> _midsR (_lq3);
      _midsR (0) = _qts (0) - _steps/2.0; 
      _midsR (_lq2) = _qts (_lq1) + _steps/2.0;
      for (int i = 0; i < _lq1; ++i){
        _midsR (i + 1) = (_qts (i) + _qts (i + 1)) / 2.0;
      }

      ublas::vector<int> _multisampR (_lq3);
      _multisampR.assign (ublas::zero_vector<int> (_lq3)); 
      for (int j = 0; j < _N; ++j){
        if (unif (j) <= _prbsR (0)){
          _multisampR (0) += 1;
        }
      }
      for (int i = 1; i < _lq2; ++i){
        for (int j = 0; j < _N; ++j){
          if (unif (j) <= _prbsR (i) && unif (j) > _prbsR (i - 1)){
            _multisampR (i) += 1;
          }
        }
      }
      for (int j = 0; j < _N; ++j){
        if (unif (j) > _prbsR (_lq1)){
          _multisampR (_lq2) += 1;
        }
      }

  //   ublas::vector<int> _multisampR (_lq1);
  //  genMultinom::ranMultinom (_mrgR, _N, _r, _multisampR);

      int K = 0;
      int incr = 0; 
      ublas::vector<double> _fnlRsamp (_N);
      for (int i = 0; i < _lq2; ++i){
        incr = K + _multisampR (i);
        for (int j = K; j < incr; ++j){
          _fnlRsamp (j) = _midsR (i); 
        }
        K += _multisampR (i);
      }
      
      meanR = 0.0;
    /* _fnlRsamp is a sample from the distribution function of R under regime (beta1, beta2) */
      for (int i = 0; i < _N; ++i){
        meanR += _fnlRsamp (i);
      }
      meanR /= static_cast<double>(_N);
      
      return true;
    }
  } // end namespace
  
  # endif
