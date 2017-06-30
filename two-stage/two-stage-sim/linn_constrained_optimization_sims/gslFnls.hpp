  /* @file gslFnls.hpp
 */
  #ifndef _GSLFNLS_HPP
  #define _GSLFNLS_HPP
  
  #include <math.h>
  #include <boost/numeric/ublas/vector.hpp>
  #include <boost/numeric/ublas/vector_proxy.hpp>
  #include <boost/numeric/ublas/matrix.hpp>
  #include <boost/math/distributions/normal.hpp>
  using boost::math::normal;
  #include <boost/range/algorithm/sort.hpp>
  #include <boost/any.hpp>
  #include <gsl/gsl_randist.h>
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_blas.h>
  #include <cmath>
  #include <algorithm>    // std::unique, std::distance
  #include <gsl/gsl_multimin.h>
  #include <iostream>
  #include <limits>

  #include "estFnls.hpp"
  #include "genDist.hpp"

  // @fn preOpt
  // @brief function that combines all constant paramerters into a
  // single vector to conform to GSL optimizer format.
  // @note this function is for optimizing the mean only, not the
  // barrier function

  namespace prep {
    namespace ublas = boost::numeric::ublas;
    bool preOpt (const ublas::vector<double>& _R, 
      const int& _J,  
      const ublas::vector<double>& _X1,
      const ublas::vector<double>& _bHat20,
      const ublas::vector<double>& _bHat21,
      const double& _sdR, 
      const ublas::vector<double>& _thetaHat,
      const ublas::vector<double>& _gammaHat,
      const ublas::vector<double>& _normSamp,
      const ublas::vector<double>& _unif,
      const int& _N,
      const double& _max, 
      const int& _lq2,
      ublas::vector<double>& mypars){

      int _n = _R.size (); 
      ublas::vector<double> _prSeq (2);
      _prSeq (0) = 0.01;
      _prSeq (1) = 0.99;
      ublas::vector<double> _init (2);
      _init.assign (ublas::zero_vector<double> (2));
    // function located in estFnls.hpp
      qnts::estQnts (_R, _prSeq, _init);
      double _dist = _init (1) - _init (0);
      double _lower = _init (0) - _dist/4.0;
      double _upper = _init (1) + _dist/4.0;

    // create equally spaced grid of 200 points from _lower to
    // _upper 
      ublas::vector<double> _grid (_lq2);
      _grid (0) = _lower;
      double _steps = (_upper - _lower)/static_cast<double>(_lq2 - 1);
      for (int i = 1; i < _lq2; ++i){
        _grid (i) = _grid (i - 1) + _steps;
      }

      mypars (0) = static_cast<double>(_n);
      mypars (1) = static_cast<double>(_lq2);
      mypars (2) = static_cast<double>(_J);
      for (int i = 0; i < _n; ++i){
        mypars (i + 3) = _X1 (i);
      }
      for (int i = 0; i < 2; ++i){
        mypars (3 + _n + i) = _bHat20 (i);
        mypars (3 + _n + 2 + i) = _bHat21 (i);
      }
      mypars (3 + _n + 4) = _sdR;
      for (int i = 0; i < 4; ++i){
        mypars (3 + _n + 4 + 1 + i) = _thetaHat (i);
        mypars (3 + _n + 4 + 1 + 4 + i) = _gammaHat (i);
      }
      for (int i = 0; i < _J; ++i){
        mypars (3 + _n + 4 + 1 + 8 + i) = _normSamp (i);
      }
      for (int i = 0; i < _lq2; ++i){
        mypars (3 + _n + 4 + 1 + 8 + _J + i) = _grid (i);
      }
      mypars (3 + _n + 4 + 1 + 8 + _J + _lq2) = 
      static_cast<double>(_N);         
      mypars (3 + _n + 4 + 1 + 8 + _J + _lq2 + 1) = _max;  
      for (int i = 0; i < _N; ++i){
        mypars (3 + _n + 4 + 1 + 8 + _J + _lq2 + 2 + i) = _unif (i);  
      }

      return true;

    }
  } // end namespace


  // @fn gslFn
  // @brief function to be passed into GSL optimizer that 
  // estimates the mean under regime \eta.

  namespace gslfnl {
    namespace ublas = boost::numeric::ublas;
    double gslFn (const gsl_vector *v, void *params){ 
      double eta10 = gsl_vector_get (v, 0);
      double eta11 = gsl_vector_get (v, 1);
      double eta20 = gsl_vector_get (v, 2);
      double eta21 = gsl_vector_get (v, 3);

      ublas::vector<double> p = *(ublas::vector<double> *)params; 
    // first param is the sample size
      int _n = static_cast<int>(p (0));
    // second param is length of quantile vector
      int _lq2 = static_cast<int>(p (1));
    // third param is number of MC samples to approx inner integral
      int _J = static_cast<int>(p (2));
    // fourth param is X1
      ublas::vector<double> X1 (_n);
      ublas::matrix<double> H1 (_n, 2);
      for (int i = 0; i < _n; ++i){
        X1 (i) = p (i + 3);
        H1 (i, 0) = 1.0;
        H1 (i, 1) = X1 (i);
      }
      ublas::vector<double> bHat20 (2);
      ublas::vector<double> bHat21 (2);
      for (int i = 0; i < 2; ++i){
        bHat20 (i) = p (3 + _n + i);
        bHat21 (i) = p (3 + _n + 2 + i);
      }
      double sdR = p (3 + _n + 2 + 2);
      ublas::vector<double> thetaHat (4);
      ublas::vector<double> gammaHat (4);
      for (int i = 0; i < 4; ++i){
        thetaHat (i) = p (3 + _n + 2 + 2 + 1 + i);
        gammaHat (i) = p (3 + _n + 2 + 2 + 1 + 4 + i);
      }
      ublas::vector<double> normSamp (_J);
      for (int i = 0; i < _J; ++i){
        normSamp (i) = p (3 + _n + 2 + 2 + 1 + 4 + 4 + i);
      }
      ublas::vector<double> qts (_lq2);
      for (int i = 0; i < _lq2; ++i){
        qts (i) = p (3 + _n + 2 + 2 + 1 + 4 + 4 + _J + i);
      }
    // Number of U[0, 1] samples to transform to multinomial
      int _N = static_cast<int>(p (3 + _n + 2 + 2 + 1 + 4 + 4 + _J + 
       _lq2)); 
      double maximize = p (3 + _n + 2 + 2 + 1 + 4 + 4 + _J + _lq2 + 1);

      ublas::vector<double> unif (_N);
      for (int i = 0; i < _N; ++i){
        unif (i) = p (3 + _n + 4 + 1 + 8 + _J + _lq2 + 2 + i);
      }

      double norm1 = sqrt (eta10*eta10 + eta11*eta11); 
      double norm2 = sqrt (eta20*eta20 + eta21*eta21); 
      eta10 /= norm1;
      eta11 /= norm1;
      eta20 /= norm2;
      eta21 /= norm2;
    /* Evaluate the distribution function of R under (eta1, eta2)
       regime */  
      ublas::vector<double> _prbsR (_lq2);
      double cdfVal = 0.0;
      for (int i = 0; i < _lq2; ++i){
      // function located in estFnls.hpp
        dist::distEval (H1, bHat20, bHat21, sdR, thetaHat, gammaHat,
          eta10, eta11, eta20, eta21, normSamp, qts (i),
          cdfVal);   
        _prbsR (i) = cdfVal;
      }

    /* Get interval probabilities and the quantile midpoints of those
       intervals */  
      int _lq1 = _lq2 - 1;
      int _lq3 = _lq2 + 1;
      double _steps = qts (1) - qts (0);
      ublas::vector<double> _midsR (_lq3);
      _midsR (0) = qts (0) - _steps/2.0; 
      _midsR (_lq2) = qts (_lq1) + _steps/2.0;
      for (int i = 0; i < _lq1; ++i){
        _midsR (i + 1) = (qts (i) + qts (i + 1)) / 2.0;
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

      //std::cout << _multisampR << "This is _multisampR \n";
    //ublas::vector<int> _multisampR (_lq1);
    //genMultinom::ranMultinom (_mrgR, _N, _multisampR);

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

    /* _fnlRsamp is a sample from the distribution function of R under regime (beta1, beta2); 
    */

      double meanR = 0.0;
      for (int i = 0; i < _N; ++i){
        meanR += _fnlRsamp (i);
      }

    // If trying to maximize, need to convert to min the negative mean.
    // Else if trying to minimize, return the mean.

      bool _max = false;
      if (maximize > 0.0){
        _max = true;
      }

      if (_max)
      { 
        meanR /= -1.0*static_cast<double>(_N);
 
        return meanR;
      }

      meanR /= static_cast<double>(_N);
 
      return meanR;
    }
  } // end namespace

  // @fn main
  // @brief calls the GSL optimizer for optimizing the mean of a
  // response under regime \eta (not the barrier function)

  namespace gslmain {
    namespace ublas = boost::numeric::ublas;
    int main(const ublas::vector<double>& mypars,
      ublas::vector<double>& etaHat1,
      ublas::vector<double>& etaHat2){

    // declare the minimizer algorithm to use
    // Nelder Mead is the only one that doesn't require gradient
      const gsl_multimin_fminimizer_type *T = 
      gsl_multimin_fminimizer_nmsimplex;

    // this is a pointer to a class that will have everything needed
    // to optimize 
      gsl_multimin_fminimizer *s = NULL;

    // ss is a pointer to a gsl_vector of step sizes
      gsl_vector *ss;

    // x is a pointer to a gsl_vector of variables to minimize the 
    // function
      gsl_vector *x;

    // the class that holds all information for the objective 
    // function
      gsl_multimin_function minex_func;

      size_t iter = 0; 
    // number of iterations
      int status; 
    // convergence status
      double size; 
    // a distance for how far the minimum is away from
    // current iter 
      int dim=4; 
    // dimension of input to objective function

    // set initial values
      x=gsl_vector_alloc(dim);
      gsl_vector_set(x,0,etaHat1 (0));
      gsl_vector_set(x,1,etaHat1 (1));
      gsl_vector_set(x,2,etaHat2 (0));
      gsl_vector_set(x,3,etaHat2 (1));

    // set initial step sizes (these adpatively change with this algorithm)
      ss=gsl_vector_alloc(dim);
      gsl_vector_set_all(ss,3.0);

      ublas::vector<double> pvec (mypars.size ());
      for (int i = 0; i < mypars.size (); ++i){
        pvec (i) = mypars (i);
      }

      minex_func.n=dim; 
    // number of function components
      minex_func.f=&gslfnl::gslFn; 
    // functor
      minex_func.params=&pvec; 
    // initial values
    
    // allocate memory for the minimizer class
      s=gsl_multimin_fminimizer_alloc(T,dim);
      gsl_multimin_fminimizer_set(s,&minex_func,x,ss);
      
    // minimize
      do{ 
        iter++;
        
      // to iterate once its a function call
      // the return value is the status of conergence
      // its up to the user to decide when to stop
      // if this function returns true, then there was an error
        status=gsl_multimin_fminimizer_iterate(s);
        if(status)
         break;
       
      // not sure what exactly this is, but its used to evaluate convergence
      // its some sort of estimate as to how far the algorithm estimates the
      // minimum to be from the current iteration
       size=gsl_multimin_fminimizer_size(s);
       
      // this actually tests for convergence.  1e-2 is the tolerance
       status=gsl_multimin_test_size(size,1e-2);
       
      // GSL_SUCCESS is a global constant within the GSL header files
  //       if(status == GSL_SUCCESS)
  //     {
  //       printf("converged to minimium at\n");
  //     }

       etaHat1 (0) = static_cast<double>(gsl_vector_get(s->x,0));
       etaHat1 (1) = static_cast<double>(gsl_vector_get(s->x,1));
       etaHat2 (0) = static_cast<double>(gsl_vector_get(s->x,2));
       etaHat2 (1) = static_cast<double>(gsl_vector_get(s->x,3));
      // minEst = s->fval;
       
     } 
    // GSL_CONTINUE is a global constant within the GSL header files
     while(status == GSL_CONTINUE && iter < 100);
     
    // free the memory allocated above
     gsl_vector_free(x);
     gsl_vector_free(ss);
     gsl_multimin_fminimizer_free(s);
     
     return 0;
   }
  } // end namespace
  

  // @fn preOptB
  // @brief combines parameters into a single vector for the barrier
  // optimization method performed using the GSL optimizer.

  namespace prepB {
    namespace ublas = boost::numeric::ublas;
    bool preOptB (const ublas::vector<double>& _Y, 
      const ublas::vector<double>& _Z, 
      const int& _J, 
      const ublas::vector<double>& _X1,
      const ublas::vector<double>& _bHatY20,
      const ublas::vector<double>& _bHatY21,
      const ublas::vector<double>& _bHatZ20,
      const ublas::vector<double>& _bHatZ21,
      const double& _sdY, 
      const double& _sdZ, 
      const ublas::vector<double>& _thetaHatY,
      const ublas::vector<double>& _gammaHatY,
      const ublas::vector<double>& _thetaHatZ,
      const ublas::vector<double>& _gammaHatZ,
      const ublas::vector<double>& _normSampY,
      const ublas::vector<double>& _normSampZ,
      const ublas::vector<double>& _unifY,
      const ublas::vector<double>& _unifZ,
      const int& _N,
      const double& _lambda, 
      const double& _kappa, 
      const int& _lq2,
      ublas::vector<double>& mypars){

      int _n = _Y.size (); 
      ublas::vector<double> _prSeq (2);
      _prSeq (0) = 0.01;
      _prSeq (1) = 0.99;
      ublas::vector<double> _initY (2);
      ublas::vector<double> _initZ (2);
      _initY.assign (ublas::zero_vector<double> (2));
      _initZ.assign (ublas::zero_vector<double> (2));
    // function located in estFnls.hpp
      qnts::estQnts (_Y, _prSeq, _initY);
      qnts::estQnts (_Z, _prSeq, _initZ);
      double _distY = _initY (1) - _initY (0);
      double _lowerY = _initY (0) - _distY/4.0;
      double _upperY = _initY (1) + _distY/4.0;
      double _distZ = _initZ (1) - _initZ (0);
      double _lowerZ = _initZ (0) - _distZ/4.0;
      double _upperZ = _initZ (1) + _distZ/4.0;

    // create equally spaced grids of 200 points from _lower to
    // _upper 
      ublas::vector<double> _gridY (_lq2);
      ublas::vector<double> _gridZ (_lq2);
      _gridY (0) = _lowerY;
      double _stepsY = (_upperY - _lowerY)/static_cast<double>(_lq2 -
        1);
      _gridZ (0) = _lowerZ;
      double _stepsZ = (_upperZ - _lowerZ)/static_cast<double>(_lq2 -
        1);
      for (int i = 1; i < _lq2; ++i){
        _gridY (i) = _gridY (i - 1) + _stepsY;
        _gridZ (i) = _gridZ (i - 1) + _stepsZ;
      }

      mypars (0) = static_cast<double>(_n);
      mypars (1) = static_cast<double>(_lq2);
      mypars (2) = static_cast<double>(_J);
      for (int i = 0; i < _n; ++i){
        mypars (i + 3) = _X1 (i);
      }
      for (int i = 0; i < 2; ++i){
        mypars (3 + _n + i) = _bHatY20 (i);
        mypars (3 + _n + 2 + i) = _bHatY21 (i);
      }
      for (int i = 0; i < 2; ++i){
        mypars (3 + _n + 4 + i) = _bHatZ20 (i);
        mypars (3 + _n + 4 + 2 + i) = _bHatZ21 (i);
      }
      mypars (3 + _n + 4 + 4) = _sdY;
      mypars (3 + _n + 4 + 4 + 1) = _sdZ;
      for (int i = 0; i < 4; ++i){
        mypars (3 + _n + 4 + 4 + 2 + i) = _thetaHatY (i);
        mypars (3 + _n + 4 + 4 + 2 + 4 + i) = _gammaHatY (i);
      }
      for (int i = 0; i < 4; ++i){
        mypars (3 + _n + 8 + 2 + 8 + i) = _thetaHatZ (i);
        mypars (3 + _n + 8 + 2 + 8 + 4 + i) = _gammaHatZ (i);
      }
      for (int i = 0; i < _J; ++i){
        mypars (3 + _n + 8 + 2 + 8 + 8 + i) = _normSampY (i);
      }
      for (int i = 0; i < _J; ++i){
        mypars (3 + _n + 8 + 2 + 8 + 8 + _J + i) = _normSampZ (i);
      }
      for (int i = 0; i < _lq2; ++i){
        mypars (3 + _n + 8 + 2 + 16 + 2*_J + i) = _gridY (i);
      }
      for (int i = 0; i < _lq2; ++i){
        mypars (3 + _n + 8 + 2 + 16 + 2*_J + _lq2 + i) = _gridZ (i);
      }
      mypars (3 + _n + 8 + 2 + 16 + 2*_J + 2*_lq2) =
      static_cast<double>(_N); 
      mypars (3 + _n + 8 + 2 + 16 + 2*_J + 2*_lq2 + 1) =
      static_cast<double>(_lambda); 
      mypars (3 + _n + 8 + 2 + 16 + 2*_J + 2*_lq2 + 2) =
      static_cast<double>(_kappa); 
      for (int i = 0; i < _N; ++i){
        mypars (3 + _n + 8 + 2 + 16 + 2*_J + 2*_lq2 + 3 + i) = _unifY (i);
      }
      for (int i = 0; i < _N; ++i){
        mypars (3 + _n + 8 + 2 + 16 + 2*_J + 2*_lq2 + 3 + _N + i) = _unifZ (i); 
      }
   
    return true;
    
  }
} // end namespace

// @fn gslFn
// @brief function to be passed into GSL optimizer that 
// estimates the mean under regime \eta.

  namespace gslfnlB {
    namespace ublas = boost::numeric::ublas;
    double gslFnB (const gsl_vector *v, void *params){ 
      double eta10 = gsl_vector_get (v, 0);
      double eta11 = gsl_vector_get (v, 1);
      double eta20 = gsl_vector_get (v, 2);
      double eta21 = gsl_vector_get (v, 3);

      ublas::vector<double> p = *(ublas::vector<double> *)params;
    // first param is the sample size
      int _n = static_cast<int>(p (0));
    // second param is length of quantile vector
      int _lq2 = static_cast<int>(p (1));
    // third param is number of MC samples to approx inner integral
      int _J = static_cast<int>(p (2));
    // fourth param is X1
      ublas::vector<double> X1 (_n);
      ublas::matrix<double> H1 (_n, 2);
      for (int i = 0; i < _n; ++i){
        X1 (i) = p (i + 3);
        H1 (i, 0) = 1.0;
        H1 (i, 1) = X1 (i);
      }
      ublas::vector<double> bHatY20 (2);
      ublas::vector<double> bHatY21 (2);
      for (int i = 0; i < 2; ++i){
        bHatY20 (i) = p (3 + _n + i);
        bHatY21 (i) = p (3 + _n + 2 + i);
      }
      ublas::vector<double> bHatZ20 (2);
      ublas::vector<double> bHatZ21 (2);
      for (int i = 0; i < 2; ++i){
        bHatZ20 (i) = p (3 + _n + 4 + i);
        bHatZ21 (i) = p (3 + _n + 4 + 2 + i);
      }
      double sdY = p (3 + _n + 4 + 4);
      double sdZ = p (3 + _n + 4 + 4 + 1);
      ublas::vector<double> thetaHatY (4);
      ublas::vector<double> gammaHatY (4);
      for (int i = 0; i < 4; ++i){
        thetaHatY (i) = p (3 + _n + 4 + 4 + 2 + i);
        gammaHatY (i) = p (3 + _n + 4 + 4 + 2 + 4 + i);
      }
      ublas::vector<double> thetaHatZ (4);
      ublas::vector<double> gammaHatZ (4);
      for (int i = 0; i < 4; ++i){
        thetaHatZ (i) = p (3 + _n + 4 + 4 + 2 + 8 + i);
        gammaHatZ (i) = p (3 + _n + 4 + 4 + 2 + 8 + 4 + i);
      }
      ublas::vector<double> normSampY (_J);
      for (int i = 0; i < _J; ++i){
        normSampY (i) = p (3 + _n + 4 + 4 + 2 + 8 + 8 + i);
      }
      ublas::vector<double> normSampZ (_J);
      for (int i = 0; i < _J; ++i){
        normSampZ (i) = p (3 + _n + 4 + 4 + 2 + 8 + 8 + _J + i);
      }
      ublas::vector<double> qtsY (_lq2);
      for (int i = 0; i < _lq2; ++i){
        qtsY (i) = p (3 + _n + 4 + 4 + 2 + 8 + 8 + 2*_J + i);
      }
      ublas::vector<double> qtsZ (_lq2);
      for (int i = 0; i < _lq2; ++i){
        qtsZ (i) = p (3 + _n + 4 + 4 + 2 + 8 + 8 + 2*_J + _lq2 + i); 
      }
      int _N = static_cast<int>(p (3 + _n + 8 + 2 + 16 + 2*_J +
       2*_lq2)); 
      double lambda = static_cast<double>(p (3 + _n + 8 + 2 + 16 + 2*_J + 2*_lq2 + 1));
      double kappa = static_cast<double>(p (3 + _n + 8 + 2 + 16 + 2*_J + 2*_lq2 + 2));
      ublas::vector<double> unifY (_N);
      for (int i = 0; i < _N; ++i){
        unifY (i) = p (3 + _n + 4 + 4 + 2 + 8 + 8 + 2*_J + 2*_lq2 + 3 + i); 
      }
      ublas::vector<double> unifZ (_N);
      for (int i = 0; i < _N; ++i){
        unifZ (i) = p (3 + _n + 4 + 4 + 2 + 8 + 8 + 2*_J + 2*_lq2 + 3 + _N + i); 
      }

      double norm1 = sqrt (eta10*eta10 + eta11*eta11); 
      double norm2 = sqrt (eta20*eta20 + eta21*eta21); 
      eta10 /= norm1;
      eta11 /= norm1;
      eta20 /= norm2;
      eta21 /= norm2;

    /* Evaluate the distribution function of Y under (eta1, eta2)
       regime */  
      ublas::vector<double> _prbsY (_lq2);
      double cdfValY = 0.0;
      for (int i = 0; i < _lq2; ++i){
        dist::distEval (H1, bHatY20, bHatY21, sdY, thetaHatY, gammaHatY, eta10, eta11, eta20, eta21, normSampY, qtsY (i),
          cdfValY);   
        _prbsY (i) = cdfValY;
      }

    /* Get interval probabilities and the quantile midpoints of those
       intervals */  
      int _lq1 = _lq2 - 1;
      int _lq3 = _lq2 + 1;
      ublas::vector<double> _midsY (_lq3);
      double _stepsY = qtsY (1) - qtsY (0);
      _midsY (0) = qtsY (0) - _stepsY/2.0; 
      _midsY (_lq2) = qtsY (_lq1) + _stepsY/2.0;
      for (int i = 0; i < _lq1; ++i){
        _midsY (i + 1) = (qtsY (i) + qtsY (i + 1)) / 2.0;
      }

      ublas::vector<int> _multisampY (_lq3);
      _multisampY.assign (ublas::zero_vector<int> (_lq3)); 
      for (int j = 0; j < _N; ++j){
        if (unifY (j) <= _prbsY (0)){
          _multisampY (0) += 1;
        }
      }
      for (int i = 1; i < _lq2; ++i){
        for (int j = 0; j < _N; ++j){
          if (unifY (j) <= _prbsY (i) && unifY (j) > _prbsY (i - 1)){
            _multisampY (i) += 1;
          }
        }
      }
      for (int j = 0; j < _N; ++j){
        if (unifY (j) > _prbsY (_lq1)){
          _multisampY (_lq2) += 1;
        }
      }

      int K = 0;
      int incr = 0; 
      ublas::vector<double> _fnlYsamp (_N);
      for (int i = 0; i < _lq2; ++i){
        incr = K + _multisampY (i);
        for (int j = K; j < incr; ++j){
          _fnlYsamp (j) = _midsY (i); 
        }
        K += _multisampY (i);
      }

    /* _fnlYsamp is a sample from the distribution function of Y under
       regime (eta1, eta2) */
      double meanY = 0.0;
      for (int i = 0; i < _N; ++i){
        meanY += _fnlYsamp (i);
      }

      meanY /= 1.0*static_cast<double>(_N);

    /* Evaluate the distribution function of Z under (eta1, eta2)
       regime */  
      ublas::vector<double> _prbsZ (_lq2);
      double cdfValZ = 0.0;
      for (int i = 0; i < _lq2; ++i){
        dist::distEval (H1, bHatZ20, bHatZ21, sdZ, thetaHatZ, gammaHatZ, eta10, eta11, eta20, eta21, normSampZ, 
          qtsZ (i), cdfValZ);   
        _prbsZ (i) = cdfValZ;
      }

    /* Get interval probabilities and the quantile midpoints of those
       intervals */  
      ublas::vector<double> _midsZ (_lq3);
      double _stepsZ = qtsZ (1) - qtsZ (0);
      _midsZ (0) = qtsZ (0) - _stepsZ/2.0; 
      _midsZ (_lq2) = qtsZ (_lq1) + _stepsZ/2.0;
      for (int i = 0; i < _lq1; ++i){
        _midsZ (i + 1) = (qtsZ (i) + qtsZ (i + 1)) / 2.0;
      }

      ublas::vector<int> _multisampZ (_lq3);
      _multisampZ.assign (ublas::zero_vector<int> (_lq3)); 
      for (int j = 0; j < _N; ++j){
        if (unifZ (j) <= _prbsZ (0)){
          _multisampZ (0) += 1;
        }
      }
      for (int i = 1; i < _lq2; ++i){
        for (int j = 0; j < _N; ++j){
          if (unifZ (j) <= _prbsZ (i) && unifZ (j) > _prbsZ (i - 1)){
            _multisampZ (i) += 1;
          }
        }
      }
      for (int j = 0; j < _N; ++j){
        if (unifZ (j) > _prbsZ (_lq1)){
          _multisampZ (_lq2) += 1;
        }
      }

      K = 0;
      incr = 0; 
      ublas::vector<double> _fnlZsamp (_N);
      for (int i = 0; i < _lq2; ++i){
        incr = K + _multisampZ (i);
        for (int j = K; j < incr; ++j){
          _fnlZsamp (j) = _midsZ (i); 
        }
        K += _multisampZ (i);
      }

    /* _fnlZsamp is a sample from the distribution function of Z under
       regime (eta1, eta2) */
      double meanZ = 0.0;
      for (int i = 0; i < _N; ++i){
        meanZ += _fnlZsamp (i);
      }

      meanZ /= 1.0*static_cast<double>(_N);

      double barrier = 0.0;
      if (kappa - meanZ <= 0)
      {
        barrier = std::numeric_limits<double>::max();
      }
      else
      { 
        barrier = -1.0*meanY - lambda*log (kappa - meanZ);
      }

      return barrier;
    }
  } // end namespace

  namespace gslmainB {
    namespace ublas = boost::numeric::ublas;
    int mainB (const ublas::vector<double>& mypars,
      ublas::vector<double>& etaHat1,
      ublas::vector<double>& etaHat2,
      double& minEst){

    // declare the minimizer algorithm to use
    // Nelder Mead is the only one that doesn't require gradient
      const gsl_multimin_fminimizer_type *T = 
      gsl_multimin_fminimizer_nmsimplex;

    // this is a pointer to a class that will have everything needed
    // to optimize 
      gsl_multimin_fminimizer *s = NULL;

    // ss is a pointer to a gsl_vector of step sizes
      gsl_vector *ss;

    // x is a pointer to a gsl_vector of variables to minimize the function
      gsl_vector *x;

    // the class that holds all information for the objective function
      gsl_multimin_function minex_func;

      size_t iter = 0; 
    // number of iterations
      int status; 
    // convergence status
      double size; 
    // a distance for how far the minimum is away from
    // current iter 
      int dim=4; 
    // dimension of input to objective function

    // set initial values
      x=gsl_vector_alloc(dim);
      gsl_vector_set(x,0,etaHat1 (0));
      gsl_vector_set(x,1,etaHat1 (1));
      gsl_vector_set(x,2,etaHat2 (0));
      gsl_vector_set(x,3,etaHat2 (1));

    // set initial step sizes (these adpatively change with this algorithm)
      ss=gsl_vector_alloc(dim);
      gsl_vector_set_all(ss,3.0);

      ublas::vector<double> pvec (mypars.size ());
      for (int i = 0; i < mypars.size (); ++i){
        pvec (i) = mypars (i);
      }

      minex_func.n=dim; 
    // number of function components
      minex_func.f=&gslfnlB::gslFnB; 
    // functor
      minex_func.params=&pvec; 
    // initial values

    // allocate memory for the minimizer class
      s=gsl_multimin_fminimizer_alloc(T,dim);
      gsl_multimin_fminimizer_set(s,&minex_func,x,ss);

    // minimize
      do {
        iter++;

      // to iterate once its a function call
      // the return value is the status of conergence
      // its up to the user to decide when to stop
      // if this function returns true, then there was an error
        status=gsl_multimin_fminimizer_iterate(s);
        if(status)
         break;

      // not sure what exactly this is, but its used to evaluate convergence
      // its some sort of estimate as to how far the algorithm estimates the
      // minimum to be from the current iteration
       size=gsl_multimin_fminimizer_size(s);

      // this actually tests for convergence.  1e-2 is the tolerance
       status=gsl_multimin_test_size(size,1e-2);

      // GSL_SUCCESS is a global constant within the GSL header files
      //     if(status == GSL_SUCCESS)
      //{
      //  printf("converged to minimium at\n");
      //}

       etaHat1 (0) = static_cast<double>(gsl_vector_get(s->x,0));
       etaHat1 (1) = static_cast<double>(gsl_vector_get(s->x,1));
       etaHat2 (0) = static_cast<double>(gsl_vector_get(s->x,2));
       etaHat2 (1) = static_cast<double>(gsl_vector_get(s->x,3));
       minEst = s->fval;

     } 

    // GSL_CONTINUE is a global constant within the GSL header files
     while(status == GSL_CONTINUE && iter < 100);

    // free the memory allocated above
     gsl_vector_free(x);
     gsl_vector_free(ss);
     gsl_multimin_fminimizer_free(s);

     return 0;
   }
  } // end namespace
  
  # endif
  
  