/* @file ols.hpp
 */

#ifndef _OLS_HPP
#define _OLS_HPP

#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

/* 
   @fn sdCmpt
   @brief Function to calculate the standard deviation of a vector 
   @param myvec numeric vector 
   @return sd standard deviation of myvec
 */
namespace sdev{
  namespace ublas = boost::numeric::ublas;
  bool sdCmpt (const ublas::vector<double>& myvec,
	       double& sd){

    int _d1 = myvec.size ();
    double _mean = sum (myvec); 
    _mean /= static_cast<double>(_d1);

    double _sdSum = 0.0;
    for (int i = 0; i < _d1; ++i){
      _sdSum += (myvec (i) - _mean) * (myvec (i) - _mean);  
    }
    _sdSum /= static_cast<double>(_d1 - 1);

    sd = sqrt (_sdSum);

    return true;
  }
} // end namespace

/* 
   @fn lm (const ublas::matrix<double>& x, const ublas::vector<double>& y,
           ublas::vector<double>& bHat)
   @brief Uses least squares to fit y ~ x^{\T}\beta + eps.  The fitting is 
          done using LU decomposition.
   @param x n x p matrix of doubles, the design matrix
   @param y n vector of doubles, responses
   @param bHat return vector to store the regression coefficients
   @return Indicator of success.
*/

namespace ols {
  namespace ublas = boost::numeric::ublas;
  using namespace std;
  bool lm (const ublas::matrix<double>& x, 
	   const ublas::vector<double>& y, 
	   ublas::vector<double>& bHat,
	   ublas::vector<double>& resids,
	   double& sdRes){
    // We will use the normal equations x'xb = x'y to obtain the 
    // ols coefficients
    ublas::matrix<double> xtx = ublas::prod (trans (x), x); // declare x'x
    ublas::matrix<double> xtxInv (x.size2(), x.size2()); 
    // will store (x'x)^-1
    xtxInv.assign (ublas::identity_matrix<double>(x.size2()));

    // To solve we write x'x = LU and solve the two triangular systems
    // Lz = x'y to obtain z then Ub = z
    ublas::permutation_matrix<std::size_t> pm (xtx.size1());
    int res = lu_factorize (xtx, pm);
    if (res != 0){std::cout << "Inversion failed.\n";}
    lu_substitute (xtx, pm, xtxInv);

    ublas::vector<double> xty = ublas::prod (trans (x), y); // x'y 
    bHat.resize (x.size2()); 
    bHat = ublas::prod (xtxInv, xty); // (x'x)^{-1}x'y 

    /* Compute fitted values and residuals from the regression */  
    ublas::vector<double> fitted = ublas::prod (x, bHat);
    resids = y - fitted;
    sdev::sdCmpt (resids, sdRes);

    return true; 
  }
}


#endif
