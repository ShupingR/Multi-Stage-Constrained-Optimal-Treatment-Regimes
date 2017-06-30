  /* @file genModel.hpp */

  #ifndef _GENMODEL_HPP
  #define _GENMODEL_HPP
  
  #include <math.h>
  #include <ctime>
  #include <boost/numeric/ublas/vector.hpp>
  #include <boost/numeric/ublas/vector_proxy.hpp>
  #include <boost/numeric/ublas/matrix.hpp>
  #include <boost/numeric/ublas/matrix_proxy.hpp>
  #include <gsl/gsl_randist.h>
  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_rng.h>
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_blas.h>
  
  /* @fn generate
   @brief Function to generate training sample of size X1.size () from
          our generative model
   @param r gsl random number generator       
   @param mu mean of X1, X1 i.i.d N (mu, 1);
   @param beta10 main effect parameter vector for the mean of X2
   @param beta11 interaction effect parameter vector for the mean of
             X2: E (X2 \mid H1, A1) = H1^{T}\beta10 + A1*H1^{T}\beta11 
   @param gamma0 main effect parameter for the variance function
   @param gamma1 interaction effect parameter for the var function:
                 var = exp (H1^{T}\gamma0 + A1*H1^{T}\gamma1);
   @param betaY20 main effect parameter for mean of Y
   @param betaY21 interaction effect parameter for mean of Y:
             E (Y \mid H2, A2) = H2^{T}\betaY20 + A2*H2^{T}\betaY21; 
   @param betaZ20 main effect parameter for mean of Z
   @param betaZ21 interaction effect parameter for mean of Z
             E (Z \mid H2, A2) = H2^{T}\betaZ20 + A2*H2^{T}\betaZ21; 
   @param sqrtSigma square root of the covariance matrix of \epsilon_Y
                    and \epsilon_Z
   @return Training data set: X1, A1, X2, A2, H1, H2, Y, Z
 */

   namespace genModel{
    namespace ublas = boost::numeric::ublas;
    bool generate (const gsl_rng* r,
      const int& _n,
      const ublas::vector<double>& mu, 
      const ublas::vector<double>& beta10, 
      const ublas::vector<double>& beta11, 
      const ublas::vector<double>& gamma0, 
      const ublas::vector<double>& gamma1, 
      const ublas::vector<double>& betaY20,  
      const ublas::vector<double>& betaY21,  
      const ublas::vector<double>& betaZ20,  
      const ublas::vector<double>& betaZ21,  
      const ublas::matrix<double>& sqrtSigma, 
      ublas::vector<double>& X1,
      ublas::vector<int>& A1,
      ublas::vector<double>& xi,
      ublas::vector<double>& X2,
      ublas::vector<int>& A2,
      ublas::matrix<double>& H1,
      ublas::matrix<double>& H2,
      ublas::vector<double>& epsY,
      ublas::vector<double>& epsZ,
      ublas::vector<double>& Y,
      ublas::vector<double>& Z){

      ublas::vector<double> epsYZ (2);
    /* create gsl_vector types to store the randomly generated numbers */
      gsl_vector *_x1Ran = gsl_vector_alloc (_n);
      gsl_vector *_a1Ran = gsl_vector_alloc (_n);
      gsl_vector *_a2Ran = gsl_vector_alloc (_n);
      gsl_vector *_xiRan = gsl_vector_alloc (_n);
      gsl_vector *_eYRan = gsl_vector_alloc (_n);
      gsl_vector *_eZRan = gsl_vector_alloc (_n);

      ublas::vector<double> eps (2);
      eps.assign (ublas::vector<double> (2));

      for (int i = 0; i < _n; ++i){
      H1 (i, 0) = 1;
      gsl_vector_set (_x1Ran, i, gsl_ran_gaussian (r, 1.0));
      H1 (i, 1) = mu[0] + static_cast<double>(gsl_vector_get (_x1Ran, i));   
      X1 (i) = H1 (i, 1);

      gsl_vector_set (_a1Ran, i, gsl_ran_gaussian (r, 1.0));
      gsl_vector_set (_a2Ran, i, gsl_ran_gaussian (r, 1.0));
      double txt1 = static_cast<double>(gsl_vector_get (_a1Ran, i));
      double txt2 = static_cast<double>(gsl_vector_get (_a2Ran, i));
      A1 (i) = static_cast<int>((txt1 >= 0) - (txt1 < 0));
      A2 (i) = static_cast<int>((txt2 >= 0) - (txt2 < 0));

      H2 (i, 0) = 1;       
      gsl_vector_set (_xiRan, i, gsl_ran_gaussian (r, 1.0));         
      xi (i) = static_cast<double>(gsl_vector_get (_xiRan, i));       
      ublas::matrix_row<ublas::matrix<double> > mr1 (H1, i);      

      double eta = sqrt (exp (inner_prod (mr1, gamma0) + 
        A1 (i)*inner_prod (mr1, gamma1)));       
      X2 (i) = inner_prod (mr1, beta10) + 
        A1 (i)*inner_prod (mr1, beta11) + eta*xi (i); 

      H2 (i, 1) = X2 (i); 

      gsl_vector_set (_eYRan, i, gsl_ran_gaussian (r, 1.0));
      gsl_vector_set (_eZRan, i, gsl_ran_gaussian (r, 1.0));
      epsYZ (0) = static_cast<double>(gsl_vector_get (_eYRan, i));
      epsYZ (1) = static_cast<double>(gsl_vector_get (_eZRan, i));
      eps = prod (sqrtSigma, epsYZ);
      epsY (i) = eps (0);
      epsZ (i) = eps (1);

      ublas::matrix_row<ublas::matrix<double> > mr2 (H2, i);
      Y (i) = inner_prod (mr2, betaY20) + 
      A2 (i)*inner_prod (mr2, betaY21) + eps (0);
      Z (i) = inner_prod (mr2, betaZ20) + 
      A2 (i)*inner_prod (mr2, betaZ21) + eps (1);
    } 

    // free the memory allocated to r and all gsl_vectors
    gsl_vector_free (_x1Ran);
    gsl_vector_free (_a1Ran); 
    gsl_vector_free (_a2Ran); 
    gsl_vector_free (_xiRan); 
    gsl_vector_free (_eYRan); 
    gsl_vector_free (_eZRan); 

    return true;
  }
} // end namespace

   namespace genEval{
    namespace ublas = boost::numeric::ublas;
    bool genOpt (const ublas::vector<double>& mu, 
      const ublas::vector<double>& beta10, 
      const ublas::vector<double>& beta11, 
      const ublas::vector<double>& gamma0, 
      const ublas::vector<double>& gamma1, 
      const ublas::vector<double>& betaY20,  
      const ublas::vector<double>& betaY21,  
      const ublas::vector<double>& betaZ20,  
      const ublas::vector<double>& betaZ21,  
      const ublas::matrix<double>& sqrtSigma, 
      const ublas::vector<double>& X1,
      const ublas::vector<double>& xi,
      const ublas::vector<double>& epsY,
      const ublas::vector<double>& epsZ,
      const ublas::vector<double>& eta1,
      const ublas::vector<double>& eta2,
      ublas::vector<double>& Yeval,
      ublas::vector<double>& Zeval){

      int _n = X1.size ();
      ublas::matrix<double> H1 (_n, 2);
      ublas::matrix<double> H2 (_n, 2);
      ublas::vector<double> A1 (_n);
      ublas::vector<double> A2 (_n);
      ublas::vector<double> X2 (_n);
      
      for (int i = 0; i < _n; ++i){
      H1 (i, 0) = 1;
      H1 (i, 1) = X1 (i);   

      double txt1 = eta1 (0) + eta1 (1)*X1 (i);
      A1 (i) = static_cast<int>((txt1 >= 0) - (txt1 < 0));

      H2 (i, 0) = 1;       
      ublas::matrix_row<ublas::matrix<double> > mr1 (H1, i);      
      double eta = sqrt (exp (inner_prod (mr1, gamma0) + 
        A1 (i)*inner_prod (mr1, gamma1)));       
      X2 (i) = inner_prod (mr1, beta10) + 
        A1 (i)*inner_prod (mr1, beta11) + eta*xi (i); 
      H2 (i, 1) = X2 (i); 
      double txt2 = eta2 (0) + eta2 (1)*X2 (i);
      A2 (i) = static_cast<int>((txt2 >= 0) - (txt2 < 0));

      ublas::matrix_row<ublas::matrix<double> > mr2 (H2, i);
      Yeval (i) = inner_prod (mr2, betaY20) + A2 (i)*inner_prod (mr2, betaY21) + epsY (i);
      Zeval (i) = inner_prod (mr2, betaZ20) + A2 (i)*inner_prod (mr2, betaZ21) + epsZ (i);
    } 

    return true;
  }
} // end namespace

#endif
