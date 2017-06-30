/* @file genDist.hpp
 */
#ifndef _GENDIST_HPP
#define _GENDIST_HPP

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

  /*
  @fn ranNorm 
  @brief Generates a vector of J iid N(0, 1) samples
  */

  namespace genNorm{
    namespace ublas = boost::numeric::ublas;
    bool ranNorm (const int& J, 
      const gsl_rng*& r,
      ublas::vector<double>& samp){

    /* create gsl_vector types to store the randomly generated numbers */
      gsl_vector *sample = gsl_vector_alloc (J);

      for (int i = 0; i < J; ++i){
        gsl_vector_set (sample, i, gsl_ran_gaussian (r, 1.0));
        samp (i) = static_cast<double>(gsl_vector_get (sample, i));
      }  

    // free the memory allocated to gsl_vectors
      gsl_vector_free (sample);

      return true;
    }
  } // end namespace


  namespace genUnif{
    namespace ublas = boost::numeric::ublas;
    bool ranUnif (const int& _N,
      const gsl_rng*& r,
      ublas::vector<double>& unif){

      for (int i = 0; i < _N; ++i){
        double u = gsl_rng_uniform (r);
        unif (i) = u;
      }  

      return true;
    }
  } // end namespace

  /*  
    @fn ranMultinom
    @brief Generates _n iid samples from multinomial distribution with probs 
  */

    namespace genMultinom{
      namespace ublas = boost::numeric::ublas;
      bool ranMultinom (const ublas::vector<double>& probs, 
        const int _n,  
        const gsl_rng*& r,
        ublas::vector<int>& samp){

    // Number of groups
        size_t _G = probs.size ();
        unsigned int _ss = _n;
        unsigned int sample[_ss];

    // Convert ublas probs to std::vector
        double reprobs[_G];
        for (int i = 0; i < _G; ++i){
          reprobs[i] = probs (i);
        }

    /* create gsl_vector types to store the randomly generated numbers */
    //    gsl_vector *_samp = gsl_vector_alloc (_G);
        gsl_ran_multinomial (r, _G, _ss, reprobs, sample);

        std::cout << sample << " This is the multinom sample \n"; 

        for (int i = 0; i < _G; ++i){
          samp (i) = sample [i];
        }

        return true;
      }
  } // end namespace




  #endif
