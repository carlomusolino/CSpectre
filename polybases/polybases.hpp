/**
 * @file polybases.hpp
 * @brief Implementation of FunctionalBase classes.
 * @author Carlo Musolino (musolino@itp.uni-frankfurt.de)
 * Implementation of abstract class FunctionalBase and subclasses. These
 * deal with the definition of collocation nodes, interpolation domain and
 * calculation of the spectral coefficients for function decomposition.
 */
#ifndef _FUNCBASES_HPP_
#define _FUNCBASES_HPP_


#include <cmath>
#include <vector>
#include <assert.h>
#include <iostream>
#include "chebyshev.hpp"

namespace FunctionalBases {

  /**
   * @brief Abstract class 
   */
  template <class T>
  class FunctionalBase {
    unsigned int N;
  public:
    //constructor
    FunctionalBase<T>(unsigned int N) : N(N) {};
    // virtual member functions 
    virtual inline T evaluate_function(const T& x, const unsigned int n){};
    virtual void calc_spectral_coeffs(const std::vector<T>& f,std::vector<T>& ftilde){};
    virtual void calc_function_values(const std::vector<T>& ftilde, std::vector<T>& f){};
    virtual inline void calc_deriv(std::vector<T>& Lij) {};
    virtual inline void calc_second_deriv(std::vector<T>& Lij) {};
    virtual inline void calc_times_x(std::vector<T>& Lij) {};
    // access 
    virtual inline void get_nodes(std::vector<T>& pts)  {};
    virtual inline void get_weights(std::vector<T>& w)  {};
    virtual inline int get_N() { return N; };
    virtual inline void print_nodes() {};
    virtual inline void print_weights() {};
    // destructor
    virtual ~FunctionalBase<T>() {} ;
  };

  /**
   * @brief Implementation of Chebyshev polynomial basis.
   * Templated subclass of abstract class FunctionalBases implementing a Chebyshev polynomial basis. 
   */
  template <class T>
  class ChebyshevBase: public FunctionalBase<T> {
    unsigned int N; //! Order
    std::vector<T> nodes; //! Gauss-Lobatto collocation points
    std::vector<T> weights; //! Gaussian quadrature weights 
  public:
    // constructor ----------------------
    /**
     * @brief Constructor
     * @param N order of the polynomial basis.
     */
    ChebyshevBase<T>(unsigned int N) : N(N), FunctionalBase<T>(N){
      calc_nodes_and_weights();
    };
    // class methods ---------------------
    /**
     * @brief Initialise value of Gauss-Lobatto nodes and Gaussian quadrature weights.
     */
    inline void calc_nodes_and_weights();
    /**
     * @brief Wrapper for nth-order polynomial evaluation. Routines implemented in chebyshev.hpp
     * @param x point of evaluation
     * @param n order of the polynomial
     */
    inline T evaluate_function(const T& x, const unsigned int n)  {
      return Chebyshev::Tn<T>(x,n);
    };
    /**
     * @brief Calculate spectral coefficients of a function given its values at collocation points 
     * @param f values of f at collocation points
     * @param ftilde output vector
     */
    inline void calc_spectral_coeffs(const std::vector<T>& f,std::vector<T>& ftilde) ;
    inline void calc_function_values(const std::vector<T>& ftilde, std::vector<T>& f) ;
    inline void calc_deriv(std::vector<T>& Lij);
    inline void calc_second_deriv(std::vector<T>& Lij);
    inline void calc_times_x(std::vector<T>& Lij);
    // access ----------------
    inline void print_nodes() {
      std::cout << "Length of nodes vector: " << nodes.size() << "\n";
      for (const auto& val : nodes) std::cout << val << " "; 
      std::cout << "\n";
    }
    inline void print_weights() {
      std::cout << "Length of weights vector: " << weights.size() << "\n";
      for (const auto& val : weights) std::cout << val << " "; 
      std::cout << "\n";
    }
    //! Return the nodes
    inline void get_nodes(std::vector<T>& pts)  {
      pts = nodes;
    };
    //! Return the weights
    inline void get_weights(std::vector<T>& w)  {
      w = weights;
    };
    inline void change_N(const unsigned int n) {
      N = n;
      calc_nodes_and_weights();
    }
    inline int get_N() {
      return N;
    }
    // destructor -----------------
    //! Default destructor
    ~ChebyshevBase<T>() {};
  };

  template <class T> inline void ChebyshevBase<T>::calc_nodes_and_weights()
    {
      weights.clear();
      nodes.clear();
      weights.push_back(static_cast<T>(M_PI/(2*N)));
	    nodes.push_back(static_cast<T>(1));
	    for (int i=1; i<N; i++){
        nodes.push_back(static_cast<T>(std::cos(M_PI*i/N)));    
        weights.push_back(static_cast<T>(M_PI/N));
	    }
	    nodes.push_back(static_cast<T>(-1));
	    weights.push_back(static_cast<T>(M_PI/(2*N)));
    }

  template <class T>  inline void ChebyshevBase<T>::calc_spectral_coeffs(const std::vector<T>& f,std::vector<T>& ftilde){
    assert(f.size()==N+1);
    ftilde.clear();
    T gamma_n, ft_n;
    for(int n=0; n<N+1; n++){
      gamma_n = static_cast<T>(0);
      ft_n = static_cast<T>(0);
      for(int i=0; i < N+1; i++){
          gamma_n += std::pow(evaluate_function(nodes[i],n),2) * weights[i];
          ft_n += f[i] * evaluate_function(nodes[i],n) * weights[i];
      }
      ftilde.push_back(ft_n / gamma_n) ;
    }  
  };

  template <class T>  inline void ChebyshevBase<T>::calc_function_values(const std::vector<T>& ftilde, std::vector<T>& f) {
    assert(ftilde.size()==N+1);
    f.clear();
    for (const auto& xval: nodes) {
        T tmp = static_cast<T>(0);
        for(int n=0; n<N+1; n++) tmp += ftilde[n] * evaluate_function(xval,n);
        f.push_back(tmp);
    }
  };

  template <class T>  inline void ChebyshevBase<T>::calc_deriv(std::vector<T>& Lij)
  {
    Lij.clear();
    for (int i=0; i<N+1; i++){
      for (int j=0; j<N+1; j++){
        T val = static_cast<T>(0);
        double delta_0n = (i==0) ? 1.0 : 0.0;
        if (((i+j)%2)&&(j>i)) val = static_cast<T>( 2.0/(1.0 + delta_0n) * j );
        Lij.push_back(val);
      }
    }
  }

  template <class T>  inline void ChebyshevBase<T>::calc_second_deriv(std::vector<T>& Lij)
  {
    Lij.clear();
    for (int i=0; i<N+1; i++){
      for (int j=0; j<N+1; j++){
        T val = static_cast<T>(0);
        double delta_0n = (i==0) ? 1.0 : 0.0;
        if ( (!((i+j)%2)) && (j>(i+1)) ) val = static_cast<T>( 1.0/(1.0 + delta_0n) * j * ( j * j - i * i));
        Lij.push_back(val);
      }
    }
  }

  template <class T>  inline void ChebyshevBase<T>::calc_times_x(std::vector<T>& Lij)
  {
    Lij.clear();
    for (int i=0; i<N+1; i++){
      for (int j=0; j<N+1; j++){
        T val = static_cast<T>(0);
        int im1 = i-1;
        double delta_0nm1 = (im1==0) ? 1.0 : 0.0;
        if (j==(i-1)){ 
          val = static_cast<T>( 0.5*(1.0+delta_0nm1) ); 
          }
        else if (j==(i+1)){
           val = static_cast<T>(0.5); 
           }
        Lij.push_back(val);
      }
    }
  }

} //namespace FunctionalBases

#endif 


