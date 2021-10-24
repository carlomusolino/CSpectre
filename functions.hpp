/**
 * @file functions.hpp
 * @brief This file contains source code for the Function class
 * @author Carlo Musolino
 * @but No known bugs
 */

#ifndef _MY_FUNCS_H
#define _MY_FUNCS_H

#include "polybases/polybases.hpp"

namespace FunctionalBases {
/** 
 * @brief namespace that contains the function class 
 * The function class is meant as a wrapper to contain
 * physical and collocation space representations of a
 * function. It contains a pointer to a FunctionalBase
 * object which is used by the decompose() method to 
 * compute the spectral coefficients.
 */
namespace Functions {
  /**
   * @brief Function class
   */
    template <class T, class FuncBase, unsigned int N>
    class Function {
      std::vector<T> f_i;   //! Physical space representation of f
        std::vector<T> ft_i;  //! Collocation space representation of f
        FuncBase* basis; //! FunctionalBase used for spectral decomposition
        public:
      // constructors -----------------
      /** 
       * @brief Constructor based on function values at grid nodes
       * @param f_i vector containing function values at grid nodes
       */
        Function<T,FuncBase,N>(const std::vector<T>& f_i): f_i(f_i) {
            assert(f_i.size()==N);
            basis = new FuncBase(N);
            decompose();
        };
      /**
       * @brief Constructor based on function pointer 
       * @param func analytic function, will be evaluated on the grid and decomposed
       */
        Function<T,FuncBase,N>(void func(const std::vector<T>&, std::vector<T>&)) {
            basis = new FuncBase(N) ;
            std::vector<T> n;
            basis->get_nodes(n);
            (*func)(n,f_i);
            decompose();
        }
      // members --------------------
      /**
       * @brief evaluate function at a point
       * @param x point at which to evaluate f
       * @param f_x reference to output value
       */
      inline void eval(const T& x, T& f_x);
      /**
       * @brief evaluate function at a vector of points
       * @param x points at which to evaluate f
       * @param f_x reference to output vector
       */
      inline void eval(const std::vector<T>& x, std::vector<T>& f_x);
      /**
       * @brief perform spectral decomposition of f
       * Use FunctionalBasis* member to compute the spectral 
       * decomposition of f. Depends on detailed implementation of 
       * FunctionalBasis abstract class.
       */
      inline void decompose(){
            basis->calc_spectral_coeffs(f_i,ft_i);
        }
        // access---------------------
      /**
       * @brief access to spectral coefficients
       * @param y reference to output vector
       */
        inline void get_spectral_coeffs(std::vector<T>& y){ y = ft_i;};
      /**
       * @brief access to function values on grid nodes
       * @param y reference to output vector
       */
        inline void get_func_vals(std::vector<T>& y){ y = f_i; };
        // destructor -------------
      /**
       * @brief Destructor. Clean up basis member.
       */
        ~Function<T,FuncBase,N>() { delete basis;};
    };


template <class T,class FuncBase, unsigned int N> inline void Function<T,FuncBase,N>::eval(const T& x, T& f_x)
{   
    for(int n=0; n<ft_i.size(); n++) f_x += ft_i[n] * basis->evaluate_function(x,n);
};

template <class T,class FuncBase, unsigned int N> inline void Function<T,FuncBase,N>::eval(const std::vector<T>& x, std::vector<T>& f_x){
    f_x.clear();
    for (const auto& xval: x) {
        T tmp = static_cast<T>(0);
        for(int n=0; n<ft_i.size(); n++) tmp += ft_i[n] * basis->evaluate_function(xval,n);
        f_x.push_back(tmp);
    }
};
            
}
}
