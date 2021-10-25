/**
 * @file linear_diff_ops.hpp
 * @author Carlo Musolino (musolino@itp.uni-frankfurt.de)
 * @brief Implementation of linear differential operator classes
 */
#ifndef _MY_LINEAR_OPERATORS_HPP
#define _MY_LINEAR_OPERATORS_HPP

#include "../polybases/polybases.hpp"
#include <algorithm>
#include <functional>
#include <iostream>

using namespace FunctionalBases;

namespace Operators {

    template<class T>
    class LinearOperator {
        unsigned int N;
        std::vector<T> L_ij;
        public:
        LinearOperator<T>(unsigned int N=0): N(N){} //! Constructor
        void print_Lij() {
            int imax = N;
            for (int i=0;i<imax+1;i++){
                for (int j=0;j<imax+1;j++){
                    int ind2d = j + (N+1)*i;
                    std::cout << L_ij[ind2d] << "\t"; 
                }
                std::cout << "\n";
            }
        }
        void print_N() { std::cout << N << "\n"; }
        void set_Lij(const std::vector<T>& Lij){
            L_ij = Lij;
        }
        void set_N(const int n){
            N = n;
        }
        LinearOperator<T>& operator+=(const LinearOperator<T>& rhs) 
        {                      
            //std::transform(L_ij.begin(),L_ij.end(),rhs.L_ij.begin(),rhs.L_ij.end(),std::back_inserter(tmp),std::plus<T>());
            for(int i=0; i<L_ij.size();i++){
                L_ij[i] += rhs.L_ij[i];
            }
            return *this; 
        }
        const LinearOperator<T> operator+(const LinearOperator<T>& rhs) {
            return LinearOperator<T>(*this) += rhs;
        }
        LinearOperator<T>& operator*=(const T& rhs) {
             for(auto& val : L_ij) val *= rhs;
             return *this;
         }

        const LinearOperator<T> operator*(const T&& alpha) {
            return LinearOperator<T>(*this) *= alpha;
        }
    };

    template<class T>
    class Derivative: public LinearOperator<T> {
        FunctionalBase<T>* basis;
        public:
        Derivative<T>(FunctionalBase<T>* b): basis(b), LinearOperator<T>() {
            std::vector<T> L_ij ;
            basis->calc_deriv(L_ij);
            LinearOperator<T>::set_N(basis->get_N());
            LinearOperator<T>::set_Lij(L_ij);
        }
        
    };

    
    template<class T>
    class SecondDerivative: public LinearOperator<T> {
        FunctionalBase<T>* basis;
        public:
        SecondDerivative<T>(FunctionalBase<T>* b): basis(b), LinearOperator<T>() {
            std::vector<T> L_ij ;
            basis->calc_second_deriv(L_ij);
            LinearOperator<T>::set_N(basis->get_N());
            LinearOperator<T>::set_Lij(L_ij);
        }
        
    }; 

    
    template<class T>
    class TimesX: public LinearOperator<T> {

        FunctionalBase<T>* basis;
        public:
        TimesX<T>(FunctionalBase<T>* b): basis(b), LinearOperator<T>() {
            std::vector<T> L_ij ;
            basis->calc_times_x(L_ij);
            LinearOperator<T>::set_N(basis->get_N());
            LinearOperator<T>::set_Lij(L_ij);
        }
    };
};

#endif