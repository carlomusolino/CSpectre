#include <vector>
#include "polybases/polybases.hpp"

using namespace FunctionalBases;

namespace Functions {

    template <class T, class FuncBase>
    class Function {
        std::vector<T> f_i;   // physical space representation of f
        std::vector<T> ft_i;  // collocation space representation of f
        unsigned int N;
        B basis; // should be a derived class of FunctionBase
        public:
        // creators
        Function<T,B>(const std::vector<T>& f_i): f_i(f_i), N(f_i.size()), basis(N) {
            interpolate();
        };
        Function<T,B>(const unsigned int N, void * func): N(N), basis(N) {
            std::vector<T> n;
            basis->get_nodes(n);
            func(n,f_i);
            decompose();
        }
        // members
        inline void eval(const T& x, T& f_x);
        inline void eval(const std::vector<T>& x, std::vector<T>& f_x);
        inline void decompose(){
            basis->calc_spectral_coeffs(f_i,ft_i);
        }
        // access
        inline void get_spectral_coeffs();
        inline void get_func_vals();
        // destructor
        ~Function<T,B>() {};
    };


template <class T,class B> inline void Function<T,B>::eval(const T& x, T& f_x)
{   
    for(int n=0; i<ft_i.size(); i++) f_x += ft_i[n] * basis->evaluate_func(x,n);
};

template <class T,class B> inline void Function<T,B>::eval(const std::vector<T>& x, std::vector<T>& f_x){
    f_x.clear();
    for (const auto& xval: x) {
        T tmp = static_cast<T>(0);
        for(int n=0; i<ft_i.size(); i++) tmp += ft_i[n] * basis->evaluate_func(xval,n);
        f_x.push_back(tmp)
    }
};
            
};
