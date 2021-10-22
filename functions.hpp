#include <vector>
#include "polybases/polybases.hpp"

using namespace FunctionalBases;

namespace Functions {

    template <class T, class FuncBase, unsigned int N>
    class Function {
        std::vector<T> f_i;   // physical space representation of f
        std::vector<T> ft_i;  // collocation space representation of f
        FuncBase* basis; // should be a derived class of FunctionBase
        public:
        // creators
        Function<T,FuncBase,N>(const std::vector<T>& f_i): f_i(f_i) {
            assert(f_i.size()==N);
            basis = new FuncBase(N);
            decompose();
        };
        Function<T,FuncBase,N>(void func(const std::vector<T>&, std::vector<T>&)) {
            basis = new FuncBase(N) ;
            std::vector<T> n;
            basis->get_nodes(n);
            (*func)(n,f_i);
            decompose();
        }
        // members
        inline void eval(const T& x, T& f_x);
        inline void eval(const std::vector<T>& x, std::vector<T>& f_x);
        inline void decompose(){
            basis->calc_spectral_coeffs(f_i,ft_i);
        }
        // access
        inline void get_spectral_coeffs(std::vector<T>& y){ y = ft_i;};
        inline void get_func_vals(std::vector<T>& y){ y = f_i; };
        // destructor
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
            
};
