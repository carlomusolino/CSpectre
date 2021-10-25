#include "../polybases/polybases.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>

using namespace FunctionalBases;
using std::cout;

inline void func(const std::vector<double> x, std::vector<double>& y);
inline void generate_linspace(const unsigned int NPOINTS, std::vector<double>& x);
inline void print_vector(const std::vector<double> v);
template <class T> inline void compute_interp_func(const std::vector<T> x, const std::vector<T>& ytilde, std::shared_ptr<FunctionalBase<T>> basis, std::vector<T>& y );

int main(){
    std::shared_ptr<FunctionalBase<double>> cheb4 {new ChebyshevBase<double>(4u)};
    std::shared_ptr<FunctionalBase<double>> cheb6 {new ChebyshevBase<double>(6u)};

    std::vector<double> ft4_i,f4_i,ft6_i,f6_i;
    std::vector<double> n6,n4;
    std::vector<double> f4,f6;
    std::vector<double> x;

    cheb4->print_weights();
    cheb6->print_weights();

    cheb4->print_nodes();
    cheb6->print_nodes();

    cheb4->get_nodes(n4);
    cheb6->get_nodes(n6);
    func(n4,f4_i);
    cout << f4_i.size() << "\n";
    func(n6,f6_i);
    cheb6->calc_spectral_coeffs(f6_i,ft6_i);
    cheb4->calc_spectral_coeffs(f4_i,ft4_i);

    generate_linspace(100,x);

    compute_interp_func<double>(x, ft4_i, cheb4, f4 );
    compute_interp_func<double>(x, ft6_i, cheb6, f6 );


    //print_vector(f4_i);
    //print_vector(f4);

    print_vector(ft6_i);
    //print_vector(f6);

    std::ofstream outfile4,outfile6;
    outfile4.open("interp_test_N4.txt");
    outfile6.open("interp_test_N6.txt");

    //for(int i=0;i<n4.size();i++) outfile4 << n4[i] << " " << f4_i[i] << " " << f4[i] << "\n";

    //for(int i=0;i<n6.size();i++) outfile6 << n6[i] << " " << f6_i[i] << " " << f6[i] << "\n";
    for(int i=0; i<x.size(); i++) outfile4 << x[i] << " " << f4[i] << "\n";

    for(int i=0; i<x.size(); i++) outfile6 << x[i] << " " << f6[i] << "\n";
    
    outfile4.close();
    outfile6.close();


};


template <class T> inline void compute_interp_func(const std::vector<T> x, const std::vector<T>& ytilde,std::shared_ptr<FunctionalBase<T>> basis, std::vector<T>& y )
{
    y.clear();
    for (auto& val : x) {
        T tmp = static_cast<T>(0.0);
        for(int N=0; N<ytilde.size();N++)
        {
            tmp += basis->evaluate_function(val,N) * ytilde[N];
        }
        y.push_back(tmp);
    }
};

inline void generate_linspace(const unsigned int NPOINTS, std::vector<double>& x){
    x.clear();
    double dx = 2.0 / NPOINTS;
    for(int i=0; i<NPOINTS; i++) x.push_back(i * dx- 1.0);
}

inline void func(const std::vector<double> x, std::vector<double>& y) {
    y.clear();
    for (auto& val: x) y.push_back( std::pow(std::cos(M_PI*val/2),3) + std::pow(val+1,3)/8.0 );
}
    
inline void print_vector(const std::vector<double> v)
{
    for (auto& val: v) cout << val << " ";
    cout << "\n";
};