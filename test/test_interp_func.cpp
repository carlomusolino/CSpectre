#include "../functions.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>

using namespace FunctionalBases;
using namespace Functions;
using std::cout;


inline void generate_linspace(const unsigned int NPOINTS, std::vector<double>& x);
inline void ffunc(const std::vector<double>& x, std::vector<double>& y);
inline void print_vector(const std::vector<double> v);

int main(){
    /*
    std::shared_ptr<FunctionalBase<double>> cheb4 {new ChebyshevBase<double>(4u)};
    std::shared_ptr<FunctionalBase<double>> cheb6 {new ChebyshevBase<double>};
    */

    /*
    std::unique_ptr<Function<double,ChebyshevBase<double>,6u>> func { 
        new Function<double,ChebyshevBase<double>,6u>(&ffunc) } 
    */
   using func4 = Function<double,ChebyshevBase<double>,4u>;
   using func6 = Function<double,ChebyshevBase<double>,6u>;

   std::unique_ptr<func4> f_4 {new func4(&ffunc) };
   std::unique_ptr<func6> f_6 {new func6(&ffunc) };
    
    std::vector<double> f4_i,f6_i;
    f_4->get_func_vals(f4_i);
    f_6->get_func_vals(f6_i);
    print_vector(f4_i);
    print_vector(f6_i);

    std::vector<double> x, y4,y6;
    generate_linspace(100,x);
    f_4->eval(x,y4);
    f_6->eval(x,y6);

    std::ofstream outfile4,outfile6;
    outfile6.open("interp_test_N6.txt");
    for(int i = 0; i<y6.size(); i++) outfile6 << x[i] << " " << y6[i] <<"\n";
    outfile4.open("interp_test_N4.txt");
    for(int i = 0; i<y4.size(); i++) outfile4 << x[i] << " " << y4[i] <<"\n";
    outfile6.close();
    outfile4.close();

}

inline void generate_linspace(const unsigned int NPOINTS, std::vector<double>& x){
    x.clear();
    double dx = 2.0 / NPOINTS;
    for(int i=0; i<NPOINTS; i++) x.push_back(i * dx- 1.0);
}

inline void ffunc(const std::vector<double>& x, std::vector<double>& y) {
    y.clear();
    for (auto& val: x) y.push_back( std::pow(std::cos(M_PI*val/2),3) + std::pow(val+1,3)/8.0 );
}

inline void print_vector(const std::vector<double> v)
{
    for (auto& val: v) cout << val << " ";
    cout << "\n";
}