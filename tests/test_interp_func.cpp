#include "../polybases/polybases.hpp"
#include "../functions.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>

using namespace FunctionalBases;
using namespace Functions;
using std::cout;



inline void func(const std::vector<double> x, std::vector<double>& y);

int main(){
    /*
    std::shared_ptr<FunctionalBase<double>> cheb4 {new ChebyshevBase<double>(4u)};
    std::shared_ptr<FunctionalBase<double>> cheb6 {new ChebyshevBase<double>};
    */

   std::unique_ptr<Function<double,ChebyshevBase<double>*>> func { new Function<double,ChebyshevBase<double>*>(6,&func) } ;

   


}

inline void func(const std::vector<double> x, std::vector<double>& y) {
    y.clear();
    for (auto& val: x) y.push_back( std::pow(std::cos(M_PI*val/2),3) + std::pow(val+1,3)/8.0 );
}