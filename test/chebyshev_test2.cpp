#include "../polybases/polybases.hpp"
#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>

using namespace FunctionalBases;
using std::cout;

int main()
{   

    const int NPOINTS = 100;
  
    std::shared_ptr<FunctionalBase<double>> cheb6 {new ChebyshevBase<double>(6u)};


    std::ofstream outfile;
    outfile.open("chebyshev_test2.txt");

    outfile << "# Chebyshev polynomial test N = 0,1,2,3,4 \n";
  
    double x = -1.0;
    double dx = 2.0 / NPOINTS;
    for(int i=0; i<NPOINTS; i++){
        outfile << x << " " ; 
        outfile << cheb6->evaluate_function(x,0) << " " << cheb6->evaluate_function(x,1) << " " << cheb6->evaluate_function(x,2) << " ";
        outfile << cheb6->evaluate_function(x,3) << " " << cheb6->evaluate_function(x,4) << "\n";
        x += dx;
    }
  outfile.close();
}