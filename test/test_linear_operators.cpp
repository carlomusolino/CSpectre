#include "../ODE/linear_diff_ops.hpp"
#include <memory>
#include <iostream>

using std::cout;
using namespace Operators;
using namespace FunctionalBases;

int main()
{   
    FunctionalBase<double>* cheb6 {new ChebyshevBase<double>(6u)};

    
    LinearOperator<double> * L1 = new Derivative<double>(cheb6) ;
    LinearOperator<double> * L2 = new TimesX<double>(cheb6);

    LinearOperator<double> L3;

    L3 = *L1 + (*L2)*2.0;

    L1->print_Lij();
    cout << "\n\n";

    L2->print_Lij();
    cout << "\n\n";
    
    L3.print_Lij();

    delete L1;
    delete L2;
    delete cheb6;
}