#include "../polybases/polybases.hpp"
#include <iostream>
#include <fstream>
using namespace FunctionalBases;

int main(){

    FunctionalBase<double> * cheb2 = new ChebyshevBase<double>(2u);

    std::vector<double> w, x;
    cheb2->get_nodes(x);
    cheb2->get_weights(w);
    std::cout << "nodes " << x.size() << "\n";
    for (const auto xi : x) std::cout << xi << "\n";
    std::cout << "weights" << "\n";
    for (const auto wi : w) std::cout << wi << "\n";

    cheb2->print_nodes();

    delete cheb2 ; 
}