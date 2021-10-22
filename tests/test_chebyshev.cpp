#include <iostream>
#include <fstream>
#include "../polybases/chebyshev.hpp"



int main()
{
  const int NPOINTS = 100;
  
  std::ofstream outfile;
  outfile.open("chebyshev_test.txt");

  outfile << "# Chebyshev polynomial test N = 0,1,2,3,4 \n";
  
  double x = -1.0;
  double dx = 2.0 / NPOINTS;
  for(int i=0; i<NPOINTS; i++){
    outfile << x << " " ; 
    outfile << Chebyshev::Tn(x,0) << " " << Chebyshev::Tn(x,1) << " " << Chebyshev::Tn(x,2) << " ";
    outfile << Chebyshev::Tn(x,3) << " " << Chebyshev::Tn(x,4) << "\n";
    x += dx;
  }
  outfile.close();
}
