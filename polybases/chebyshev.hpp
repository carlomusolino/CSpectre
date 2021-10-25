/**
 * @brief Somewhat efficient implementation of Chebyshev polynomial evaluation
 * @author Carlo Musolino (musolino@itp.uni-frankfurt.de)
 */
#ifndef _MY_CHEBYSHEV_HPP
#define _MY_CHEBYSHEV_HPP
#endif 

#include <cmath>

namespace Chebyshev {


  /*
  template <unsigned int N> inline auto Tn(const double& x)
  {
    return ( (2*x*Tn<N-1>(x)) - Tn<N-2>(x)  );
  }

  template <> inline auto Tn<0>(const double& x)
  {
    return 1.0;
  }

  template <> inline auto Tn<1>(const double& x)
  {
    return x; 
  }
  */

  template <class C> inline auto Tn(const C& x, unsigned int N)
  {
    switch(N){
    case 0: {
      return static_cast<C>(1);
    };
    case 1: {
      return x;
    };
    default:
    {
    return ( static_cast<C>(2)*x*Tn(x,N-1) - Tn(x,N-2) ) ;
    };
    }
  }


}
