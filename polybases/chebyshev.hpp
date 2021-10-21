#include <cmath>

namespace Chebyshev {


  template <unsigned int N> inline auto T(const double& x)
  {
    return ( (2*x*T<N-1>(x)) - T<N-2>(x)  );
  }

  template <> inline auto T<0>(const double& x)
  {
    return 1.0;
  }

  template <> inline auto T<1>(const double& x)
  {
    return x; 
  }

}
