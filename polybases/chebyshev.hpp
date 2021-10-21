#include <cmath>
#include <stdio>
#include <stdlib>

namespace Chebyshev {


  template <class T, unsigned int N> inline auto T(const T& x)
  {
    return ( (static_cast<T>(2)*x*<T,N-1>T(x)) - <T,N-2>T(x)  );
  }

  template <class T,0> inline auto T(const T& x)
  {
    return static_cast<T>(1);
  }

  template <class T,1> inline auto T(const T& x)
  {
    return x; 
  }

}
