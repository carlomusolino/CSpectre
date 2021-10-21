#include <cmath>
#include <vector>
#include <iostream>
#include "chebyshev.hpp"

namespace FunctionalBases {

  template <class T>
  class FunctionalBase {
    unsigned int N;
    std::vector<T> nodes;
    std::vector<T> weights;
  public:
    template <unsigned int n>
    virtual inline auto evaluate_function(const T& x);
    virtual inline void get_nodes(std::vector<T>* pts);
    virtual inline void get_weights(std::vector<T>* w);
  };

  template <class T>
  class ChebyshevBase public: FunctionalBase {
    unsigned int N;
    std::vector<T> nodes;
    std::vector<T> weights;
  public:
    ChebyshevBase(unsigned int N) : N(N)
      {
	weights.push_back(static_cast<T>(M_PI/(2*(N+1))));
	nodes.push_back(static_cast<T>(1))
	for (int i=1, i<N, i++){
	  nodes.push_back(static_cast<T>(std::cos(M_PI*i/N)));
	  weights.push_back(static_cast<T>(M_PI/N));
	}
	nodes.push_back(static_cast<T>(-1));
	weights.push_back(static_cast<T>(M_PI/(2*(N+1))));
      }
    template <unsigned int n> inline auto evaluate_function(const T& x){
      return Chebyshev::T<n>(x);
    }
    inline void get_nodes(std::vector<T>* pts){
      pts = &nodes;
    }
    inline void get_weights(std::vector<T>* w){
      w = &weights;
    }
  };


}
