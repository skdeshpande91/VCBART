//
//  draw_multinomial.cpp
//  
//
//  Created by Sameer Deshpande on 4/6/20.
//

#include "draw_multinomial.h"
size_t draw_multinomial(size_t &R, std::vector<double> &probs, RNG &gen)
{
  size_t x = 0;
  double cumsum = 0.0;
  double unif = gen.uniform();
  for(size_t r = 0; r < R; r++){
    cumsum += probs[r];
    if(unif < cumsum){
      x = r;
      break;
    }
  }
  return(x);
}
