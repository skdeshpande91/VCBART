//
//  draw_multinomial.h
//    A quick multinomial sampler
//    Probably faster to use Rcpp's existing one
//    but let's roll our own for the moment
//  Created by Sameer Deshpande on 4/6/20.
//

#ifndef _GUARD_draw_multinomial_h
#define _GUARD_draw_multinomial_h
#include<RcppArmadillo.h>
#include "rng.h"


size_t draw_multinomial(size_t &R, std::vector<double> &probs, RNG &gen);

#endif /* draw_multinomial_hpp */
