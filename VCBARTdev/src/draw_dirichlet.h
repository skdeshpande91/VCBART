//
//  draw_dirichlet.h
//  
//
//  Created by Sameer Deshpande on 4/7/20.
//

#ifndef _GUARD_draw_dirichlet_h
#define _GUARD_draw_dirichlet_h

#include "rng.h"
#include <stdio.h>


void draw_dirichlet(std::vector<double> &theta, double &alpha, size_t &R, RNG &gen); // for dirichlet(alpha/R, ..., alpha/R)

void draw_dirichlet(std::vector<double> &theta, std::vector<double> &alpha, size_t &R, RNG &gen);


#endif /* draw_dirichlet_h */
