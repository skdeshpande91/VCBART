//
//  update_rho.hpp
//  
//
//  Created by Cecilia Balocchi on 2/22/21.
//

#ifndef GUARD_update_rho_h
#define GUARD_update_rho_h

#include <RcppArmadillo.h>
#include <cmath>
#include "tree.h"
#include "info.h"
#include "funs.h"
#include <stdio.h>

#endif /* update_rho_hpp */

// update rho when we have global sigma
void update_rho_cs(double &rho_eps, double &sigma, data_info &di, RNG &gen);
void update_rho_cs(double &rho_eps, double &sigma, double &xi_sd, double &xi_mean, int iter, data_info &di, RNG &gen);