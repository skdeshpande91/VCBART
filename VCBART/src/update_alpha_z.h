//
//  update_alpha_z.h
//  
//
//

#ifndef _GUARD_update_alpha_z_h
#define _GUARD_update_alpha_z_h
#include "rng.h"
#include "draw_multinomial.h"
#include <stdio.h>

void update_alpha_z(double &alpha_z, double &rho_alpha, std::vector<double> &theta, size_t &R, size_t &N_u, double &a, double &b, RNG &gen);

#endif /* update_alpha_z */
