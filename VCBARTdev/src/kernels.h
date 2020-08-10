//
//  kernels.hpp
//  
//
//  Created by Sameer Deshpande on 8/31/19.
//

#ifndef _GUARD_kernels_h
#define _GUARD_kernels_h


#include <RcppArmadillo.h>
#include <cmath>
#include <stdio.h>

double ar_kernel(double z1, double z2, double rho);
double cs_kernel(double z1, double z2, double rho);



#endif /* kernels_hpp */
