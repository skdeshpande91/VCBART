//
//  kernels.cpp
//  
//
//  Created by Sameer Deshpande on 8/31/19.
//

#include "kernels.h"

double ar_kernel(double z1, double z2, double rho)
{
  return(pow(rho, abs(z1 - z2)));
}



double cs_kernel(double z1, double z2, double rho)
{
  if(z1 == z2) return(1.0);
  else return(rho);
}
