#ifdef RAYTRACE_TG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allvars.h"

#define pi 3.1415927
#define c 2.99792458e10
#define h_nu 6.6262e-27
#define k_B 1.3806e-16
#define pc 3.085678e18

double heat_ion_rates(int rad_type, double L3, double T3)
  {
    int i = 0;
    int N_steps = 10000;
    int flag_sun = 0;
    double A_0 = 6.3e-18;
    double sigma = 2.0*pow(pi,5.0)*pow(k_B,4.0)/15.0/pow(h_nu,3.0)/pow(c,2.0);

    double L0 = pow(10.0, 5.568)*3.827e33;
    double L1 = pow(10.0, 6.095)*3.827e33;
    double L2 = pow(10.0, 6.574)*3.827e33;
    double T0 = pow(10.0, 4.922);
    double T1 = pow(10.0, 4.975);
    double T2 = pow(10.0, 4.999);

    double L = 0;
    double T_eff = 0;
    double prefactor = 0;
    double Z_HI = 1.0;
    double Z_HeI = 0.89;
    double Z_HeII = 2.0;
    double k_LW = 1.1e8;
    double nu_L_LW = 2.7e15;
    double nu_L_HI = 3.3e15;
    double nu_L_HeI = 5.95e15;
    double nu_L_HeII = 1.32e16;
    double nu_max = 1.0e20;
    double heat_HI = 0.0;
    double heat_HeI = 0.0;
    double heat_HeII = 0.0;
    double ion_LW = 0.0;
    double ion_HI = 0.0;
    double ion_HeI = 0.0;
    double ion_HeII = 0.0;
    double nu_LW[N_steps];
    double nu_HI[N_steps];
    double nu_HeI[N_steps];
    double nu_HeII[N_steps];
    double dnu_LW[N_steps];
    double dnu_HI[N_steps];
    double dnu_HeI[N_steps];
    double dnu_HeII[N_steps];
    double epsilon_HI[N_steps];
    double epsilon_HeI[N_steps];
    double epsilon_HeII[N_steps];
    double sigma_HI[N_steps];
    double sigma_HeI[N_steps];
    double sigma_HeII[N_steps];
    double B_LW[N_steps];
    double B_HI[N_steps];
    double B_HeI[N_steps];
    double B_HeII[N_steps];
    double x, y;

    int num=7;
    double heat_ion_vec[num];
 
    if(All.ray_flag_sun == 0)
      {
        L = L0;
        T_eff = T0;
      }

    if(All.ray_flag_sun == 1)
      {
        L = L1;
        T_eff = T1;
      }

    if(All.ray_flag_sun == 2)
      {
        L = L2;
        T_eff = T2;
      }

    if(All.ray_flag_sun == 3)
      {
        L = L3;
        T_eff = T3;
      }

    if(L > 0)
     prefactor = L/4.0/sigma/pow(T_eff,4.0)/pow(pc,2.0);
    else
     prefactor = 0;    

    for(i = 0; i < N_steps; i++)
      {
        nu_LW[i] = nu_L_LW*pow(nu_L_HI / nu_L_LW, (i + 0.5) / (double)N_steps);
        nu_HI[i] = nu_L_HI*pow(nu_L_HeII / nu_L_HI, (i + 0.5) / (double)N_steps);
        nu_HeI[i] = nu_L_HeI*pow(nu_L_HeII / nu_L_HeI, (i + 0.5) / (double)N_steps);
        nu_HeII[i] = nu_L_HeII*pow(nu_max / nu_L_HeII, (i + 0.5) / (double)N_steps);

        dnu_LW[i] = nu_L_LW*pow(nu_L_HI / nu_L_LW, (i + 1) / (double)N_steps) - nu_L_LW*pow(nu_L_HI / nu_L_LW, i / (double)N_steps);
        dnu_HI[i] = nu_L_HI*pow(nu_L_HeII / nu_L_HI, (i + 1) / (double)N_steps) - nu_L_HI*pow(nu_L_HeII / nu_L_HI, i / (double)N_steps);
        dnu_HeI[i] = nu_L_HeI*pow(nu_L_HeII / nu_L_HeI, (i + 1) / (double)N_steps) - nu_L_HeI*pow(nu_L_HeII / nu_L_HeI, i / (double)N_steps);
        dnu_HeII[i] = nu_L_HeII*pow(nu_max / nu_L_HeII, (i + 1) / (double)N_steps) - nu_L_HeII*pow(nu_max / nu_L_HeII, i / (double)N_steps);

        epsilon_HI[i] = sqrt(nu_HI[i]/nu_L_HI-1.0);
        epsilon_HeI[i] = sqrt(nu_HeI[i]/nu_L_HeI-1.0);
        epsilon_HeII[i] = sqrt(nu_HeII[i]/nu_L_HeII-1.0);

        x= (nu_HeI[i]/3.286e15) - .4434;
        y= pow(x,2) + 4.563;

        sigma_HI[i] = A_0/pow(Z_HI,2.0)*pow(nu_L_HI/nu_HI[i],4.0)*exp(4.0-4.0*atan(epsilon_HI[i])/epsilon_HI[i])/(1.0-exp(-2.0*pi/epsilon_HI[i]));
        //sigma_HeI[i] = A_0/pow(Z_HeI,2.0)*pow(nu_L_HeI/nu_HeI[i],4.0)*exp(4.0-4.0*atan(epsilon_HeI[i])/epsilon_HeI[i])/(1.0-exp(-2.0*pi/epsilon_HeI[i]));
        sigma_HeI[i] = 9.492e-16*(pow(x-1,2) + 4.158)*pow(y,-1.953)*pow(1 + .825*pow(y,0.25),-3.188); 
        sigma_HeII[i] = A_0/pow(Z_HeII,2.0)*pow(nu_L_HeII/nu_HeII[i],4.0)*exp(4.0-4.0*atan(epsilon_HeII[i])/epsilon_HeII[i])/(1.0-exp(-2.0*pi/epsilon_HeII[i]));

        B_LW[i] = 2.0*h_nu*pow(nu_LW[i],3.0)/pow(c,2.0)/(exp(h_nu*nu_LW[i]/k_B/T_eff)-1.0);
        B_HI[i] = 2.0*h_nu*pow(nu_HI[i],3.0)/pow(c,2.0)/(exp(h_nu*nu_HI[i]/k_B/T_eff)-1.0);
        B_HeI[i] = 2.0*h_nu*pow(nu_HeI[i],3.0)/pow(c,2.0)/(exp(h_nu*nu_HeI[i]/k_B/T_eff)-1.0);
        B_HeII[i] = 2.0*h_nu*pow(nu_HeII[i],3.0)/pow(c,2.0)/(exp(h_nu*nu_HeII[i]/k_B/T_eff)-1.0);

        heat_HI += prefactor*B_HI[i]*sigma_HI[i]*(1.0-nu_L_HI/nu_HI[i])*dnu_HI[i];
        heat_HeI += prefactor*B_HeI[i]*sigma_HeI[i]*(1.0-nu_L_HeI/nu_HeI[i])*dnu_HeI[i];
        heat_HeII += prefactor*B_HeII[i]*sigma_HeII[i]*(1.0-nu_L_HeII/nu_HeII[i])*dnu_HeII[i];

        //ion_LW += prefactor*k_LW*B_LW[i]*dnu_LW[i];
        ion_HI += prefactor*B_HI[i]*sigma_HI[i]/h_nu/nu_HI[i]*dnu_HI[i];
        ion_HeI += prefactor*B_HeI[i]*sigma_HeI[i]/h_nu/nu_HeI[i]*dnu_HeI[i];
        ion_HeII += prefactor*B_HeII[i]*sigma_HeII[i]/h_nu/nu_HeII[i]*dnu_HeII[i];
      }
   
    ion_LW = prefactor*k_LW*B_LW[i-4830];

    if(All.NumCurrentTiStep == 1)
      printf("sigma_H0 = %lg, sigma_He0 = %lg \n", sigma_HI[i-100], sigma_HeI[i-100]);

/*
    printf("heat_HI = %g\n", heat_HI);
    printf("heat_HeI = %g\n", heat_HeI);
    printf("heat_HeII = %g\n", heat_HeII);

    printf("ion_LW = %g\n", ion_LW);
    printf("ion_HI = %g\n", ion_HI);
    printf("ion_HeI = %g\n", ion_HeI);
    printf("ion_HeII = %g\n", ion_HeII);
*/
    heat_ion_vec[0] = heat_HI;
    heat_ion_vec[1] = heat_HeI;
    heat_ion_vec[2] = heat_HeII;
    heat_ion_vec[3] = ion_HI;
    heat_ion_vec[4] = ion_HeI;
    heat_ion_vec[5] = ion_HeII;
    heat_ion_vec[6] = ion_LW;

   if(rad_type == 0) 
     return heat_HI;

   else if(rad_type == 1)
     return heat_HeI;
    
   else if(rad_type == 2)
     return heat_HeII;

   else if(rad_type == 3)
     return ion_HI;

   else if(rad_type == 4)
     return ion_HeI;

   else if(rad_type == 5)
     return ion_HeII;

   else if(rad_type == 6)
     return ion_LW;

   else
   {
    printf("rad_type not properly set!\n");
    exit(0);
   }
}

#endif
