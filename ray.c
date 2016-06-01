#ifdef RAYTRACE_TG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"

void raytrace_TG(double dt_raytrace)
  {
    int i, j, k, l, m, n, flag_bug, flag_bug_tot, flag_theta, flag_phi, flag_r, N_H_ionized_parts_local, N_He_ionized_parts_local, N_H_ionized_parts_tot, N_He_ionized_parts_tot;

    double a, a3, hubble, hubble2, hubble3, elec_here, mu_here, gamma_here, r_here, nh_here, entropy_here, r_k_lower, r_k_upper, dr_k, N_photons_per_ray, d_N, N_atoms, N_recombinations, sim_theta, sim_theta_left, sim_theta_right, sim_phi, sim_phi_left, sim_phi_right, sim_r, sim_r2, sim_r_lower, sim_r_upper, t_start, t_end, HI, HeI, rate_HII, rate_HeII, rate_HeIII, rate_LW, eq_HII, eq_HeII, eq_HeIII, eq_LW, t_HII, t_HeII, t_HeIII, t_LW;
    double x_shift, y_shift, z_shift, theta_shift, phi_shift;
    double nu_min_H = 3.3e15;
    double nu_min_He = 1.32e16;
    double star_mass, n_max= 1.e5, n_cut = 1.e-2;
    double prad_dir, sigma_ion = 6.e-18;
    double lum_tot, Teff, r1, r2, min_fac=1.0;
    double sim_NH2, e_fac, fshield, NH2_fac, h2_formrate, temp_here, muh2, muh2in, h2frac;

    char buf[500];

    FILE *infile, *outfile, *outfile2;

    flag_bug_tot = flag_bug = N_H_ionized_parts_tot = N_He_ionized_parts_tot = N_H_ionized_parts_local = N_He_ionized_parts_local = 0;

    ray.r_min_local = ray.r_min = ray.center_z_local = ray.center_z = ray.center_y_local = ray.center_y = ray.center_x_local = ray.center_x = -1;
    ray.r_min_local = 100.0;
    r1 = 5.e-5; 
 
    if(All.ComovingIntegrationOn)
      {
        a = All.Time;
        a3 = a*a*a;
        hubble = All.HubbleParam;
        hubble2 = hubble*hubble;
        hubble3 = hubble*hubble2;
      }
    else
      a = a3 = hubble = hubble2 = hubble3 = 1.0;

    if(All.NumCurrentTiStep == 0)
      {
        ray.flag_switch = 0;
        ray.flag_continue_start = 0;
        ray.flag_continue_end = 1;
        ray.alpha_B_H = 2.6e-13;
        ray.alpha_B_He = 1.3e-12;
        ray.temp_HII_estimate = 1.0e4;
        ray.temp_HeIII_estimate = 2.0e4;

        if(All.ray_flag_sun == 0)
          {
            ray.lifetime = 3.7e6*SEC_PER_YEAR;
            ray.r_max = 500.0;
            ray.HI_ion_rate = 4.55e-7*1.e-20;
            ray.HeI_ion_rate = 4.18e-7*1.e-20;
            ray.HeII_ion_rate = 6.73e-9*1.e-20;
            ray.LW_rate = 1.28e-7*1.e-20;
          }
        else if(All.ray_flag_sun == 1)
          {
            ray.lifetime = 2.7e6*SEC_PER_YEAR;
            ray.r_max = 150.0;
            ray.HI_ion_rate = 1.32e-6;
            ray.HeI_ion_rate = 1.43e-6;
            ray.LW_rate = 3.38e-7;
          }
        else if(All.ray_flag_sun == 2)
          {
            ray.lifetime = 2.2e6*SEC_PER_YEAR;
            ray.r_max = 200.0;
            ray.HI_ion_rate = 3.69e-6;
            ray.HeI_ion_rate = 4.29e-6;
            ray.LW_rate = 9.07e-7;
          }
        else if(All.ray_flag_sun == 3)            //Get mass from sink particle          
          {
            ray.lifetime = 2.2e6*SEC_PER_YEAR;
            //ray.r_max = 200.0;
            ray.r_max = 1.5e-1;    

            ray.HI_ion_rate = heat_ion_rates(3, All.lum_tot, All.Teff);
            ray.HeI_ion_rate = heat_ion_rates(4, All.lum_tot, All.Teff);           
            ray.LW_rate = heat_ion_rates(6, All.lum_tot, All.Teff);
       
            if(ThisTask == 0)
              printf("HI_ion_rate = %lg, HeI_ion_rate = %lg, LW_rate = %lg \n", ray.HI_ion_rate, ray.HeI_ion_rate, ray.LW_rate);
          }

        else
          {
            printf("Flag_sun not set! Aborting...\n");

            exit(0);
          }

        ray.N_rays = 0;

        for(i = 0; i < N_theta; i++)
          {
            ray.theta[i] = (2.0*i+1.0)*PI/2.0/(double)N_theta;

            ray.N_phi[i] = 2.0*N_theta*sqrt(pow(sin(ray.theta[i]),2.0));

            if(ray.N_phi[i] == 0)
              ray.N_phi[i] = 1;

            for(j = 0; j < ray.N_phi[i]; j++)
              {
                ray.phi[i][j] = (2.0*j+1.0)*PI/(double)ray.N_phi[i];

                ray.N_rays++;

                ray.H_front[i][j] = ray.He_front[i][j] = 0;
              }

            for(j = ray.N_phi[i]; j < 2*N_theta; j++)
              ray.phi[i][j] = ray.H_front[i][j] = ray.He_front[i][j] = 0;
          }

        if(ThisTask == 0)
          printf("Found %d rays and %d bins\n", ray.N_rays, ray.N_rays*N_shells);
      }

    if(All.NumCurrentTiStep == 0 && ray.flag_continue_start == 1)
      {
        sprintf(buf, "%s/%s.front2", All.OutputDir, All.SnapshotFileBase);

        infile = fopen(buf, "r");

        for(i = 0; i < N_theta; i++)
          for(j = 0; j < ray.N_phi[i]; j++)
            {
              //fscanf(infile, "%f", &ray.H_front[i][j]);
              //fscanf(infile, "%f", &ray.He_front[i][j]);
 
              fscanf(infile, "%lg %lg\n", &ray.H_front[i][j], &ray.He_front[i][j]);

              if(ThisTask == 0 && i%10 == 0 && j%10 == 0)
                printf("read the i-front file! ray.H_front = %lg\n", ray.H_front[i][j]);
            }

        fclose(infile);
      }

    MPI_Barrier(MPI_COMM_WORLD);

    if(ray.flag_continue == 1 && All.Time > ray.time_end)
      {
        for(n = 0; n < N_gas; n++)
          {
          SphP[n].Ray_H_coeff = SphP[n].Ray_He_coeff = SphP[n].Ray_LW_coeff = -1.0;
          SphP[n].Prad = 0;
          for(m = 0; m < 3; m++)
           SphP[n].Prad_dir[m] = 0;
          }

        for(i = 0; i < N_theta; i++)
          for(j = 0; j < 2 * N_theta; j++)
            ray.H_front[i][j] = ray.He_front[i][j] = 0;

        ray.flag_continue = ray.flag_start = 0;

        ray.time_start = ray.time_end = All.TimeMax;

        if(ThisTask == 0)
          printf("Raytracing ended!\n");
      }

    if(ray.flag_start == 1 && ray.flag_continue == 0)
      {
        ray.flag_continue = 1;
        ray.time_start = All.Time;
        ray.time_end = pow(pow(All.Time, 1.5) + 1.5 * HUBBLE * hubble * sqrt(All.Omega0) * ray.lifetime, 2.0/3.0);

        if(ThisTask == 0)
          {
            printf("Raytracing start: %g\n", ray.time_start);
            printf("Raytracing end: %g\n", ray.time_end);
          }
      }

    if(ray.flag_continue == 0)  
      return;

    if(ThisTask == 0)
      printf("\nTime Step %d: doing raytracing for %g years...\n", All.NumCurrentTiStep, dt_raytrace/SEC_PER_YEAR);

    t_start = second();

    if(All.ray_flag_sun == 0)
      {
        ray.N_H_photons = 2.80e49*dt_raytrace*1.e-50;
        ray.N_He_photons = 7.19e47*dt_raytrace*1.e-50;
      }
    else if(All.ray_flag_sun == 1)
      {
        ray.N_H_photons = 9.14e49*dt_raytrace;
        ray.N_He_photons = 4.14e48*dt_raytrace;
      }
    else if(All.ray_flag_sun == 2)
      {
        ray.N_H_photons = 2.70e50*dt_raytrace;
        ray.N_He_photons = 1.54e49*dt_raytrace;
      }
    else if(All.ray_flag_sun == 3)               //evolving ionization rate based on sink mass
      { 

        //ray.Q_H_ion = lum_calc(0, All.star_mass, All.mdot, nu_min_H, dt_raytrace);
        //ray.Q_He_ion = lum_calc(0, All.star_mass, All.mdot, nu_min_He, dt_raytrace);
        //All.lum_tot = lum_calc(1, All.star_mass, All.mdot, nu_min_H, dt_raytrace);
        //All.Teff = lum_calc(2, All.star_mass, All.mdot, nu_min_H, dt_raytrace);

        ray.N_H_photons = ray.Q_H_ion*dt_raytrace;        
        ray.N_He_photons = ray.Q_He_ion*dt_raytrace;

        ray.HI_ion_rate = heat_ion_rates(3, All.lum_tot, All.Teff);
        ray.HeI_ion_rate = heat_ion_rates(4, All.lum_tot, All.Teff);
        ray.LW_rate = heat_ion_rates(6, All.lum_tot, All.Teff);

        if(ThisTask == 0)
         printf("dt_raytrace = %lg N_H_photons = %lg, N_He_photons = %lg, star_mass = %lg, mdot = %lg All.star_rad = %lg\n", dt_raytrace, ray.N_H_photons, ray.N_He_photons, All.star_mass, All.mdot, All.star_rad);

        if(ThisTask == 0)
          printf("HI_ion_rate = %lg, HeI_ion_rate = %lg, LW_rate = %lg \n", ray.HI_ion_rate, ray.HeI_ion_rate, ray.LW_rate);
      }
    else
      {
        printf("Flag_sun not set! Aborting...\n");

        exit(0);
      }

    for(i = 0; i < N_theta; i++)
      for(j = 0; j < 2*N_theta; j++)
        for(k = 0; k < N_shells; k++)
          {
          ray.count_tot[i][j][k] = ray.count_local[i][j][k] = 0;
          ray.nh_tot[i][j][k] = ray.nh_local[i][j][k] = 0;
          //ray.ne_local[i][j][k] = ray.ne_tot[i][j][k] = 0;
          ray.NH2_local[i][j][k] = ray.NH2_tot[i][j][k] = 0;
          }

    for(n = 0; n < N_gas; n++)
      {
      if(P[n].ID == All.ray_center_ID)
        {
          ray.center_x_local = P[n].Pos[0];
          ray.center_y_local = P[n].Pos[1];
          ray.center_z_local = P[n].Pos[2];

          //if(All.Teff > 1.e4)
          //  ray.r_min_local = SphP[n].Hsml;
          //else
          //  ray.r_min_local = 2.0e0*All.SofteningGas;
        }
      if(P[n].ID ==  3385271)
        {
          ray.center_x_local2 = P[n].Pos[0];
          ray.center_y_local2 = P[n].Pos[1];
          ray.center_z_local2 = P[n].Pos[2];
        }
      }

    MPI_Allreduce(&ray.center_x_local, &ray.center_x, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&ray.center_x_local2, &ray.center_x2, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if(ray.r_min < 0)
      {
        printf("Warning! Source ID not found! Aborting...\n");

        exit(0);
      }


//ARS trying new algorithm for minimum i-front distance

    for(n = 0; n < N_gas; n++)
      {
       if(P[n].ID > 0)
         {
         sim_r = sqrt((P[n].Pos[0]-ray.center_x)*(P[n].Pos[0]-ray.center_x)+(P[n].Pos[1]-ray.center_y)*(P[n].Pos[1]-ray.center_y)+(P[n].Pos[2]-ray.center_z)*(P[n].Pos[2]-ray.center_z));
 
        if(min_fac*sim_r < ray.r_min_local && sim_r > 0.0)
          ray.r_min_local = min_fac*sim_r;
          }
       }

//end new algorithm
    MPI_Allreduce(&ray.r_min_local, &ray.r_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);


    if(ThisTask == 0)
      {
      printf("ID = %d, x = %15.11g, y = %15.11g, z = %15.11g, r_min = %g\n", All.ray_center_ID, ray.center_x, ray.center_y, ray.center_z, ray.r_min);
      printf("ID = %d, x2 = %15.11g, y2 = %15.11g, z2 = %15.11g, r_min = %g\n", All.ray_center_ID, ray.center_x2, ray.center_y2, ray.center_z2, ray.r_min);
      }

   for(k=0; k<100; k++)
     ray.r[k] = ray.r_min + k*r1;

   l = (int) (2.0*N_shells*log(ray.r[99]/ray.r_min)/log(ray.r_max/ray.r_min) - 1)/2.0;

   for(k=100; k<N_shells; k++)
     {
       //ray.r[k] = ray.r[49] * pow(ray.r_max / ray.r[49], (2.0 * (k-50) + 1.0) / 2.0 / N_shells);
     ray.r[k] = ray.r_min * pow(ray.r_max / ray.r_min, (2.0 * (k - 99 + l) + 1.0) / 2.0 / N_shells);
     }


//    for(k = 0; k < N_shells; k++)
//      ray.r[k] = ray.r_min * pow(ray.r_max / ray.r_min, (2.0 * k + 1.0) / 2.0 / N_shells);


    MPI_Barrier(MPI_COMM_WORLD);

   if(ThisTask == 0 && All.NumCurrentTiStep < 3 /*All.NumCurrentTiStep % 1000 == 0*/)
     {
     for(k = 0; k < N_shells; k++)
       printf("k = %d, r = %lg\n", k, ray.r[k]);
     }

    for(n = 0; n < N_gas; n++)
      if(P[n].ID > 0 && SphP[n].Density*All.UnitDensity_in_cgs*hubble2/a3*HYDROGEN_MASSFRAC/PROTONMASS > 1.e0)  //skip the low-density far-out gas
        {
          i = j = k = flag_theta = flag_phi = flag_r = 0;

          sim_r2 = sqrt((P[n].Pos[0]-ray.center_x2)*(P[n].Pos[0]-ray.center_x2)+(P[n].Pos[1]-ray.center_y2)*(P[n].Pos[1]-ray.center_y2)+(P[n].Pos[2]-ray.center_z2)*(P[n].Pos[2]-ray.center_z2));

          nh_here = SphP[n].Density*All.UnitDensity_in_cgs*hubble2/a3*HYDROGEN_MASSFRAC/PROTONMASS;

          if(sim_r2 < 1.e-5)
             {
             //printf("Upping the density! nh_old = %lg nh_new = %lg\n", nh_here, 1.e12);
             //nh_here = 1.e13;
             }

          if(SphP[n].sink > 0)
          {
           //nh_here = 16.*nh_here;
           //printf("Increasing density of sinks, ID = %d, nh = %lg\n", P[n].ID, nh_here);
          }

          if(P[n].ID == All.ray_center_ID)
          {
           //nh_here = 0;
           nh_here = fmax(All.n_sink, n_max);
           printf("Setting star_mass density to zero, ID = %d, nh = %lg\n", P[n].ID, nh_here);
          }

          sim_r = sqrt((P[n].Pos[0]-ray.center_x)*(P[n].Pos[0]-ray.center_x)+(P[n].Pos[1]-ray.center_y)*(P[n].Pos[1]-ray.center_y)+(P[n].Pos[2]-ray.center_z)*(P[n].Pos[2]-ray.center_z));

          sim_r_lower = sim_r-1.0/2.0*SphP[n].Hsml;

          if(sim_r_lower < 0.0)
            sim_r_lower = 0.0;

          sim_r_upper = sim_r+1.0/2.0*SphP[n].Hsml;

          if(sim_r == 0.0)
            sim_theta = PI/2.0;
          else
            sim_theta = acos((P[n].Pos[2]-ray.center_z)/sim_r);

          if(sim_r == 0.0)
            sim_theta_left = 0.0;
          else
            sim_theta_left = sim_theta-1.0/2.0*SphP[n].Hsml/sim_r;

          if(sim_theta_left < 0.0)
            sim_theta_left = 0.0;

          if(sim_r == 0.0)
            sim_theta_right = PI;
          else
            sim_theta_right = sim_theta+1.0/2.0*SphP[n].Hsml/sim_r;

          if(sim_theta_right > PI)
            sim_theta_right = PI;

          if(P[n].Pos[0]-ray.center_x > 0.0 && P[n].Pos[1]-ray.center_y >= 0.0)
            sim_phi = atan((P[n].Pos[1]-ray.center_y)/(P[n].Pos[0]-ray.center_x));
          else if(P[n].Pos[0]-ray.center_x < 0.0 && P[n].Pos[1]-ray.center_y >= 0.0)
            sim_phi = PI+atan((P[n].Pos[1]-ray.center_y)/(P[n].Pos[0]-ray.center_x));
          else if(P[n].Pos[0]-ray.center_x < 0.0 && P[n].Pos[1]-ray.center_y < 0.0)
            sim_phi = PI+atan((P[n].Pos[1]-ray.center_y)/(P[n].Pos[0]-ray.center_x));
          else if(P[n].Pos[0]-ray.center_x > 0.0 && P[n].Pos[1]-ray.center_y < 0.0)
            sim_phi = 2.0*PI+atan((P[n].Pos[1]-ray.center_y)/(P[n].Pos[0]-ray.center_x));
          else if(P[n].Pos[0]-ray.center_x == 0.0 && P[n].Pos[1]-ray.center_y > 0.0)
            sim_phi = PI/2.0;
          else if(P[n].Pos[0]-ray.center_x == 0.0 && P[n].Pos[1]-ray.center_y < 0.0)
            sim_phi = 3.0/2.0*PI;
          else
            sim_phi = PI;

          if(sim_r == 0.0 || sin(sim_theta) == 0.0)
            sim_phi_left = 0.0;
          else
            sim_phi_left = sim_phi-1.0/2.0*SphP[n].Hsml/sim_r/sin(sim_theta);

          if(sim_phi_left < 0.0)
            sim_phi_left = 0.0;

          if(sim_r == 0.0 || sin(sim_theta) == 0.0)
            sim_phi_right = 2.0*PI;
          else
            sim_phi_right = sim_phi+1.0/2.0*SphP[n].Hsml/sim_r/sin(sim_theta);

          if(sim_phi_right > 2.0*PI)
            sim_phi_right = 2.0*PI;

          do
            {
              if(ray.theta[i] >= sim_theta_left && ray.theta[i] <= sim_theta_right)
                {
                  flag_theta = 1;

                  do
                    {
                      if(ray.phi[i][j] >= sim_phi_left && ray.phi[i][j] <= sim_phi_right)
                        {
                          flag_phi = 1;

                          do
                            {
                              if(ray.r[k] >= sim_r_lower && ray.r[k] <= sim_r_upper)
                                {
                                  flag_r = 1;

                                  ray.nh_local[i][j][k] += nh_here*nh_here*nh_here/*(1.0 - SphP[n].TracAbund[IHP])*/;
                                  //ray.ne_local[i][j][k] += /*nh_here*nh_here**/nh_here*(SphP[n].TracAbund[IHP] + SphP[n].TracAbund[IDP] + SphP[n].TracAbund[IHEP]);                                  
                                  if(k == 0 && ray.NH2_local[i][j][k] <  1.e16*1.e14*1.e14*1.0*(ray.r[k] - 0))
                                  //if(SphP[n].sink > 0)
                                    ray.NH2_local[i][j][k] +=  1.e16*1.e16*1.e16*1.0*(ray.r[k] - 0);
 
                                  if(k > 0)  //ARS adding up H2 column densities
                                    {
                                    ray.NH2_local[i][j][k] +=  nh_here*nh_here*nh_here*SphP[n].TracAbund[IH2]*(ray.r[k] - ray.r[k-1]);
                                    }
                                  
                                  ray.count_local[i][j][k] += nh_here*nh_here;
                                }
                              else if(ray.r[k] > sim_r_upper)
                                flag_r = 2;
                              else
                                flag_r = 0;

                              k++;
                            }
                          while((flag_r == 0 || flag_r == 1) && k < N_shells);

                          k = 0;
                        }
                      else if(ray.phi[i][j] > sim_phi_right)
                        flag_phi = 2;
                      else flag_phi = 0;

                      j++;
                    }
                  while((flag_phi == 0 || flag_phi == 1) && j < ray.N_phi[i]);

                  j = 0;
                }
              else if(ray.theta[i] > sim_theta_right)
                flag_theta = 2;
              else
                flag_theta = 0;

              i++;
            }
          while((flag_theta == 0 || flag_theta == 1) && i < N_theta);

          i = 0;
        }

    MPI_Allreduce(&ray.count_local[0][0][0], &ray.count_tot[0][0][0], 2*N_theta*N_theta*N_shells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(&ray.nh_local[0][0][0], &ray.nh_tot[0][0][0], 2*N_theta*N_theta*N_shells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //MPI_Allreduce(&ray.ne_local[0][0][0], &ray.ne_tot[0][0][0], 2*N_theta*N_theta*N_shells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(&ray.NH2_local[0][0][0], &ray.NH2_tot[0][0][0], 2*N_theta*N_theta*N_shells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    for(i = 0; i < N_theta; i++)
      for(j = 0; j < ray.N_phi[i]; j++)
        for(k = 0; k < 1; k++)
          {
          if(ray.nh_tot[i][j][0] <  fmax(All.n_sink, n_max)*fmax(All.n_sink, n_max)*fmax(All.n_sink, n_max))
            {
            ray.nh_tot[i][j][0] = fmax(All.n_sink, n_max)*fmax(All.n_sink, n_max)*fmax(All.n_sink, n_max);
            ray.count_tot[i][j][0] =  fmax(All.n_sink, n_max)*fmax(All.n_sink, n_max);
            }
          }


    for(i = 0; i < N_theta; i++)
      for(j = 0; j < ray.N_phi[i]; j++)
        for(k = 0; k < N_shells; k++)
          if(ray.count_tot[i][j][k] > 0)
            {
            ray.nh_tot[i][j][k] = ray.nh_tot[i][j][k] / ray.count_tot[i][j][k];
            //ray.ne_tot[i][j][k] = ray.ne_tot[i][j][k] / ray.count_tot[i][j][k];
            ray.NH2_tot[i][j][k] = ray.NH2_tot[i][j][k] / ray.count_tot[i][j][k];
            }
          else
            {
            ray.nh_tot[i][j][k] = ray.nh_tot[i][j][k-1];  //ARS makin' some changes
            //ray.ne_tot[i][j][k] = ray.ne_tot[i][j][k-1];
            ray.NH2_tot[i][j][k] = ray.NH2_tot[i][j][k-1];
            //ray.nh_tot[i][j][k] = 0;
            //ray.ne_tot[i][j][k] = 0;
            //ray.NH2_tot[i][j][k] = 0;
            }


    for(i = 0; i < 25; i++)
      for(j = 0; j < ray.N_phi[i]; j++)
        for(k = 0; k < 1; k++)
          {
            //printf("central dens = %lg\n", ray.nh_tot[i][j][0]);
          }


    for(i = 0; i < N_theta; i++)
      for(j = 0; j < ray.N_phi[i]; j++)
        {
          k = 0;

          N_photons_per_ray = ray.N_H_photons / (double)ray.N_rays;

          do
            {

              r_k_lower = ray.r_min * pow(ray.r_max / ray.r_min, k / (double)N_shells);

              r_k_upper = ray.r_min * pow(ray.r_max / ray.r_min, (k + 1) / (double)N_shells);

              dr_k = r_k_upper - r_k_lower;

              N_atoms = ABHE * 4.0*PI / (double)ray.N_rays * ray.r[k] * ray.r[k] * dr_k * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm * a3 / hubble3 * ray.nh_tot[i][j][k];

              N_recombinations = ray.alpha_B_He * dt_raytrace * (1.0 + 2.0 * ABHE) * ABHE * ray.nh_tot[i][j][k] * ray.nh_tot[i][j][k] * 4.0 * PI / (double)ray.N_rays * ray.r[k] * ray.r[k] * dr_k * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm * a3 / hubble3;

              d_N = N_recombinations;

              if(ray.He_front[i][j] > r_k_lower && ray.He_front[i][j] < r_k_upper)
                d_N += (r_k_upper - ray.He_front[i][j]) / dr_k * N_atoms;

              if(ray.He_front[i][j] <= r_k_lower)
                d_N += N_atoms;

              if(d_N > N_photons_per_ray)
                {
                  if(ray.He_front[i][j] > r_k_lower && ray.He_front[i][j] < r_k_upper)
                    {
                      if(N_photons_per_ray <= (ray.He_front[i][j] - r_k_lower) / dr_k * N_recombinations)
                        ray.He_front[i][j] = r_k_lower + N_photons_per_ray / ((ray.He_front[i][j] - r_k_lower) / dr_k * N_recombinations) * (ray.He_front[i][j] - r_k_lower);
                      else
                        ray.He_front[i][j] += (N_photons_per_ray - (ray.He_front[i][j] - r_k_lower) / dr_k * N_recombinations) / ((r_k_upper - ray.He_front[i][j]) / dr_k * (N_atoms + N_recombinations)) * (r_k_upper - ray.He_front[i][j]);
                    }
                  else
                    ray.He_front[i][j] = r_k_lower + N_photons_per_ray / d_N * dr_k;  //<-- This is why ifront is always just whatever sink smoothing link is!!!!!!!  (r_k_lower ~ sink_hsml)
                }

              N_photons_per_ray -= d_N;

              k++;
///////////////////////////////////////////////////////
//Athena's pretend shielding

              if(ray.nh_tot[i][j][k] > n_cut)
                {
                ray.He_front[i][j] = 0;
                N_photons_per_ray=0;
                }

///////////////////////////////////////////////////////
            }
          while(N_photons_per_ray > 0 && k < N_shells);

          if(ray.He_front[i][j] < ray.r_min)
            k = 0;
          else
            k = (int)(N_shells * log(ray.He_front[i][j] / ray.r_min) / log(ray.r_max / ray.r_min));

          N_photons_per_ray = ray.N_H_photons / (double)ray.N_rays;

          do
            {

              r_k_lower = ray.r_min * pow(ray.r_max / ray.r_min, k / (double)N_shells);

              r_k_upper = ray.r_min * pow(ray.r_max / ray.r_min, (k + 1) / (double)N_shells);

              dr_k = r_k_upper - r_k_lower;

              N_atoms = (1.0 + ABHE) * 4.0 * PI / (double)ray.N_rays * ray.r[k] * ray.r[k] * dr_k * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm * a3 / hubble3 * ray.nh_tot[i][j][k];

              N_recombinations = ray.alpha_B_H * dt_raytrace * (1.0 + ABHE) * (1.0 + ABHE) * ray.nh_tot[i][j][k] * ray.nh_tot[i][j][k] * 4.0 * PI / (double)ray.N_rays * ray.r[k] * ray.r[k] * dr_k * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm * a3 / hubble3;

              d_N = N_recombinations;

              if(ray.H_front[i][j] > r_k_lower && ray.H_front[i][j] < r_k_upper)
                d_N += (r_k_upper - ray.H_front[i][j]) / dr_k * N_atoms;

              if(ray.H_front[i][j] <= r_k_lower)
                d_N += N_atoms;

              if(d_N > N_photons_per_ray)
                {
                  if(ray.H_front[i][j] > r_k_lower && ray.H_front[i][j] < r_k_upper)
                    {
                      if(N_photons_per_ray <= (ray.H_front[i][j] - r_k_lower) / dr_k * N_recombinations)
                        ray.H_front[i][j] = r_k_lower + N_photons_per_ray / ((ray.H_front[i][j] - r_k_lower) / dr_k * N_recombinations) * (ray.H_front[i][j] - r_k_lower);
                      else
                        ray.H_front[i][j] += (N_photons_per_ray - (ray.H_front[i][j] - r_k_lower) / dr_k * N_recombinations) / ((r_k_upper - ray.H_front[i][j]) / dr_k * (N_atoms + N_recombinations)) * (r_k_upper - ray.H_front[i][j]);
                    }
                  else
                    ray.H_front[i][j] = r_k_lower + N_photons_per_ray / d_N * dr_k;
                }

              N_photons_per_ray -= d_N;

              k++;
///////////////////////////////////////////////////////
//Athena's pretend shielding

              if(ray.nh_tot[i][j][k] > n_cut)
                {
                ray.H_front[i][j] = 0;
                N_photons_per_ray=0;
                }

///////////////////////////////////////////////////////
            }
          while(N_photons_per_ray > 0 && k < N_shells);
        }

    outfile2=fopen("NH2.dat", "a");

    for(n = 0; n < N_gas; n++)
      if(P[n].ID > 0)
        {
          nh_here = SphP[n].Density*All.UnitDensity_in_cgs*hubble2/a3*HYDROGEN_MASSFRAC/PROTONMASS;

          sim_r = sqrt((P[n].Pos[0]-ray.center_x)*(P[n].Pos[0]-ray.center_x)+(P[n].Pos[1]-ray.center_y)*(P[n].Pos[1]-ray.center_y)+(P[n].Pos[2]-ray.center_z)*(P[n].Pos[2]-ray.center_z));

          if(sim_r == 0.0)
            sim_theta = PI/2.0;
          else
            sim_theta = acos((P[n].Pos[2]-ray.center_z)/sim_r);

          if(P[n].Pos[0]-ray.center_x > 0.0 && P[n].Pos[1]-ray.center_y >= 0.0)
            sim_phi = atan((P[n].Pos[1]-ray.center_y)/(P[n].Pos[0]-ray.center_x));
          else if(P[n].Pos[0]-ray.center_x < 0.0 && P[n].Pos[1]-ray.center_y >= 0.0)
            sim_phi = PI+atan((P[n].Pos[1]-ray.center_y)/(P[n].Pos[0]-ray.center_x));
          else if(P[n].Pos[0]-ray.center_x < 0.0 && P[n].Pos[1]-ray.center_y < 0.0)
            sim_phi = PI+atan((P[n].Pos[1]-ray.center_y)/(P[n].Pos[0]-ray.center_x));
          else if(P[n].Pos[0]-ray.center_x > 0.0 && P[n].Pos[1]-ray.center_y < 0.0)
            sim_phi = 2.0*PI+atan((P[n].Pos[1]-ray.center_y)/(P[n].Pos[0]-ray.center_x));
          else if(P[n].Pos[0]-ray.center_x == 0.0 && P[n].Pos[1]-ray.center_y > 0.0)
            sim_phi = PI/2.0;
          else if(P[n].Pos[0]-ray.center_x == 0.0 && P[n].Pos[1]-ray.center_y < 0.0)
            sim_phi = 3.0/2.0*PI;
          else
            sim_phi = PI;

          i = sim_theta/PI*N_theta;

          j = sim_phi/2.0/PI*ray.N_phi[i];

          if(sim_r <= ray.H_front[i][j] && nh_here < n_cut)
            {
              if(sim_r < ray.r_min)
                sim_r = ray.r_min;

#ifdef CHEMCOOL
              if(SphP[n].Ray_H_coeff < 0.0)
                {
                  r_here = sim_r * a / hubble * 1.0e3;

                  SphP[n].TracAbund[IHP] = 1.0 / (1.0 + ray.alpha_B_H / ray.HI_ion_rate * r_here * r_here * nh_here);
                  SphP[n].TracAbund[IDP] = All.DeutAbund * SphP[n].TracAbund[IHP];
                  SphP[n].TracAbund[IHEP] = ABHE / (1.0 + ray.alpha_B_H / ray.HeI_ion_rate * r_here * r_here * nh_here);
                  SphP[n].TracAbund[IHEPP] = 0.0;

                  mu_here = 1.0 / (HYDROGEN_MASSFRAC * (1.0 + ABHE + SphP[n].TracAbund[IHP] + SphP[n].TracAbund[IHEP] + 2.0 * SphP[n].TracAbund[IHEPP]));

                  gamma_here = 5.0/3.0;

                  entropy_here = BOLTZMANN / mu_here / PROTONMASS * All.UnitMass_in_g / All.UnitEnergy_in_cgs * ray.temp_HII_estimate / pow(SphP[n].Density / a3, gamma_here - 1.0);

                  SphP[n].Entropy = entropy_here;
                }
#endif
              SphP[n].Ray_H_coeff = sim_r * a / hubble * 1.0e3;

              //SphP[n].Prad = (1.6e-13*nh_here/sigma_ion)*(nu_min_H*6.626e-27/3.e10);
              //SphP[n].Prad = 1.6e-13*pow(nh_here*SphP[n].TracAbund[IHP],2)*
              //               (SphP[i].Hsml*a/hubble*1.0e3*3.0857e18)*(nu_min_H*6.626e-27/3.e10);
              SphP[n].Prad = 1.6e-13*pow(nh_here*SphP[n].TracAbund[IHP],1)*(P[n].Mass*1.e10*2e33/.7/1.67e-24)*
                             pow((SphP[i].Hsml*a/hubble*1.0e3*3.0857e18),-2)*(nu_min_H*6.626e-27/3.e10);
              SphP[n].Prad = SphP[n].Prad * pow(a3, SphP[n].Gamma) / All.UnitPressure_in_cgs / (All.HubbleParam*All.HubbleParam);

              x_shift = P[n].Pos[0]-ray.center_x;
              y_shift = P[n].Pos[1]-ray.center_y;
              z_shift = P[n].Pos[2]-ray.center_z;
              phi_shift = atan(y_shift/x_shift);
              theta_shift = acos(z_shift/sim_r);

              //prad_dir = (1.-SphP[n].TracAbund[IHP])*ray.Q_H_ion * (nu_min_H*6.626e-27/3.e10) / (4.*3.14159) * pow(sim_r*a/hubble*1.0e3*3.0857e18,-2);
              prad_dir = nh_here*(1.-SphP[n].TracAbund[IHP])*ray.HI_ion_rate/SphP[n].Ray_H_coeff/SphP[n].Ray_H_coeff*
                         (SphP[i].Hsml*a/hubble*1.0e3*3.0857e18)*(nu_min_H*6.626e-27/3.e10);

              SphP[n].Prad_dir[0] = cos(phi_shift)*sin(theta_shift)*prad_dir;
              SphP[n].Prad_dir[1] = sin(phi_shift)*sin(theta_shift)*prad_dir;
              SphP[n].Prad_dir[2] = cos(theta_shift)*prad_dir;

              for(m = 0; m < 3; m++)
                SphP[n].Prad_dir[m] = SphP[n].Prad_dir[m] * pow(a3, SphP[n].Gamma) / All.UnitPressure_in_cgs / (All.HubbleParam*All.HubbleParam); 

              //SphP[n].Prad = 0;
              printf("raytracing Prad = %lg Prad_dir[0]=%lg Prad_dir[1]=%lg Prad_dir[2]=%lg, phi_shift = %lg, theta_shift = %lg\n", SphP[n].Prad, SphP[n].Prad_dir[0], SphP[n].Prad_dir[1], SphP[n].Prad_dir[2], phi_shift, theta_shift);

              N_H_ionized_parts_local++;
            }
          else
            {
            SphP[n].Ray_H_coeff = -1.0;
            SphP[n].Prad = 0;
            for(m = 0; m < 3; m++)
              SphP[n].Prad_dir[m] = 0;
            }

          if(sim_r <= ray.He_front[i][j] && nh_here < n_cut)  //just add an outright density constraint
            {
              if(sim_r < ray.r_min)
                sim_r = ray.r_min;
#ifdef CHEMCOOL
              if(SphP[n].Ray_He_coeff < 0.0)  
                {
                  r_here = sim_r * a / hubble * 1.0e3;

                  SphP[n].TracAbund[IHP] = 1.0 / (1.0 + ray.alpha_B_H / ray.HI_ion_rate * r_here * r_here * nh_here);
                  SphP[n].TracAbund[IDP] = All.DeutAbund * SphP[n].TracAbund[IHP];
                  SphP[n].TracAbund[IHEPP] = ABHE / (1.0 + ray.alpha_B_He / ray.HeII_ion_rate * r_here * r_here * nh_here);
                  SphP[n].TracAbund[IHEP] = ABHE - SphP[n].TracAbund[IHEPP];

                  mu_here = 1.0 / (HYDROGEN_MASSFRAC * (1.0 + ABHE + SphP[n].TracAbund[IHP] + SphP[n].TracAbund[IHEP] + 2.0 * SphP[n].TracAbund[IHEPP]));

                  gamma_here = 5.0/3.0;

                  entropy_here = BOLTZMANN / mu_here / PROTONMASS * All.UnitMass_in_g / All.UnitEnergy_in_cgs * ray.temp_HeIII_estimate / pow(SphP[n].Density / a3, gamma_here - 1.0);

                  SphP[n].Entropy = entropy_here;
                }
#endif
              SphP[n].Ray_He_coeff = sim_r * a / hubble * 1.0e3;

              N_He_ionized_parts_local++;
            }
          else
            SphP[n].Ray_He_coeff = -1.0;

          if(sim_r < ray.r_min)
            sim_r = ray.r_min;
#ifdef CHEMCOOL
          if(SphP[n].Ray_LW_coeff < 0.0)
            {
              r_here = sim_r * a / hubble * 1.0e3;

              if(sim_r <= ray.H_front[i][j])
                SphP[n].TracAbund[IH2] = 0.0;
              else
                {

                 mu_here = 1.0 / (HYDROGEN_MASSFRAC * (1.0 + ABHE + SphP[n].TracAbund[IHP] + SphP[n].TracAbund[IHEP] + 2.0 * SphP[n].TracAbund[IHEPP]));
                 HI = fmax(1.0 - 2.0 * SphP[n].TracAbund[IH2] - SphP[n].TracAbund[IHP] - SphP[n].TracAbund[IHD], 0.0);      

                 if(SphP[n].TracAbund[IH2] > 0.1)
                   {
                   mu_here=1.2195;
                   h2frac=2.0*SphP[n].TracAbund[IH2];

                   muh2in=(0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0);
                   muh2=pow(muh2in, -1.0);

                   if(muh2 >= 1.22)
                     {
                     mu_here=muh2;
                     }
                   }

                  //ARS accounting for 3-body rates here
                  temp_here = (SphP[n].Gamma-1.0)/BOLTZMANN * (SphP[n].Entropy/(SphP[n].Gamma-1.0))*pow((SphP[n].Density/a3),(SphP[n].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs) * PROTONMASS * mu_here;
                  
                  h2_formrate =  HI * nh_here * (1.3e-9 * SphP[n].HM + 6.4e-10 * SphP[n].H2II)
                                 + (5.5e-29/temp_here)*HI*HI*HI*nh_here*nh_here 
                                 + (0.125e0*5.5e-29/temp_here)*HI*HI*SphP[n].TracAbund[IH2]*nh_here*nh_here; 

                
                   rate_LW = h2_formrate - SphP[n].TracAbund[IH2] * ray.LW_rate / r_here / r_here;

                   eq_LW = h2_formrate / (ray.LW_rate / r_here / r_here);

                   t_LW = fmin(fabs((eq_LW - SphP[n].TracAbund[IH2]) / rate_LW), 1.0e30);

                   //SphP[n].TracAbund[IH2] = SphP[n].TracAbund[IH2] * exp(- dt_raytrace / t_LW) + eq_LW * (1.0 - exp(- dt_raytrace / t_LW));
                }

              SphP[n].TracAbund[IHD] = 3.0 * All.DeutAbund * SphP[n].TracAbund[IH2];
            }
#endif
           SphP[n].Ray_LW_coeff = sim_r * a / hubble * 1.0e3;


	   /////////////////////////////////////////////////////////////
           //add up H2 column density here
           /////////////////////////////////////////////////////////////
           //sim_NH2 = 1.e0;
           sim_NH2 = ray.NH2_tot[i][j][0]* a / hubble * 1.0e3 *3.0857e18;
           //if(nh_here > 1.e2)
            for(k = 0; k < N_shells; k++)
               if(sim_r > ray.r[k])
                 sim_NH2 += ray.NH2_tot[i][j][k]* a / hubble * 1.0e3 *3.0857e18;

           SphP[n].Ray_NH2 = sim_NH2;
           if(/*SphP[n].Ray_NH2 < 1.e0 ||*/ SphP[n].Ray_NH2 != SphP[n].Ray_NH2)
             {
             printf("low H2! Ray_NH2 = %lg, dens = %lg, sim_r = %lg\n", SphP[n].Ray_NH2, nh_here, sim_r);
             sim_NH2 = 1.e0;
             }

           if(nh_here > 1.e14)
             {
             printf("high H2! Ray_NH2 = %lg, dens = %lg, sim_r = %lg\n", SphP[n].Ray_NH2, nh_here, sim_r);
             }


           if(nh_here > 1.e4 && nh_here < 1.e6 && All.NumCurrentTiStep == 30)
             {
             NH2_fac = SphP[n].Ray_NH2/5.e14;
             e_fac = 8.5e-4*pow((1.0 + NH2_fac),0.5);
             fshield = 0.965/(1.0 +pow((NH2_fac/3.0),2));
             fshield = fshield + 0.035/(pow((1.0 + NH2_fac),0.5))/exp(e_fac);

             //fprintf(outfile2, "%8d %5d %5d %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %8d\n", All.NumCurrentTiStep, i, j, nh_here, sim_r, sim_NH2, e_fac, fshield,  P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], SphP[n].Hsml, P[n].ID);  
             }
             
           if(nh_here > 1.e6 && nh_here < 1.e12 && All.NumCurrentTiStep == 40)
             {
             NH2_fac = SphP[n].Ray_NH2/5.e14;
             e_fac = 8.5e-4*pow((1.0 + NH2_fac),0.5);
             fshield = 0.965/(1.0 +pow((NH2_fac/3.0),2));
             fshield = fshield + 0.035/(pow((1.0 + NH2_fac),0.5))/exp(e_fac);
    

            // if(SphP[n].Ray_NH2 < 2.e0)
            //   printf("low H2!  nh = %lg, sim_r = %lg, sim_NH2 = %lg, e_fac = %lg, fshield = %lg\n", nh_here, sim_r, sim_NH2, e_fac, fshield);
 
             //fprintf(outfile2, "%8d %5d %5d %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %8d\n", All.NumCurrentTiStep, i, j, nh_here, sim_r, sim_NH2, e_fac, fshield, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], SphP[n].Hsml, P[n].ID);
             }
          
	  //////////////////////////////////////////////////////////          

        }
    
    fclose(outfile2);

    if(ray.flag_continue_end == 1 && ThisTask == 0)
 //     if(All.NumCurrentTiStep > 1000)
      {
        if(ray.flag_switch == 0)
          {
            sprintf(buf, "%s/%s.front1", All.OutputDir, All.SnapshotFileBase);

            ray.flag_switch = 1;
          }
        else
          {
            sprintf(buf, "%s/%s.front2", All.OutputDir, All.SnapshotFileBase);

            ray.flag_switch = 0;
          }

        outfile = fopen(buf, "w");

        for(i = 0; i < N_theta; i++)
          for(j = 0; j < ray.N_phi[i]; j++)
            fprintf(outfile, "%f %f\n", ray.H_front[i][j], ray.He_front[i][j]);

        fclose(outfile);
      }

    MPI_Allreduce(&N_H_ionized_parts_local, &N_H_ionized_parts_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&N_He_ionized_parts_local, &N_He_ionized_parts_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(&flag_bug, &flag_bug_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    t_end = second();

    if(ThisTask == 0)
      {
        printf("HII particles = %d\n", N_H_ionized_parts_tot);
        printf("HeIII particles = %d\n", N_He_ionized_parts_tot);
        printf("Time consumed = %g\n", t_end-t_start);
        printf("Raytracing done!\n\n");
      }

    if(flag_bug_tot != 0)
      {
        savepositions(All.SnapshotFileCount++);

        exit(0);
      }
  }
	#endif
