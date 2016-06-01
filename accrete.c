#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#define workdir /nobackupp1/astacy/

accrete(void)
{
 MPI_Barrier(MPI_COMM_WORLD);
 int j=0, sink_tot, i, k, rank_here, num, num_tot, my_rank, local_accID_arr[10], local_sink_num, sink_num_max;

 for(i = 0; i < N_gas; i++)
   if(SphP[i].sink > 0.5)
     {
     j++;
     }

 MPI_Allreduce(&j, &sink_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 if(ThisTask == 0)
   printf("sink_tot = %d\n", sink_tot);

 MPI_Barrier(MPI_COMM_WORLD);

 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 num = 1;
 MPI_Allreduce(&num, &num_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

 All.sink_num=0;
 local_sink_num=0;
 
 MPI_Barrier(MPI_COMM_WORLD);

 //printf("my_rank1 = %d\n", my_rank);
 for(i = 0; i < N_gas; i++) 
   {
    if(SphP[i].sink > 0.5 /*&& my_rank == rank_here*/)
      {
      local_sink_num++;
      printf("Found a sinkie!, ID = %d\n", P[i].ID);
      }      
   }
 
 MPI_Barrier(MPI_COMM_WORLD);


 MPI_Allreduce(&local_sink_num, &sink_num_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

 if(my_rank == 0)
  {
   for(k=0; k < sink_tot; k++)
     {
     //printf("All.sink_num = %d, All.accID_arr = %d\n", All.sink_num, All.accID_arr[k]);
     //printf("sink_num_max = %d\n", sink_num_max);
     }
  } 

for(i = 0; i < N_gas; i++)
   if(SphP[i].sink > 0.5)
      {
      printf("new sinkval = %g \n", SphP[i].sink);
      }

 //printf("my_rank2 = %d\n", my_rank);
 for(k=0; k < sink_num_max; k++)
  {
  MPI_Barrier(MPI_COMM_WORLD);
  accrete_gas_AS(0);
  }
 MPI_Barrier(MPI_COMM_WORLD);

}




accrete_gas_AS(int accID)
{
  FILE *sinkmergers;
  FILE *sinkdat;
  FILE *sinkangmom;
  int i, j, k, h, g, l, num_sink_part, local_sink_count_this_timestep, PCUsink_rank, tot_sink_ngb, merger_number;
  int ncount, my_rank, tot_sink_count, sum_of_ranks, tot_PCUsink_rank, local_sink_ngb, y, y_tot, sink_ID, sink_ID_exchange, tot_sink_ID_exchange,m;
  double dt, hubble_a = 0, a3inv = 0, ne = 1, a3, num_dens_sink;
  double deltat, unew, ffH2, yH, mass_sink, TT_sink; 
  double x_sink, y_sink, z_sink, r_sink, d_sink, sink_Mass, ent_sink, pres_sink;
  double numpos_x, numpos_y, numpos_z, numvel_x, numvel_y, numvel_z, mu, u, numTT;
  double mass_ngb, numpos_ngb_x, numpos_ngb_y, numpos_ngb_z, numvel_ngb_x, numvel_ngb_y, numvel_ngb_z, numTT_ngb;
  double tot_mass_ngb, tot_numpos_ngb_x, tot_numpos_ngb_y, tot_numpos_ngb_z, tot_numvel_ngb_x, tot_numvel_ngb_y, tot_numvel_ngb_z, tot_numTT_ngb;
  double r_sink_phys, mu_sink, Temp, Temp_part, u_part, mu_part, kin_engy, pot_engy, vy_sink, vx_sink, vz_sink;
  double Temp_tot=0.0, Temp_avg;
  int n_acc=0;
  double  SinkCriticalDensity,  OuterAccrRadius, hubble_param, hubble_param2, a;

  double UnitEnergy_in_cgs, UnitTime_in_s, UnitMass_in_g, UnitVelocity_in_cm_per_s, UnitLength_in_cm;
  double jx, jy, jz, jtot, d_sink_phys, vx, vy, vz, vx_phys, vy_phys, vz_phys, jcent_sink;
  double muh2in, muh2, h2frac;
  int number;
  double cgs_fac;
  char fsinkmergers[200], fsinkdat[200], fsinkangmom[200];

  sprintf(fsinkdat,"/nobackupp1/astacy/sinkdat_bin");
  sprintf(fsinkangmom,"/nobackupp1/astacy/sinkangmom_bin");
  sprintf(fsinkmergers,"/nobackupp1/astacy/sinkmergers_bin");


UnitLength_in_cm= 3.085678e21;
UnitVelocity_in_cm_per_s= 1.0e5;
UnitMass_in_g= 1.989e43;
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

  MPI_Status status;


  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a = All.Time;
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = All.Hubble * sqrt(All.Omega0 / (All.Time * All.Time * All.Time)
				   + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) +
				   All.OmegaLambda);
      hubble_param = All.HubbleParam;
      hubble_param2 = hubble_param*hubble_param;
      
    }
  else
    a = a3inv = hubble_param = hubble_param2  = 1.0;

  a3 = 1.e0/a3inv;
  cgs_fac = 1.e5*3.086e18*1.e3*All.Time/(hubble_param); 

 
  SinkCriticalDensity = All.SinkCriticalDens / All.UnitDensity_in_cgs * PROTONMASS / HYDROGEN_MASSFRAC * a3 / hubble_param2;
  //SinkCriticalDensity = All.SinkCriticalDens / All.UnitDensity_in_cgs * PROTONMASS * 2.27 * a3 / hubble_param2;

  //OuterAccrRadius = All.ROuter / a * hubble_param;
  OuterAccrRadius = 2.0e0*All.SofteningGas;


  /* JJ -- Do a loop over all processors, to allow for multiple sinks to be created. */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  y = 1;
  MPI_Allreduce(&y, &y_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  for(l=0; l < y_tot; l++)
    {
      /* JJ -- In order to be sure that multiple processors are not making sinks at one time, an all_reduce statement. */
      MPI_Allreduce(&y, &y_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      

	  /* JJ - Here we set up the sink particles.  The criterion presently is just n > 4*10^7 cm^-3.  Only one sink particle can be created per timestep.  */
	  /* JJ - Set these to zero, in order that the Bcast command has something to send in the event that no sink is created. */
      
      tot_sink_ngb   = 0;
      local_sink_ngb = 0;
      local_sink_count_this_timestep = 0;
      tot_sink_count = 0;
      tot_PCUsink_rank = 0;
      PCUsink_rank = 0;
      sink_ID_exchange = 0;
      tot_sink_ID_exchange = 0;
      
      x_sink       = 0.e0;
      y_sink       = 0.e0;
      z_sink       = 0.e0;
      r_sink       = 0.e0;
      
      mass_ngb     = 0.e0;
      
      numpos_ngb_x = 0.e0;
      numpos_ngb_y = 0.e0;
      numpos_ngb_z = 0.e0;
      
      numvel_ngb_x = 0.e0;
      numvel_ngb_y = 0.e0;
      numvel_ngb_z = 0.e0;
      
      numTT_ngb = 0.e0;
      numTT     = 0.e0;
      Temp_tot = 0.0;      
      
      //printf("accrete.c line 190, N_gas = %d\n", N_gas);

      for(i = 1; i <= N_gas; i++)
	{      
          
          n_acc=0;

	  //if(SphP[i].sink > 0.5e0 && my_rank == l && local_sink_count_this_timestep < 0.5e0)
          //if(P[i].ID == accID  && my_rank == l && local_sink_count_this_timestep < 0.5e0)
          if(SphP[i].sink > 0.5e0 && SphP[i].sink < 1.5e0 && my_rank == l && local_sink_count_this_timestep < 0.5e0)
	    {
	     
             printf("ID OF SINK UNDERGOING ACCRETION: %d\n", P[i].ID);
 
	      num_sink_part = i;
	      
	      local_sink_count_this_timestep += 1;

              Temp_tot = 0.0; 

              SphP[i].Density = SinkCriticalDensity;

              mu=1.2195;
              h2frac=2.0*SphP[i].TracAbund[IH2];

              muh2in=(0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0);
              muh2=pow(muh2in, -1.0);

              if(muh2 >= 1.22){
               mu=muh2;
              }


              SphP[i].Entropy = 2000.0*BOLTZMANN / (pow(SinkCriticalDensity*a3inv,(SphP[i].Gamma - 1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs) * PROTONMASS * mu);
              SphP[i].Hsml = All.SofteningGas;
              SphP[i].Pressure = SphP[i].Entropy * pow(SphP[i].Density, SphP[i].Gamma);
 
	      SphP[i].sink  = 2.0e0; 

	      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	      
	      PCUsink_rank = my_rank;
	      
	      num_dens_sink = (SphP[i].Density * a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam)/(1.22 * PROTONMASS);
	      

	      	      /* These two are to be broadcast so that mergers can be tracked. */
              sink_ID = P[num_sink_part].ID;
              sink_Mass     = P[i].Mass;

	      mass_sink     = P[i].Mass;
	      
	      numpos_x  = P[i].Pos[0]*P[i].Mass;
	      numpos_y  = P[i].Pos[1]*P[i].Mass;
	      numpos_z  = P[i].Pos[2]*P[i].Mass;
	      
	      numvel_x  = P[i].Vel[0]*P[i].Mass; 
	      numvel_y  = P[i].Vel[1]*P[i].Mass;
	      numvel_z  = P[i].Vel[2]*P[i].Mass;

	      mass_ngb     = 0.e0;
	      
	      numpos_ngb_x = 0.e0;
	      numpos_ngb_y = 0.e0;
	      numpos_ngb_z = 0.e0;
	      
	      numvel_ngb_x = 0.e0;
	      numvel_ngb_y = 0.e0;
	      numvel_ngb_z = 0.e0;
	      
	      numTT_ngb = 0.e0;
	      
	      u=(SphP[i].Entropy/(SphP[i].Gamma-1.0))*pow((SphP[i].Density/a3),(SphP[i].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs);
	      

	      numTT = ((SphP[i].Gamma-1.0)/BOLTZMANN) * u * PROTONMASS * 1.22 * P[i].Mass;
	      
	      x_sink  = P[i].Pos[0];
	      y_sink  = P[i].Pos[1];
	      z_sink  = P[i].Pos[2];
              vx_sink = P[i].Vel[0];
              vy_sink = P[i].Vel[1];
              vz_sink = P[i].Vel[2];

              r_sink = OuterAccrRadius;

	      /* JJ -- Now calculate the Bondi accretion radius, to be used as r_sink. */

              Temp = (SphP[i].Gamma-1.0)/BOLTZMANN * u * PROTONMASS * 1.22;

	      r_sink_phys = r_sink * All.Time / All.HubbleParam;

	      
	      for(k = 1; k <= i-1; k++)
		{
		  d_sink = sqrt((P[k].Pos[0] - x_sink)*(P[k].Pos[0] - x_sink) + (P[k].Pos[1] - y_sink)*(P[k].Pos[1] - y_sink) + (P[k].Pos[2] - z_sink)*(P[k].Pos[2] - z_sink));

                  d_sink_phys = d_sink * 3.086e18*1.e3*All.Time/(hubble_param);  //convert to cgs!

                  vx = P[k].Vel[0] - vx_sink;
                  vy = P[k].Vel[1] - vy_sink;
                  vz = P[k].Vel[2] - vz_sink;

                  vx_phys = vx*pow(All.Time, 0.5)*1.e5;
                  vy_phys = vy*pow(All.Time, 0.5)*1.e5;
                  vz_phys = vz*pow(All.Time, 0.5)*1.e5;

                  jx = (P[k].Pos[1]-y_sink)*vz  - (P[k].Pos[2]-z_sink)*vy;
                  jy = (P[k].Pos[2]-z_sink)*vx  - (P[k].Pos[0]-x_sink)*vz;
                  jz = (P[k].Pos[0]-x_sink)*vy  - (P[k].Pos[1]-y_sink)*vx;

                  cgs_fac = 1.e5*3.086e18*1.e3*All.Time/(hubble_param);
                    
                  jx =  jx*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);  //convert to cgs units
                  jy =  jy*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);
                  jz =  jz*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);

                  jtot = pow(jx*jx + jy*jy + jz*jz, 0.5);
                  jcent_sink = pow((6.67e-8*sink_Mass*(1.e10/hubble_param)*1.98892e33*d_sink_phys),0.5);

	          u_part=(SphP[k].Entropy*pow((SphP[k].Density/a3),(SphP[k].Gamma-1.0))/(SphP[k].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs);

   
	          Temp_part=(SphP[k].Gamma-1.0)/BOLTZMANN * u_part * PROTONMASS * 1.22;
 
	          pot_engy =  GRAVITY * mass_sink * All.UnitMass_in_g * All.HubbleParam / (d_sink * All.UnitLength_in_cm * All.Time);
	          kin_engy =  0.5e0 * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s * ((P[k].Vel[0] - vx_sink)*(P[k].Vel[0] - vx_sink) + (P[k].Vel[1] - vy_sink)*(P[k].Vel[1] - vy_sink) + (P[k].Vel[2] - vz_sink)*(P[k].Vel[2] - vz_sink));

  	          if(d_sink < r_sink /*&& P[k].Ti_endstep == All.Ti_Current*/ && jtot < jcent_sink /*&& SphP[k].sink < 0.5*/ && SphP[k].sink != -1)
		    {
		      
                      n_acc=n_acc+1;
		      Temp_tot=Temp_tot+ (SphP[k].Entropy/(SphP[k].Gamma-1.0))*pow(SphP[k].Density*a3inv,(SphP[k].Gamma-1.0));

                      printf("before file print %15.6g %8d %15.6g %15.6g %15.11g \n", All.Time, n_acc, r_sink_phys, Temp_avg, jcent_sink);

                      Temp_avg=(SphP[k].Entropy/(SphP[k].Gamma-1.0))*pow(SphP[k].Density*a3inv,(SphP[k].Gamma-1.0));
                      sinkdat=fopen(fsinkdat,"a");                             
                      fprintf(sinkdat,"%17.13g %8d %15.6g %15.6g %15.6g %10d %10d %15.6g\n", All.Time, n_acc, r_sink_phys, (SphP[k].Entropy/(SphP[k].Gamma-1.0))*pow(SphP[k].Density*a3inv,(SphP[k].Gamma-1.0)), SphP[num_sink_part].Entropy, P[k].ID, P[num_sink_part].ID, SphP[num_sink_part].Pressure);
                      fclose(sinkdat);

                      sinkangmom=fopen(fsinkangmom,"a");                             
                      fprintf(sinkangmom,"%17.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %10d\n", All.Time, d_sink_phys, vx_phys, vy_phys, vz_phys, jx, jy, jz, jtot, sink_Mass, sink_ID);
                      fclose(sinkangmom);


                      printf("after file print %15.6g %8d %15.6g %15.6g \n", All.Time, n_acc, r_sink_phys, Temp_avg);
 
		      mass_sink += P[k].Mass;
		      
		      numpos_x  += P[k].Pos[0]*P[k].Mass;
		      numpos_y  += P[k].Pos[1]*P[k].Mass;
		      numpos_z  += P[k].Pos[2]*P[k].Mass;
		      
		      numvel_x  += P[k].Vel[0]*P[k].Mass; 
		      numvel_y  += P[k].Vel[1]*P[k].Mass;
		      numvel_z  += P[k].Vel[2]*P[k].Mass;
		      
		      u=(SphP[k].Entropy/(SphP[k].Gamma-1.0))*pow((SphP[k].Density/a3),(SphP[k].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs);


		      numTT += ((SphP[k].Gamma-1.0)/BOLTZMANN) * u * PROTONMASS * 1.22 * P[k].Mass;
		      

                      if(SphP[k].sink > 0.5e0)
			{

			  sinkmergers=fopen(fsinkmergers,"a");
			  fprintf(sinkmergers,"%17.13g %8d %8d %15.6g %15.6g \n", All.Time, P[k].ID, P[num_sink_part].ID, P[k].Mass, P[num_sink_part].Mass);
			  /*	  fflush(Fdsinkmergers);*/
			  fclose(sinkmergers);
			  
			  if(P[k].Mass > P[num_sink_part].Mass)
			    P[num_sink_part].ID = P[k].ID;
			  
                          printf("line 226 of accretemulti.c -- merger has occured. P[k].ID = %d.  P[num_sink_part].ID = %d.", P[k].ID,  P[num_sink_part].ID);

			}

		      
		      /* JJ -- The below is taken from Allvars.h.  We shift all of the Sph data, so that the particles that go into the sink are excluded.  */
		      
		      local_sink_ngb += 1;
		      
		      N_gas   -= 1;
		      NumPart -= 1;
		      All.TotN_gas   -= 1;
		      All.TotNumPart -= 1;
                      header.npart[0] -= 1;
                      header.npartTotal[0] -= 1;
                      All.MassTable[0] = 0.e0;		     
 
		      for(g = k; g <= NumPart; g++) 
			{

			  P[g].Pos[0] = P[g+1].Pos[0];
			  P[g].Pos[1] = P[g+1].Pos[1];           
			  P[g].Pos[2] = P[g+1].Pos[2];           

			  P[g].Mass   = P[g+1].Mass; 
			  
			  P[g].Vel[0] = P[g+1].Vel[0];   
			  P[g].Vel[1] = P[g+1].Vel[1]; 
			  P[g].Vel[2] = P[g+1].Vel[2]; 

                          P[g].GravAccel[0] = P[g+1].GravAccel[0];
                          P[g].GravAccel[1] = P[g+1].GravAccel[1];
                          P[g].GravAccel[2] = P[g+1].GravAccel[2];			  
          
#ifdef PMGRID
			  P[g].GravPM[0] = P[g+1].GravPM[0];
			  P[g].GravPM[1] = P[g+1].GravPM[1];           
			  P[g].GravPM[2] = P[g+1].GravPM[2];    
#endif
#ifdef FORCETEST     
			  P[g].GravAccelDirect[0] = P[g+1].GravAccelDirect[0];
			  P[g].GravAccelDirect[1] = P[g+1].GravAccelDirect[1];           
			  P[g].GravAccelDirect[2] = P[g+1].GravAccelDirect[2]; 
#endif
          
                          P[g].Potential = P[g+1].Potential;
                          P[g].OldAcc = P[g+1].OldAcc;
                          P[g].ID = P[g+1].ID;
                          P[g].Type = P[g+1].Type;
                          P[g].Ti_endstep = P[g+1].Ti_endstep;
                          P[g].Ti_begstep = P[g+1].Ti_begstep;
#ifdef FLEXSTEPS 
                          P[g].FlexStepGrp = P[g+1].FlexStepGrp;
#endif
                          P[g].GravCost = P[g+1].GravCost;
#ifdef PSEUDOSYMMETRIC
                          P[g].AphysOld = P[g+1].AphysOld;
#endif 
			}
		      
		      for(g = k; g <= N_gas; g++) 
			{	       
#ifndef POLYTROPE
                          SphP[g].Entropy = SphP[g+1].Entropy;
#endif
                          SphP[g].Density = SphP[g+1].Density;
                          SphP[g].Hsml = SphP[g+1].Hsml;
                          SphP[g].Left = SphP[g+1].Left;                        
                          SphP[g].Right = SphP[g+1].Right;
                          SphP[g].NumNgb = SphP[g+1].NumNgb;	
                          SphP[g].Pressure = SphP[g+1].Pressure;
                          SphP[g].Prad = SphP[g+1].Prad;
                          for(number = 0; number < 3; number++)
                            SphP[g].Prad_dir[number] = SphP[g+1].Prad_dir[number];

#ifndef POLYTROPE
                          SphP[g].DtEntropy = SphP[g+1].DtEntropy; 
#endif
                          for(number = 0; number < 3; number++)
			    SphP[g].HydroAccel[number] = SphP[g+1].HydroAccel[number];

                          for(number = 0; number < 3; number++)
			    SphP[g].VelPred[number] = SphP[g+1].VelPred[number];

                          SphP[g].DivVel = SphP[g+1].DivVel;
                          SphP[g].CurlVel = SphP[g+1].CurlVel;

                          for(number = 0; number < 3; number++)
			    SphP[g].Rot[number] = SphP[g+1].Rot[number];

                          SphP[g].DhsmlDensityFactor = SphP[g+1].DhsmlDensityFactor;
                          SphP[g].MaxSignalVel = SphP[g+1].MaxSignalVel;
#ifdef CHEMCOOL
                          SphP[g].DtEntropyVisc = SphP[g+1].DtEntropyVisc;
                          SphP[g].Gamma = SphP[g+1].Gamma;

                          for(number = 0; number < TRAC_NUM; number++)
			    SphP[g].TracAbund[number] = SphP[g+1].TracAbund[number];

                          SphP[g].EntropyOut = SphP[g+1].EntropyOut;

                          for(number = 0; number < TRAC_NUM; number++)
	 		    SphP[g].TracAbundOut[number] = SphP[g+1].TracAbundOut[number];

                          SphP[g].DustTemp = SphP[g+1].DustTemp;
                          SphP[g].HM = SphP[g+1].HM;
                          SphP[g].H2II = SphP[g+1].H2II;
#ifdef RAYTRACE
                          for(number = 0; number < 6; number++)
	 		    SphP[g].TotalColumnDensity[number] = SphP[g+1].TotalColumnDensity[number];

                          for(number = 0; number < 6; number++)
	 		    SphP[g].H2ColumnDensity[number] = SphP[g+1].H2ColumnDensity[number];

                          for(number = 0; number < 6; number++)
	 		    SphP[g].COColumnDensity[number] = SphP[g+1].COColumnDensity[number];
#endif /* RAYTRACE */
#endif /* CHEMCOOL */
#ifdef RAYTRACE_TG
                          SphP[g].Ray_H_coeff = SphP[g+1].Ray_H_coeff;
                          SphP[g].Ray_He_coeff = SphP[g+1].Ray_He_coeff;
                          SphP[g].Ray_LW_coeff = SphP[g+1].Ray_LW_coeff;
#endif
#ifdef METALS_TG
                          SphP[g].Sigma = SphP[g+1].Sigma;
                          SphP[g].const_A = SphP[g+1].const_A;
                          SphP[g].const_B = SphP[g+1].const_B;
                          SphP[g].Metallicity = SphP[g+1].Metallicity;
#endif
                          SphP[g].sink = SphP[g+1].sink;
			  
			}
		      /*		  printf("k before = %d.  i before = %d.", k, i);*/
		      
		      k -= 1;
		      i -= 1;
		      
		      /*		  printf("k after = %d. i after = %d", k, i);*/
		    }
		}

	      
	      /* JJ -- update the number of the sink particle, as it may have changed in the shifting of the struct already carried out. */
	      num_sink_part = num_sink_part - local_sink_ngb;
	      
	      for(k = i+1; k <= N_gas; k++)
		{
	      d_sink = sqrt((P[k].Pos[0] - x_sink)*(P[k].Pos[0] - x_sink) + (P[k].Pos[1] - y_sink)*(P[k].Pos[1] - y_sink) + (P[k].Pos[2] - z_sink)*(P[k].Pos[2] - z_sink));

              d_sink_phys = d_sink * 3.086e18*1.e3*All.Time/(hubble_param);

              vx = P[k].Vel[0] - vx_sink;
              vy = P[k].Vel[1] - vy_sink;
              vz = P[k].Vel[2] - vz_sink;

              vx_phys = vx*pow(All.Time, 0.5)*1.e5;
              vy_phys = vy*pow(All.Time, 0.5)*1.e5;
              vz_phys = vz*pow(All.Time, 0.5)*1.e5;

              jx = (P[k].Pos[1]-y_sink)*vz  - (P[k].Pos[2]-z_sink)*vy;
              jy = (P[k].Pos[2]-z_sink)*vx  - (P[k].Pos[0]-x_sink)*vz;
              jz = (P[k].Pos[0]-x_sink)*vy  - (P[k].Pos[1]-y_sink)*vx;

              jx =  jx*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);  //convert to cgs units
              jy =  jy*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);
              jz =  jz*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);

              jtot = pow(jx*jx + jy*jy + jz*jz, 0.5);
              jcent_sink = pow((6.67e-8*sink_Mass*(1.e10/hubble_param)*1.98892e33*d_sink_phys),0.5);

	      u_part=(SphP[k].Entropy*pow((SphP[k].Density/a3),(SphP[k].Gamma-1.0))/(SphP[k].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs);
   

	      Temp_part=(SphP[k].Gamma-1.0)/BOLTZMANN * u_part * PROTONMASS * 1.22;

	      pot_engy =  GRAVITY * mass_sink * All.UnitMass_in_g * All.HubbleParam / (d_sink * All.UnitLength_in_cm * All.Time);
	      kin_engy =  0.5e0 * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s * ((P[k].Vel[0] - vx_sink)*(P[k].Vel[0] - vx_sink) + (P[k].Vel[1] - vy_sink)*(P[k].Vel[1] - vy_sink) + (P[k].Vel[2] - vz_sink)*(P[k].Vel[2] - vz_sink));

	      if(d_sink < r_sink /*&& P[k].Ti_endstep == All.Ti_Current*/ && jtot < jcent_sink /*&& SphP[k].sink < 0.5*/ && SphP[k].sink != -1)
		  { 

                  n_acc=n_acc+1;
        	  Temp_tot=Temp_tot+ (SphP[k].Entropy/(SphP[k].Gamma-1.0))*pow(SphP[k].Density*a3inv,(SphP[k].Gamma-1.0));
 
                  printf("before file print %15.6g %8d %15.6g %15.6g %15.11g \n", All.Time, n_acc, r_sink_phys, Temp_avg, jcent_sink);

                  Temp_avg=(SphP[k].Entropy/(SphP[k].Gamma-1.0))*pow(SphP[k].Density*a3inv,(SphP[k].Gamma-1.0));
                  sinkdat=fopen(fsinkdat,"a");
                  fprintf(sinkdat,"%17.13g %8d %15.6g %15.6g %15.6g %10d %10d %15.6g\n", All.Time, n_acc, r_sink_phys, (SphP[k].Entropy/(SphP[k].Gamma-1.0))*pow(SphP[k].Density*a3inv,(SphP[k].Gamma-1.0)), SphP[num_sink_part].Entropy, P[k].ID, P[num_sink_part].ID, SphP[num_sink_part].Pressure);
                  fclose(sinkdat);



                  sinkangmom=fopen(fsinkangmom,"a");
                  fprintf(sinkangmom,"%17.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %10d\n", All.Time, d_sink_phys, vx_phys, vy_phys, vz_phys, jx, jy, jz, jtot, sink_Mass, sink_ID);
                  fclose(sinkangmom);


                  printf("after file print %15.6g %8d %15.6g %15.6g \n", All.Time, n_acc, r_sink_phys, Temp_avg); 

                  if(SphP[k].sink > 0.5e0)
			{
				  sinkmergers=fopen(fsinkmergers,"a");
				  fprintf(sinkmergers,"%17.13g %8d %8d %15.6g %15.6g \n", All.Time, P[k].ID, P[num_sink_part].ID, P[k].Mass, P[num_sink_part].Mass);
				  /*  fflush(Fdsinkmergers);*/
				  fclose(sinkmergers);
				  
				  if(P[k].Mass > P[num_sink_part].Mass)
				    P[num_sink_part].ID = P[k].ID;
				  
				  merger_number = All.merger_number;
				  All.merger_historyTIME[merger_number] = All.Time;
				  All.merger_historyIDGONE[merger_number] = P[k].ID;
				  All.merger_historyIDGROW[merger_number] = P[num_sink_part].ID;
				  All.merger_historyMASSGONE[merger_number] = P[k].Mass;
				  All.merger_historyMASSGROW[merger_number] = P[num_sink_part].Mass;
				  All.merger_number += 1;

				  printf("line 450 of accretemulti.c -- merger has occured. P[k].ID = %d.  P[num_sink_part].ID = %d.", P[k].ID,  P[num_sink_part].ID);

				}


			      
			      mass_sink += P[k].Mass;
			      
			      numpos_x  += P[k].Pos[0]*P[k].Mass;
			      numpos_y  += P[k].Pos[1]*P[k].Mass;
			      numpos_z  += P[k].Pos[2]*P[k].Mass;
			      
			      numvel_x  += P[k].Vel[0]*P[k].Mass; 
			      numvel_y  += P[k].Vel[1]*P[k].Mass;
			      numvel_z  += P[k].Vel[2]*P[k].Mass;
			      
			      u=(SphP[k].Entropy/(SphP[k].Gamma-1.0))*pow((SphP[k].Density/a3),(SphP[k].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs);

			      numTT += ((SphP[k].Gamma-1.0)/BOLTZMANN) * u * PROTONMASS * 1.22 * P[k].Mass;
			      
			      
			      /* JJ -- The below is taken from Allvars.h.  We shift all of the Sph data, so that the particles that go into the sink are excluded.  */
			      
			      local_sink_ngb += 1;
			      
			      N_gas   -= 1;
			      NumPart -= 1;
			      All.TotN_gas   -= 1;
			      All.TotNumPart -= 1;
			      header.npart[0] -= 1;
			      header.npartTotal[0] -= 1;		     
			      All.MassTable[0] = 0.e0;
	 
			      for(g = k; g <= NumPart; g++) 
				{
				  
				  P[g].Pos[0] = P[g+1].Pos[0];
				  P[g].Pos[1] = P[g+1].Pos[1];           
				  P[g].Pos[2] = P[g+1].Pos[2];           

				  P[g].Mass   = P[g+1].Mass; 
				  
				  P[g].Vel[0] = P[g+1].Vel[0];   
				  P[g].Vel[1] = P[g+1].Vel[1]; 
				  P[g].Vel[2] = P[g+1].Vel[2]; 

				  P[g].GravAccel[0] = P[g+1].GravAccel[0];
				  P[g].GravAccel[1] = P[g+1].GravAccel[1];
				  P[g].GravAccel[2] = P[g+1].GravAccel[2];			  
		  
	#ifdef PMGRID
				  P[g].GravPM[0] = P[g+1].GravPM[0];
				  P[g].GravPM[1] = P[g+1].GravPM[1];           
				  P[g].GravPM[2] = P[g+1].GravPM[2];    
	#endif
	#ifdef FORCETEST     
				  P[g].GravAccelDirect[0] = P[g+1].GravAccelDirect[0];
				  P[g].GravAccelDirect[1] = P[g+1].GravAccelDirect[1];           
				  P[g].GravAccelDirect[2] = P[g+1].GravAccelDirect[2]; 
	#endif
		  
				  P[g].Potential = P[g+1].Potential;
				  P[g].OldAcc = P[g+1].OldAcc;
				  P[g].ID = P[g+1].ID;
				  P[g].Type = P[g+1].Type;
				  P[g].Ti_endstep = P[g+1].Ti_endstep;
				  P[g].Ti_begstep = P[g+1].Ti_begstep;
	#ifdef FLEXSTEPS 
				  P[g].FlexStepGrp = P[g+1].FlexStepGrp;
	#endif
				  P[g].GravCost = P[g+1].GravCost;
	#ifdef PSEUDOSYMMETRIC
				  P[g].AphysOld = P[g+1].AphysOld;
	#endif 
				}
			      
			      for(g = k; g <= N_gas; g++) 
				{	       
	#ifndef POLYTROPE
				  SphP[g].Entropy = SphP[g+1].Entropy;
	#endif
				  SphP[g].Density = SphP[g+1].Density;
				  SphP[g].Hsml = SphP[g+1].Hsml;
				  SphP[g].Left = SphP[g+1].Left;                        
				  SphP[g].Right = SphP[g+1].Right;
				  SphP[g].NumNgb = SphP[g+1].NumNgb;	
				  SphP[g].Pressure = SphP[g+1].Pressure;
                                  SphP[g].Prad = SphP[g+1].Prad;
                                  for(number = 0; number < 3; number++)
                                    SphP[g].Prad_dir[number] = SphP[g+1].Prad_dir[number];

	#ifndef POLYTROPE
				  SphP[g].DtEntropy = SphP[g+1].DtEntropy; 
	#endif
				  for(number = 0; number < 3; number++)
				    SphP[g].HydroAccel[number] = SphP[g+1].HydroAccel[number];

				  for(number = 0; number < 3; number++)
				    SphP[g].VelPred[number] = SphP[g+1].VelPred[number];

				  SphP[g].DivVel = SphP[g+1].DivVel;
				  SphP[g].CurlVel = SphP[g+1].CurlVel;

				  for(number = 0; number < 3; number++)
				    SphP[g].Rot[number] = SphP[g+1].Rot[number];

				  SphP[g].DhsmlDensityFactor = SphP[g+1].DhsmlDensityFactor;
				  SphP[g].MaxSignalVel = SphP[g+1].MaxSignalVel;
	#ifdef CHEMCOOL
				  SphP[g].DtEntropyVisc = SphP[g+1].DtEntropyVisc;
				  SphP[g].Gamma = SphP[g+1].Gamma;

				  for(number = 0; number < TRAC_NUM; number++)
				    SphP[g].TracAbund[number] = SphP[g+1].TracAbund[number];

				  SphP[g].EntropyOut = SphP[g+1].EntropyOut;

				  for(number = 0; number < TRAC_NUM; number++)
				    SphP[g].TracAbundOut[number] = SphP[g+1].TracAbundOut[number];

				  SphP[g].DustTemp = SphP[g+1].DustTemp;
				  SphP[g].HM = SphP[g+1].HM;
				  SphP[g].H2II = SphP[g+1].H2II;
	#ifdef RAYTRACE
				  for(number = 0; number < 6; number++)
				    SphP[g].TotalColumnDensity[number] = SphP[g+1].TotalColumnDensity[number];

				  for(number = 0; number < 6; number++)
				    SphP[g].H2ColumnDensity[number] = SphP[g+1].H2ColumnDensity[number];

				  for(number = 0; number < 6; number++)
				    SphP[g].COColumnDensity[number] = SphP[g+1].COColumnDensity[number];
	#endif /* RAYTRACE */
	#endif /* CHEMCOOL */
	#ifdef RAYTRACE_TG
				  SphP[g].Ray_H_coeff = SphP[g+1].Ray_H_coeff;
				  SphP[g].Ray_He_coeff = SphP[g+1].Ray_He_coeff;
				  SphP[g].Ray_LW_coeff = SphP[g+1].Ray_LW_coeff;
	#endif
	#ifdef METALS_TG
				  SphP[g].Sigma = SphP[g+1].Sigma;
				  SphP[g].const_A = SphP[g+1].const_A;
				  SphP[g].const_B = SphP[g+1].const_B;
				  SphP[g].Metallicity = SphP[g+1].Metallicity;
	#endif
				  SphP[g].sink = SphP[g+1].sink;
				}	  
			      k -= 1;
			    }
			}
		      
		    }
		}
	      
	      
	      /* Check that only one sink particle was created.  If not, then the necessary MPI gets too complicated for this scheme, and the code needs to stop. */
	      
	      MPI_Allreduce(&local_sink_count_this_timestep, &tot_sink_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	      
	      
      if(tot_sink_count > 0.5e0)
	{
	  
	  /* Check other processors for neighboring particles that must be added to the sink particle, if a sink is created. */
	  
	  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
	  
	  /* Tell the other processors where the new sink particle is and within what distance from it particles must incorporated into it.*/
	  
	  /* JJ - This line will allow all processors to find out which processors has the sink. */
	  MPI_Allreduce(&PCUsink_rank, &tot_PCUsink_rank, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);   
	  PCUsink_rank = tot_PCUsink_rank;       
	  
	  MPI_Bcast(&x_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
	  MPI_Bcast(&y_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
	  MPI_Bcast(&z_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
          MPI_Bcast(&vx_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
          MPI_Bcast(&vy_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
          MPI_Bcast(&vz_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
          MPI_Bcast(&pres_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
          MPI_Bcast(&ent_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
	  MPI_Bcast(&r_sink, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
	  MPI_Bcast(&sink_ID, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
	  MPI_Bcast(&sink_Mass, 1, MPI_DOUBLE, PCUsink_rank, MPI_COMM_WORLD);
	  
	  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
	  
	  if(my_rank != PCUsink_rank)
	    {
	      for(h = 1; h <= N_gas; h++)
		{
		  d_sink = sqrt((P[h].Pos[0] - x_sink)*(P[h].Pos[0] - x_sink) + (P[h].Pos[1] - y_sink)*(P[h].Pos[1] - y_sink) + (P[h].Pos[2] - z_sink)*(P[h].Pos[2] - z_sink));

                  d_sink_phys = d_sink * 3.086e18*1.e3*All.Time/(hubble_param);

                  vx = P[h].Vel[0] - vx_sink;
                  vy = P[h].Vel[1] - vy_sink;
                  vz = P[h].Vel[2] - vz_sink;

                  vx_phys = vx*pow(All.Time, 0.5)*1.e5;
                  vy_phys = vy*pow(All.Time, 0.5)*1.e5;
                  vz_phys = vz*pow(All.Time, 0.5)*1.e5;

                  jx = (P[h].Pos[1]-y_sink)*vz  - (P[h].Pos[2]-z_sink)*vy;
                  jy = (P[h].Pos[2]-z_sink)*vx  - (P[h].Pos[0]-x_sink)*vz;
                  jz = (P[h].Pos[0]-x_sink)*vy  - (P[h].Pos[1]-y_sink)*vx;

                  jx =  jx*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);  //convert to cgs units
                  jy =  jy*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);
                  jz =  jz*pow(All.Time, 0.5)*1.e5*3.086e18*1.e3*All.Time/(hubble_param);

                  jtot = pow(jx*jx + jy*jy + jz*jz, 0.5);
                  jcent_sink = pow((6.67e-8*sink_Mass*(1.e10/hubble_param)*1.98892e33*d_sink_phys),0.5);

	          u_part=(SphP[h].Entropy*pow((SphP[h].Density/a3),(SphP[h].Gamma-1.0))/(SphP[h].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs);
   
	          Temp_part=(SphP[h].Gamma-1.0)/BOLTZMANN * u_part * PROTONMASS * 1.22;

	          pot_engy =  GRAVITY * mass_sink * All.UnitMass_in_g * All.HubbleParam / (d_sink * All.UnitLength_in_cm * All.Time);
	          kin_engy =  0.5e0 * All.UnitVelocity_in_cm_per_s * All.UnitVelocity_in_cm_per_s * ((P[h].Vel[0] - vx_sink)*(P[h].Vel[0] - vx_sink) + (P[h].Vel[1] - vy_sink)*(P[h].Vel[1] - vy_sink) + (P[h].Vel[2] - vz_sink)*(P[h].Vel[2] - vz_sink));

	          if(d_sink < r_sink /*&& P[h].Ti_endstep == All.Ti_Current*/ && jtot < jcent_sink /*&& SphP[h].sink < 0.5*/ && SphP[h].sink != -1)
		    {
	
                      n_acc=n_acc+1;
		      Temp_tot=Temp_tot+ (SphP[h].Entropy/(SphP[h].Gamma-1.0))*pow(SphP[h].Density*a3inv,(SphP[h].Gamma-1.0));

                      printf("before file print %15.6g %8d %15.6g %15.6g %15.11g \n", All.Time, n_acc, r_sink_phys, Temp_avg, jcent_sink);
 
                      Temp_avg=(SphP[h].Entropy/(SphP[h].Gamma-1.0))*pow(SphP[h].Density*a3inv,(SphP[h].Gamma-1.0));
                      sinkdat=fopen(fsinkdat,"a");
                      fprintf(sinkdat,"%17.13g %8d %15.6g %15.6g %15.6g %10d %10d %15.6g\n", All.Time, n_acc, r_sink_phys, (SphP[h].Entropy/(SphP[h].Gamma-1.0))*pow(SphP[h].Density*a3inv,(SphP[h].Gamma-1.0)), ent_sink, P[h].ID, sink_ID, pres_sink);
                      fclose(sinkdat);
 
                      sinkangmom=fopen(fsinkangmom,"a");
                      fprintf(sinkangmom,"%17.13g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %15.6g %10d\n", All.Time, d_sink_phys, vx_phys, vy_phys, vz_phys, jx, jy, jz, jtot, sink_Mass, sink_ID);

                      fclose(sinkangmom);

                      printf("after file print %15.6g %8d %15.6g %15.6g \n", All.Time, n_acc, r_sink_phys, Temp_avg); 

		      mass_ngb      += P[h].Mass;    
		      
		      numpos_ngb_x  += P[h].Pos[0]*P[h].Mass;
		      numpos_ngb_y  += P[h].Pos[1]*P[h].Mass;
		      numpos_ngb_z  += P[h].Pos[2]*P[h].Mass;
		      
		      numvel_ngb_x  += P[h].Vel[0]*P[h].Mass; 
		      numvel_ngb_y  += P[h].Vel[1]*P[h].Mass;
		      numvel_ngb_z  += P[h].Vel[2]*P[h].Mass;
		      
		      u=(SphP[h].Entropy/(SphP[h].Gamma-1.0))*pow((SphP[h].Density/a3),(SphP[h].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs);
		      
		      
		      numTT_ngb += ((SphP[h].Gamma-1.0)/BOLTZMANN) * u * PROTONMASS * 1.22 * P[h].Mass;
		      
		      
                     if(SphP[h].sink > 0.5e0)
			{

			  sinkmergers=fopen(fsinkmergers,"a");
			  fprintf(sinkmergers,"%17.13g %8d %8d %15.6g %15.6g \n", All.Time, P[h].ID, sink_ID, P[h].Mass, sink_Mass);
			  /*	  fflush(Fdsinkmergers);*/
			  fclose(sinkmergers);
			  
			  if(P[h].Mass > sink_Mass)
			    sink_ID_exchange = P[h].ID;
			  
			  merger_number = All.merger_number;
			  All.merger_historyTIME[merger_number] = All.Time;
			  All.merger_historyIDGONE[merger_number] = P[h].ID;
			  All.merger_historyIDGROW[merger_number] = sink_ID;
			  All.merger_historyMASSGONE[merger_number] = P[h].Mass;
			  All.merger_historyMASSGROW[merger_number] = sink_Mass;
			  All.merger_number += 1;
			  
			printf("line 753 of accretemulti.c -- merger has occured. P[h].ID = %d.  sink_ID = %d.", P[h].ID, sink_ID);  
			}



		      /* JJ -- The below is taken from Allvars.h.  We shift all of the Sph data, so that the particles that go into the sink are excluded.  */
		      
		      local_sink_ngb += 1;
		      
		      N_gas   -= 1;
		      NumPart -= 1;
		      All.TotN_gas   -= 1;
		      All.TotNumPart -= 1;
                      header.npart[0] -= 1;
                      header.npartTotal[0] -= 1;
                      All.MassTable[0] = 0.e0;		     
 
		      for(g = h; g <= NumPart; g++) 
			{		      
			  
			  P[g].Pos[0] = P[g+1].Pos[0];
			  P[g].Pos[1] = P[g+1].Pos[1];           
			  P[g].Pos[2] = P[g+1].Pos[2];           

			  P[g].Mass   = P[g+1].Mass; 
			  
			  P[g].Vel[0] = P[g+1].Vel[0];   
			  P[g].Vel[1] = P[g+1].Vel[1]; 
			  P[g].Vel[2] = P[g+1].Vel[2]; 

                          P[g].GravAccel[0] = P[g+1].GravAccel[0];
                          P[g].GravAccel[1] = P[g+1].GravAccel[1];
                          P[g].GravAccel[2] = P[g+1].GravAccel[2];			  
          
#ifdef PMGRID
			  P[g].GravPM[0] = P[g+1].GravPM[0];
			  P[g].GravPM[1] = P[g+1].GravPM[1];           
			  P[g].GravPM[2] = P[g+1].GravPM[2];    
#endif
#ifdef FORCETEST     
			  P[g].GravAccelDirect[0] = P[g+1].GravAccelDirect[0];
			  P[g].GravAccelDirect[1] = P[g+1].GravAccelDirect[1];           
			  P[g].GravAccelDirect[2] = P[g+1].GravAccelDirect[2]; 
#endif
          
                          P[g].Potential = P[g+1].Potential;
                          P[g].OldAcc = P[g+1].OldAcc;
                          P[g].ID = P[g+1].ID;
                          P[g].Type = P[g+1].Type;
                          P[g].Ti_endstep = P[g+1].Ti_endstep;
                          P[g].Ti_begstep = P[g+1].Ti_begstep;
#ifdef FLEXSTEPS 
                          P[g].FlexStepGrp = P[g+1].FlexStepGrp;
#endif
                          P[g].GravCost = P[g+1].GravCost;
#ifdef PSEUDOSYMMETRIC
                          P[g].AphysOld = P[g+1].AphysOld;
#endif 
			}
		      
		      for(g = h; g <= N_gas; g++) 
			{	       
#ifndef POLYTROPE
                          SphP[g].Entropy = SphP[g+1].Entropy;
#endif
                          SphP[g].Density = SphP[g+1].Density;
                          SphP[g].Hsml = SphP[g+1].Hsml;
                          SphP[g].Left = SphP[g+1].Left;                        
                          SphP[g].Right = SphP[g+1].Right;
                          SphP[g].NumNgb = SphP[g+1].NumNgb;	
                          SphP[g].Pressure = SphP[g+1].Pressure;
                          SphP[g].Prad = SphP[g+1].Prad;
                          for(number = 0; number < 3; number++)
                            SphP[g].Prad_dir[number] = SphP[g+1].Prad_dir[number];

#ifndef POLYTROPE
                          SphP[g].DtEntropy = SphP[g+1].DtEntropy; 
#endif
                          for(number = 0; number < 3; number++)
			    SphP[g].HydroAccel[number] = SphP[g+1].HydroAccel[number];

                          for(number = 0; number < 3; number++)
			    SphP[g].VelPred[number] = SphP[g+1].VelPred[number];

                          SphP[g].DivVel = SphP[g+1].DivVel;
                          SphP[g].CurlVel = SphP[g+1].CurlVel;

                          for(number = 0; number < 3; number++)
			    SphP[g].Rot[number] = SphP[g+1].Rot[number];

                          SphP[g].DhsmlDensityFactor = SphP[g+1].DhsmlDensityFactor;
                          SphP[g].MaxSignalVel = SphP[g+1].MaxSignalVel;
#ifdef CHEMCOOL
                          SphP[g].DtEntropyVisc = SphP[g+1].DtEntropyVisc;
                          SphP[g].Gamma = SphP[g+1].Gamma;

                          for(number = 0; number < TRAC_NUM; number++)
			    SphP[g].TracAbund[number] = SphP[g+1].TracAbund[number];

                          SphP[g].EntropyOut = SphP[g+1].EntropyOut;

                          for(number = 0; number < TRAC_NUM; number++)
	 		    SphP[g].TracAbundOut[number] = SphP[g+1].TracAbundOut[number];

                          SphP[g].DustTemp = SphP[g+1].DustTemp;
                          SphP[g].HM = SphP[g+1].HM;
                          SphP[g].H2II = SphP[g+1].H2II;
#ifdef RAYTRACE
                          for(number = 0; number < 6; number++)
	 		    SphP[g].TotalColumnDensity[number] = SphP[g+1].TotalColumnDensity[number];

                          for(number = 0; number < 6; number++)
	 		    SphP[g].H2ColumnDensity[number] = SphP[g+1].H2ColumnDensity[number];

                          for(number = 0; number < 6; number++)
	 		    SphP[g].COColumnDensity[number] = SphP[g+1].COColumnDensity[number];
#endif /* RAYTRACE */
#endif /* CHEMCOOL */
#ifdef RAYTRACE_TG
                          SphP[g].Ray_H_coeff = SphP[g+1].Ray_H_coeff;
                          SphP[g].Ray_He_coeff = SphP[g+1].Ray_He_coeff;
                          SphP[g].Ray_LW_coeff = SphP[g+1].Ray_LW_coeff;
#endif
#ifdef METALS_TG
                          SphP[g].Sigma = SphP[g+1].Sigma;
                          SphP[g].const_A = SphP[g+1].const_A;
                          SphP[g].const_B = SphP[g+1].const_B;
                          SphP[g].Metallicity = SphP[g+1].Metallicity;
#endif
                          SphP[g].sink = SphP[g+1].sink;	      

			}
		      h -= 1;
		    }
		}
	      
	    }
	  
	  
	  MPI_Allreduce(&local_sink_ngb, &tot_sink_ngb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
	  MPI_Allreduce(&mass_ngb, &tot_mass_ngb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	  MPI_Allreduce(&numpos_ngb_x, &tot_numpos_ngb_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&numpos_ngb_y, &tot_numpos_ngb_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&numpos_ngb_z, &tot_numpos_ngb_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&numvel_ngb_x, &tot_numvel_ngb_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&numvel_ngb_y, &tot_numvel_ngb_y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&numvel_ngb_z, &tot_numvel_ngb_z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Allreduce(&numTT_ngb, &tot_numTT_ngb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
          MPI_Allreduce(&sink_ID_exchange, &tot_sink_ID_exchange, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
      
	  if(my_rank == PCUsink_rank)
	    {
	      
	      mass_sink += tot_mass_ngb;
	      
	      numpos_x  += tot_numpos_ngb_x;
	      numpos_y  += tot_numpos_ngb_y;
	      numpos_z  += tot_numpos_ngb_z;
	      
	      numvel_x  += tot_numvel_ngb_x; 
	      numvel_y  += tot_numvel_ngb_y;
	      numvel_z  += tot_numvel_ngb_z;
	      
	      numTT     += tot_numTT_ngb;

	      if(tot_sink_ID_exchange > 0.5e0)
		P[num_sink_part].ID = tot_sink_ID_exchange;	      

	      P[num_sink_part].Mass = mass_sink;
	      
	      P[num_sink_part].Pos[0] =   numpos_x/mass_sink;
	      P[num_sink_part].Pos[1] =   numpos_y/mass_sink;
	      P[num_sink_part].Pos[2] =   numpos_z/mass_sink;
	      
	      P[num_sink_part].Vel[0] =   numvel_x/mass_sink;
	      P[num_sink_part].Vel[1] =   numvel_y/mass_sink;	  
	      P[num_sink_part].Vel[2] =   numvel_z/mass_sink;
	      

	      TT_sink          =   numTT/mass_sink;
	      u                =   TT_sink * BOLTZMANN/((SphP[num_sink_part].Gamma-1.0)*PROTONMASS*1.22);


            }
	}
    }
     
}


