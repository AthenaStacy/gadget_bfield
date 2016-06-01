#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include <limits.h>

#include "allvars.h"
#include "proto.h"

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over single
 *  timesteps. The loop terminates when the cpu-time limit is reached, when a
 *  `stop' file is found in the output directory, or when the simulation ends
 *  because we arrived at TimeMax.
 */
void run(void)
{
  FILE *fd;
  int stopflag = 0;
  char stopfname[200], contfname[200];
  double t0, t1, tstart, tend, nh_local, nh_max, tot_nh_max, nh_max_nosink, tot_nh_max_nosink;
  double dt_part, ch2_local, ch2_max, tot_ch2_max, fac3, nsmooth=0, nsmooth_tot=0, bmag, errDivB; 
  double vx_tot, vy_tot, vz_tot, vx_tot_all, vy_tot_all, vz_tot_all;
  int nsinks, i, j, k; /*SINKS*/
  double res_mass, prad_avg, pres_avg, prad_tot, pres_tot, sinkmass_sum, sinkmass_sum_tot;
  double tstep2, tstep1;
  int ion, ion_tot, sink_tot_acc;
  double a3, a3inv, hubble_param, hubble_param2, Temp, SinkCriticalDensity;

  sprintf(stopfname, "%sstop", All.OutputDir);
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  do				/* main loop */
    {
      t0 = second();

      j=0;
      k=0;
      ion = 0;
      sink_tot_acc = sinkmass_sum = 0;   
      prad_avg = pres_avg = 0;
      //res_mass=.035;
      //res_mass=0.30;
      //res_mass = 6000.0;
      //res_mass = 100.0;  //res_mass for bin_zoom!
      //res_mass = 12.5;  //res_mass for bin_zoom_ref!
      //res_mass = 1.75;  //res_mass for bin_zoom_ref2 and ref3!
      res_mass = 0.9;  //res_mass for bin_zoom10_ref4
      //res_mass = 0.22; //res_mass for bin_zoom_ref4

      a3=All.Time*All.Time*All.Time;
      a3inv=1.e0/a3;

      hubble_param = All.HubbleParam;
      hubble_param2 = hubble_param*hubble_param;


      if(All.max_dens > 0)
        {
          tot_nh_max = tot_nh_max_nosink = nh_max = nh_max_nosink = nh_local = 0;
          ch2_max = tot_ch2_max = vx_tot = vy_tot = vz_tot = 0;
          SinkCriticalDensity = All.SinkCriticalDens / All.UnitDensity_in_cgs * PROTONMASS / HYDROGEN_MASSFRAC * a3 / hubble_param2;

          for(i = 0; i < N_gas; i++)
            if(P[i].ID > 0)
              {
                nh_local = SphP[i].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/All.Time/All.Time/All.Time*HYDROGEN_MASSFRAC/PROTONMASS;

                Temp =  (SphP[i].Gamma-1.0)/BOLTZMANN * (SphP[i].Entropy/(SphP[i].Gamma-1.0))*pow((SphP[i].Density/a3),(SphP[i].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs) * PROTONMASS * 1.22;

                //if(SphP[i].Density > (1.22 * PROTONMASS * 1.e3)/(a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam) && All.NumCurrentTiStep < 10)
                //  printf("highT, ID = %d, Temp = %lg, elec = %lg nh = %g  Pres = %lg, Dens = %lg, Gam = %lg\n", P[i].ID, Temp, SphP[i].TracAbund[IHP], nh_local, SphP[i].Pressure, SphP[i].Density, SphP[i].Gamma);

                if(SphP[i].Density > (1.22 * PROTONMASS * 1.e3)/(a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam) && Temp > 10000.0 && All.NumCurrentTiStep % 10 == 0)
                //if(SphP[i].Density > (1.22 * PROTONMASS * 10)/(a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam))
                  printf("highT, ID = %d, Temp = %lg, elec = %lg nh = %g  Pres = %lg, Dens = %lg, Gam = %lg\n", P[i].ID, Temp, SphP[i].TracAbund[IHP], nh_local, SphP[i].Pressure, SphP[i].Density, SphP[i].Gamma);


                if(SphP[i].sink > 0)
                  {
                  printf("sink = %lg, ID = %d, Temp = %lg, gamma = %lg\n", SphP[i].sink, P[i].ID, Temp, SphP[i].Gamma);
                  SphP[i].Density = SinkCriticalDensity;
                  SphP[i].Entropy = 2000.0*BOLTZMANN / (pow(SinkCriticalDensity*a3inv,(SphP[i].Gamma - 1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs) * PROTONMASS * 2.27);
                  SphP[i].DtEntropy = 0;
                  SphP[i].Hsml = All.SofteningGas;
                  SphP[i].Pressure = SphP[i].Entropy * pow(SphP[i].Density, SphP[i].Gamma);
                  sink_tot_acc++;
                  }

                if(P[i].Ti_endstep > P[i].Ti_begstep || SphP[i].MaxSignalVel > 0)
                  {
                  //dt_part = (double)(time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;
                  //dt_part = get_drift_factor(P[i].Ti_begstep, P[i].Ti_endstep);
                  //ch2_local = pow(10.0 * SphP[i].Hsml / dt_part, 2);

                  fac3 = pow(All.Time, 3 * (1 - SphP[i].Gamma) / 2.0);
                  ch2_local = pow(SphP[i].MaxSignalVel * fac3 * pow(All.Time,0.5) * All.UnitVelocity_in_cm_per_s, 2); 
                  ch2_local = pow(5.e5, 2);
                  }

                if(nh_local > nh_max)
                  nh_max = nh_local;

                if(ch2_local > ch2_max)
                  ch2_max = ch2_local;

                vx_tot = vx_tot + P[i].Vel[0];
                vy_tot = vy_tot + P[i].Vel[1];
                vz_tot = vz_tot + P[i].Vel[2];
              }


         if(sink_tot_acc > 0)
           nh_max = All.SinkCriticalDens;

         MPI_Allreduce(&nh_max, &tot_nh_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         MPI_Allreduce(&ch2_max, &tot_ch2_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         MPI_Allreduce(&vx_tot, &vx_tot_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(&vy_tot, &vy_tot_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(&vz_tot, &vz_tot_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         All.ch2 = tot_ch2_max;  All.ch = pow(All.ch2, 0.5);
         All.v_bulk[0] = vx_tot_all/ ((double) (All.TotN_gas)); 
         All.v_bulk[1] = vy_tot_all/ ((double) (All.TotN_gas));
         All.v_bulk[2] = vz_tot_all/ ((double) (All.TotN_gas));

         if(tot_nh_max >= 0.9*All.SinkCriticalDens)
            All.MassTable[0]=0.0;


          if(tot_nh_max > All.max_dens)
            {
              savepositions(All.SnapshotFileCount++);
              printf("We've gone past max_dens!\n");
              exit(0);
            }
        }

     if(.5*pow((res_mass/1.0e10*All.HubbleParam)/(tot_nh_max/All.UnitDensity_in_cgs/All.HubbleParam/All.HubbleParam*All.Time*All.Time*All.Time/HYDROGEN_MASSFRAC*PROTONMASS),1.0/3.0) <  0.001)
        All.SofteningGas = 0.5*pow((res_mass/1.0e10*All.HubbleParam)/(tot_nh_max/All.UnitDensity_in_cgs/All.HubbleParam/All.HubbleParam*All.Time*All.Time*All.Time/HYDROGEN_MASSFRAC*PROTONMASS),1.0/3.0);

     if(ThisTask == 0)
        {
          printf("dens_max = %g\n", tot_nh_max);
          printf("ch = %lg, ch2 = %g, vx = %lg, vy = %lg, vz = %lg\n", 
                  All.ch, All.ch2, All.v_bulk[0], All.v_bulk[1], All.v_bulk[2]);
          printf("next_time = %15.11g, ti_next = %d\n", All.TimeBegin * exp(All.Ti_nextoutput * All.Timebase_interval), All.Ti_nextoutput);
          printf("softening = %g\n", All.SofteningGas);
          printf("soft_table = %g\n", All.SofteningTable[0]);
          printf("min_soft = %g\n", All.MinGasHsml);
          fflush(stdout);
        }


      //particle_shift2();

      find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
					 * If needed, this function will also write an output file
					 * at the desired time.
					 */

      every_timestep_stuff();	/* write some info to log-files */

      domain_Decomposition();	/* do domain decomposition if needed */

      N_sinks = 0;
      for(i = 0; i < NumPart; i++)
       {
         if(P[i].Type == 5)
           N_sinks++;
       }

      particle_shift3();

      compute_accelerations(0);	/* compute accelerations for 
				 * the particles that are to be advanced  
				 */

      if(All.NumCurrentTiStep > 0)
         {
         sci_calc(0,0);
         bfield_evol(0,0);
         }

      nsmooth=0;
      for(i = 0; i < N_gas; i++)
        {
        bmag = pow(SphP[i].bfield[0]*SphP[i].bfield[0] + SphP[i].bfield[1]*SphP[i].bfield[1] + SphP[i].bfield[2]*SphP[i].bfield[2], 0.5);
        errDivB = SphP[i].DivB * SphP[i].Hsml / bmag;

        if(All.NumCurrentTiStep > 2 && errDivB > 0.1) nsmooth++;
        }

      MPI_Allreduce(&nsmooth, &nsmooth_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if(ThisTask == 0) printf("nsmooth_tot = %lg\n", nsmooth_tot);

      //particle_shift2();

      /* check whether we want a full energy statistics */
      if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
	{
#ifdef COMPUTE_POTENTIAL_ENERGY
	  compute_potential();
#endif
	  energy_statistics();	/* compute and output energy statistics */
	  All.TimeLastStatistics += All.TimeBetStatistics;
	}

/*SINKS*/
      MPI_Barrier(MPI_COMM_WORLD);

      tstart = second();

      if(ThisTask == 0)
         printf("line 144 of run.c - before sink()\n");

      if(All.NumCurrentTiStep % 10 == 0 || All.NumCurrentTiStep == 0 || All.flag_sink == 1)
       {
       sink();
       }

      for(i = 0; i < N_gas; i++)
        if(SphP[i].sink > 0.5)
             {    
             printf("sinkval = %g \n", SphP[i].sink);
             }

       if(ThisTask == 0)
         printf("line 144 of run.c - after sink()\n");

       if(All.NumCurrentTiStep % 10 == 0 || All.NumCurrentTiStep == 0 || All.flag_sink == 1)
       {
       accrete();
       }

      for(i = 0; i < N_gas; i++)
        if(SphP[i].sink > 0.5)
             {    
             printf("new new sinkval = %g \n", SphP[i].sink);
             }

      tend = second();

      All.CPU_Sinks += timediff(tstart,tend);

      MPI_Barrier(MPI_COMM_WORLD);

      if(nsinks)
        {
          All.NumForcesSinceLastDomainDecomp = All.TotNumPart * All.TreeDomainUpdateFrequency + 1;
        }
/*SINKS*/

      advance_and_find_timesteps();	/* 'kick' active particles in
					 * momentum space and compute new
					 * timesteps for them
					 */

#ifdef TURBULENCE 
      if(N_gas>0)
        {
          tstart = second();
 
          rsk_turbdriving();
 
          tend = second();
          All.CPU_Turbulence+= timediff(tstart,tend);
        }
#endif 

      All.NumCurrentTiStep++;

      /* Check whether we need to interrupt the run */
      if(ThisTask == 0)
	{
	  /* Is the stop-file present? If yes, interrupt the run. */
	  if((fd = fopen(stopfname, "r")))
	    {
	      fclose(fd);
	      stopflag = 1;
	      unlink(stopfname);
	    }

	  /* are we running out of CPU-time ? If yes, interrupt run. */
	  if(CPUThisRun > 0.85 * All.TimeLimitCPU)
	    {
	      printf("reaching time-limit. stopping.\n");
	      stopflag = 2;
	    }
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag)
	{
	  restart(0);		/* write restart file */
	  MPI_Barrier(MPI_COMM_WORLD);

	  if(stopflag == 2 && ThisTask == 0)
	    {
	      if((fd = fopen(contfname, "w")))
		fclose(fd);
	    }

	  if(stopflag == 2 && All.ResubmitOn && ThisTask == 0)
	    {
	      close_outputfiles();
	      system(All.ResubmitCommand);
	    }
	  return;
	}

      /* is it time to write a regular restart-file? (for security) */
      if(ThisTask == 0)
	{
	  if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
	    {
	      All.TimeLastRestartFile = CPUThisRun;
	      stopflag = 3;
	    }
	  else
	    stopflag = 0;
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag == 3)
	{
	  restart(0);		/* write an occasional restart file */
	  stopflag = 0;
	}

      t1 = second();

      All.CPU_Total += timediff(t0, t1);
      CPUThisRun += timediff(t0, t1);
    }
  while(All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax);

  restart(0);
/*
  savepositions(All.SnapshotFileCount++);*/	/* write a last snapshot
						 * file at final time (will
						 * be overwritten if
						 * All.TimeMax is increased
						 * and the run is continued)
						 */
}


/*! This function finds the next synchronization point of the system (i.e. the
 *  earliest point of time any of the particles needs a force computation),
 *  and drifts the system to this point of time.  If the system drifts over
 *  the desired time of a snapshot file, the function will drift to this
 *  moment, generate an output, and then resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, flag, *temp, i, nskip=10;
  long long int min_glob, min;
  double timeold;
  double t0, t1, hubble_a;
#ifdef RAYTRACE_TG
  int task_max, loc_max, tot_loc_max, list_loc_max[NTask], n_check;
  double dt_raytrace, nh_local, nh_max, tot_nh_max, mass_max, tot_mass_max, ray_dist2, list_nh_max[NTask], list_mass_max[NTask];
  double nu_min_H = 3.3e15;
  double nu_min_He = 1.32e16;
  double c_s = 16.5967; //soundspeed of 20,000K gas in km/s
#endif

  if (All.ComovingIntegrationOn) { /* comoving variables */
    hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) 
        + All.OmegaLambda;
    hubble_a = All.Hubble * All.HubbleParam * sqrt(hubble_a);
  }
  else hubble_a = 1.0;

  All.HubbleA = hubble_a;

  t0 = second();

  timeold = All.Time;

  /*SINK - must skip any accreted particles at the beginning of the SPH particle list*/
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 5 || P[i].ID >= 0)
	break;
    }


  if(i == NumPart)
    min = INT_MAX;
  else
    {
      for(n = i+1, min = P[i].Ti_endstep; n < NumPart; n++)
	{
	  if(P[n].Type == 0 && P[n].ID < 0) /*SINK*/
	    continue;

          if(min > P[n].Ti_endstep)
            min = P[n].Ti_endstep;
	}
    } 

  MPI_Allreduce(&min, &min_glob, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);

  /* We check whether this is a full step where all particles are synchronized */
  flag = 1;
  for(n = 0; n < NumPart; n++)
    {
      if(P[n].Type == 0 && P[n].ID < 0) /*SINK*/
	continue;
      if(P[n].Ti_endstep > min_glob)
        flag = 0;
    }

  MPI_Allreduce(&flag, &Flag_FullStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

#ifdef PMGRID
  if(min_glob >= All.PM_Ti_endstep)
    {
      min_glob = All.PM_Ti_endstep;
      Flag_FullStep = 1;
    }
#endif

  /* Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to be active */
  for(n = 0, NumForceUpdate = 0; n < NumPart; n++)
    {
      if(P[n].Type == 0 && P[n].ID < 0) /*SINK*/
	continue;
      if(P[n].Ti_endstep == min_glob)
#ifdef SELECTIVE_NO_GRAVITY
        if(!((1 << P[n].Type) & (SELECTIVE_NO_GRAVITY)))
#endif
          NumForceUpdate++;
    }

  /* note: NumForcesSinceLastDomainDecomp has type "long long" */
  temp = malloc(NTask * sizeof(int));
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);
  for(n = 0; n < NTask; n++)
    All.NumForcesSinceLastDomainDecomp += temp[n];
  free(temp);



  t1 = second();

  All.CPU_Predict += timediff(t0, t1);

#ifdef RAYTRACE_TG
  tot_nh_max = nh_max = nh_local = tot_mass_max = mass_max = tot_loc_max = loc_max = task_max = 0;

  for(n = 0; n < NTask; n++)
    list_nh_max[n] = list_mass_max[n] = list_loc_max[n] = 0;

  for(i = 0; i < N_gas; i++)
    if(P[i].ID > 0)
      {
        nh_local = SphP[i].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/All.Time/All.Time/All.Time*HYDROGEN_MASSFRAC/PROTONMASS;

        if(nh_local > nh_max)
        //if(P[i].ID == 52576121)
          {
            nh_max = nh_local;
            mass_max = P[i].Mass;
            loc_max = i;
          }
      }

  MPI_Allgather(&nh_max, 1, MPI_DOUBLE, &list_nh_max, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  MPI_Allgather(&mass_max, 1, MPI_DOUBLE, &list_mass_max, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  MPI_Allgather(&loc_max, 1, MPI_INT, &list_loc_max, 1, MPI_INT, MPI_COMM_WORLD);

  for(n = 0; n < NTask; n++)
    if(list_nh_max[n] > tot_nh_max)
      {
        tot_nh_max = list_nh_max[n];
        tot_mass_max = list_mass_max[n];
        tot_loc_max = list_loc_max[n];
        task_max = n;
      }

  if(ThisTask == task_max)
    ray_dist2 = (P[tot_loc_max].Pos[0] - All.BoxSize / 2.0) * (P[tot_loc_max].Pos[0] - All.BoxSize / 2.0) + (P[tot_loc_max].Pos[1] - All.BoxSize / 2.0) * (P[tot_loc_max].Pos[1] - All.BoxSize / 2.0) + (P[tot_loc_max].Pos[2] - All.BoxSize / 2.0) * (P[tot_loc_max].Pos[2] - All.BoxSize / 2.0);

  MPI_Bcast(&ray_dist2, 1, MPI_DOUBLE, task_max, MPI_COMM_WORLD);

  if(ThisTask == task_max)
    {
    printf("Densest particle (ID %d) with nh = %g at x = %g, y = %g, z = %g and mass %g\n", P[tot_loc_max].ID, tot_nh_max, P[tot_loc_max].Pos[0], P[tot_loc_max].Pos[1], P[tot_loc_max].Pos[2], P[tot_loc_max].Mass / All.HubbleParam * 1.0e10);
    for(n = 0; n < N_gas; n++)
      if(P[n].ID == P[tot_loc_max].ID)          //center ray at most dense particle
          {
           All.star_mass = P[n].Mass/All.HubbleParam/1.e-10;           //sink mass in solar masses
           printf("star_mass =%lg\n", All.star_mass);
           printf("ray_dist2 = %lg, ray_r_max_sink =%lg\n", ray_dist2, All.ray_r_max_sink);
           }

     printf("All.t_s = %lg, All.t_s0 = %lg\n", All.t_s, All.t_s0);
     All.flag_sink = 0;
     if(All.NumCurrentTiStep == 10 || All.t_s - All.t_s0 > 1.e10)
         {
         All.t_s0 = All.t_s;
         All.flag_sink = 1;
         printf("All.numtot = %d, All.t_s0 = %lg\n", All.numtot, All.t_s0);
         }
      }
      
  MPI_Bcast(&All.t_s0, 1, MPI_DOUBLE, task_max, MPI_COMM_WORLD);
  MPI_Bcast(&All.flag_sink, 1, MPI_INT, task_max, MPI_COMM_WORLD);

  if(tot_nh_max >= 0.99*All.ray_crit_dens && ray.flag_continue == 0 && ray_dist2 < All.ray_r_max_sink * All.ray_r_max_sink)
    {
      ray.flag_start = 1;

      if(ThisTask == task_max)
        All.ray_center_ID = P[tot_loc_max].ID;

      MPI_Bcast(&All.ray_center_ID, 1, MPI_INT, task_max, MPI_COMM_WORLD);

      if(ThisTask == task_max)
        printf("Found starp (ID %d) at x = %g, y = %g, z = %g with mass %g\n", P[tot_loc_max].ID, P[tot_loc_max].Pos[0], P[tot_loc_max].Pos[1], P[tot_loc_max].Pos[2], P[tot_loc_max].Mass / All.HubbleParam * 1.0e10);

    }

  // ARS asks why we cannot have a SINK be a star particle (starp)?
  if(tot_nh_max >  0.99*All.SinkCriticalDens && ray_dist2 < All.ray_r_max_sink * All.ray_r_max_sink)
    {
      if(ThisTask == task_max)
        printf("Problem! A sink instead of a starp will form (ID %d) at x = %g, y = %g, z = %g with mass %g! Aborting...\n", P[tot_loc_max].ID, P[tot_loc_max].Pos[0], P[tot_loc_max].Pos[1], P[tot_loc_max].Pos[2], P[tot_loc_max].Mass / All.HubbleParam * 1.0e10);

      //exit(0);
    }


  if(ray.flag_start == 1 && ray.flag_continue == 0 || (ray.flag_start == 1 && All.NumCurrentTiStep % nskip == 0))  //ARS asks: Where should the parentheses go?
    {
      if(ray.flag_start == 1 && ray.flag_continue == 0)
        {
          All.Time_last_raytrace = All.Time;

          dt_raytrace = fmax(All.Timebase_interval, All.MinSizeTimestep) / hubble_a * All.UnitTime_in_s;
        }
      else
        dt_raytrace = fmax(fmax((min_glob - All.Time_last_raytrace) * All.Timebase_interval, All.MinSizeTimestep), All.Timebase_interval) / hubble_a * All.UnitTime_in_s;

      All.t_s = All.t_s + dt_raytrace;   //time in seconds

      raytrace_TG(dt_raytrace);

      All.Time_last_raytrace = min_glob;
    }
#endif

//ARS NOTE!! - if savepositions function is calle OUTSIDE of "min_glob >= ..." conditional, and if not careful, the Entropy values written out to the snapshot will be WRONG!!!

  if(All.NumCurrentTiStep > 2)
        {
        All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
        domain_Decomposition();
        compute_potential();
        //savepositions(All.SnapshotFileCount++); /* write snapshot file */
        }

  while(min_glob >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
#ifdef CHEMCOOL
      All.NeedAbundancesForOutput = 1;
#endif
      move_particles(All.Ti_Current, All.Ti_nextoutput, 0);

      All.Ti_Prev = All.Ti_Current;
      All.Ti_Current = All.Ti_nextoutput;

      if(All.ComovingIntegrationOn)
	All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
      else
	All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

#ifdef OUTPUTPOTENTIAL
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
      domain_Decomposition();
      compute_potential();
#endif

      if(All.NumCurrentTiStep > 2)
        savepositions(All.SnapshotFileCount++);	/* write snapshot file */

#ifdef CHEMCOOL
      All.NeedAbundancesForOutput = 0;
#endif

      All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
#ifdef CHEMCOOL
      All.Ti_nextnextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
#endif
    }

  move_particles(All.Ti_Current, min_glob, 1);

  All.Ti_Prev = All.Ti_Current;
  All.Ti_Current = min_glob;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  All.TimeStep = All.Time - timeold;
}



/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
long long int find_next_outputtime(long long int ti_curr)
{
  int i, iter = 0, n;
  double next, time;
  long long int ti, ti_next;

  double t_dyn = 0.0;
  double res_mass = 0.0;
  double dyn_fraction = 0.0;
  double dens_dyn, next_time;
  double nh_local, nh_max, tot_nh_max;
  double hsm_local, hsm_min, tot_hsm_min;

  //res_mass = 0.035;
  //res_mass = 0.30;
  //res_mass = 6000.0;
  //res_mass = 100.0;
  //res_mass = 12.5;  //res_mass for bin_zoom_ref
  //res_mass = 1.75;  //res_mass for bin_zoom_ref2 and bin_zoom_ref3
  res_mass = 0.9;  //res_mass for bin_zoom10_ref4
  //res_mass = 0.22; //res_mass for bin_zoom_ref4

  dyn_fraction = 1.e0;
  //dyn_fraction = 1.e-1;
  //dyn_fraction = 4.e-2;
  tot_nh_max = nh_max = nh_local = 0;

  ti_next = -1;

  if(All.OutputListOn)
    {

       tot_nh_max = nh_max = nh_local = 0;
       tot_hsm_min = hsm_min = hsm_local = 1.e30;
 
       for(n = 0; n < N_gas; n++)
         if(P[n].ID > 0)
           {
             nh_local = SphP[n].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/All.Time/All.Time/All.Time*HYDROGEN_MASSFRAC/PROTONMASS;
             hsm_local = 3.085677e16*SphP[n].Hsml*All.Time/All.HubbleParam; //smoothing length in physical km             

            // if(All.NumCurrentTiStep == 0)
            //   hsm_local = pow(2.e33*res_mass*nh_local/HYDROGEN_MASSFRAC*PROTONMASS, -0.33333);
  
             if(nh_local > nh_max)
               nh_max = nh_local;

             if(hsm_local < hsm_min)
                hsm_min = hsm_local;
           }

       //nh_max = All.SinkCriticalDens;

       MPI_Allreduce(&nh_max, &tot_nh_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
       MPI_Allreduce(&hsm_min, &tot_hsm_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

       dens_dyn = tot_nh_max/HYDROGEN_MASSFRAC*PROTONMASS;
       t_dyn = 1.0/sqrt(GRAVITY*dens_dyn);
       if(ThisTask == 0)
         printf("tot_nh_max = %lg, t_dyn = %lg\n",tot_nh_max,  t_dyn);

       //ARS using Courant condition instead of t_ff condition:
       t_dyn = tot_hsm_min/3.0;   //(smoothing length divided by speed of 3 km/s)
       if(ThisTask == 0)
         printf("tot_hsm_min = %lg, t_cross = %lg\n", tot_hsm_min, t_dyn);

       next_time = pow(3.0/2.0*HUBBLE*All.HubbleParam*sqrt(All.Omega0)*dyn_fraction*t_dyn+pow(All.Time,3.0/2.0),2.0/3.0);

       if(ThisTask ==0)
         printf("next_dyn_time = %15.11g tot_nh_max %lg\n", next_time, tot_nh_max);

       for(i = 0; i < All.OutputListLength; i++)
         {
          time = All.OutputListTimes[i];

          if(next_time < time)
            time = next_time;

	  if(time >= All.TimeBegin && time <= All.TimeMax)
	    {
	      if(All.ComovingIntegrationOn)
		ti = log(time / All.TimeBegin) / All.Timebase_interval;
	      else
		ti = (time - All.TimeBegin) / All.Timebase_interval;

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    ti_next = ti;

		  if(ti_next > ti)
		    ti_next = ti;
		}
	    }
	}
    }
  else
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      else
	{
	  if(All.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
/*
      time = All.TimeOfFirstSnapshot;
*/
      time = All.TimeBegin;

      iter = 0;
/*
      while(time < All.TimeBegin)
	{
	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      printf("%e %e %e \n",time, All.TimeBegin, All.TimeBetSnapshot);
	      endrun(110);
	    }
	}
*/
      while(time <= All.TimeMax)
	{
	  if(All.ComovingIntegrationOn)
	    ti = log(time / All.TimeBegin) / All.Timebase_interval;
	  else
	    ti = (time - All.TimeBegin) / All.Timebase_interval;

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(111);
	    }
	}
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g, ti_next= %d\n\n", next, ti_next);
    }

  return ti_next;
}




/*! This routine writes one line for every timestep to two log-files.  In
 *  FdInfo, we just list the timesteps that have been done, while in FdCPU the
 *  cumulative cpu-time consumption in various parts of the code is stored.
 */
void every_timestep_stuff(void)
{
  double z;

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
	{
	  z = 1.0 / (All.Time) - 1;
	  fprintf(FdInfo, "\nBegin Step %d, Time: %15.11g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z, All.TimeStep,
		  log(All.Time) - log(All.Time - All.TimeStep));
	  printf("\nBegin Step %d, Time: %15.11g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
		 All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo, "\nBegin Step %d, Time: %15.11g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time,
		  All.TimeStep);
	  printf("\nBegin Step %d, Time: %15.11g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
	  fflush(FdInfo);
	}

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);
      fprintf(FdCPU,
             "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f ",
             All.CPU_Total, All.CPU_Gravity, All.CPU_Hydro, All.CPU_Domain, All.CPU_Potential,
             All.CPU_Predict, All.CPU_TimeLine, All.CPU_Snapshot, All.CPU_TreeWalk, All.CPU_TreeConstruction,
             All.CPU_CommSum, All.CPU_Imbalance, All.CPU_HydCompWalk, All.CPU_HydCommSumm,
             All.CPU_HydImbalance, All.CPU_EnsureNgb, All.CPU_PM, All.CPU_Peano, All.CPU_Sinks);
#ifdef TURBULENCE
      fprintf(FdCPU, "%10.2f ", All.CPU_Turbulence);
#endif
#ifdef CHEMCOOL
      fprintf(FdCPU, "%10.2f ", All.CPU_Chemcool);
#ifdef RAYTRACE
      fprintf(FdCPU, "%10.2f ", All.CPU_Raytrace);
#endif /* RAYTRACE */
#endif /* CHEMCOOL */
      fprintf(FdCPU, "\n");
      fflush(FdCPU);
    }

  set_random_numbers();
}


/*! This routine first calls a computation of various global quantities of the
 *  particle distribution, and then writes some statistics about the energies
 *  in the various particle components to the file FdEnergy.
 */
void energy_statistics(void)
{
  int i;

  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ",
	      All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);
#ifdef CHEMCOOL
      for (i=0; i < TRAC_NUM; i++) {
        fprintf(FdEnergy, "%g ", SysState.MolAbund[i]);
      }
#endif
      fprintf(FdEnergy, "\n");
      fflush(FdEnergy);
    }
}

void particle_shift2()
{  
   int i;
   for(i = 0; i < N_gas; i++)
      {
      
      SphP[i].VelPred[0] = SphP[i].VelPred[1] = SphP[i].VelPred[2] = 0;
      P[i].Vel[0] =  P[i].Vel[1] = P[i].Vel[2] = 0;
      
      }
}

void particle_shift3()
{
   int i, p, q;
   for(i = 0; i < N_gas; i++)
      {

      //SphP[i].VelPred[0] = SphP[i].VelPred[1] = SphP[i].VelPred[2] = 0;
      //P[i].Vel[0] =  P[i].Vel[1] = P[i].Vel[2] = 0;

      for(p=0; p<3; p++)
         for(q=0; q<3; q++)
            SphP[i].jacob[p][q] = SphP[i].jacob2[p][q] = 0;

      }
}

