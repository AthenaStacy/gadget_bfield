#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"

/*! \file timestep.c 
 *  \brief routines for 'kicking' particles in momentum space and assigning new timesteps
 */

static double fac1, fac2, fac3, hubble_a, atime, a3inv;
static double dt_displacement = 0;


/*! This function advances the system in momentum space, i.e. it does apply
 *  the 'kick' operation after the forces have been computed. Additionally, it
 *  assigns new timesteps to particles. At start-up, a half-timestep is
 *  carried out, as well as at the end of the simulation. In between, the
 *  half-step kick that ends the previous timestep and the half-step kick for
 *  the new timestep are combined into one operation.
 */
void advance_and_find_timesteps(void)
{
  int i, j, no;
  long long int ti_step, ti_min, tend, tstart;
  double dt_entr, dt_entr2, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2, t0, t1;
  double t2, t3;
  double aphys;
  double gam_fac;
#if !defined(CHEMCOOL) && !defined(POLYTROPE)
  double minentropy;
#endif
  FLOAT dv[3];

#ifdef FLEXSTEPS
  long long int ti_grp;
#endif
#if defined(PSEUDOSYMMETRIC) && !defined(FLEXSTEPS)
  double apred, prob;
  long long int ti_step2;
#endif
#ifdef PMGRID
  double dt_gravkickA, dt_gravkickB;
#endif
#ifdef MAKEGLASS
  double disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;
#endif
#ifdef CHEMCOOL
  double old_entropy, dt_entr_cc;
#endif

  t0 = second();

  if(All.ComovingIntegrationOn)
    {
      fac1 = 1 / (All.Time * All.Time);
#ifndef CHEMCOOL
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
#endif
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

      hubble_a = All.Hubble * sqrt(hubble_a);
      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    fac1 = fac2 = fac3 = hubble_a = a3inv = atime = 1;

#ifdef NOPMSTEPADJUSTMENT
  dt_displacement = All.MaxSizeTimestep;
#else
  if(Flag_FullStep || dt_displacement == 0)
    find_dt_displacement_constraint(hubble_a * atime * atime);
#endif

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;

  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      /* make sure that we reconstruct the domain/tree next time because we don't kick the tree nodes in this case */
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
    }
#endif


#ifdef MAKEGLASS
  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
	continue;
      for(j = 0; j < 3; j++)
	{
	  P[i].GravPM[j] *= -1;
	  P[i].GravAccel[j] *= -1;
	  P[i].GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
	}

      disp = sqrt(P[i].GravAccel[0] * P[i].GravAccel[0] +
		  P[i].GravAccel[1] * P[i].GravAccel[1] + P[i].GravAccel[2] * P[i].GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
	     dmean, globmax, sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0, dispmax = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
	continue;
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].GravAccel[j] = 0;
	}
    }
#endif

  /* Update column densities before doing chemistry update below.
   * NB. We only need to do this when using the raytracing approx.
   */
#if defined RAYTRACE && defined CHEMCOOL
  t2 = second();

  if (All.PhotochemApprox == 2) {
    raytrace();
  }

  t3 = second();
  All.CPU_Raytrace += timediff(t2, t3);
#endif /* RAYTRACE && CHEMCOOL */

  /* Now assign new timesteps and kick */

  if((All.Ti_Current % (4 * All.PresentMinStep)) == 0)
    if(All.PresentMinStep < TIMEBASE)
      All.PresentMinStep *= 2;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].ID < 0) /*SINK*/
	continue;
      //if(P[i].Type == 0 && SphP[i].sink > 0.5)  /*other SINK*/
      //  continue;
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  ti_step = get_timestep(i, &aphys, 0);
	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

	  if(ti_step < All.PresentMinStep)
	    All.PresentMinStep = ti_step;
	}
    }

  ti_step = All.PresentMinStep;
  MPI_Allreduce(&ti_step, &All.PresentMinStep, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);

  if(dt_displacement < All.MaxSizeTimestep)
    ti_step = (long long int) (dt_displacement / All.Timebase_interval);
  else
    ti_step = (long long int) (All.MaxSizeTimestep / All.Timebase_interval);

  /* make it a power 2 subdivision */
  ti_min = TIMEBASE;
  while(ti_min > ti_step)
    ti_min >>= 1;
  All.PresentMaxStep = ti_min;


  if(ThisTask == 0)
    printf("Syn Range = %g  PresentMinStep = %lld  PresentMaxStep = %lld \n",
	   (double) All.PresentMaxStep / All.PresentMinStep, All.PresentMinStep, All.PresentMaxStep);


  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
	continue;
      if(P[i].Ti_endstep == All.Ti_Current)
	{
          if(P[i].Type == 5) /*SINK*/
	    ti_step = All.PresentMinStep;
          //else if(P[i].Type == 0 && SphP[i].sink > 0.5)  /*other SINK*/
          //  ti_step = All.PresentMinStep;
	  else
  	    ti_step = get_timestep(i, &aphys, 0);

	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

#ifdef FLEXSTEPS
	  ti_grp = P[i].FlexStepGrp % All.PresentMaxStep;
	  ti_grp = (ti_grp / All.PresentMinStep) * All.PresentMinStep;
	  ti_step = ((P[i].Ti_endstep + ti_grp + ti_step) / ti_step) * ti_step - (P[i].Ti_endstep + ti_grp);
#else

#ifdef PSEUDOSYMMETRIC
	  if(P[i].Type != 0)
	    {
	      if(P[i].Ti_endstep > P[i].Ti_begstep)
		{
		  apred = aphys + ((aphys - P[i].AphysOld) / (P[i].Ti_endstep - P[i].Ti_begstep)) * ti_step;
		  if(fabs(apred - aphys) < 0.5 * aphys)
		    {
		      ti_step2 = get_timestep(i, &apred, -1);
		      ti_min = TIMEBASE;
		      while(ti_min > ti_step2)
			ti_min >>= 1;
		      ti_step2 = ti_min;

		      if(ti_step2 < ti_step)
			{
			  get_timestep(i, &apred, ti_step);
			  prob =
			    ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
									  P[i].Ti_begstep)) / ti_step;
			  if(prob < get_random_number(P[i].ID))
			    ti_step /= 2;
			}
		      else if(ti_step2 > ti_step)
			{
			  get_timestep(i, &apred, 2 * ti_step);
			  prob =
			    ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
									  P[i].Ti_begstep)) / ti_step;
			  if(prob < get_random_number(P[i].ID + 1))
			    ti_step *= 2;
			}
		    }
		}
	      P[i].AphysOld = aphys;
	    }
#endif

#ifdef SYNCHRONIZATION
	  if(ti_step > (P[i].Ti_endstep - P[i].Ti_begstep))	/* timestep wants to increase */
	    {
	      if(((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
		ti_step = P[i].Ti_endstep - P[i].Ti_begstep;	/* leave at old step */
	    }
#endif
#endif /* end of FLEXSTEPS */

	  if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	    ti_step = 0;

	  if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	    ti_step = TIMEBASE - All.Ti_Current;

	  tstart = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;	/* midpoint of old step */
	  tend = P[i].Ti_endstep + ti_step / 2;	/* midpoint of new step */

	  if(All.ComovingIntegrationOn)
	    {
	      dt_entr = (tend - tstart) * All.Timebase_interval;
	      dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
	      dt_gravkick = get_gravkick_factor(tstart, tend);
	      dt_hydrokick = get_hydrokick_factor(tstart, tend);
	      dt_gravkick2 = get_gravkick_factor(P[i].Ti_endstep, tend);
	      dt_hydrokick2 = get_hydrokick_factor(P[i].Ti_endstep, tend);
	    }
	  else
	    {
	      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
	      dt_gravkick2 = dt_hydrokick2 = dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
	    }

	  P[i].Ti_begstep = P[i].Ti_endstep;
	  P[i].Ti_endstep = P[i].Ti_begstep + ti_step;

	  /* do the kick */

	  for(j = 0; j < 3; j++)
	    {
	      dv[j] = P[i].GravAccel[j] * dt_gravkick;
	      P[i].Vel[j] += dv[j];
	    }

	  if(P[i].Type == 0 && SphP[i].sink > -4 && SphP[i].sink < 0.5)	/* SPH stuff */ //ARS added "&& SphP[i].sink < 0.5"
	    {
              gam_fac = pow(All.Time,3*(5.0/3.0) -2)/pow(All.Time,3*(SphP[i].Gamma) - 2);
	      for(j = 0; j < 3; j++)
		{
		  dv[j] += SphP[i].HydroAccel[j] * dt_hydrokick * gam_fac;
		  P[i].Vel[j] += SphP[i].HydroAccel[j] * dt_hydrokick * gam_fac;

		  SphP[i].VelPred[j] =
		    P[i].Vel[j] - dt_gravkick2 * P[i].GravAccel[j] - dt_hydrokick2 * SphP[i].HydroAccel[j] * gam_fac;
#ifdef PMGRID
		  SphP[i].VelPred[j] += P[i].GravPM[j] * dt_gravkickB;
#endif
		}

#ifdef CHEMCOOL
              t2 = second();
	      old_entropy = SphP[i].Entropy;

              if (tend > All.Ti_nextoutput) {
		/* We're evolving past an output time, so we need to do this in two parts */
		dt_entr_cc = (double)(All.Ti_nextoutput - tstart) * All.Timebase_interval;

                if(SphP[i].sink < 0.5 && SphP[i].sink > -4)
 		  do_chemcool_step(dt_entr_cc, &SphP[i], 1, All.NumCurrentTiStep, ThisTask, P[i].ID);
  
		/* Save a copy of the abundances and entropy at the output time */
                for(j = 0; j < TRAC_NUM; j++) {
                  SphP[i].TracAbundOut[j] = SphP[i].TracAbund[j];
  	        }
		SphP[i].EntropyOut = SphP[i].Entropy;

		dt_entr_cc = (double)(tend - All.Ti_nextoutput) * All.Timebase_interval;

                if(SphP[i].sink < 0.5 && SphP[i].sink > -4)
		  do_chemcool_step(dt_entr_cc, &SphP[i], 1, All.NumCurrentTiStep, ThisTask, P[i].ID);
	      }
	      else {
                if(SphP[i].sink < 0.5 && SphP[i].sink > -4)
 		  do_chemcool_step(dt_entr, &SphP[i], 1, All.NumCurrentTiStep, ThisTask, P[i].ID);
                if (ti_step == TIMEBASE - All.Ti_Current) {
		  /* Last step... */
                  for(j = 0; j < TRAC_NUM; j++) {
                    SphP[i].TracAbundOut[j] = SphP[i].TracAbund[j];
  	          }
		  SphP[i].EntropyOut = SphP[i].Entropy;
		}
	      }

 	      if (dt_entr > 0 && SphP[i].sink < 0.5) {
  	        SphP[i].DtEntropy = (SphP[i].Entropy - old_entropy) / dt_entr;
	      }
	      else {
		/* Actually, it's undefined, but this ensures that we store _some_ value
		 * here, rather than just leaving whatever junk was there before
		 */
		SphP[i].DtEntropy = 0.0;
	      }
              t3 = second();
              All.CPU_Chemcool += timediff(t2, t3);
#else /* CHEMCOOL */
#ifndef POLYTROPE
	      /* In case of cooling, we prevent that the entropy (and
	         hence temperature decreases by more than a factor 0.5 */

	      if(SphP[i].DtEntropy * dt_entr > -0.5 * SphP[i].Entropy && SphP[i].sink < 0.5)
		SphP[i].Entropy += SphP[i].DtEntropy * dt_entr;
	      else
		SphP[i].Entropy *= 0.5;

	      if(All.MinEgySpec)
		{
		  minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
		  if(SphP[i].Entropy < minentropy)
		    {
		      SphP[i].Entropy = minentropy;
		      SphP[i].DtEntropy = 0;
		    }
		}

	      /* In case the timestep increases in the new step, we
	         make sure that we do not 'overcool' when deriving
	         predicted temperatures. The maximum timespan over
	         which prediction can occur is ti_step/2, i.e. from
	         the middle to the end of the current step */

	      dt_entr = ti_step / 2 * All.Timebase_interval;
	      if(SphP[i].Entropy + SphP[i].DtEntropy * dt_entr < 0.5 * SphP[i].Entropy && SphP[i].sink < 0.5)
		SphP[i].DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
#endif /* POLYTROPE */
#endif /* CHEMCOOL */
	    }


	  /* if tree is not going to be reconstructed, kick parent nodes dynamically.
	   */
	  if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
	    {
	      no = Father[i];
	      while(no >= 0)
		{
		  for(j = 0; j < 3; j++)
		    Extnodes[no].vs[j] += dv[j] * P[i].Mass / Nodes[no].u.d.mass;

		  no = Nodes[no].u.d.father;
		}
	    }
	}
    }

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
	    ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      tstart = (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
      tend = All.PM_Ti_endstep + ti_step / 2;

      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(tstart, tend);
      else
	dt_gravkick = (tend - tstart) * All.Timebase_interval;

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      if(All.ComovingIntegrationOn)
	dt_gravkickB = -get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkickB =
	  -((All.PM_Ti_begstep + All.PM_Ti_endstep) / 2 - All.PM_Ti_begstep) * All.Timebase_interval;

      for(i = 0; i < NumPart; i++)
	{
	  if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
	    continue;
	  for(j = 0; j < 3; j++)	/* do the kick */
	    P[i].Vel[j] += P[i].GravPM[j] * dt_gravkick;

	  if(P[i].Type == 0  && SphP[i].sink > -4 && SphP[i].sink < 0.5)   //ARS added "&& SphP[i].sink < 0.5"
	    {
	      if(All.ComovingIntegrationOn)
		{
		  dt_gravkickA = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_gravkick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
		  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_hydrokick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
		}
	      else
		dt_gravkickA = dt_hydrokick =
		  (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

              gam_fac = pow(All.Time,3*(5.0/3.0) -2)/pow(All.Time,3*(SphP[i].Gamma) - 2);

	      for(j = 0; j < 3; j++)
		SphP[i].VelPred[j] = P[i].Vel[j]
		  + P[i].GravAccel[j] * dt_gravkickA
		  + SphP[i].HydroAccel[j] * dt_hydrokick * gam_fac + P[i].GravPM[j] * dt_gravkickB;
	    }
	}
    }
#endif

  t1 = second();
  All.CPU_TimeLine += timediff(t0, t1);
}




/*! This function normally (for flag==0) returns the maximum allowed timestep
 *  of a particle, expressed in terms of the integer mapping that is used to
 *  represent the total simulated timespan. The physical acceleration is
 *  returned in `aphys'. The latter is used in conjunction with the
 *  PSEUDOSYMMETRIC integration option, which also makes of the second
 *  function of get_timestep. When it is called with a finite timestep for
 *  flag, it returns the physical acceleration that would lead to this
 *  timestep, assuming timestep criterion 0.
 */
long long int get_timestep(int p,		/*!< particle index */
		 double *aphys,	/*!< acceleration (physical units) */
		 long long int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding
				   aphys */ )
{

  int printfile = 0, q;
  FILE *CoolTime;
  char fCoolTime[200];
  double ax, ay, az, ac, csnd;
  double dt = 0, dt_courant = 0, dt_ch = 0, dt_accel = 0;
  long long int ti_step;
  double entropy_here, gamma_here, a3, nh_here, Temp, dt_ff, dt_bfd[3], dt_sci;

  sprintf(fCoolTime,"%s/CoolTime_hiacc", All.OutputDir);

#ifdef CHEMCOOL
  double dt_cooling = -1, dt_entr = -1;
  long long int tend, time_after_next;
#endif

#ifdef CONDUCTION
  double dt_cond;
#endif

  a3=All.Time*All.Time*All.Time;

  if(flag == 0)
    {
      ax = fac1 * P[p].GravAccel[0];
      ay = fac1 * P[p].GravAccel[1];
      az = fac1 * P[p].GravAccel[2];

#ifdef PMGRID
      ax += fac1 * P[p].GravPM[0];
      ay += fac1 * P[p].GravPM[1];
      az += fac1 * P[p].GravPM[2];
#endif

      if(P[p].Type == 0)
	{
        if(SphP[p].sink > -4 && SphP[p].sink < 0.5)  //ARS added "&& SphP[i].sink < 0.5"
         {
#ifdef CHEMCOOL
    	    if(All.ComovingIntegrationOn)
              {
                 fac2 = 1 / pow(All.Time, 3 * SphP[p].Gamma - 2);
	      }
#endif
  	    ax += fac2 * SphP[p].HydroAccel[0];
	    ay += fac2 * SphP[p].HydroAccel[1];
	    az += fac2 * SphP[p].HydroAccel[2];
          }
	}

      ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
      *aphys = ac;
    }
  else
    ac = *aphys;

  if(ac == 0)
    ac = 1.0e-30;

  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
	  dt = flag * All.Timebase_interval;
	  dt /= hubble_a;	/* convert dloga to physical timestep  */

/* XXX correction for variable softening length is missing below */

	  ac = 2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / (dt * dt);

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
          if(P[p].Type == 0)
            ac = 2 * All.ErrTolIntAccuracy * atime * SphP[p].Hsml / (dt * dt);
#endif
	  *aphys = ac;
	  return flag;
	}
      dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(P[p].Type == 0)
	dt = dt_accel = sqrt(2 * All.ErrTolIntAccuracy * atime * SphP[p].Hsml / 2.8 / ac);
#endif
      break;
    default:
      endrun(888);
      break;
    }

  if(P[p].Type == 0 && SphP[p].sink < 0.5 && SphP[p].sink > -4)
    {
#ifdef CHEMCOOL
      csnd = sqrt(SphP[p].Gamma * SphP[p].Pressure / SphP[p].Density);
#else
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].Density);
#endif

      if(All.ComovingIntegrationOn) {
#ifdef CHEMCOOL
        fac3 = pow(All.Time, 3 * (1 - SphP[p].Gamma) / 2.0);

        //ARS speeding things up by trying to account for variable gamma values of neighbor particles and... yeah...
        nh_here = SphP[p].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/a3*HYDROGEN_MASSFRAC/PROTONMASS;
        //if(nh_here > 1.e10 && SphP[p].Gamma > 1.4) fac3 = pow(All.Time, 3 * (1 - 1.4) / 2.0);
#endif
	dt_courant = 2 * All.CourantFac * All.Time * SphP[p].Hsml / (fac3 * SphP[p].MaxSignalVel);

        if(All.NumCurrentTiStep > 0)
          dt_ch = 2 * All.CourantFac * All.Time * SphP[p].Hsml / 
                     (All.ch / pow(All.Time,0.5) / All.UnitVelocity_in_cm_per_s);

        if(dt_ch < dt_courant && All.NumCurrentTiStep > 0) dt_courant = dt_ch;
      }
      else {
	dt_courant = 2 * All.CourantFac * SphP[p].Hsml / SphP[p].MaxSignalVel;
      }

      if(dt_courant < dt)
	dt = dt_courant;
    }

  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     hubble_a=1.
   */
  dt *= hubble_a;

#ifdef BFF 
if(P[p].Type == 0)
   {  
   for(q=0; q<3; q++)
      {
      dt_bfd[q] = 1.0 * dmax(1.e-20, fabs(SphP[p].bfield[q])) / SphP[p].dBdt[q]; 
      if(fabs(dt_bfd[q]) < dt) 
         {
         //printf("nh_here = %lg, bfd = %lg, dBdt = %lg, dtB = %lg, dt = %lg\n", 
         //       nh_here, SphP[p].bfield[q], SphP[p].dBdt[q], dt_bfd[q], dt);
         //dt = fabs(dt_bfd[q]);
         }
      }

   dt_sci = 0.2 * dmax(1.e-20, fabs(SphP[p].Sci)) / SphP[p].dScidt;
   if(fabs(dt_sci) < dt) 
       {
        //printf("nh_here = %lg, sci = %lg, dScidt = %lg, dt_sci = %lg, dt = %lg\n", 
        //        nh_here, SphP[p].Sci, SphP[p].dScidt, dt_sci, dt);
       //dt = fabs(dt_sci);
       }
  }
#endif

#ifdef CHEMCOOL
  if(P[p].Type == 0 && SphP[p].sink < 0.5 && SphP[p].sink > -4)
    {
      dt_cooling = GetCoolTime(&SphP[p], p);

      if(SphP[p].sink > 0.5)
        dt_cooling = 1.e20;


      if(dt_cooling < dt)
        dt = dt_cooling;

      /* Additional timestep constraint to prevent large entropy variations in a single step */
      if (SphP[p].DtEntropy > 0.0 || dt_cooling*All.UnitTime_in_s/hubble_a < 1.1*dt_cool_min) {
        dt_entr = (fabs(ENTROPY_TIMESTEP_FACTOR * SphP[p].Entropy / SphP[p].DtEntropy));

      if(SphP[p].sink > 0.5)
        {
        dt_entr = 1.e20;
        }

        if(dt_entr < dt) {
          dt = dt_entr;
	}
      }
    }


///////////////////////////////////////////////////////////////

  nh_here = 0;
  //CoolTime=fopen(fCoolTime,"a");
  if(P[p].Type == 0)
     {
     nh_here = SphP[p].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/a3*HYDROGEN_MASSFRAC/PROTONMASS;
    
     Temp =  (SphP[p].Gamma-1.0)/BOLTZMANN * (SphP[p].Entropy/(SphP[p].Gamma-1.0))*pow((SphP[p].Density/a3),(SphP[p].Gamma-1.0))*(All.UnitPressure_in_cgs/All.UnitDensity_in_cgs) * PROTONMASS * 1.22;

     dt_ff = ENTROPY_TIMESTEP_FACTOR * pow(GRAVITY*PROTONMASS*nh_here ,-0.5);
     dt_ff /= All.UnitTime_in_s;
     dt_ff *= hubble_a;
     if((nh_here > 1.e13 && All.NumCurrentTiStep > 0 && All.NumCurrentTiStep % 1000 == 0) /*|| (All.NumCurrentTiStep > 10 && dt*All.UnitTime_in_s/hubble_a < dt_cool_min)*/)
        {
        printf("Part-ID=%d  dt=%g  dt_ff=%g  dt_cour=%g  dt_cool=%g  dt_entr=%g, acx = %g, acy = %g, acz = %g\n",
           (int) P[p].ID, dt, dt_ff, dt_courant, dt_cooling, dt_entr, SphP[p].HydroAccel[0], SphP[p].HydroAccel[1], SphP[p].HydroAccel[2]);
        printf("Part-ID=%d  dt=%g  dt_ff=%g  dt_cour=%g  dt_cool=%g  dt_entr=%g, dt_accel = %g, nh_here=%g \n", (int) P[p].ID, dt*All.UnitTime_in_s/hubble_a, 
           dt_ff*All.UnitTime_in_s/hubble_a, dt_courant*All.UnitTime_in_s/hubble_a, dt_cooling*All.UnitTime_in_s/hubble_a, dt_entr*All.UnitTime_in_s/hubble_a, dt_accel*All.UnitTime_in_s/hubble_a, nh_here);
        }
/*
   if(All.NumCurrentTiStep < 2 && nh_here > 1.e9)
        {
        CoolTime=fopen(fCoolTime,"a");
        fprintf(CoolTime, "%15.6g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.6g %15.6g %15.6g, %15.6g, %15.6g, %15.6g, %15.6g, %15.6g\n", nh_here, Temp, SphP[p].TracAbund[IH2], P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], dt_ff*All.UnitTime_in_s/hubble_a, dt_courant*All.UnitTime_in_s/hubble_a, dt_accel*All.UnitTime_in_s/hubble_a,  dt_cooling*All.UnitTime_in_s/hubble_a, dt_entr*All.UnitTime_in_s/hubble_a, SphP[p].Pressure, SphP[p].Prad, SphP[p].Prad_dir[0], SphP[p].Prad_dir[1], SphP[p].Prad_dir[2]);
        fclose(CoolTime);
        }
*/
     }
   //fclose(CoolTime);

//////////////////////////////////////////////////////////////


  /* Limit timestep to prevent us from evolving past multiple output times. 
   * We can (and often must) overshoot one, but overshooting several causes
   * more trouble than it is worth
   */
  tend = P[p].Ti_endstep + (dt / All.Timebase_interval) / 2;
  if (tend > All.Ti_nextoutput) {
    time_after_next = All.Ti_nextnextoutput;
    if (time_after_next > 0 && time_after_next <= TIMEBASE && tend >= time_after_next) {
      /* The first check here is necessary in case find_next_outputtime returns a negative
       * value (which it will do if there's no further valid output time). We use a factor
       * of 1.95 rather than 2 in the line below to ensure that we fall just before the
       * output time 
       */
      dt = 1.95 * All.Timebase_interval * (time_after_next - P[p].Ti_endstep);
    }
  }
#endif /* CHEMCOOL */

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

  if(dt >= dt_displacement)
    dt = dt_displacement;

/* 
  nh_here = SphP[p].Density*All.UnitDensity_in_cgs*All.HubbleParam*All.HubbleParam/a3*HYDROGEN_MASSFRAC/PROTONMASS;
  if(dt < 1.e-11 && nh_here > 5.e13)
       {
       printf("Part-ID=%d  dt=%g  dt_ff=%g  dt_cour=%g  dt_cool=%g  dt_entr=%g, acx = %g, acy = %g, acz = %g\n",
           (int) P[p].ID, dt, dt_ff, dt_courant, dt_cooling, dt_entr, SphP[p].HydroAccel[0], SphP[p].HydroAccel[1], SphP[p].HydroAccel[2]);
       printf("Part-ID=%d  dt=%g  dt_ff=%g  dt_cour=%g  dt_cool=%g  dt_entr=%g, nh = %g\n", (int) P[p].ID, dt*All.UnitTime_in_s/hubble_a,
           dt_ff*All.UnitTime_in_s/hubble_a, dt_courant*All.UnitTime_in_s/hubble_a, dt_cooling*All.UnitTime_in_s/hubble_a, dt_entr*All.UnitTime_in_s/hubble_a, nh_here);
       }
*/

  //if(All.NumCurrentTiStep < 5 /*&& nh_here > 1.e1*/)
  //  dt = All.MinSizeTimestep;

  if(dt < All.MinSizeTimestep || dt != dt)
    {
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");
     if(P[p].Type == 0)
        {
          printf
            ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxsignalvel=%g dt0=%g eps=%g\n",
             (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
             SphP[p].Hsml, SphP[p].MaxSignalVel,
             sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
             All.SofteningTable[P[p].Type]);
         }
      else
        {
          printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
                 P[p].Pos[2]);
        }
      fflush(stdout);

#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

      if(P[p].Type == 0)
	{
	  printf
	    ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxsignalvel=%g dt0=%g eps=%g\n",
	     (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
	     SphP[p].Hsml, SphP[p].MaxSignalVel,
	     sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
	     All.SofteningTable[P[p].Type]);
#ifdef CHEMCOOL
          printf("Part-ID=%d  dt=%g  dt_cool=%g  dt_entr=%g\n", (int) P[p].ID, dt, dt_cooling,
	         dt_entr);
#endif
	}
      else
	{
	  printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
		 P[p].Pos[2]);
	}
      fflush(stdout);
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(P[p].Type == 0 && SphP[p].Hsml < 1.e-10)
       printf("low Hsml! Hsml = %lg dt = %lg\n", SphP[p].Hsml, dt); 
#endif

  ti_step = dt / All.Timebase_interval;

 if(dt < 1.1*All.MinSizeTimestep && All.NumCurrentTiStep > 5)
    {
     printf("\n low timestep!\n");
      if(P[p].Type == 0)
        {
          printf
            ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxsignalvel=%g dt0=%g eps=%g\n",
             (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
             SphP[p].Hsml, SphP[p].MaxSignalVel,
             sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
             All.SofteningTable[P[p].Type]);
#ifdef CHEMCOOL
          printf("Part-ID=%d  dt=%g  dt_cool=%g  dt_entr=%g sink = %g\n", (int) P[p].ID, dt, dt_cooling,
                 dt_entr, SphP[p].sink);
#endif
        }
      else
        {
          printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
                 P[p].Pos[2]);
        }
      fflush(stdout);
    }

  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g tibase=%g ti_step=%lld ac=%g xyz=(%g|%g|%g) tree=(%g|%g|%g)\n\n",
	     ThisTask, (int) P[p].ID, dt, All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].GravAccel[0], P[p].GravAccel[1], P[p].GravAccel[2]);
#ifdef CHEMCOOL
      printf("Part-ID=%d  dt=%g  dt_cool=%g\n", (int) P[p].ID, dt, dt_cooling);
#endif
#ifdef PMGRID
      printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
      if(P[p].Type == 0)
	printf("hydro-frc=(%g|%g|%g)\n", SphP[p].HydroAccel[0], SphP[p].HydroAccel[1], SphP[p].HydroAccel[2]);

      fflush(stdout);
      endrun(818);
    }

  return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global
 *  timestep of the system based on the rms velocities of particles. For
 *  cosmological simulations, the criterion used is that the rms displacement
 *  should be at most a fraction MaxRMSDisplacementFac of the mean particle
 *  separation. Note that the latter is estimated using the assigned particle
 *  masses, separately for each particle type. If comoving integration is not
 *  used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, j, type, *temp;
  int count[6];
  long long count_sum[6];
  double v[6], v_sum[6], mim[6], min_mass[6];
  double dt, dmean, asmth = 0;

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(i = 0; i < NumPart; i++)
	{
	  if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
	    continue;
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(mim[P[i].Type] > P[i].Mass)
	    mim[P[i].Type] = P[i].Mass;
	  count[P[i].Type]++;
	}

      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      temp = malloc(NTask * 6 * sizeof(int));
      MPI_Allgather(count, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
      for(i = 0; i < 6; i++)
	{
	  count_sum[i] = 0;
	  for(j = 0; j < NTask; j++)
	    count_sum[i] += temp[j * 6 + i];
	}
      free(temp);

      for(type = 0; type < 6; type++)
	{
	  if(count_sum[type] > 0)
	    {
	      if(type == 0)
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
	      else
		dmean =
		  pow(min_mass[type] /
		      ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);

	      dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		asmth = All.Asmth[1];
#endif
	      if(asmth < dmean)
		dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

	      if(ThisTask == 0)
		printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		       type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);

	      if(dt < dt_displacement)
		dt_displacement = dt;
	    }
	}

      if(ThisTask == 0)
	printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}
