#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and
 *  some auxiliary quantities are computed.  If the number of neighbours
 *  obtained falls outside the target range, the correct smoothing
 *  length is determined iteratively, if needed.
 */


#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif


/*! This function computes the local density for each active SPH particle,
 *  the number of neighbours in the current smoothing radius, and the
 *  divergence and curl of the velocity field.  The pressure is updated as
 *  well.  If a particle with its smoothing region is fully inside the
 *  local domain, it is not exported to the other processors. The function
 *  also detects particles that have a number of neighbours outside the
 *  allowed tolerance range. For these particles, the smoothing length is
 *  adjusted accordingly, and the density computation is executed again.
 *  Note that the smoothing length is not allowed to fall below the lower
 *  bound set by MinGasHsml.
 */
void density(void)
{
  long long ntot, ntotleft;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n, ndone, npleft, maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask, place, nexport;
  double dt_entr, tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm, timengb, sumtimengb;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  MPI_Status status;
  double a3;
#ifdef METALS_TG
  int metal_disperse;
  long long ntotsave;
  double a, hubble_param, hubble_a, dt, exp_func, old_met, old_met_tot, new_met, new_met_tot, M_metals_local, M_metals_tot;

  if(All.ComovingIntegrationOn)
    {
      a = All.Time;

      hubble_param = All.HubbleParam;

      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time) + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda;

      hubble_a = All.Hubble * sqrt(hubble_a);
    }
  else
    a = hubble_a = hubble_param = 1.0;
#endif

  if(All.ComovingIntegrationOn)
    {
       a3 = All.Time * All.Time * All.Time;
    }
  else
    a3 = 1.0;


#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(P[n].ID < 0 || SphP[n].sink > 0.5) /*SINK*/
      	continue;
      SphP[n].Left = SphP[n].Right = 0;

      if(P[n].Ti_endstep == All.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);

#ifdef METALS_TG
  ntotsave = ntot;

  M_metals_tot = M_metals_local = 0.0;

  for(metal_disperse = 0; metal_disperse < 2; metal_disperse++) {

  ntot = ntotsave;
#endif

  /* we will repeat the whole thing for those particles where we didn't
   * find enough neighbours
   */
  do
    {
      i = 0;			/* begin with this index */
      ntotleft = ntot;		/* particles left for all tasks together */
      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeDensity - NTask; i++) {
	    if(P[i].ID < 0 || SphP[i].sink > 0.5) /*SINK*/
	     continue;

	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
		ndone++;

		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;

#ifdef METALS_TG
		density_evaluate(i, 0, metal_disperse);
#else
		density_evaluate(i, 0);
#endif

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			DensDataIn[nexport].Pos[0] = P[i].Pos[0];
			DensDataIn[nexport].Pos[1] = P[i].Pos[1];
			DensDataIn[nexport].Pos[2] = P[i].Pos[2];
			DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
			DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
			DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];

                        DensDataIn[nexport].bfield[0] = SphP[i].bfield[0];
                        DensDataIn[nexport].bfield[1] = SphP[i].bfield[1];
                        DensDataIn[nexport].bfield[2] = SphP[i].bfield[2];
                        
			DensDataIn[nexport].Hsml = SphP[i].Hsml;
                        DensDataIn[nexport].Sci = SphP[i].Sci;
#ifdef METALS_TG
                        DensDataIn[nexport].Sigma = SphP[i].Sigma;
#endif
			DensDataIn[nexport].Index = i;
			DensDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }
	  }
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

	  for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	  tstart = second();

	  MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  tend = second();
	  timeimbalance += timediff(tstart, tend);


	  /* now do the particles that need to be exported */

	  for(level = 1; level < (1 << PTask); level++)
	    {
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				       recvTask, TAG_DENS_A,
				       &DensDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);


	      tstart = second();
	      for(j = 0; j < nbuffer[ThisTask]; j++)
#ifdef METALS_TG
		density_evaluate(j, 1, metal_disperse);
#else
		density_evaluate(j, 1);
#endif
	      tend = second();
	      timecomp += timediff(tstart, tend);

	      /* do a block to explicitly measure imbalance */
	      tstart = second();
	      MPI_Barrier(MPI_COMM_WORLD);
	      tend = second();
	      timeimbalance += timediff(tstart, tend);

	      /* get the result */
	      tstart = second();
	      for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	      for(ngrp = level; ngrp < (1 << PTask); ngrp++)
		{
		  maxfill = 0;
		  for(j = 0; j < NTask; j++)
		    {
		      if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			  maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		    }
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* send the results */
			  MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B,
				       &DensDataPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
			      source = j + noffset[recvTask];
			      place = DensDataIn[source].Index;
#ifdef METALS_TG
			      if(metal_disperse == 1) {
		              SphP[place].const_A += DensDataPartialResult[source].const_A;
		              SphP[place].const_B += DensDataPartialResult[source].const_B;
			      }
			      else {
#endif
			      SphP[place].NumNgb += DensDataPartialResult[source].Ngb;
			      SphP[place].Density += DensDataPartialResult[source].Rho;
			      SphP[place].DivVel += DensDataPartialResult[source].Div;
                              SphP[place].DivB += DensDataPartialResult[source].Div_B;

			      SphP[place].DhsmlDensityFactor += DensDataPartialResult[source].DhsmlDensity;

			      SphP[place].Rot[0] += DensDataPartialResult[source].Rot[0];
			      SphP[place].Rot[1] += DensDataPartialResult[source].Rot[1];
			      SphP[place].Rot[2] += DensDataPartialResult[source].Rot[2];

                              SphP[place].GradSci[0] += DensDataPartialResult[source].GSci[0];
                              SphP[place].GradSci[1] += DensDataPartialResult[source].GSci[1];
                              SphP[place].GradSci[2] += DensDataPartialResult[source].GSci[2];

                              SphP[place].VelRel[0] += DensDataPartialResult[source].VRel[0];
                              SphP[place].VelRel[1] += DensDataPartialResult[source].VRel[1];
                              SphP[place].VelRel[2] += DensDataPartialResult[source].VRel[2];
#ifdef METALS_TG
			      SphP[place].Sigma += DensDataPartialResult[source].Sigma;
			      }
#endif
			    }
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);

	      level = ngrp - 1;
	    }

	  MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
	  for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}


      /* do final operations on results */
      tstart = second();
#ifdef METALS_TG
      if(metal_disperse == 0) {
#endif
      for(i = 0, npleft = 0; i < N_gas; i++)
	{
	  if(P[i].ID < 0 || SphP[i].sink > 0.5) /*SINK*/
	    continue;
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      {
		SphP[i].DhsmlDensityFactor =
		  1.0 / (1.0 + SphP[i].Hsml * SphP[i].DhsmlDensityFactor / ((double)(NUMDIMS) * SphP[i].Density));

		SphP[i].CurlVel = sqrt(SphP[i].Rot[0] * SphP[i].Rot[0] +
				       SphP[i].Rot[1] * SphP[i].Rot[1] +
				       SphP[i].Rot[2] * SphP[i].Rot[2]) / SphP[i].Density;

		SphP[i].DivVel /= SphP[i].Density;
                SphP[i].DivB /= SphP[i].Density;
                SphP[i].GradSci[0] /= SphP[i].Density;
                SphP[i].GradSci[1] /= SphP[i].Density;
                SphP[i].GradSci[2] /= SphP[i].Density;
                SphP[i].VelRel[0] /= SphP[i].Density;
                SphP[i].VelRel[1] /= SphP[i].Density;
                SphP[i].VelRel[2] /= SphP[i].Density;
#ifdef POLYTROPE
                SphP[i].Pressure = get_pressure(SphP[i].Density);
#else /* POLYTROPE */
		dt_entr = (double)((All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2.0)) * All.Timebase_interval;
#ifdef CHEMCOOL
		SphP[i].Pressure =
		  (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, SphP[i].Gamma);
#else /* CHEMCOOL */
		SphP[i].Pressure =
		  (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
#endif /* CHEMCOOL */
#endif /* POLYTROPE */
	      }

	      /* now check whether we had enough neighbours */

	      if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation) ||
		 (SphP[i].NumNgb > (All.DesNumNgb + All.MaxNumNgbDeviation)
		  && SphP[i].Hsml > (1.01 * All.MinGasHsml)))
		{
		  /* need to redo this particle */
		  npleft++;

		  if(SphP[i].Left > 0 && SphP[i].Right > 0)
		    if((SphP[i].Right - SphP[i].Left) < 1.0e-3 * SphP[i].Left)
		      {
			/* this one should be ok */
			npleft--;
			P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
			continue;
		      }

		  if(SphP[i].NumNgb < (All.DesNumNgb - All.MaxNumNgbDeviation))
		    SphP[i].Left = dmax(SphP[i].Hsml, SphP[i].Left);
		  else
		    {
		      if(SphP[i].Right != 0)
			{
			  if(SphP[i].Hsml < SphP[i].Right)
			    SphP[i].Right = SphP[i].Hsml;
			}
		      else
			SphP[i].Right = SphP[i].Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (int) P[i].ID, SphP[i].Hsml, SphP[i].Left, SphP[i].Right,
			 (float) SphP[i].NumNgb, SphP[i].Right - SphP[i].Left, P[i].Pos[0], P[i].Pos[1],
			 P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(SphP[i].Right > 0 && SphP[i].Left > 0)
		    SphP[i].Hsml = pow(0.5 * (pow(SphP[i].Left, 3) + pow(SphP[i].Right, 3)), 1.0 / 3.0);
		  else
		    {
		      if(SphP[i].Right == 0 && SphP[i].Left == 0)
			endrun(8188);	/* can't occur */

		      if(SphP[i].Right == 0 && SphP[i].Left > 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1.0 - (SphP[i].NumNgb -
				     All.DesNumNgb) / ((double)(NUMDIMS) * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml *= 1.26;
			}

		      if(SphP[i].Right > 0 && SphP[i].Left == 0)
			{
			  if(P[i].Type == 0 && fabs(SphP[i].NumNgb - All.DesNumNgb) < 0.5 * All.DesNumNgb)
			    {
			      SphP[i].Hsml *=
				1.0 - (SphP[i].NumNgb -
				     All.DesNumNgb) / ((double)(NUMDIMS) * SphP[i].NumNgb) * SphP[i].DhsmlDensityFactor;
			    }
			  else
			    SphP[i].Hsml /= 1.26;
			}
		    }

		  if(SphP[i].Hsml < All.MinGasHsml)
		    SphP[i].Hsml = All.MinGasHsml;
		}
	      else
		P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
	    }
	}
#ifdef METALS_TG
      }
      else {
 	ntot = 0;
      }
#endif
      tend = second();
      timecomp += timediff(tstart, tend);


#ifdef METALS_TG
      if(metal_disperse == 0) {
#endif
      numlist = malloc(NTask * sizeof(int) * NTask);
      MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0, ntot = 0; i < NTask; i++)
	ntot += numlist[i];
      free(numlist);

      if(ntot > 0)
	{
	  if(iter == 0)
	    tstart_ngb = second();

	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
      else
	tend_ngb = second();
#ifdef METALS_TG
      }
#endif
    }
  while(ntot > 0);

  /* mark as active again */
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0 && P[i].ID < 0) /*ACCRETED*/
        continue;

      if(P[i].Ti_endstep < 0)
        P[i].Ti_endstep = -P[i].Ti_endstep - 1;
    }
#ifdef METALS_TG
  if(metal_disperse == 0)
    {
      for(n = 0; n < N_gas; n++)
        {
          if(P[n].ID < 0 || SphP[n].sink > 0.5) /*SINK*/
            continue;

          if(P[n].Ti_endstep == All.Ti_Current)
            {
              SphP[n].Sigma = sqrt(a*SphP[n].Sigma/SphP[n].NumNgb);
            }
        }
    }
  }

  old_met = old_met_tot = new_met = new_met_tot = 0.0;

  for(n = 0; n < N_gas; n++)
    {
      if(P[n].ID < 0 || SphP[n].sink > 0.5) /*SINK*/
        continue;

      if(P[n].Ti_endstep == All.Ti_Current)
        {
          dt = (All.Ti_Current-P[n].Ti_begstep)*All.Timebase_interval/hubble_a;

          old_met += SphP[n].Metallicity*P[n].Mass;

          exp_func = exp(-SphP[n].const_A*dt);

          SphP[n].Metallicity *= exp_func;

          SphP[n].Metallicity += SphP[n].const_B/SphP[n].const_A*(1.0-exp_func);

          new_met += SphP[n].Metallicity*P[n].Mass;
        }
    }

  MPI_Allreduce(&old_met, &old_met_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&new_met, &new_met_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(n = 0; n < N_gas; n++)
    {
      if(P[n].ID < 0 || SphP[n].sink > 0.5)  /*SINK*/
        continue;

      if(P[n].Ti_endstep == All.Ti_Current && new_met_tot > 0.0)
        {
          SphP[n].Metallicity *= old_met_tot/new_met_tot;
        }

      M_metals_local += SphP[n].Metallicity*P[n].Mass;
    }

  MPI_Allreduce(&M_metals_local, &M_metals_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("Metal mass = %g\n", M_metals_tot*All.UnitMass_in_g/SOLAR_MASS*Z_SOLAR/hubble_param);
#endif

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);


  /* collect some timing information */
  if(iter > 0)
    timengb = timediff(tstart_ngb, tend_ngb);
  else
    timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
      All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
#ifdef METALS_TG
void density_evaluate(int target, int mode, int metal_disperse)
#else
void density_evaluate(int target, int mode)
#endif
{
  int ii, j, n, startnode, numngb, numngb_inbox;
  double h, h2, fac, hinv, hinv3, hinv4;
  double rho, divv, wk, dwk, divb;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz, rotv[3], gsci[3], vrel[3], sci, dBx, dBy, dBz;
  double weighted_numngb, dhsmlrho;
  FLOAT *pos, *vel, *bfd;
#ifdef METALS_TG
  double a, hubble_param, dens, D_i, D_j, D_avg, dwk_r, K_ab;
  double sigma = 0.0;
  double const_A = 0.0;
  double const_B = 0.0;

  if(All.ComovingIntegrationOn)
    {
      a = All.Time;

      hubble_param = All.HubbleParam;
    }
  else
    a = hubble_param = 1.0;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      bfd = SphP[target].bfield;
      sci = SphP[target].Sci;
      h = SphP[target].Hsml;
#ifdef METALS_TG
      D_i = 2.0*h*SphP[target].Sigma;

      if(D_i < 1.0e-20)
        D_i = 1.0e-20;
#endif
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      bfd = DensDataGet[target].bfield;
      sci = DensDataGet[target].Sci;
      h = DensDataGet[target].Hsml;
#ifdef METALS_TG
      D_i = 2.0*h*DensDataGet[target].Sigma;

      if(D_i < 1.0e-20)
        D_i = 1.0e-20;
#endif
    }

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  rho = divv = rotv[0] = rotv[1] = rotv[2] = vrel[0] = vrel[1] = vrel[2] = 0;
  gsci[0] = gsci[1] = gsci[2] = sci = divb = 0;
  weighted_numngb = 0;
  dhsmlrho = 0;

  startnode = All.MaxPart;
  numngb = 0;
  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];


	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;
          r = sqrt(r2);

	  if(r2 < h2)
	    {
	      numngb++;

	      u = r * hinv;


	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		  dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		}

	      mass_j = P[j].Mass;

              if(P[j].ID < 0 || SphP[j].sink > 0.5) /*SINK*/
               {
               printf("Skipping sink! sink = %g, ID = %d\n", SphP[j].sink, P[j].ID);
               mass_j = 0;
               }


#ifdef METALS_TG
              if(metal_disperse == 1) {
              D_j = 2.0*SphP[j].Hsml*SphP[j].Sigma;

              if(D_j < 1.0e-20)
                D_j = 1.0e-20;

              D_avg = 4.0*D_i*D_j/(D_i+D_j);

              if(r == 0.0)
                {
                  dwk_r = (double) fabs(hinv4/h*KERNEL_COEFF_4);
                }
              else
                dwk_r = (double) fabs(dwk/r);

              K_ab = mass_j/SphP[j].Density*D_avg*dwk_r*hubble_param/a;

              const_A += K_ab;
              const_B += K_ab*SphP[j].Metallicity;
              }
              else {
              sigma += (P[j].Vel[0]-vel[0])*(P[j].Vel[0]-vel[0])+(P[j].Vel[1]-vel[1])*(P[j].Vel[1]-vel[1])+(P[j].Vel[2]-vel[2])*(P[j].Vel[2]-vel[2]);
#endif
	      rho += mass_j * wk;

	      weighted_numngb += NORM_COEFF * wk / hinv3;

	      dhsmlrho += -mass_j * ((double)(NUMDIMS) * hinv * wk + u * dwk);

	      if(r > 0)
		{
		  fac = mass_j * dwk / r;

		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];

                  dBx = bfd[0] - SphP[j].bfield[0];
                  dBy = bfd[1] - SphP[j].bfield[1];
                  dBz = bfd[2] - SphP[j].bfield[2];

		  divv -= fac * (dx * dvx + dy * dvy + dz * dvz);
                  divb -= fac * (dx * dBx + dy * dBy + dz * dBz);

                  gsci[0] -= fac * dx * (sci - SphP[j].Sci);
                  gsci[1] -= fac * dy * (sci - SphP[j].Sci); 
                  gsci[2] -= fac * dz * (sci - SphP[j].Sci);   

		  rotv[0] += fac * (dz * dvy - dy * dvz);
		  rotv[1] += fac * (dx * dvz - dz * dvx);
		  rotv[2] += fac * (dy * dvx - dx * dvy);

                  vrel[0] += mass_j * wk * dvx; 
                  vrel[1] += mass_j * wk * dvy;
                  vrel[2] += mass_j * wk * dvz; 
		}
#ifdef METALS_TG
	        }
#endif
	    }
	}
    }
  while(startnode >= 0);

#ifdef METALS_TG
  if(metal_disperse == 0) {
#endif
  if(mode == 0)
    {
      SphP[target].NumNgb = weighted_numngb;
      SphP[target].Density = rho;
      SphP[target].DivVel = divv;
      SphP[target].DivB = divb;
      SphP[target].DhsmlDensityFactor = dhsmlrho;
      SphP[target].Rot[0] = rotv[0];
      SphP[target].Rot[1] = rotv[1];
      SphP[target].Rot[2] = rotv[2];
      SphP[target].GradSci[0] = gsci[0];
      SphP[target].GradSci[1] = gsci[1];
      SphP[target].GradSci[2] = gsci[2];
      SphP[target].VelRel[0] = vrel[0];
      SphP[target].VelRel[1] = vrel[1];
      SphP[target].VelRel[2] = vrel[2];
#ifdef METALS_TG
      SphP[target].Sigma = sigma;
#endif
    }
  else
    {
      DensDataResult[target].Rho = rho;
      DensDataResult[target].Div = divv;
      DensDataResult[target].Div_B = divb;
      DensDataResult[target].Ngb = weighted_numngb;
      DensDataResult[target].DhsmlDensity = dhsmlrho;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];
      DensDataResult[target].GSci[0] = gsci[0];
      DensDataResult[target].GSci[1] = gsci[1];
      DensDataResult[target].GSci[2] = gsci[2];
      DensDataResult[target].GSci[2] = gsci[2];
      DensDataResult[target].VRel[0] = vrel[0];
      DensDataResult[target].VRel[1] = vrel[1];
      DensDataResult[target].VRel[2] = vrel[2];
#ifdef METALS_TG
      DensDataResult[target].Sigma = sigma;
#endif
    }
#ifdef METALS_TG
  }
  else {
    if(mode == 0)
      {
        SphP[target].const_A = const_A;
        SphP[target].const_B = const_B;
      }
    else
      {
        DensDataResult[target].const_A = const_A;
        DensDataResult[target].const_B = const_B;
      }
  }
#endif
}




/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void *a, const void *b)
{
  if(((struct densdata_in *) a)->Task < (((struct densdata_in *) b)->Task))
    return -1;

  if(((struct densdata_in *) a)->Task > (((struct densdata_in *) b)->Task))
    return +1;

  return 0;
}
