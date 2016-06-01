#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


/*! \file predict.c
 *  \brief drift particles by a small time interval
 *
 *  This function contains code to implement a drift operation on all the
 *  particles, which represents one part of the leapfrog integration scheme.
 */


/*! This function drifts all particles from the current time to the future:
 *  time0 - > time1
 *
 *  If there is no explicit tree construction in the following timestep, the
 *  tree nodes are also drifted and updated accordingly. Note: For periodic
 *  boundary conditions, the mapping of coordinates onto the interval
 *  [0,All.BoxSize] is only done before the domain decomposition, or for
 *  outputs to snapshot files.  This simplifies dynamic tree updates, and
 *  allows the domain decomposition to be carried out only every once in a
 *  while.
 */
void move_particles(long long int time0, long long int time1, int mode)
{
  int i, j;
  double dt_drift, dt_gravkick, dt_hydrokick, dt_entr, dt_part, fac3, hubble_a;
  double t0, t1;
#ifdef CHEMCOOL
  double t2, t3;
#endif
  double a3;
  double gam_fac;

  t0 = second();

  if(All.ComovingIntegrationOn)
    {
      dt_drift = get_drift_factor(time0, time1);
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
      a3 = All.Time * All.Time * All.Time;
      hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
          + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time)
          + All.OmegaLambda;
      hubble_a = All.Hubble * All.HubbleParam * sqrt(hubble_a);
      if(ThisTask == 0)
        printf("dt_drift = %lg, dt_grav = %lg, dt_hydro = %lg\n", dt_drift, dt_gravkick, dt_hydrokick);
    }
  else
    {
      dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
      a3 = 1.0;
    }

#ifdef CHEMCOOL
  t2 = second();
#endif

  for(i = 0; i < NumPart; i++)
    {

      if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
	continue;
      for(j = 0; j < 3; j++)
	P[i].Pos[j] += P[i].Vel[j] * dt_drift;



      if(P[i].Type == 0 && SphP[i].sink >0)
        continue;

      if(P[i].Type == 0)
	{
          gam_fac = pow(All.Time,3*(5.0/3.0) -2)/pow(All.Time,3*(SphP[i].Gamma) - 2);
#ifdef PMGRID
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] +=
	      (P[i].GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick * gam_fac;
#else
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick * gam_fac;
#endif
	  SphP[i].Density *= exp(-SphP[i].DivVel * dt_drift);
	  SphP[i].Hsml *= exp(0.333333333333 * SphP[i].DivVel * dt_drift);

	  if(SphP[i].Hsml < All.MinGasHsml)
	    SphP[i].Hsml = All.MinGasHsml;


          if(SphP[i].Prad > 0)
           printf("predict acc_dir[0] = %lg, acc_dir[1] = %lg, acc_dir[2] = %lg, acc[0] = %lg, acc[1] = %lg, acc[2] = %lg\n", SphP[i].HydroAccel[0], SphP[i].HydroAccel[1], SphP[i].HydroAccel[2], P[i].GravAccel[0], P[i].GravAccel[1], P[i].GravAccel[2]); 

#ifdef POLYTROPE
          SphP[i].Pressure = get_pressure(SphP[i].Density);
#else
	  dt_entr = (double)(time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;
#ifdef CHEMCOOL
          if (All.NeedAbundancesForOutput == 1) {
           /* If dt_entr is positive, we need to evolve the chemical network forward in time
            * to get the required values. On the other hand, if dt_entr is negative (i.e. this
            * particle has already evolved past the output time), then we need to go back in
            * time. Since we can't evolve the network backwards, we instead simply use the
            * values at the output time that we computed & stored earlier.
            * 
            * Note that the timestep is constrained such that we can only have overshot one
            * output time, so the values we stored are certain to be the values for that 
            * output time.
            */
           if (dt_entr > 0 && SphP[i].sink < 0.5) { /* Output is in this particle's future */
             do_chemcool_step(dt_entr, &SphP[i], 0, All.NumCurrentTiStep, ThisTask, P[i].ID);
           }
            SphP[i].Pressure = SphP[i].EntropyOut * pow(SphP[i].Density, SphP[i].Gamma);
         }
         else {
           SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, SphP[i].Gamma);
         }
#else
	  SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
#endif /* CHEMCOOL */
#endif /* POLYTROPE */

#ifdef BFF   
#endif

	}
    }

#ifdef CHEMCOOL
  t3 = second();
#endif 

  /* if domain-decomp and tree are not going to be reconstructed, update dynamically.  */
  if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
      for(i = 0; i < Numnodestree; i++)
	for(j = 0; j < 3; j++)
	  Nodes[All.MaxPart + i].u.d.s[j] += Extnodes[All.MaxPart + i].vs[j] * dt_drift;

      force_update_len();

      force_update_pseudoparticles();
    }

  t1 = second();

#ifdef CHEMCOOL
  All.CPU_Chemcool += timediff(t2, t3);
#endif
  All.CPU_Predict += timediff(t0, t1);
}



/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  for(i = 0; i < NumPart; i++) 
    {
      if(P[i].Type == 0 && P[i].ID < 0) /*SINK*/
      continue;
    for(j = 0; j < 3; j++)
      {
	while(P[i].Pos[j] < 0)
	  P[i].Pos[j] += boxsize[j];

	while(P[i].Pos[j] >= boxsize[j])
	  P[i].Pos[j] -= boxsize[j];
      }
    }
}
#endif

#ifdef POLYTROPE
FLOAT get_pressure(FLOAT density) 
{
  FLOAT p, p1, p2, d1, d2, dd, logdens;
  int indx;

  /* Rescale from code units to CGS number density of H nuclei */
  density *= All.UnitDensity_in_cgs / ((1.0 + 4.0 *ABHE)*PROTONMASS);
  logdens = log10(density);
  if (logdens <= All.MinTabulatedDensity) {
    /* At densities below the minimum tabulated density, pressure scales as p ~ rho^All.PolyIndexLowDensity */
    p  = All.EOSPressure[0];
    p *= pow((pow(10, logdens) / pow(10, All.MinTabulatedDensity)), All.PolyIndexLowDensity);
  }
  else if (logdens >= All.MaxTabulatedDensity) {
    /* At densities above the maximum tabulated density, pressure scales as p ~ rho^All.PolyIndexHighDensity */
    p  = All.EOSPressure[All.EOSFullTableSize - 1];
    p *= pow((pow(10, logdens) / pow(10, All.MaxTabulatedDensity)), All.PolyIndexHighDensity);
  }
  else {
    indx = floor((logdens - All.MinTabulatedDensity) / All.EOSDensDel);
    d1 = All.EOSDensity[indx];
    d2 = All.EOSDensity[indx+1];
    p1 = All.EOSPressure[indx];
    p2 = All.EOSPressure[indx+1];
    dd = (logdens - d1) / (d2 - d1);
    p  = pow(10, log10(p1) + dd * (log10(p2) - log10(p1)));
  }

  /* Convert p back into code units */
  p /= All.UnitPressure_in_cgs;
  return p;
}

FLOAT get_energy(FLOAT density) 
{
  FLOAT e, e1, e2, d1, d2, dd, logdens;
  int indx;

  /* Rescale from code units to CGS number density of H nuclei */
  density *= All.UnitDensity_in_cgs / ((1.0 + 4.0 *ABHE)*PROTONMASS);
  logdens = log10(density);
  if (logdens <= All.MinTabulatedDensity) {
    /* At densities below the minimum tabulated density, energy scales as e ~ rho^All.PolyIndexLowDensity */
    e  = All.EOSEnergy[0];
    e *= pow((pow(10, logdens) / pow(10, All.MinTabulatedDensity)), All.PolyIndexLowDensity);
  }
  else if (logdens >= All.MaxTabulatedDensity) {
    /* At densities above the maximum tabulated density, energy scales as e ~ rho^All.PolyIndexHighDensity */
    e  = All.EOSEnergy[All.EOSFullTableSize - 1];
    e *= pow((pow(10, logdens) / pow(10, All.MaxTabulatedDensity)), All.PolyIndexHighDensity);
  }
  else {
    indx = floor((logdens - All.MinTabulatedDensity) / All.EOSDensDel);
    d1 = All.EOSDensity[indx];
    d2 = All.EOSDensity[indx+1];
    e1 = All.EOSEnergy[indx];
    e2 = All.EOSEnergy[indx+1];
    dd = (logdens - d1) / (d2 - d1);
    e  = pow(10, log10(e1) + dd * (log10(e2) - log10(e1)));
  }

  /* Convert e back into code units. Note that conversion factor
   * for the internal energy density is the same as that for the
   * pressure
   */
  e /= All.UnitPressure_in_cgs;
  return e;
}


#endif /* POLYTROPE */

double calc_det(double matrix00, double matrix01, double matrix02, double matrix10, double matrix11, double matrix12, double matrix20, double matrix21, double matrix22)
{
   int p, q, r;
   double def_fac = 0;
   double determinant, trA1, trA2, trA3, A2[3][3], A3[3][3], matrix[3][3];

   matrix[0][0] = matrix00 + def_fac;
   matrix[0][1] = matrix01;
   matrix[0][2] = matrix02;
   matrix[1][0] = matrix10;
   matrix[1][1] = matrix11 + def_fac;
   matrix[1][2] = matrix12;
   matrix[2][0] = matrix20;
   matrix[2][1] = matrix21;
   matrix[2][2] = matrix22 + def_fac;

   determinant =    matrix[0][0]*matrix[1][1]*matrix[2][2]
                  + matrix[0][1]*matrix[1][2]*matrix[2][0]
                  + matrix[0][2]*matrix[1][0]*matrix[2][1]
                  - matrix[0][2]*matrix[1][1]*matrix[2][0]
                  - matrix[0][1]*matrix[1][0]*matrix[2][2]
                  - matrix[0][0]*matrix[1][2]*matrix[2][1];

   return(determinant);
}


void bfield_evol(int i, double dt_part)
{
int twoD = 0;
int p, q, numsmooth=0;
double a1, a2, bmag, errDivB, dfac, vfac, vexp, dis[3], vel[3], dt_drift;
double jacob[3][3], jacob2[3][3], det_jacob, jacob_inv[3][3], jacob_fin[3][3];               
double bfield[3], bfieldy, bfieldz, dBdt[3], gsci[3], powell[3]; 

a1 = a2 = 1.0; vexp = 0;

 for(i = 0; i < N_gas; i++)
  {

  dt_part = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
  dt_part = dt_part / All.HubbleA * All.UnitTime_in_s;

  if(All.ComovingIntegrationOn)
    {
    a1 = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
    a2 = All.TimeBegin * exp(P[i].Ti_endstep * All.Timebase_interval);
    }

  if(fabs(SphP[i].jacob[0][0]) <= 0 || dt_part <= 0) continue;

  dfac = a2 / All.HubbleParam *  All.UnitLength_in_cm;
  vfac = pow(a2, -1) * All.UnitVelocity_in_cm_per_s;

  dt_drift = get_drift_factor(P[i].Ti_begstep, P[i].Ti_endstep);
  for(p=0; p<3; p++)
    {
    //dis[p] = (P[i].Vel[p] - All.v_bulk[p]) * dt_drift * dfac;
    dis[p] = SphP[i].VelRel[p] * dt_drift * dfac;

    vel[p] = dis[p] / dt_part; 

    //vel[p] = SphP[i].VelRel[p] * vfac + vexp;
    //vel[p] = (SphP[i].VelPred[p] - All.v_bulk[p]) * vfac
    }

  bmag = pow(SphP[i].bfield[0]*SphP[i].bfield[0] + SphP[i].bfield[1]*SphP[i].bfield[1] + SphP[i].bfield[2]*SphP[i].bfield[2], 0.5);
  errDivB = SphP[i].DivB * SphP[i].Hsml / bmag;


  //if(All.NumCurrentTiStep % 250 == 0 && All.NumCurrentTiStep > 2)
  //if(All.NumCurrentTiStep % 100 == 0 && All.NumCurrentTiStep > 2)
  //if(All.NumCurrentTiStep % 50 == 0 && All.NumCurrentTiStep > 2)
  //if(All.NumCurrentTiStep % 20 == 0 && All.NumCurrentTiStep > 2)

  if(All.NumCurrentTiStep > 2 && errDivB > 0.5)
    for(p=0; p<3; p++)
        {
        SphP[i].bfield[p] =SphP[i].Bsmooth[p]/SphP[i].ksum;
        //SphP[i].Sci = SphP[i].GradSci[0] = SphP[i].GradSci[1] = SphP[i].GradSci[2] = 0; 
        }


  for(p=0; p<3; p++)
    bfield[p] = SphP[i].bfield[p];

  for(p=0; p<3; p++)
    for(q=0; q<3; q++)  
        {
        jacob_inv[p][q] = jacob_fin[p][q] = jacob[p][q] = jacob2[p][q] = 0;
        jacob[p][q] = SphP[i].jacob[p][q];
        jacob2[p][q] = SphP[i].jacob2[p][q];
        }

  if(twoD == 1)
    det_jacob = jacob[0][0]*jacob[1][1] - jacob[0][1]*jacob[1][0];
  else
    det_jacob = calc_det(jacob[0][0], jacob[0][1], jacob[0][2], jacob[1][0], jacob[1][1], jacob[1][2], jacob[2][0], jacob[2][1], jacob[2][2]);

  if(twoD == 1)
    {
    jacob_inv[0][0] =  jacob[1][1];
    jacob_inv[0][1] = -jacob[0][1];
    jacob_inv[1][0] = -jacob[1][0];
    jacob_inv[1][1] =  jacob[0][0];
    }
  else
    {
    jacob_inv[0][0] = jacob[1][1]*jacob[2][2] - jacob[1][2]*jacob[2][1];
    jacob_inv[0][1] = jacob[0][2]*jacob[2][1] - jacob[0][1]*jacob[2][2];
    jacob_inv[0][2] = jacob[0][1]*jacob[1][2] - jacob[0][2]*jacob[1][1];

    jacob_inv[1][0] = jacob[1][2]*jacob[2][0] - jacob[1][0]*jacob[2][2];
    jacob_inv[1][1] = jacob[0][0]*jacob[2][2] - jacob[0][2]*jacob[2][0];
    jacob_inv[1][2] = jacob[0][2]*jacob[1][0] - jacob[0][0]*jacob[1][2];

    jacob_inv[2][0] = jacob[1][0]*jacob[2][1] - jacob[1][1]*jacob[2][0];
    jacob_inv[2][1] = jacob[0][1]*jacob[2][0] - jacob[0][0]*jacob[2][1];
    jacob_inv[2][2] = jacob[0][0]*jacob[1][1] - jacob[0][1]*jacob[1][0];
    }

  for(p=0; p<3; p++)
    for(q=0; q<3; q++)
       jacob_inv[p][q] = jacob_inv[p][q] / det_jacob;


  jacob_fin[0][0] = jacob_inv[0][0]*jacob2[0][0] + jacob_inv[0][1]*jacob2[0][1] + jacob_inv[0][2]*jacob2[0][2];
  jacob_fin[0][1] = jacob_inv[1][0]*jacob2[0][0] + jacob_inv[1][1]*jacob2[0][1] + jacob_inv[1][2]*jacob2[0][2];
  jacob_fin[0][2] = jacob_inv[2][0]*jacob2[0][0] + jacob_inv[2][1]*jacob2[0][1] + jacob_inv[2][2]*jacob2[0][2];

  jacob_fin[1][0] = jacob_inv[0][0]*jacob2[1][0] + jacob_inv[0][1]*jacob2[1][1] + jacob_inv[0][2]*jacob2[1][2];
  jacob_fin[1][1] = jacob_inv[1][0]*jacob2[1][0] + jacob_inv[1][1]*jacob2[1][1] + jacob_inv[1][2]*jacob2[1][2];
  jacob_fin[1][2] = jacob_inv[2][0]*jacob2[1][0] + jacob_inv[2][1]*jacob2[1][1] + jacob_inv[2][2]*jacob2[1][2];

  jacob_fin[2][0] = jacob_inv[0][0]*jacob2[2][0] + jacob_inv[0][1]*jacob2[2][1] + jacob_inv[0][2]*jacob2[2][2];
  jacob_fin[2][1] = jacob_inv[1][0]*jacob2[2][0] + jacob_inv[1][1]*jacob2[2][1] + jacob_inv[1][2]*jacob2[2][2];
  jacob_fin[2][2] = jacob_inv[2][0]*jacob2[2][0] + jacob_inv[2][1]*jacob2[2][1] + jacob_inv[2][2]*jacob2[2][2];

  if(twoD == 1)
    det_jacob = jacob_fin[0][0]*jacob_fin[1][1] - jacob_fin[0][1]*jacob_fin[1][0];
  else
    det_jacob = calc_det(jacob_fin[0][0], jacob_fin[0][1], jacob_fin[0][2], jacob_fin[1][0], jacob_fin[1][1], jacob_fin[1][2], jacob_fin[2][0], jacob_fin[2][1], jacob_fin[2][2]);


  SphP[i].bfield[0] = bfield[0] * jacob_fin[0][0] + bfield[1] * jacob_fin[0][1] + bfield[2] * jacob_fin[0][2];
  SphP[i].bfield[1] = bfield[1] * jacob_fin[1][1] + bfield[0] * jacob_fin[1][0] + bfield[2] * jacob_fin[1][2];
  SphP[i].bfield[2] = bfield[2] * jacob_fin[2][2] + bfield[0] * jacob_fin[2][0] + bfield[1] * jacob_fin[2][1];

  for(p=0; p<3; p++)
     SphP[i].bfield[p] = SphP[i].bfield[p] / det_jacob;

  for(p=0; p<3; p++)
     gsci[p] = SphP[i].GradSci[p] / dfac;
     //gsci[p] = 0;

  for(p=0; p<3; p++)
     powell[p] = vel[p] * (SphP[i].DivB / dfac);     

  for(p=0; p<3; p++)
    dBdt[p] = ((SphP[i].bfield[p] - bfield[p]) / dt_part) - gsci[p] - powell[p];

  for(p=0; p<3; p++)
    SphP[i].bfield[p] = bfield[p] + (dBdt[p] * dt_part); 

  for(p=0; p<3; p++)
    if(fabs(SphP[i].bfield[p]) > 1.1*fabs(bfield[p]) && fabs(bfield[p]) > 1.e-18 && All.NumCurrentTiStep > 100)
      {
      //printf("bx_old = %lg, by_old = %lg, bz_old = %lg bx = %lg by = %lg bz = %lg\n",       
      //     bfield[0], bfield[1], bfield[2], SphP[i].bfield[0], SphP[i].bfield[1],  SphP[i].bfield[2]);
    
      //SphP[i].bfield[p] = 1.1*bfield[p];
      }

  for(p=0; p<3; p++)
    SphP[i].dBdt[p] = (SphP[i].bfield[p] - bfield[p]) / (dt_part * All.HubbleA);

  for(p=0; p<3; p++)
     for(q=0; q<3; q++)
         SphP[i].jacob[p][q] = SphP[i].jacob2[p][q] = 0;


  if(i%100000 == 0)
     printf("Bx1 = %lg, By1 = %lg Bz1 = %lg, Bx2 = %lg, By2 = %lg, Bz2 = %lg,  gscix = %lg gsciy = %lg gsciz = %lg, powx = %lg, powy = %lg, powz = %lg, vx = %lg, vy = %lg, vz = %lg, dt_part = %lg\n", 
           bfield[0], bfield[1], bfield[2], SphP[i].bfield[0], SphP[i].bfield[1], SphP[i].bfield[2], 
           gsci[0]*dt_part, gsci[1]*dt_part, gsci[2]*dt_part, 
           powell[0]*dt_part, powell[1]* dt_part, powell[2]*dt_part, vel[0], vel[1], vel[2], dt_part);

  }
}


void sci_calc(int i, double dt_part)
{
double tau, c_tau, fac3, vfac, dfac, ch_local, ch2_local, csnd, sci_old, sci_new, dSci_dt; 
double sigma_ch=1.0, sigma_p=1.0, t1, t2, a1, a2;
int p;

a1 = a2 = 1.0;

 for(i = 0; i < N_gas; i++)
  {

  dt_part = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
  dt_part = dt_part / All.HubbleA * All.UnitTime_in_s;

  if(All.ComovingIntegrationOn)
    {
    a1 = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
    a2 = All.TimeBegin * exp(P[i].Ti_endstep * All.Timebase_interval);
    }

  if(fabs(SphP[i].jacob[0][0]) <= 0 || dt_part <= 0) continue;

  fac3 = pow(a2, 3 * (1 - SphP[i].Gamma) / 2.0);
  vfac = fac3 * sqrt(a2) * All.UnitVelocity_in_cm_per_s;
  dfac = a2 / All.HubbleParam *  All.UnitLength_in_cm;

  ch_local = pow(sigma_ch,0.5) * SphP[i].MaxSignalVel / 2 * vfac;  //velocity in physical km/s
  ch2_local = ch_local * ch_local; 

  //c_tau = SphP[i].MaxSignalVel / 2 * vfac;
  c_tau = All.ch;  

  csnd = sqrt(SphP[i].Gamma * SphP[i].Pressure / SphP[i].Density) * vfac;
  if(csnd > c_tau) c_tau = csnd;

  tau = sigma_p * c_tau / (SphP[i].Hsml * dfac);
  tau = 1./tau;

//////////////////////////////////////initialize Sci
 if(All.NumCurrentTiStep < 5)
    SphP[i].Sci = 0;
///////////////////////////////////////////////////////


  //t1 = - (ch2_local * (SphP[i].DivB / dfac)); 
  t1 = - (All.ch2 * (SphP[i].DivB / dfac)); 

  t2 = - (SphP[i].Sci / tau); 
  //t2 = 0;

  dSci_dt = t1 + t2;

  sci_old = SphP[i].Sci; 
  sci_new = sci_old + dSci_dt * dt_part;

  if(i%100000 == 0)  
  //if(fabs(sci_new) > 1.1*fabs(sci_old) && fabs(sci_new) > 1.e-18 && All.NumCurrentTiStep > 100)
    {
    printf("sci_old = %lg, sci_new = %lg, bx = %lg by = %lg bz = %lg t1 = %lg, t2 = %lg, tau = %lg, dt_part = %lg\n",
           sci_old, sci_new, SphP[i].bfield[0], SphP[i].bfield[1], SphP[i].bfield[2], t1, t2, tau, dt_part);

    printf("ch = %lg, c_tau = %lg, csnd = %lg\n", ch_local, c_tau, csnd);
    //sci_new = 1.1*sci_old;
    }

  SphP[i].Sci = sci_new;

  SphP[i].dScidt = (SphP[i].Sci - sci_old) / (dt_part * All.HubbleA); 

  }
}


