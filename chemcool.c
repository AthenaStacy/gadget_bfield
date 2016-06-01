#ifdef CHEMCOOL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "f2c.h"

/*initialization of chemcool*/
void chemcool_init(void)
{
  if(ThisTask==0)
    {
      printf("initialize cooling functions...\n");
      fflush(stdout);
    }

  COOLINMO();
  CHEMINMO();
  INIT_TOLERANCES();
  LOAD_H2_TABLE();

  if(ThisTask==0)
    {
      printf("initialization of cooling functions finished.\n");
      fflush(stdout);
    }

}

/* Compute new entropy and abundances at end of timestep dt. 
 * Mode = 0 ==> update TracAbundOut, EntropyOut
 * Mode = 1 ==> update TracAbund, Entropy, Gamma
 */
void do_chemcool_step(double dt, struct sph_particle_data *current_particle, int mode, int cur_ti_step, int this_task, int part_id)
{
  double entropy_init, entropy_after_visc_update, entropy_new;
  double rho, dl, timestep, divv, energy, gamma, ekn;
  double yn, a3inv, hubble_a, abh2, abe;
  double abundances[TRAC_NUM];
#ifdef RAYTRACE
  double col_tot[6], col_H2[6], col_CO[6];
#endif
  int i;

  if (All.ComovingIntegrationOn) { /* comoving variables */
    a3inv    =  1 / (All.Time*All.Time*All.Time); 
    hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) 
        + All.OmegaLambda;
    hubble_a = All.Hubble * All.HubbleParam * sqrt(hubble_a);
    rho      = current_particle->Density * a3inv;
    dl       = All.Time * current_particle->Hsml / All.HubbleParam;
    timestep = dt / hubble_a;
    divv     = 3.0*hubble_a + current_particle->DivVel / sqrt(All.Time); /* XXX: check this */
    COOLR.redshift = (1.0 / All.Time) - 1.0;
  }
  else {
    rho      = current_particle->Density;
    dl       = current_particle->Hsml;
    /* We assume that All.Time = 0 at All.InitRedshift, and that
     * all redshifts of interest are >> 1, so we can use a simple
     * approximation for the lookback times; if this isn't the case,
     * then we should probably be using comoving coordinates anyway.
     * Note that we assume that the redshift is fixed for the duration 
     * of the timestep (or in other words that dt << t_Hubble)
     */
     COOLR.redshift = pow((3. * All.Hubble * All.HubbleParam *
			    sqrt(All.Omega0) * All.Time / 2.)  + 
			    pow(1 + All.InitRedshift, -1.5), -2./3.) - 1;
    timestep = dt;
    divv     = current_particle->DivVel;
  }

  /* Set correct dust temperature in coolr common block */
  COOLR.tdust = current_particle->DustTemp;
  //COOLR.Prad = current_particle->Prad;

#ifdef RAYTRACE_TG
  COOLR.ray_H_coeff = current_particle->Ray_H_coeff;
  COOLR.ray_He_coeff = current_particle->Ray_He_coeff;
  COOLR.ray_LW_coeff = current_particle->Ray_LW_coeff;
  COOLR.ray_NH2 = current_particle->Ray_NH2;
#endif

#ifdef METALS_TG
  /* Set elemental abundances */
  COOLR.abundo  = current_particle->Metallicity * COOLR.abratio_o;
  COOLR.abundc  = current_particle->Metallicity * COOLR.abratio_c;
  COOLR.abundsi = current_particle->Metallicity * COOLR.abratio_si;
#endif

  COOLI.ipart_id = part_id;

  /* At the moment, the entropy update due to viscous forces is operator
   * split from the entropy update due to cooling & chemistry
   */

  entropy_init              = current_particle->Entropy;
  entropy_after_visc_update = entropy_init + dt * current_particle->DtEntropyVisc;
  gamma  = current_particle->Gamma;

  /* 'energy' is internal energy density [in code units] */ 
  energy = entropy_after_visc_update * pow(rho, gamma) / (gamma - 1.0);

  if (energy < rho * All.MinEgySpec) {
      energy = rho * All.MinEgySpec;
  }

  /* Convert to cgs units */
  rho      *= All.UnitDensity_in_cgs;
  timestep *= All.UnitTime_in_s;
  energy   *= All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
  dl       *= All.UnitLength_in_cm;
  divv     *= All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm;
  for (i=0; i<TRAC_NUM; i++) {
    abundances[i] = current_particle->TracAbund[i];
  }
  yn        = rho / ((1.0 + 4.0 * ABHE) * PROTONMASS);

  if (All.ComovingIntegrationOn)
    {
      energy *= All.HubbleParam*All.HubbleParam;
      yn *=  All.HubbleParam*All.HubbleParam;
    }

#ifdef RAYTRACE
  for (i=0; i<6; i++) {
    col_tot[i] = current_particle->TotalColumnDensity[i] / (All.UnitLength_in_cm * All.UnitLength_in_cm);
    col_H2[i]  = current_particle->H2ColumnDensity[i]    / (All.UnitLength_in_cm * All.UnitLength_in_cm);
#ifdef CO_SHIELDING
    col_CO[i]  =  current_particle->COColumnDensity[i]    / (All.UnitLength_in_cm * All.UnitLength_in_cm);
#endif
  }
#ifdef CO_SHIELDING
  EVOLVE_ABUNDANCES(&timestep, &dl, &yn, &divv, &energy, abundances, &cur_ti_step, &this_task, &part_id, col_tot, col_H2, col_CO);
#else
  EVOLVE_ABUNDANCES(&timestep, &dl, &yn, &divv, &energy, abundances, &cur_ti_step, &this_task, &part_id, col_tot, col_H2);
#endif /* CO_SHIELDING */
#else
  EVOLVE_ABUNDANCES(&timestep, &dl, &yn, &divv, &energy, abundances, &cur_ti_step, &this_task, &part_id);
#endif /* RAYTRACE */

#ifdef CHEMCOOL
  current_particle->HM = COOLR.HM;
  current_particle->H2II = COOLR.H2II;
#endif

  /* Compute our new value of gamma */
  abe  = compute_electron_fraction((FLOAT*)abundances);
  abh2 = abundances[IH2];
  ekn  = energy / (BOLTZMANN * (1.0 + ABHE + abe - abh2) * yn);
  COMPUTE_GAMMA(&abh2, &ekn, &gamma);
  current_particle->Gamma = gamma;

  /* Convert back to code units from cgs */
  energy *= pow(All.UnitLength_in_cm, 3) / All.UnitEnergy_in_cgs;
  rho    /= All.UnitDensity_in_cgs;

  if (All.ComovingIntegrationOn)
    energy = energy/All.HubbleParam/All.HubbleParam;

 //ARS trying to adjust for too-low heating rates within the i-front
  //#ifdef RAYTRACE_TG
    //if(COOLR.ray_H_coeff > 0 && COOLR.ray_H_coeff < 1.e-1)
  //energy = energy*pow(1.e-1/COOLR.ray_H_coeff, 2);
  //#endif

  entropy_new = (gamma - 1.0) * energy / pow(rho, gamma);

  /* Mode 0: store new abundances, new entropy in TracAbundOut[], EntropyOut.
   * Don't update actual values.
   *
   * Mode 1: store new abundances, new entropy in TracAbund[], Entropy.
   * Store updated value of Gamma.
   */

  if (mode == 0) {
    for (i=0; i<TRAC_NUM; i++) {
      current_particle->TracAbundOut[i] = abundances[i];
    }
    current_particle->EntropyOut = entropy_new;
  }
  else if (mode == 1) {
    for (i=0; i<TRAC_NUM; i++) {
      current_particle->TracAbund[i] = abundances[i];
    }
    current_particle->Gamma   = gamma;
    current_particle->Entropy = entropy_new;
  }
  else {
    printf("Unknown mode: %d!\n", mode);
    endrun(101);
  }
}

/* Calculate the cooling time for use in timestep.c */

double GetCoolTime(struct sph_particle_data *current_particle, int p)
{
  FILE *CoolTime;
  double entropy, energy, rho, divv, dl, yn, gamma, abe, abh2, ekn, temp;
  int    i, nsp, ipar[NIPAR];
  double dtcool, a3inv, t_start, hubble_a;
  double rpar[NRPAR], y[NSPEC], ydot[NSPEC];
  double dtvisc;

  entropy         = current_particle->Entropy;
  rho             = current_particle->Density;
  divv            = current_particle->DivVel;
  dl              = current_particle->Hsml;
 
  /* Multiply by appropriate scale-factors if using comoving coords */

  if (All.ComovingIntegrationOn) {
    a3inv    = 1/(All.Time*All.Time*All.Time); 
    hubble_a = All.Hubble * All.HubbleParam * sqrt(All.Omega0 / 
               (All.Time*All.Time*All.Time) + (1 - All.Omega0 - 
                All.OmegaLambda) / (All.Time*All.Time) + All.OmegaLambda);
    rho *= a3inv;
    dl  *= All.Time / All.HubbleParam;
    divv = 3.0*hubble_a + divv / sqrt(All.Time); /* XXX: check this */
  }

  //COOLR.Prad = current_particle->Prad;
#ifdef RAYTRACE_TG
  COOLR.ray_H_coeff = current_particle->Ray_H_coeff;
  COOLR.ray_He_coeff = current_particle->Ray_He_coeff;
  COOLR.ray_LW_coeff = current_particle->Ray_LW_coeff;
  COOLR.ray_NH2 = current_particle->Ray_NH2;
#endif

#ifdef METALS_TG
  /* Set elemental abundances */
  COOLR.abundo  = current_particle->Metallicity * COOLR.abratio_o;
  COOLR.abundc  = current_particle->Metallicity * COOLR.abratio_c;
  COOLR.abundsi = current_particle->Metallicity * COOLR.abratio_si;
#endif

  for (i=0; i<TRAC_NUM; i++) { 
    y[i] = current_particle->TracAbund[i];
  }

  gamma  = current_particle->Gamma;
  /* 'energy' is internal energy density [in code units] */ 
  energy = entropy * pow(rho, gamma) / (gamma - 1.0);

  if (energy < rho * All.MinEgySpec) {
      energy = rho * All.MinEgySpec;
  }

  /* Rescale to cgs for FORTRAN routines */

  energy *= All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
  rho    *= All.UnitDensity_in_cgs;
  divv   *= All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm;
  dl     *= All.UnitLength_in_cm;

  /* Set redshift */
  if (All.ComovingIntegrationOn) {
    COOLR.redshift = (1.0 / All.Time) - 1.0;
  }
  else {
    /* See note in do_chemcool_step() above */
     COOLR.redshift = pow((3. * All.Hubble * All.HubbleParam *
			    sqrt(All.Omega0) * All.Time / 2.)  + 
			    pow(1 + All.InitRedshift, -1.5), -2./3.) - 1;
  }

  yn = rho / ((1.0 + 4.0 * ABHE) * PROTONMASS); /* 1.4 if 10% He by number */

  if (All.ComovingIntegrationOn)
    {
      energy *= All.HubbleParam*All.HubbleParam;
      yn *=  All.HubbleParam*All.HubbleParam;
    }

//#ifdef METALS_TG
  abe  = compute_electron_fraction((FLOAT*)y);
  abh2 = y[IH2];
  ekn  = energy / (BOLTZMANN * (1.0 + ABHE + abe - abh2) * yn);
  gamma  = current_particle->Gamma;
  temp = (gamma - 1.0) * ekn;
#ifdef METALS_TG
  if (temp > 2e4) {
    COOLI.iflag_highT = 1;
  }
  else {
    COOLI.iflag_highT = 0;
  }
#else
  COOLI.iflag_highT = 0;
#endif

  nsp = NSPEC;
  t_start = -1.; /* Nothing in rate_eq depends on t_start, so it doesn't
                    matter what value we give it */
  rpar[0] = yn;
  rpar[1] = dl;
  rpar[2] = divv;
#ifdef RAYTRACE
  for (i=0; i<6; i++) {
    rpar[i+3] = current_particle->TotalColumnDensity[i] / (All.UnitLength_in_cm * All.UnitLength_in_cm);
    rpar[i+9] = current_particle->H2ColumnDensity[i]    / (All.UnitLength_in_cm * All.UnitLength_in_cm);
#ifdef CO_SHIELDING
    rpar[i+15] = current_particle->COColumnDensity[i]   / (All.UnitLength_in_cm * All.UnitLength_in_cm);
#endif
  }
#endif
  ipar[0] = 0;

  y[ITMP] = energy;

  RATE_EQ(&nsp, &t_start, y, ydot, rpar, ipar);

#ifdef CHEMCOOL
  current_particle->HM = COOLR.HM;
  current_particle->H2II = COOLR.H2II;
#endif

  if (ydot[ITMP] == 0.0) {
    /* Cooling time is formally infinite. Since we can't return infinity,
       however, we make do with a very big number: 10^20 seconds. */
    dtcool = 1e20;
  }
  else {
    /* We assume that the energy is non-zero */
    dtcool = DTCOOL_SCALE_FACTOR * y[ITMP] / ydot[ITMP];
  }

  dtcool = sqrt(dtcool * dtcool); /* make sure timestep is not negative */  

/*
  if(All.NumCurrentTiStep == 0 && rho > 1.e-16)
   {
   CoolTime=fopen("/nobackupp1/astacy/CoolTime","a");
   fprintf(CoolTime,"%15.11g %8d %15.6g %15.11g %15.11g\n", All.Time, P[p].ID, rho, temp, dtcool);
   fclose(CoolTime);
   } 
*/

  dtcool /= All.UnitTime_in_s; 
  /* Rescale to cosmological units if necessary */
  if (All.ComovingIntegrationOn) {
    dtcool *= hubble_a;
  }

  /* Also compute timescale on which Entropy varies due to viscous 
   * dissipation. If this is shorter than the radiative cooling 
   * timescale, then we use this as the cooling time.
   *
   * NB This is in CODE UNITS, not CGS.
   */
  dtvisc = 1e20;
  if (current_particle->DtEntropyVisc != 0.0) {
    dtvisc = fabs(ENTROPY_TIMESTEP_FACTOR * entropy / current_particle->DtEntropyVisc);
  }

  if (dtvisc < dtcool) {
    dtcool = dtvisc;
  }

//if(COOLR.ray_H_coeff > 0 && COOLR.ray_H_coeff < 1.e1 && dtcool < All.MinSizeTimestep)
//  dtcool = fmin(dtcool*pow(1.e1/COOLR.ray_H_coeff,2), All.MinSizeTimestep);

  return dtcool;
}

double initial_electron_fraction(void)
{
  double abe, abhp, abcp, absip, abop, abdp, abhep, abhcop;
  double abhepp, abmgp, absipp, abch3p;

  abe = 0.0;

  switch(All.ChemistryNetwork) {
  case 1:
    abhp   = All.InitHPlusAbund;
    abdp   = All.InitDIIAbund;
    abhep  = All.InitHeIIAbund;
    abhepp = All.InitHeIIIAbund;
    abe    = abhp + abdp + abhep + 2.0 * abhepp;
    break;
  case 2:
    abhp   = All.InitHPlusAbund;
    abdp   = All.InitDIIAbund;
    abhep  = All.InitHeIIAbund;
    abhepp = All.InitHeIIIAbund;
    abcp   = All.InitCIIAbund;
    absip  = All.InitSiIIAbund;
    abop   = All.InitOIIAbund;
    abe    = abhp + abdp + abhep + 2.0 * abhepp;
    abe    += abcp + absip + abop;
    break;
  case 3:
    abhp   = All.InitHPlusAbund;
    abdp   = All.InitDIIAbund;
    abhep  = All.InitHeIIAbund;
    abhepp = All.InitHeIIIAbund;
    abcp   = All.InitCIIAbund;
    absip  = All.InitSiIIAbund;
    abop   = All.InitOIIAbund;
    abhcop = All.InitHCOPlusAbund;
    absipp = All.InitSiIIIAbund;
    abch3p = All.InitCH3PlusAbund;
    abmgp  = All.InitMgPlusAbund; /* Set to zero currently */
    abe    = abhp + abdp + abhep + 2.0 * abhepp;
    abe    += abcp + absip + abop;
    abe    += abhcop + 2.0 * absipp + abch3p + abmgp;
    break;
 case 4:
 case 5:
    abhp  = All.InitHPlusAbund;
    abe   = abhp;
    break;
 case 7:
    abhp   = All.InitHPlusAbund;
    abhep  = All.InitHeIIAbund;
    abcp   = All.InitCIIAbund;
    abop   = All.InitOIIAbund;
    abhcop = All.InitHCOPlusAbund;
    abch3p = All.InitCH3PlusAbund;
    abe    = abhp + abhep + abcp + abop + abhcop + abch3p;
    break;
 case 8:
 case 10:
   abe = 0.0;
   break;
 default:
   break;
  }
  return abe;
}

double compute_electron_fraction(FLOAT abundances[NSPEC])
{
  double abe;

  switch(All.ChemistryNetwork) {
  case 1:
    abe = abundances[IHP] + abundances[IDP] + abundances[IHEP] 
        + 2.0 * abundances[IHEPP];
    break;
  case 2:
    abe = abundances[IHP] + abundances[IDP] + abundances[IHEP]
        + 2.0 * abundances[IHEPP] + abundances[IC] + abundances[IO]  
        + abundances[ISi];
    break;
  case 3:
    abe = abundances[IHP] + abundances[IDP] + abundances[IHEP] 
        + 2.0 * abundances[IHEPP] +  abundances[IC] + abundances[IO]  
        + abundances[ISi]  + abundances[IHCOP]
        + 2.0 * abundances[ISIPP] + abundances[ICH3P];
    break;
  case 4:
  case 5:
    abe = abundances[IHP];
    break;
  case 7:
    abe = abundances[IHP] + abundances[IHEP]
        +  abundances[IC] + abundances[IO]  + abundances[IHCOP]
        + abundances[ICH3P];
    break;
  case 8:
  case 10:
    abe = 0.0;
    break;
  default:
    abe = 0.0; 
    break;
  }
  return abe;
}

double compute_initial_gamma(void)
{
  double abh2, abe, gamma;

  /* Simple estimate of starting gamma; this should be sufficiently
   * accurate to get us going, provided that the initial H2 fraction
   * is small
   */
  abh2 = All.InitMolHydroAbund;
  abe  = initial_electron_fraction();
  gamma = (5.0 + 5.0 * ABHE - 3.0 * abh2 + 5.0 * abe) / 
          (3.0 + 3.0 * ABHE - abh2 + 3.0 * abe);
  return gamma;
}

/* Computes initial molecular weight in units of PROTONMASS */
double compute_initial_molecular_weight(void)
{
  double abe, abh2, abheI, abheII, abheIII, abhI, abhp, mu;

  /* Ignore minor corrections due to heavy elements, deuterium */
  switch(All.ChemistryNetwork) {
  case 1:
  case 2:
  case 3:
    abhp = All.InitHPlusAbund;
    abh2 = All.InitMolHydroAbund;
    abhI = 1.0 - abhp - 2.0 * abh2;

    abheIII = All.InitHeIIIAbund;
    abheII  = All.InitHeIIAbund;
    abheI   = ABHE - abheII - abheIII;

    abe = abhp + abheII + 2.0 * abheIII;
    break;
  case 7:
    abhp = All.InitHPlusAbund;
    abh2 = All.InitMolHydroAbund;
    abhI = 1.0 - abhp - 2.0 * abh2;

    abheIII = 0.0;
    abheII  = All.InitHeIIAbund;
    abheI   = ABHE - abheII;

    abe = abhp + abheII;
    break;
 case 4:
 case 5:
    abhp = All.InitHPlusAbund;
    abh2 = All.InitMolHydroAbund;
    abhI = 1.0 - abhp - 2.0 * abh2;

    abheIII = 0.0;
    abheII  = 0.0;
    abheI   = ABHE;

    abe = abhp;
    break;
  case 8:
  case 10:
    abhp = 0.0;
    abe  = 0.0;
    abh2 = All.InitMolHydroAbund;
    abhI = 1.0 - 2.0 * abh2;
    abheIII = 0.0;
    abheII  = 0.0;
    abheI   = ABHE;
    break;
  default:
    abe = abhp = abhI = abh2 = abheI = abheII = abheIII = 0.0;
    break;
  }

  mu = abhI + abhp + 2.0 * abh2 + 4.0 * (abheI + abheII + abheIII);
  mu /= abhI + abhp + abh2 + abheI + abheII + abheIII + abe;

  return mu;
}
#endif /* CHEMCOOL */
