#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file init.c
 *  \brief Code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree. Various variables of the particle data are initialised and An intial
 *  domain decomposition is performed. If SPH particles are present, the inial
 *  SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j, m, p, q;
#ifdef CHEMCOOL
  int invalid_abundance_flag;
  double max_abundance, min_abundance;
  double gamm1;
#endif
  double a3, SinkCriticalDensity,hubble_param, hubble_param2;

  All.Time = All.TimeBegin;

  hubble_param = All.HubbleParam;
  hubble_param2 = hubble_param*hubble_param;

  switch (All.ICFormat)
    {
    case 1:
#if (MAKEGLASS > 1)
      seed_glass();
#else
      read_ic(All.InitCondFile);
#endif
      break;
    case 2:
    case 3:
      read_ic(All.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;
  All.Ti_Current = 0;

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      a3 = 1;
    }

  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;

//ARS: Go here to change number of digits in snapshot number if we have over 1000 files!
  if(RestartFlag == 2)
    //All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 4) + 1;
    //All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 5) + 1;

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

#ifdef RAYTRACE_TG
  All.ray_center_ID = 0;
  All.Time_last_raytrace = 0;
  ray.flag_start = ray.flag_continue = 0;
  ray.time_start = ray.time_end = All.TimeMax;
  All.t_s = 0;
  All.t_s0 = 0;
#endif
  All.stod  = 0;         //ARS setting up protostellar evolution!

  All.alpha = 1;
  All.star_rad = 0;
  All.Tint = 0;
  All.r0 = 0; 
  All.m0 = 0; 
  All.r1 = 0; 
  All.m1 = 0; 
  All.mdot1 = 0; 
  All.r2 = 0; 
  All.m2 = 0;  
  All.e2 = 0;

  All.trans1  = 0; 
  All.trans1a = 0; 
  All.trans2  = 0; 
  All.set1    = 0;
  All.star_read = 0;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;
      P[i].OldAcc = 0;
      P[i].GravCost = 1;
      P[i].Potential = 0;
    }

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

  All.PresentMinStep = TIMEBASE;
#ifdef FLEXSTEPS
  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      P[i].FlexStepGrp = (long long int) (TIMEBASE * get_random_number(P[i].ID));
    }
#endif

/* Sanity checks for initial values */ 
#if defined CHEMCOOL && !defined METALS_TG
  invalid_abundance_flag = 0;
  if (ThisTask == 0) {
    /* H2 */
    if (All.InitMolHydroAbund > 0.5 || All.InitMolHydroAbund < 0.0) {
      invalid_abundance_flag = 1;
    }
    /* H+ */
    if (All.InitHPlusAbund > 1.0 || All.InitHPlusAbund < 0.0) {
      invalid_abundance_flag = 1;
    }
    /* HD, D+ */
    max_abundance = GADGET2_MAX(All.InitDIIAbund, All.InitHDAbund);
    min_abundance = GADGET2_MIN(All.InitDIIAbund, All.InitHDAbund);
    if (max_abundance > All.DeutAbund || min_abundance < 0.0) {
      invalid_abundance_flag = 1;
    }
    /* He+, He++ */
    max_abundance = GADGET2_MAX(All.InitHeIIAbund, All.InitHeIIIAbund);
    min_abundance = GADGET2_MIN(All.InitHeIIAbund, All.InitHeIIIAbund);
    if (max_abundance > ABHE || min_abundance < 0.0) {
      invalid_abundance_flag = 1;
    }
    /* C+, CO, HCO+, CH, CH2, CH3+ */
    max_abundance = GADGET2_MAX(All.InitCIIAbund, All.InitCOAbund);
    min_abundance = GADGET2_MIN(All.InitCIIAbund, All.InitCOAbund);
    max_abundance = GADGET2_MAX(All.InitHCOPlusAbund, max_abundance);
    min_abundance = GADGET2_MIN(All.InitHCOPlusAbund, min_abundance);
    max_abundance = GADGET2_MAX(All.InitCHAbund, max_abundance);
    min_abundance = GADGET2_MIN(All.InitCHAbund, min_abundance);
    max_abundance = GADGET2_MAX(All.InitCH2Abund, max_abundance);
    min_abundance = GADGET2_MIN(All.InitCH2Abund, min_abundance);
    max_abundance = GADGET2_MAX(All.InitCH3PlusAbund, max_abundance);
    min_abundance = GADGET2_MIN(All.InitCH3PlusAbund, min_abundance);
    if (max_abundance > All.CarbAbund || min_abundance < 0.0) {
      invalid_abundance_flag = 1;
    }
    /* O+, OH, H2O, CO */
    max_abundance = GADGET2_MAX(All.InitOIIAbund, All.InitOHAbund);
    min_abundance = GADGET2_MIN(All.InitOIIAbund, All.InitOHAbund);
    max_abundance = GADGET2_MAX(All.InitH2OAbund, max_abundance);
    min_abundance = GADGET2_MIN(All.InitH2OAbund, min_abundance);
    max_abundance = GADGET2_MAX(All.InitCOAbund, max_abundance);
    min_abundance = GADGET2_MIN(All.InitCOAbund, min_abundance);
    if (max_abundance > All.OxyAbund || min_abundance < 0.0) {
      invalid_abundance_flag = 1;
    } 
    /* C2 */
    if (All.InitC2Abund > 0.5 * All.CarbAbund || All.InitC2Abund < 0.0) {
      invalid_abundance_flag = 1;
    }

    /* O2 */
    if (All.InitO2Abund > 0.5 * All.OxyAbund || All.InitO2Abund < 0.0) {
      invalid_abundance_flag = 1;
    }
    /* Si+, Si++ */
    max_abundance = GADGET2_MAX(All.InitSiIIAbund, All.InitSiIIIAbund);
    min_abundance = GADGET2_MIN(All.InitSiIIAbund, All.InitSiIIIAbund);
    if (max_abundance > All.SiAbund || min_abundance < 0.0) {
      invalid_abundance_flag = 1;
    }
    /* Mg+ */
    if (All.InitMgPlusAbund > All.MgAbund || All.InitMgPlusAbund < 0.0) {
      invalid_abundance_flag = 1;
    }

    if (invalid_abundance_flag) {
      fprintf(stderr, "Initial abundances are invalid!\n");
      endrun(0);
    }
  }
#endif /* CHEMCOOL */

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].HydroAccel[j] = SphP[i].GradSci[j] = 0;
	}
      
      SphP[i].Sci = 0;

      for(p=0; p<3; p++)
          for(q=0; q<3; q++)
              SphP[i].jacob[p][q] = SphP[i].jacob2[p][q] = 0;

      SphP[i].dScidt = 1.e-20;
      for(p=0; p<3; p++)
        SphP[i].dBdt[p] = 1.e-20;

#ifndef POLYTROPE
      SphP[i].DtEntropy = 0;
#endif
#ifdef CHEMCOOL
      SphP[i].DtEntropyVisc = 0;
      SphP[i].EntropyOut = SphP[i].Entropy;

#endif
  
#ifdef RAYTRACE_TG
      SphP[i].Ray_H_coeff = SphP[i].Ray_He_coeff = SphP[i].Ray_LW_coeff = -1.0;
#endif
      SphP[i].Prad = 0;
      for(m = 0; m < 3; m++)
         SphP[i].Prad_dir[m] = 0;


      if(RestartFlag == 0)
	{
#ifndef CHEMCOOL
	  SphP[i].Hsml = 0;
	  SphP[i].Density = -1;
#endif
#ifdef CHEMCOOL
          SphP[i].DustTemp = All.InitDustTemp;

	  switch (All.ChemistryNetwork) {
	  case 1:
	    SphP[i].TracAbund[IH2]   = All.InitMolHydroAbund;
	    SphP[i].TracAbund[IHP]   = All.InitHPlusAbund;
	    SphP[i].TracAbund[IDP]   = All.InitDIIAbund;
	    SphP[i].TracAbund[IHD]   = All.InitHDAbund;
	    SphP[i].TracAbund[IHEP]  = All.InitHeIIAbund;
	    SphP[i].TracAbund[IHEPP] = All.InitHeIIIAbund;
	    break;
	  case 2:
	    SphP[i].TracAbund[IH2]   = All.InitMolHydroAbund;
	    SphP[i].TracAbund[IHP]   = All.InitHPlusAbund;
	    SphP[i].TracAbund[IDP]   = All.InitDIIAbund;
	    SphP[i].TracAbund[IHD]   = All.InitHDAbund;
	    SphP[i].TracAbund[IHEP]  = All.InitHeIIAbund;
	    SphP[i].TracAbund[IHEPP] = All.InitHeIIIAbund;
            SphP[i].TracAbund[IC]    = All.InitCIIAbund;
            SphP[i].TracAbund[ISi]   = All.InitSiIIAbund;
            SphP[i].TracAbund[ISIPP] = All.InitSiIIIAbund;
            SphP[i].TracAbund[IO]    = All.InitOIIAbund;
	    break;
	  case 3:
	    SphP[i].TracAbund[IH2]   = All.InitMolHydroAbund;
	    SphP[i].TracAbund[IHP]   = All.InitHPlusAbund;
	    SphP[i].TracAbund[IDP]   = All.InitDIIAbund;
	    SphP[i].TracAbund[IHD]   = All.InitHDAbund;
	    SphP[i].TracAbund[IHEP]  = All.InitHeIIAbund;
	    SphP[i].TracAbund[IHEPP] = All.InitHeIIIAbund;
            SphP[i].TracAbund[IC]    = All.InitCIIAbund;
            SphP[i].TracAbund[ISi]   = All.InitSiIIAbund;
            SphP[i].TracAbund[IO]    = All.InitOIIAbund;
	    SphP[i].TracAbund[ICO]   = All.InitCOAbund;
	    SphP[i].TracAbund[IC2]   = All.InitC2Abund;
	    SphP[i].TracAbund[IOH]   = All.InitOHAbund;
	    SphP[i].TracAbund[IH2O]  = All.InitH2OAbund;
	    SphP[i].TracAbund[IO2]   = All.InitO2Abund;
	    SphP[i].TracAbund[IHCOP] = All.InitHCOPlusAbund;
	    SphP[i].TracAbund[ICH]   = All.InitCHAbund;
	    SphP[i].TracAbund[ICH2]  = All.InitCH2Abund;
	    SphP[i].TracAbund[ISIPP] = All.InitSiIIIAbund;
	    SphP[i].TracAbund[ICH3P] = All.InitCH3PlusAbund;
	    break;
	  case 4:
	    SphP[i].TracAbund[IH2] = All.InitMolHydroAbund;
	    SphP[i].TracAbund[IHP] = All.InitHPlusAbund;
	    break;
	  case 5:
	    SphP[i].TracAbund[IH2] = All.InitMolHydroAbund;
	    SphP[i].TracAbund[IHP] = All.InitHPlusAbund;
	    SphP[i].TracAbund[ICO] = All.InitCOAbund;
	    break;
	  case 7:
	    SphP[i].TracAbund[IH2]   = All.InitMolHydroAbund;
	    SphP[i].TracAbund[IHP]   = All.InitHPlusAbund;
	    SphP[i].TracAbund[IHEP]  = All.InitHeIIAbund;
            SphP[i].TracAbund[IC]    = All.InitCIIAbund;
            SphP[i].TracAbund[IO]    = All.InitOIIAbund;
	    SphP[i].TracAbund[ICO]   = All.InitCOAbund;
	    SphP[i].TracAbund[IC2]   = All.InitC2Abund;
	    SphP[i].TracAbund[IOH]   = All.InitOHAbund;
	    SphP[i].TracAbund[IH2O]  = All.InitH2OAbund;
	    SphP[i].TracAbund[IO2]   = All.InitO2Abund;
	    SphP[i].TracAbund[IHCOP] = All.InitHCOPlusAbund;
	    SphP[i].TracAbund[ICH]   = All.InitCHAbund;
	    SphP[i].TracAbund[ICH2]  = All.InitCH2Abund;
	    SphP[i].TracAbund[ICH3P] = All.InitCH3PlusAbund;
	    break;
	  case 8:
	    SphP[i].TracAbund[IH2]   = All.InitMolHydroAbund;
	    SphP[i].TracAbund[IHD]   = All.InitHDAbund;
	    break;
          case 10:
            SphP[i].TracAbund[IH2]   = All.InitMolHydroAbund;
	    SphP[i].TracAbund[IHD]   = All.InitHDAbund;
	    SphP[i].TracAbund[IOH]   = All.InitOHAbund;
	    SphP[i].TracAbund[IH2O]  = All.InitH2OAbund;
	    SphP[i].TracAbund[IO2]   = All.InitO2Abund;
	    SphP[i].TracAbund[ICO]   = All.InitCOAbund;
	    SphP[i].TracAbund[ITD]   = All.InitDustTemp;
	  default:
	    break;
	  }

          for(j=0; j<TRAC_NUM; j++) {
	    SphP[i].TracAbundOut[j] = SphP[i].TracAbund[j];
	  }

          SphP[i].Gamma = compute_initial_gamma();
#endif /* CHEMCOOL */
#ifdef METALS_TG
          SphP[i].Metallicity = 0.0;
#endif

#ifdef SINKVAL
          SphP[i].sink = 0.0;
#endif
	}
        else {
#ifdef CHEMCOOL
	/* Starting from a restart dump, so only need to initialize 
	 * TracAbundOut and EntropyOut here (as they aren't in the 
	 * Gadget snapshot file).
	 */
	  for(j=0; j<TRAC_NUM; j++) {
	      SphP[i].TracAbundOut[j] = SphP[i].TracAbund[j];
	  } 
	  SphP[i].EntropyOut = SphP[i].Entropy;
#endif
	}
    }

  ngb_treeallocate(MAX_NGB);

  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  Flag_FullStep = 1;		/* to ensure that Peano-Hilbert order is done */

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  ngb_treebuild();		/* will build tree */

  setup_smoothinglengths();

  TreeReconstructFlag = 1;

  /* At this point, the entropy variable normally contains the 
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly. 
   * Once the density has been computed, we can convert thermal energy to entropy.
   *
   * N.B. density() is called by setup_smoothinglengths(), so by this point particle
   * densities are known
   */
#ifndef ISOTHERM_EQS
  /* With the polytopic EOS code, we directly specify pressure as a function of density, and
   * so there is no Entropy slot in the IC file
   */
#ifndef POLYTROPE
  if(header.flag_entropy_instead_u == 0)
    for(i = 0; i < N_gas; i++) {
#ifdef CHEMCOOL
      gamm1 = SphP[i].Gamma - 1.0;
      SphP[i].Entropy = gamm1 * SphP[i].Entropy / pow(SphP[i].Density / a3, gamm1);
#else
      SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif /* CHEMCOOL */
    }
#endif /* POLYTROPE */
#endif /* ISOTHERM_EQS */

#if defined RAYTRACE && defined CHEMCOOL
  if (All.PhotochemApprox == 2) {
    raytrace();
  }
#endif /* RAYTRACE && CHEMCOOL */
}


/*! This routine computes the mass content of the box and compares it to the
 *  specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e4)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {

      for(i = 0; i < N_gas; i++)
	{
	  if(P[i].ID < 0) /*SINK*/
	    continue;
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }
#ifndef TWODIMS
	  SphP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  SphP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
	}
    }

  density();
}


/*! If the code is run in glass-making mode, this function populates the
 *  simulation box with a Poisson sample of particles.
 */
#if (MAKEGLASS > 1)
void seed_glass(void)
{
  int i, k, n_for_this_task;
  double Range[3], LowerBound[3];
  double drandom, partmass;
  long long IDstart;

  All.TotNumPart = MAKEGLASS;
  partmass = All.Omega0 * (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G))
    * (All.BoxSize * All.BoxSize * All.BoxSize) / All.TotNumPart;

  All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);	/* sets the maximum number of particles that may */

  allocate_memory();

  header.npartTotal[1] = All.TotNumPart;
  header.mass[1] = partmass;

  if(ThisTask == 0)
    {
      printf("\nGlass initialising\nPartMass= %g\n", partmass);
      printf("TotNumPart= %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
    }

  /* set the number of particles assigned locally to this task */
  n_for_this_task = All.TotNumPart / NTask;

  if(ThisTask == NTask - 1)
    n_for_this_task = All.TotNumPart - (NTask - 1) * n_for_this_task;

  NumPart = 0;
  IDstart = 1 + (All.TotNumPart / NTask) * ThisTask;

  /* split the temporal domain into Ntask slabs in z-direction */

  Range[0] = Range[1] = All.BoxSize;
  Range[2] = All.BoxSize / NTask;
  LowerBound[0] = LowerBound[1] = 0;
  LowerBound[2] = ThisTask * Range[2];

  srand48(ThisTask);

  for(i = 0; i < n_for_this_task; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  drandom = drand48();

	  P[i].Pos[k] = LowerBound[k] + Range[k] * drandom;
	  P[i].Vel[k] = 0;
	}

      P[i].Mass = partmass;
      P[i].Type = 1;
      P[i].ID = IDstart + i;

      NumPart++;
    }
}
#endif
