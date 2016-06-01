#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"
#include "f2c.h"


/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are
 *  initialized to their proper values.
 */


/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
#ifdef TURBULENCE
   double tstart, tend;
   int bytes;
#endif 

  struct global_data_all_processes all;

  if(ThisTask == 0)
    {
      printf("\nThis is Gadget, version `%s'.\n", GADGETVERSION);
      printf("\nRunning on %d processors.\n", NTask);
    }

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

  allocate_commbuffers();	/* ... allocate buffer-memory for particle 
				   exchange during force computation */
  set_units();

#ifdef POLYTROPE
  define_EOS_table(All.WhichEOS);
#endif /* POLYTROPE */

#ifdef TURBULENCE
   if(!(xmatrix = malloc(bytes=sizeof(double)*(IFIELDSIZE*IFIELDSIZE*IFIELDSIZE))))
    {
      printf("failed to allocate memory for xmatrix: %d.\n",bytes);
      endrun(2);
    }
   if(!(ymatrix = malloc(bytes=sizeof(double)*(IFIELDSIZE*IFIELDSIZE*IFIELDSIZE))))
    {
      printf("failed to allocate memory for ymatrix: %d.\n",bytes);
      endrun(2);
    }
    if(!(zmatrix = malloc(bytes=sizeof(double)*(IFIELDSIZE*IFIELDSIZE*IFIELDSIZE))))
    {
      printf("failed to allocate memory for zmatrix: %d.\n",bytes);
      endrun(2);
    }
 
  tstart = second();
  icomp = 1;
  iiseed = All.Seed0;
  rsk_turbdriving_field();
  icomp = 2;
  iiseed = All.Seed1;
  rsk_turbdriving_field();
  icomp = 3;
  iiseed = All.Seed2;
  rsk_turbdriving_field();
  tend = second();
  All.CPU_Turbulence+= timediff(tstart,tend);
  if(ThisTask==0)
    {
      printf("initialization of turbulence field finished.\n");
      fflush(stdout);
    }
#endif 

#ifdef CHEMCOOL
  COOLR.deff    = All.H2RefDustEff;
#ifdef METALS_TG
  COOLR.abratio_o  = All.OxyAbund;
  COOLR.abratio_c  = All.CarbAbund;
  COOLR.abratio_si = All.SiAbund;
#endif
  COOLR.abundD  = All.DeutAbund;
  COOLR.abundmg = All.MgAbund; 
  COOLR.G0      = All.UVField;
  COOLR.phi_pah = All.PhiPAH;
  COOLR.dust_to_gas_ratio = All.DustToGasRatio;
  COOLR.AV_conversion_factor = All.AVConversionFactor;
  COOLR.cosmic_ray_ion_rate  = All.CosmicRayIonRate;
  COOLR.redshift   = All.InitRedshift;
  COOLR.AV_ext     = All.ExternalDustExtinction;
  COOLI.iphoto     = All.PhotochemApprox;
  COOLI.iflag_mn   = All.MNRateFlag;
  COOLI.iflag_ad   = All.ADRateFlag;
  COOLI.iflag_atom = All.AtomicFlag;
  COOLI.iflag_3bh2a = All.ThreeBodyFlagA;
  COOLI.iflag_3bh2b = All.ThreeBodyFlagB;
  COOLI.iflag_h3pra = All.H3PlusRateFlag;
#ifdef METALS_TG
  COOLI.iflag_fixed_ion = 1;
#else
  COOLI.iflag_fixed_ion = 0;
#endif

  chemcool_init();
#endif /* CHEMCOOL */

#ifdef RAYTRACE_TG
  COOLI.ray_flag_sun  = All.ray_flag_sun;
#endif

#if defined(PERIODIC) && (!defined(PMGRID) || defined(FORCETEST))
  ewald_init();
#endif

  open_outputfiles();

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	/* start-up seed */

#ifdef PMGRID
  long_range_init();
#endif

  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag == 2)
    {
      set_random_numbers();

      init();			/* ... read in initial model */
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets 
				   all variables in the struct `All'. 
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to 
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.  
				 */

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.BunchSizeForce = all.BunchSizeForce;
      All.BunchSizeDensity = all.BunchSizeDensity;
      All.BunchSizeHydro = all.BunchSizeHydro;
      All.BunchSizeDomain = all.BunchSizeDomain;

      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

      All.ErrTolForceAcc = all.ErrTolForceAcc;

      All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

      All.SnapFormat = all.SnapFormat;
      All.NumFilesPerSnapshot = all.NumFilesPerSnapshot;
      All.MaxNumNgbDeviation = all.MaxNumNgbDeviation;
      All.ArtBulkViscConst = all.ArtBulkViscConst;


      All.OutputListOn = all.OutputListOn;
      All.CourantFac = all.CourantFac;

      All.OutputListLength = all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);


      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);
    }

#ifdef PMGRID
  long_range_init_regionsize();
#endif

  if(All.ComovingIntegrationOn)
    init_drift_table();
/*
  if(RestartFlag == 2)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);
*/
  All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);

#ifdef CHEMCOOL
  All.Ti_nextnextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
#endif

  All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
  double meanweight;
#ifdef CHEMCOOL
  double gamm1;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
      printf("\n");
    }

#ifdef CHEMCOOL
  meanweight = compute_initial_molecular_weight();
  gamm1      = compute_initial_gamma() - 1.0;
  All.MinEgySpec = 1 / meanweight * (1.0 / gamm1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#else /* CHEMCOOL */
  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: we assume neutral gas here */
#ifdef ISOTHERM_EQS
  All.MinEgySpec = 0;
#else /* ISOTHERM_EQS */
#ifdef POLYTROPE
  All.MinEgySpec = 0;
#else /* POLYTROPE */
  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif /* POLYTROPE */
#endif /* ISOTHERM_EQS */
#endif /* CHEMCOOL */

if(RestartFlag == 2)
  All.MinEgySpec = 0;
if(ThisTask == 0)
  printf("All.MinEgySpec = %lg\n", All.MinEgySpec);

}



/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(ThisTask != 0)		/* only the root processor writes to the log files */
    return;

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");


  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
/*SINK*/
  sprintf(buf, "%s%s", All.OutputDir, All.SinkFile);
  if(!(FdSink = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
/*SINK*/
#ifdef FORCETEST
  if(RestartFlag == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      fclose(FdForceTest);
    }
#endif
}


/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
  if(ThisTask != 0)		/* only the root processor writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
/*SINK*/
  fclose(FdSink);
#ifdef FORCETEST
  fclose(FdForceTest);
#endif
}




/*! This function parses the parameterfile in a simple way.  Each paramater
 *  is defined by a keyword (`tag'), and can be either of type double, int,
 *  or character string.  The routine makes sure that each parameter
 *  appears exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag = 0;


  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;
      /*SINK*/
      strcpy(tag[nt], "SinkFile");
      addr[nt] = All.SinkFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeAllocFactor");
      addr[nt] = &All.TreeAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = DOUBLE;

#ifdef CHEMCOOL
      strcpy(tag[nt],"H2RefDustEff"); 
      addr[nt]=&All.H2RefDustEff;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"OxyAbund"); 
      addr[nt]=&All.OxyAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"CarbAbund"); 
      addr[nt]=&All.CarbAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"SiAbund"); 
      addr[nt]=&All.SiAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"DeutAbund"); 
      addr[nt]=&All.DeutAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"MgAbund"); 
      addr[nt]=&All.MgAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"UVField"); 
      addr[nt]=&All.UVField;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"PhiPAH"); 
      addr[nt]=&All.PhiPAH;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitDustTemp"); 
      addr[nt]=&All.InitDustTemp;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"DustToGasRatio"); 
      addr[nt]=&All.DustToGasRatio;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"AVConversionFactor"); 
      addr[nt]=&All.AVConversionFactor;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"CosmicRayIonRate"); 
      addr[nt]=&All.CosmicRayIonRate;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitRedshift"); 
      addr[nt]=&All.InitRedshift;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"ExternalDustExtinction"); 
      addr[nt]=&All.ExternalDustExtinction;
      id[nt++]=DOUBLE;

      strcpy(tag[nt], "PhotochemApprox");
      addr[nt] = &All.PhotochemApprox;
      id[nt++] = INT;

      strcpy(tag[nt], "ChemistryNetwork");
      addr[nt] = &All.ChemistryNetwork;
      id[nt++] = INT;

      strcpy(tag[nt], "ADRateFlag");
      addr[nt] = &All.ADRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "MNRateFlag");
      addr[nt] = &All.MNRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "AtomicFlag");
      addr[nt] = &All.AtomicFlag;
      id[nt++] = INT;

      strcpy(tag[nt], "ThreeBodyFlagA");
      addr[nt] = &All.ThreeBodyFlagA;
      id[nt++] = INT;

      strcpy(tag[nt], "ThreeBodyFlagB");
      addr[nt] = &All.ThreeBodyFlagB;
      id[nt++] = INT;

      strcpy(tag[nt], "H3PlusRateFlag");
      addr[nt] = &All.H3PlusRateFlag;
      id[nt++] = INT;

      strcpy(tag[nt],"InitMolHydroAbund"); 
      addr[nt]=&All.InitMolHydroAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHPlusAbund"); 
      addr[nt]=&All.InitHPlusAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitDIIAbund"); 
      addr[nt]=&All.InitDIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHDAbund"); 
      addr[nt]=&All.InitHDAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHeIIAbund"); 
      addr[nt]=&All.InitHeIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHeIIIAbund"); 
      addr[nt]=&All.InitHeIIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCIIAbund"); 
      addr[nt]=&All.InitCIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitSiIIAbund"); 
      addr[nt]=&All.InitSiIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitOIIAbund"); 
      addr[nt]=&All.InitOIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCOAbund"); 
      addr[nt]=&All.InitCOAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitC2Abund"); 
      addr[nt]=&All.InitC2Abund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitOHAbund"); 
      addr[nt]=&All.InitOHAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitH2OAbund"); 
      addr[nt]=&All.InitH2OAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitO2Abund"); 
      addr[nt]=&All.InitO2Abund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitHCOPlusAbund"); 
      addr[nt]=&All.InitHCOPlusAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCHAbund"); 
      addr[nt]=&All.InitCHAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCH2Abund"); 
      addr[nt]=&All.InitCH2Abund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitSiIIIAbund"); 
      addr[nt]=&All.InitSiIIIAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitCH3PlusAbund"); 
      addr[nt]=&All.InitCH3PlusAbund;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"InitMgPlusAbund"); 
      addr[nt]=&All.InitMgPlusAbund;
      id[nt++]=DOUBLE;
#endif /* CHEMCOOL */

#ifdef POLYTROPE
      strcpy(tag[nt], "WhichEOS");
      addr[nt] = &All.WhichEOS;
      id[nt++] = INT;

      strcpy(tag[nt], "EOSFullTableSize");
      addr[nt] = &All.EOSFullTableSize;
      id[nt++] = INT;
#endif /* POLYTROPE */

      /*akj*/

      strcpy(tag[nt],"TurbulenceOn");
      addr[nt]=&All.TurbulenceOn;
      id[nt++]=INT;

#ifdef TURBULENCE

      strcpy(tag[nt],"DrvTimestep");
      addr[nt]=&All.DrvTimestep;
      id[nt++]=DOUBLE;

      strcpy(tag[nt],"VelDispConst");
      addr[nt]=&All.VelDispConst;
      id[nt++]=INT;
 
      strcpy(tag[nt],"LDrv");
      addr[nt]=&All.LDrv;
      id[nt++]=DOUBLE;
 
      strcpy(tag[nt],"MachNumber");
      addr[nt]=&All.MachNumber;
      id[nt++]=DOUBLE;
 
      strcpy(tag[nt],"StartDriving");
      addr[nt]=&All.StartDriving;
      id[nt++]=DOUBLE;
 
      strcpy(tag[nt],"DrvIndx");
      addr[nt]=&All.DrvIndx;
      id[nt++]=INT;
 
      strcpy(tag[nt],"kMin");
      addr[nt]=&All.kMin;
      id[nt++]=INT;
 
      strcpy(tag[nt],"kMax");
      addr[nt]=&All.kMax;
      id[nt++]=INT;
 
      strcpy(tag[nt],"Seed0");
      addr[nt]=&All.Seed0;
      id[nt++]=INT;
 
      strcpy(tag[nt],"Seed1");
      addr[nt]=&All.Seed1;
      id[nt++]=INT;
 
      strcpy(tag[nt],"Seed2");
      addr[nt]=&All.Seed2;
      id[nt++]=INT;

      strcpy(tag[nt],"MWeight");
      addr[nt]=&All.MWeight;
      id[nt++]=DOUBLE;
#endif 

       /* SINK: add parameters for sinks */
      strcpy(tag[nt], "HSinkCreate");
      addr[nt] = &All.HSinkCreate;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "RInner");
      addr[nt] = &All.RInner;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ROuter");
      addr[nt] = &All.ROuter;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SinkCriticalDens");
      addr[nt] = &All.SinkCriticalDens;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SinkCriticalRedshift");
      addr[nt] = &All.SinkCriticalRedshift;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxNumSinks");
      addr[nt] = &All.MaxNumSinks;
      id[nt++] = INT;

      strcpy(tag[nt], "RefinementMass");
      addr[nt] = &All.RefinementMass;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "max_dens");
      addr[nt] = &All.max_dens;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ray_crit_dens");
      addr[nt] = &All.ray_crit_dens;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ray_r_max_sink");
      addr[nt] = &All.ray_r_max_sink;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ray_flag_sun");
      addr[nt] = &All.ray_flag_sun;
      id[nt++] = INT;

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case DOUBLE:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy(addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
	      system(buf3);
	    }
	}
      else
	{
	  printf("\nParameter file %s not found.\n\n", fname);
	  errorFlag = 2;
	}

      if(errorFlag != 2)
	for(i = 0; i < nt; i++)
	  {
	    if(*tag[i])
	      {
		printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
		errorFlag = 1;
	      }
	  }

      if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);


  if(All.NumFilesWrittenInParallel < 1)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be at least 1\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

#ifdef TURBULENCE
  if(All.TurbulenceOn==0)
    {
      if(ThisTask==0)
        {
          fprintf(stdout,"Code was compiled with turbulence switched on.\n");
          fprintf(stdout,"You must set `TurbulenceOn=1', or recompile the code.\n");
        }
          endrun(0);
    }
#else
  if(All.TurbulenceOn==1)
    {
      if(ThisTask==0)
        {
          fprintf(stdout,"Code was compiled with turbulence switched off.\n");
          fprintf(stdout,"You must set `TurbulenceOn=0', or recompile the code.\n");
        }
      endrun(0);
    }
#endif 


#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif


  if(All.TypeOfTimestepCriterion >= 1)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      endrun(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  endrun(0);
#endif
#endif

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}


/*! this function reads a table with a list of desired output times. The
 *  table does not have to be ordered in any way, but may not contain more
 *  than MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;
  do
    {
      if(fscanf(fd, " %lg ", &All.OutputListTimes[All.OutputListLength]) == 1)
	All.OutputListLength++;
      else
	break;
    }
  while(All.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax
 *  variable is increased, then the integer timeline needs to be
 *  adjusted. The approach taken here is to reduce the resolution of the
 *  integer timeline by factors of 2 until the new final time can be
 *  reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      endrun(556);
    }

  if(All.ComovingIntegrationOn)
    ti_end = log(TimeMax_new / All.TimeBegin) / All.Timebase_interval;
  else
    ti_end = (TimeMax_new - All.TimeBegin) / All.Timebase_interval;

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].Ti_endstep /= 2;
	}
    }

  All.TimeMax = TimeMax_new;
}

#ifdef POLYTROPE
/* Tabulate gas pressure as a function of number density of H nuclei.
   Densities are given in terms of log10 units, pressures in standard CGS */
void define_EOS_table(int option)
{
  double xi, dpres, ddens;
  double n0, n1, n2, n3, p0, p1, p2, p3;
  double g0, g1, gm, dergy;
  double density_orig[MAX_SIZE_EOS_TABLE];
  double pressure_orig[MAX_SIZE_EOS_TABLE];
  double energy_orig[MAX_SIZE_EOS_TABLE];
  int i, j;

  switch(option) {
  case 0: /* Test - isothermal gas, pure H2 (so gamma=7/5), T=20K */
    All.PolyIndexLowDensity  = 1.0;
    All.PolyIndexHighDensity = 1.0;

    All.EOSInitTableSize = 201;

    for(i=0; i<All.EOSInitTableSize; i++) {
      density_orig[i]  = 0.05 * i;
      pressure_orig[i] = 2.76e-15 * pow(10, density_orig[i]);
      energy_orig[i]   = pressure_orig[i] / 0.4;
    }

    break;
  case 1: /* Tsuribe & Omukai piecewise polytrope, Z = 1e-6 Z_sun */
    All.PolyIndexLowDensity  = 1.0;
    All.PolyIndexHighDensity = 5./3.;

    All.EOSInitTableSize = 241;

    /* For simplicity, assume gamma = 7/5 at n > 1e12 cm-3,
     *                        gamma = 5/3 at n < 1e8  cm-3,
     * and that gamma varies smoothly with log(n) in between
     */
    g0 = 5.0/3.0;
    g1 = 7.0/5.0;

    n0 = 4.0;
    p0 = 3.62e-10; /* 200K gas at n_H = 1e4 cm^-3 */

    n1 = 12.0;
    p1 = p0 * pow(pow(10, n1-n0), 1.1);

    n2 = 13.5;
    p2 = p1 * pow(pow(10, n2-n1), 0.5);
    
    n3 = 15.77;
    p3 = p2 * pow(pow(10, n3-n2), 1.1);

    for(i=0; i<All.EOSInitTableSize; i++) {
      density_orig[i]  = 4.0 + 0.05 * i;
      if (density_orig[i] < n1) {
	pressure_orig[i] = p0 * pow(pow(10, (density_orig[i] - n0)), 1.1);
      }
      else if (density_orig[i] < n2) {
	pressure_orig[i] = p1 * pow(pow(10, (density_orig[i] - n1)), 0.5);
      }
      else if (density_orig[i] < n3) {
	pressure_orig[i] = p2 * pow(pow(10, (density_orig[i] - n2)), 1.1);
      }
      else {
	pressure_orig[i] = p3 * pow(pow(10, (density_orig[i] - n3)), 1.4);
      }
      /* Internal energy */
      if (density_orig[i] < 8.0) {
        gm = g0;
      }
      else if (density_orig[i] > 12.0) {
	gm = g1;
      }
      else {
	gm = g0 + (g1 - g0) * (density_orig[i] - 8.0) / 4.0;
      }
      energy_orig[i] = pressure_orig[i] / (gm - 1.0);
    }
    break;
  case 2: /* Tsuribe & Omukai piecewise polytrope, Z = 1e-5 Z_sun */
    All.PolyIndexLowDensity  = 1.0;
    All.PolyIndexHighDensity = 5./3.;

    All.EOSInitTableSize = 241;

    /* See note on gamma for case 1 */
    g0 = 5.0/3.0;
    g1 = 7.0/5.0;

    n0 = 4.0;
    p0 = 3.62e-10; /* 200K gas at n_H = 1e4 cm^-3 */

    n1 = 11.0;
    p1 = p0 * pow(pow(10, n1-n0), 1.1);

    n2 = 13.1;
    p2 = p1 * pow(pow(10, n2-n1), 0.5);
    
    n3 = 14.47;
    p3 = p2 * pow(pow(10, n3-n2), 1.1);

    for(i=0; i<All.EOSInitTableSize; i++) {
      density_orig[i]  = 4.0 + 0.05 * i;
      if (density_orig[i] < n1) {
	pressure_orig[i] = p0 * pow(pow(10, (density_orig[i] - n0)), 1.1);
      }
      else if (density_orig[i] < n2) {
	pressure_orig[i] = p1 * pow(pow(10, (density_orig[i] - n1)), 0.5);
      }
      else if (density_orig[i] < n3) {
	pressure_orig[i] = p2 * pow(pow(10, (density_orig[i] - n2)), 1.1);
      }
      else {
	pressure_orig[i] = p3 * pow(pow(10, (density_orig[i] - n3)), 1.4);
      }
      /* Internal energy */
      if (density_orig[i] < 8.0) {
        gm = g0;
      }
      else if (density_orig[i] > 12.0) {
	gm = g1;
      }
      else {
	gm = g0 + (g1 - g0) * (density_orig[i] - 8.0) / 4.0;
      }
      energy_orig[i] = pressure_orig[i] / (gm - 1.0);
    }
    break;
  case 3: /* Tsuribe & Omukai piecewise polytrope, Z = 1e-4 Z_sun */
    All.PolyIndexLowDensity  = 1.0;
    All.PolyIndexHighDensity = 5./3.;

    All.EOSInitTableSize = 241;

    /* See note on gamma for case 1 */
    g0 = 5.0/3.0;
    g1 = 7.0/5.0;

    n0 = 4.0;
    p0 = 3.62e-10; /* 200K gas at n_H = 1e4 cm^-3 */

    n1 = 10.0;
    p1 = p0 * pow(pow(10, n1-n0), 1.1);

    n2 = 12.7;
    p2 = p1 * pow(pow(10, n2-n1), 0.5);
    
    n3 = 13.27;
    p3 = p2 * pow(pow(10, n3-n2), 1.1);

    for(i=0; i<All.EOSInitTableSize; i++) {
      density_orig[i]  = 4.0 + 0.05 * i;
      if (density_orig[i] < n1) {
	pressure_orig[i] = p0 * pow(pow(10, (density_orig[i] - n0)), 1.1);
      }
      else if (density_orig[i] < n2) {
	pressure_orig[i] = p1 * pow(pow(10, (density_orig[i] - n1)), 0.5);
      }
      else if (density_orig[i] < n3) {
	pressure_orig[i] = p2 * pow(pow(10, (density_orig[i] - n2)), 1.1);
      }
      else {
	pressure_orig[i] = p3 * pow(pow(10, (density_orig[i] - n3)), 1.4);
      }
      /* Internal energy */
      if (density_orig[i] < 8.0) {
        gm = g0;
      }
      else if (density_orig[i] > 12.0) {
	gm = g1;
      }
      else {
	gm = g0 + (g1 - g0) * (density_orig[i] - 8.0) / 4.0;
      }
      energy_orig[i] = pressure_orig[i] / (gm - 1.0);
    }
    break;
  case 4: /* Omukai et al, 2005, Z=0 */
    All.PolyIndexLowDensity  = 1.0;
    All.PolyIndexHighDensity = 5./3.;

    All.EOSInitTableSize = 17;
    density_orig[0] = 3.93; pressure_orig[0] = 2.68e-10; energy_orig[0] = 4.02e-10;
    density_orig[1] = 4.63; pressure_orig[1] = 1.34e-09; energy_orig[1] = 2.01e-09;
    density_orig[2] = 5.54; pressure_orig[2] = 1.35e-08; energy_orig[2] = 2.02e-08;
    density_orig[3] = 6.24; pressure_orig[3] = 8.81e-08; energy_orig[3] = 1.32e-07;
    density_orig[4] = 7.01; pressure_orig[4] = 6.77e-07; energy_orig[4] = 1.02e-06;
    density_orig[5] = 8.06; pressure_orig[5] = 1.10e-05; energy_orig[5] = 1.65e-05;
    density_orig[6] = 8.83; pressure_orig[6] = 7.19e-05; energy_orig[6] = 1.08e-04;
    density_orig[7] = 9.59; pressure_orig[7] = 4.39e-04; energy_orig[7] = 6.67e-04;
    density_orig[8] = 10.29; pressure_orig[8] = 2.24e-03; energy_orig[8] = 3.75e-03;
    density_orig[9] = 11.20; pressure_orig[9] = 1.71e-02; energy_orig[9] = 3.64e-02;
    density_orig[10] = 12.18; pressure_orig[10] = 1.58e-01; energy_orig[10] = 3.75e-01;
    density_orig[11] = 13.51; pressure_orig[11] = 3.97e+00; energy_orig[11] = 9.38e+00;
    density_orig[12] = 14.49; pressure_orig[12] = 4.44e+01; energy_orig[12] = 1.05e+02;
    density_orig[13] = 15.26; pressure_orig[13] = 2.75e+02; energy_orig[13] = 6.51e+02;
    density_orig[14] = 16.10; pressure_orig[14] = 1.71e+03; energy_orig[14] = 4.04e+03;
    density_orig[15] = 16.80; pressure_orig[15] = 1.00e+04; energy_orig[15] = 2.37e+04;
    density_orig[16] = 18.48; pressure_orig[16] = 7.74e+05; energy_orig[16] = 1.83e+06;
    break;
  case 5: /* Omukai et al, 2005, Z = 1e-6 Z_sun */
    All.PolyIndexLowDensity  = 1.0;
    All.PolyIndexHighDensity = 5./3.;

    All.EOSInitTableSize = 18;
    density_orig[0] = 3.93; pressure_orig[0] = 2.68e-10; energy_orig[0] = 4.02e-10;
    density_orig[1] = 4.63; pressure_orig[1] = 1.34e-09; energy_orig[1] = 2.01e-09;
    density_orig[2] = 5.54; pressure_orig[2] = 1.35e-08; energy_orig[2] = 2.02e-08;
    density_orig[3] = 6.24; pressure_orig[3] = 8.81e-08; energy_orig[3] = 1.32e-07;
    density_orig[4] = 7.01; pressure_orig[4] = 6.76e-07; energy_orig[4] = 1.02e-06;
    density_orig[5] = 8.41; pressure_orig[5] = 2.33e-05; energy_orig[5] = 3.50e-05;
    density_orig[6] = 9.24; pressure_orig[6] = 1.68e-04; energy_orig[6] = 2.54e-04;
    density_orig[7] = 9.94; pressure_orig[7] = 9.27e-04; energy_orig[7] = 1.46e-03;
    density_orig[8] = 10.50; pressure_orig[8] = 3.42e-03; energy_orig[8] = 5.97e-03;
    density_orig[9] = 11.27; pressure_orig[9] = 1.75e-02; energy_orig[9] = 3.98e-02;
    density_orig[10] = 12.18; pressure_orig[10] = 1.35e-01; energy_orig[10] = 3.19e-01;
    density_orig[11] = 13.09; pressure_orig[11] = 1.04e+00; energy_orig[11] = 2.45e+00;
    density_orig[12] = 14.00; pressure_orig[12] = 1.04e+01; energy_orig[12] = 2.47e+01;
    density_orig[13] = 14.98; pressure_orig[13] = 1.17e+02; energy_orig[13] = 2.76e+02;
    density_orig[14] = 15.54; pressure_orig[14] = 4.46e+02; energy_orig[14] = 1.06e+03;
    density_orig[15] = 15.89; pressure_orig[15] = 9.46e+02; energy_orig[15] = 2.24e+03;
    density_orig[16] = 16.52; pressure_orig[16] = 4.49e+03; energy_orig[16] = 1.06e+04;
    density_orig[17] = 18.27; pressure_orig[17] = 4.29e+05; energy_orig[17] = 1.02e+06;
    break;
  case 6:  /* Omukai et al, 2005, Z = 1e-5 Z_sun */
    All.PolyIndexLowDensity  = 1.0;
    All.PolyIndexHighDensity = 5./3.;

    All.EOSInitTableSize = 22;
    density_orig[0] = 3.93; pressure_orig[0] = 2.68e-10; energy_orig[0] = 4.02e-10;
    density_orig[1] = 4.84; pressure_orig[1] = 2.06e-09; energy_orig[1] = 3.09e-09;
    density_orig[2] = 5.96; pressure_orig[2] = 3.35e-08; energy_orig[2] = 5.03e-08;
    density_orig[3] = 6.80; pressure_orig[3] = 2.86e-07; energy_orig[3] = 4.30e-07;
    density_orig[4] = 7.22; pressure_orig[4] = 7.51e-07; energy_orig[4] = 1.13e-06;
    density_orig[5] = 8.13; pressure_orig[5] = 5.75e-06; energy_orig[5] = 8.67e-06;
    density_orig[6] = 8.76; pressure_orig[6] = 2.55e-05; energy_orig[6] = 3.88e-05;
    density_orig[7] = 9.38; pressure_orig[7] = 1.43e-04; energy_orig[7] = 2.25e-04;
    density_orig[8] = 10.01; pressure_orig[8] = 7.33e-04; energy_orig[8] = 1.27e-03;
    density_orig[9] = 10.64; pressure_orig[9] = 3.54e-03; energy_orig[9] = 6.78e-03;
    density_orig[10] = 11.34; pressure_orig[10] = 1.58e-02; energy_orig[10] = 3.73e-02;
    density_orig[11] = 11.83; pressure_orig[11] = 4.62e-02; energy_orig[11] = 1.09e-01;
    density_orig[12] = 12.32; pressure_orig[12] = 1.03e-01; energy_orig[12] = 2.44e-01;
    density_orig[13] = 12.74; pressure_orig[13] = 1.87e-01; energy_orig[13] = 4.41e-01;
    density_orig[14] = 13.30; pressure_orig[14] = 4.65e-01; energy_orig[14] = 1.10e+00;
    density_orig[15] = 14.14; pressure_orig[15] = 2.89e+00; energy_orig[15] = 6.83e+00;
    density_orig[16] = 14.91; pressure_orig[16] = 1.99e+01; energy_orig[16] = 4.71e+01;
    density_orig[17] = 15.75; pressure_orig[17] = 2.91e+02; energy_orig[17] = 6.88e+02;
    density_orig[18] = 16.52; pressure_orig[18] = 3.43e+03; energy_orig[18] = 8.12e+03;
    density_orig[19] = 17.22; pressure_orig[19] = 2.64e+04; energy_orig[19] = 6.23e+04;
    density_orig[20] = 17.92; pressure_orig[20] = 1.72e+05; energy_orig[20] = 4.08e+05;
    density_orig[21] = 19.03; pressure_orig[21] = 2.96e+06; energy_orig[21] = 7.01e+06;
    break;
  default: 
    printf("Invalid EOS selected!\n");
    endrun(1);
    break;
  }

  All.MinTabulatedDensity = density_orig[0];
  All.MaxTabulatedDensity = density_orig[All.EOSInitTableSize-1];
  All.EOSDensDel = (All.MaxTabulatedDensity - All.MinTabulatedDensity) /  (All.EOSFullTableSize - 1);

  for (i=0 ; i<All.EOSFullTableSize; ++i){
    xi = All.MinTabulatedDensity + i * All.EOSDensDel;
    for (j=All.EOSInitTableSize-1 ; j>=0 ; --j) {
      if (xi >= density_orig[j]) break;
    }
    ddens = density_orig[j+1]  - density_orig[j];
    dpres = log10(pressure_orig[j+1]) - log10(pressure_orig[j]);
    dergy = log10(energy_orig[j+1]) - log10(energy_orig[j]);
    All.EOSDensity[i]  = xi;
    All.EOSPressure[i] = pow(10, log10(pressure_orig[j]) + (xi - density_orig[j]) * dpres / ddens);
    All.EOSEnergy[i]   = pow(10, log10(energy_orig[j]) + (xi - density_orig[j]) * dergy / ddens);
  }
  return;
}
#endif /* POLYTROPE */
