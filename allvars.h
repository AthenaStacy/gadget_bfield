/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define's, typedef's, and enum's
 *     - add #include "allvars.h", delete the #ifndef ALLVARS_H conditional
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
//#include </share/apps/gsl/1.13/include/gsl/gsl_rng.h>
#include <gsl/gsl_rng.h>
#include "tags.h"
#include "chemcool_consts.h"

#define  GADGETVERSION   "2.0"   /*!< code version string */

//#define  TIMEBASE        (1<<32)
#define TIMEBASE (long long int)(pow(2,50)) 
//#define TIMEBASE (long long int)(pow(2,44))
                                 /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                  *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds to 2^29
                                  */

#define  MAXTOPNODES     500000   /*!< Maximum number of nodes in the top-level tree used for domain decomposition */


typedef  long long  peanokey;    /*!< defines the variable type used for Peano-Hilbert keys */

#define  BITS_PER_DIMENSION 18	 /*!< Bits per dimension available for Peano-Hilbert order. 
				      Note: If peanokey is defined as type int, the allowed maximum is 10.
				      If 64-bit integers are used, the maximum is 21 */

#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))  /*!< The number of different Peano-Hilbert cells */


#define  RNDTABLE         3000   /*!< gives the length of a table with random numbers, refreshed at every timestep.
				      This is used to allow application of random numbers to a specific particle
				      in a way that is independent of the number of processors used. */
#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  MAXLEN_FILENAME  100    /*!< Maximum number of characters for filenames (including the full path) */

#define dt_cool_min 3.e5

#ifndef CHEMCOOL
#ifdef   ISOTHERM_EQS
#define  GAMMA         (1.0)     /*!< index for isothermal gas */
#else
#define  GAMMA         (5.0/3)   /*!< adiabatic index of simulated gas */
#endif /* ISOTHERM_EQS */

#define  GAMMA_MINUS1  (GAMMA-1)
#endif  /* CHEMCOOL */

#define  HYDROGEN_MASSFRAC 0.76  /*!< mass fraction of hydrogen, relevant only for radiative cooling */

#define ABHE ((1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC)

#define Z_SOLAR 0.02

/* Some physical constants in cgs units */

#define  GRAVITY           6.672e-8   /*!< Gravitational constant (in cgs units) */
#define  SOLAR_MASS        1.989e33
#define  SOLAR_LUM         3.826e33
#define  RAD_CONST         7.565e-15
#define  AVOGADRO          6.0222e23
#define  GAS_CONST         8.31425e7
#define  C                 2.9979e10
#define  PLANCK            6.6262e-27
#define  CM_PER_MPC        3.085678e24
#define  ELECTRONMASS      9.10953e-28
#define  THOMPSON          6.65245e-25
#define  ELECTRONCHARGE    4.8032e-10
#define  HUBBLE            3.2407789e-18	/* in h/sec */

#ifndef CHEMCOOL
#define  BOLTZMANN         1.3806e-16
#define  PROTONMASS        1.6726e-24
#endif

/* Some conversion factors */

#define  SEC_PER_MEGAYEAR  3.155e13
#define  SEC_PER_YEAR      3.155e7

#ifdef TURBULENCE
#define  THIRD            (1.0/3.0)
#define  PI               3.14159265358979323846
#define  PI_INV           (1.0/PI)
#define  LN2              0.69314718
#define  DOTS 64 /*akj*/
#define  IFIELDSIZE 64  /*akj*/
#define  TDRVSTEPFACTOR 0.05  /*akj*/
#endif

#ifndef ASMTH
#define ASMTH 1.25  /*!< ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#endif

#ifndef RCUT
#define RCUT  4.5   /*!< RCUT gives the maximum distance (in units of the scale used for the force split) out to 
                         which short-range forces are evaluated in the short-range tree walk. */
#endif

#define MAX_NGB             20000  /*!< defines maximum length of neighbour list */

#define MAXLEN_OUTPUTLIST   500	   /*!< maxmimum number of entries in list of snapshot output times */

#define DRIFT_TABLE_LENGTH  1000   /*!< length of the lookup table used to hold the drift and kick factors */ 

#define MAXITER             1500    /*!< maxmimum number of steps for SPH neighbour iteration */


#define numtype             20 //total number of particle types in the simulation

#ifdef DOUBLEPRECISION             /*!< If defined, the variable type FLOAT is set to "double", otherwise to FLOAT */
#define FLOAT double
#else
#define FLOAT float
#endif


#ifndef  TWODIMS
#define  NUMDIMS 3                                      /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470                 /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786                 /*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */ 
#else
#define  NUMDIMS 2                                      /*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)         /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI                           /*!< Coefficient for kernel normalization. */
#endif

#ifdef TURBULENCE
extern int    icomp, iiseed;          /*akj*/   
extern double *xmatrix;
extern double *ymatrix;
extern double *zmatrix;
#endif

extern int ThisTask;		/*!< the rank of the local processor */
extern int NTask;               /*!< number of processors */
extern int PTask;	        /*!< smallest integer such that NTask <= 2^PTask */

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
extern int N_sinks;             /*!< number of sink particles on the LOCAL processor */
extern long long Ntype[numtype];      /*!< total number of particles of each type */
extern int NtypeLocal[numtype];       /*!< local number of particles of each type */

extern int NumForceUpdate;      /*!< number of active particles on local processor in current timestep  */
extern int NumSphUpdate;        /*!< number of active SPH particles on local processor in current timestep  */

extern double CPUThisRun;	/*!< Sums the CPU time for the process (current submission only) */


extern int RestartFlag;         /*!< taken from command line used to start code. 0 is normal start-up from
                                     initial conditions, 1 is resuming a run from a set of restart files, while 2
                                     marks a restart from a snapshot file. */

extern char *Exportflag;        /*!< Buffer used for flagging whether a particle needs to be exported to another process */

extern int  *Ngblist;           /*!< Buffer to hold indices of neighbours retrieved by the neighbour search routines */

extern int TreeReconstructFlag; /*!< Signals that a new tree needs to be constructed */

extern int Flag_FullStep;       /*!< This flag signals that the current step involves all particles */


extern gsl_rng *random_generator; /*!< the employed random number generator of the GSL library */

extern double RndTable[RNDTABLE]; /*!< Hold a table with random numbers, refreshed every timestep */


extern double DomainCorner[3];    /*!< gives the lower left corner of simulation volume */
extern double DomainCenter[3];    /*!< gives the center of simulation volume */
extern double DomainLen;          /*!< gives the (maximum) side-length of simulation volume */
extern double DomainFac;          /*!< factor used for converting particle coordinates to a Peano-Hilbert mesh covering the simulation volume */
extern int    DomainMyStart;      /*!< first domain mesh cell that resides on the local processor */
extern int    DomainMyLast;       /*!< last domain mesh cell that resides on the local processor */
extern int    *DomainStartList;   /*!< a table that lists the first domain mesh cell for all processors */
extern int    *DomainEndList;     /*!< a table that lists the last domain mesh cell for all processors */
extern double *DomainWork;        /*!< a table that gives the total "work" due to the particles stored by each processor */
extern int    *DomainCount;       /*!< a table that gives the total number of particles held by each processor */
extern int    *DomainCountSph;    /*!< a table that gives the total number of SPH particles held by each processor */

extern int    *DomainTask;        /*!< this table gives for each leaf of the top-level tree the processor it was assigned to */
extern int    *DomainNodeIndex;   /*!< this table gives for each leaf of the top-level tree the corresponding node of the gravitational tree */
extern FLOAT  *DomainTreeNodeLen; /*!< this table gives for each leaf of the top-level tree the side-length of the corresponding node of the gravitational tree */
extern FLOAT  *DomainHmax;        /*!< this table gives for each leaf of the top-level tree the maximum SPH smoothing length among the particles of the corresponding node of the gravitational tree */

extern struct DomainNODE
{
  FLOAT s[3];                     /*!< center-of-mass coordinates */
  FLOAT vs[3];                    /*!< center-of-mass velocities */
  FLOAT mass;                     /*!< mass of node */
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int   bitflags;                 /*!< this bit-field encodes the particle type with the largest softening among the particles of the nodes, and whether there are particles with different softening in the node */
#else
  FLOAT maxsoft;                  /*!< hold the maximum gravitational softening of particles in the 
                                       node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
#endif
#endif
}
 *DomainMoment;                   /*!< this table stores for each node of the top-level tree corresponding node data from the gravitational tree */

extern peanokey *DomainKeyBuf;    /*!< this points to a buffer used during the exchange of particle data */

extern peanokey *Key;             /*!< a table used for storing Peano-Hilbert keys for particles */
extern peanokey *KeySorted;       /*!< holds a sorted table of Peano-Hilbert keys for all particles, used to construct top-level tree */


extern int NTopnodes;             /*!< total number of nodes in top-level tree */
extern int NTopleaves;            /*!< number of leaves in top-level tree. Each leaf can be assigned to a different processor */

extern struct topnode_data
{
  int Daughter;                   /*!< index of first daughter cell (out of 8) of top-level node */
  int Pstart;                     /*!< for the present top-level node, this gives the index of the first node in the concatenated list of topnodes collected from all processors */
  int Blocks;                     /*!< for the present top-level node, this gives the number of corresponding nodes in the concatenated list of topnodes collected from all processors */
  int Leaf;                       /*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
  peanokey Size;                  /*!< number of Peano-Hilbert mesh-cells represented by top-level node */
  peanokey StartKey;              /*!< first Peano-Hilbert key in top-level node */
  long long Count;                /*!< counts the number of particles in this top-level node */
}
 *TopNodes;                       /*!< points to the root node of the top-level tree */


extern double TimeOfLastTreeConstruction; /*!< holds what it says, only used in connection with FORCETEST */



/* variables for input/output, usually only used on process 0 */

extern char ParameterFile[MAXLEN_FILENAME];  /*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo;       /*!< file handle for info.txt log-file. */
extern FILE *FdEnergy;     /*!< file handle for energy.txt log-file. */
extern FILE *FdTimings;    /*!< file handle for timings.txt log-file. */
extern FILE *FdCPU;        /*!< file handle for cpu.txt log-file. */
/*Sink*/
extern FILE *FdSink;        /*!< file handle for sink.txt log-file. */


#ifdef FORCETEST
extern FILE *FdForceTest;  /*!< file handle for forcetest.txt log-file. */
#endif


extern double DriftTable[DRIFT_TABLE_LENGTH];      /*!< table for the cosmological drift factors */
extern double GravKickTable[DRIFT_TABLE_LENGTH];   /*!< table for the cosmological kick factor for gravitational forces */
extern double HydroKickTable[DRIFT_TABLE_LENGTH];  /*!< table for the cosmological kick factor for hydrodynmical forces */

extern void *CommBuffer;   /*!< points to communication buffer, which is used in the domain decomposition, the
                                parallel tree-force computation, the SPH routines, etc. */



/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!< total particle numbers (global value) */
  long long TotN_gas;		/*!< total gas particle number (global value) */
  long long TotN_sinks;         /*!< total sink particle number (global value) */


  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one processor. */

  double BoxSize;               /*!< Boxsize in case periodic boundary conditions are used */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;/*!< maximum number of files that may be written simultaneously when
                                     writing/reading restart-files, or when writing snapshot files */ 

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force algorithm  */
  int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density computation */
  int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH hydrodynamical force computation */
  int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain decomposition */
  /*SINK*/
  int    BunchSizeSink;         /*!< number of particles fitting into the communication buffer in the sink computation */


  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				     NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				     the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				     the maximum(!) number of particles.  Note: A typical local tree for N
				     particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  double DesNumNgb;             /*!< Desired number of SPH neighbours */
  double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */


  /* some force counters  */

  long long TotNumOfForces;	             /*!< counts total number of force computations  */
  long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */


  /* system of units  */

  double G;                        /*!< Gravity-constant in internal units */
  double UnitTime_in_s;   	   /*!< factor to convert internal time unit to seconds/h */
  double UnitMass_in_g;            /*!< factor to convert internal mass unit to grams/h */
  double UnitVelocity_in_cm_per_s; /*!< factor to convert intqernal velocity unit to cm/sec */
  double UnitLength_in_cm;         /*!< factor to convert internal length unit to cm/h */
  double UnitPressure_in_cgs;      /*!< factor to convert internal pressure unit to cgs units (little 'h' still around!) */
  double UnitDensity_in_cgs;       /*!< factor to convert internal length unit to g/cm^3*h^2 */
  double UnitCoolingRate_in_cgs;   /*!< factor to convert internal cooling rate to cgs units */
  double UnitEnergy_in_cgs;        /*!< factor to convert internal energy to cgs units */
  double UnitTime_in_Megayears;    /*!< factor to convert internal time to megayears/h */
  double GravityConstantInternal;  /*!< If set to zero in the parameterfile, the internal value of the
                                        gravitational constant is set to the Newtonian value based on the system of
                                        units specified. Otherwise the value provided is taken as internal gravity constant G. */


  /* Cosmological parameters */

  double Hubble;       /*!< Hubble-constant in internal units */
  double Omega0;       /*!< matter density in units of the critical density (at z=0)*/
  double OmegaLambda;  /*!< vaccum energy density relative to crictical density (at z=0) */
  double OmegaBaryon;  /*!< baryon density in units of the critical density (at z=0)*/
  double HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute physical values for cooling physics */
  double HubbleA; 
 

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
  int ResubmitOn;               /*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative criterion */
  int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported right now - unlike gadget-1.1) */
  int OutputListOn;             /*!< flags that output times are listed in a specified file */

  int TurbulenceOn;
#ifdef TURBULENCE
  double DrvTimestep;       /* timestep between driving events (parameter) */
  double DeltaTimeDrv;     /* time between driving events -- variable - not parameter  */
  int VelDispConst;         /* switch: 1 = const vel disp, 0 = LDrv given */
  double LDrv;              /* driving luminosity */
  double MachNumber;            /* Mach number */
  double StartDriving;    /* time to turn on driver */
  int    DrvIndx;         /* power spectrum of driving field */
  int    kMin;
  int    kMax;
  int    Seed0;           /* random number */
  int    Seed1;           /* random number */
  int    Seed2;           /* random number */
  double MWeight;         /* molecular weight - in Gadget2 = 1.0 not yet implemented everywhere*/
#endif

  /* Parameters determining output frequency */

  int SnapshotFileCount;        /*!< number of snapshot that is written next */
  double TimeBetSnapshot;       /*!< simulation time interval between snapshot files */
  double TimeOfFirstSnapshot;   /*!< simulation time of first snapshot files */
  double CpuTimeBetRestartFile; /*!< cpu-time between regularly generated restart files */
  double TimeLastRestartFile;   /*!< cpu-time when last restart-file was written */
  double TimeBetStatistics;     /*!< simulation time interval between computations of energy statistics */
  double TimeLastStatistics;    /*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;         /*!< counts the number of system steps taken up to this point */


  /* Current time of the simulation, global step, and end of simulation */

  double Time;                  /*!< current time of the simulation */
  double TimeBegin;             /*!< time of initial conditions of the simulation */
  double TimeStep;              /*!< difference between current times of previous and current timestep */
  double TimeMax;	        /*!< marks the point of time until the simulation is to be evolved */


  /* variables for organizing discrete timeline */

  double Timebase_interval;     /*!< factor to convert from floating point time interval to integer timeline */
  long long int Ti_Current;               /*!< current time on integer timeline */ 
  long long int Ti_Prev;
  long long int Ti_nextoutput;            /*!< next output time on integer timeline */
#ifdef CHEMCOOL
  long long int Ti_nextnextoutput;        /*!< next+1 output time on integer timeline */
#endif
  long long int PresentMinStep;           /*!< If FLEXSTEPS is used, particle timesteps are chosen as multiples of the present minimum timestep. */
  long long int PresentMaxStep;		/*!< If FLEXSTEPS is used, this is the maximum timestep in timeline units, rounded down to the next power 2 division */
#ifdef PMGRID
  long long int PM_Ti_endstep;            /*!< begin of present long-range timestep */
  long long int PM_Ti_begstep;            /*!< end of present long-range timestep */
#endif


  /* Placement of PM grids */

#ifdef PMGRID
  double Asmth[2];              /*!< Gives the scale of the long-range/short-range split (in mesh-cells), both for the coarse and the high-res mesh */
  double Rcut[2];               /*!< Gives the maximum radius for which the short-range force is evaluated with the tree (in mesh-cells), both for the coarse and the high-res mesh */
  double Corner[2][3];          /*!< lower left corner of coarse and high-res PM-mesh */
  double UpperCorner[2][3];     /*!< upper right corner of coarse and high-res PM-mesh */
  double Xmintot[2][3];         /*!< minimum particle coordinates both for coarse and high-res PM-mesh */
  double Xmaxtot[2][3];         /*!< maximum particle coordinates both for coarse and high-res PM-mesh */
  double TotalMeshSize[2];      /*!< total extension of coarse and high-res PM-mesh */
#endif


  /* Variables that keep track of cumulative CPU consumption */

#ifdef TURBULENCE
  double  CPU_Turbulence;
#endif
  double TimeLimitCPU;          /*!< CPU time limit as defined in parameterfile */
  double CPU_TreeConstruction;  /*!< time spent for constructing the gravitational tree */
  double CPU_TreeWalk;          /*!< actual time spent for pure tree-walks */
  double CPU_Gravity;           /*!< cumulative time used for gravity computation (tree-algorithm only) */
  double CPU_Potential;         /*!< time used for computing gravitational potentials */
  double CPU_Domain;            /*!< cumulative time spent for domain decomposition */
  double CPU_Snapshot;          /*!< time used for writing snapshot files */
  double CPU_Total;             /*!< cumulative time spent for domain decomposition */
  double CPU_CommSum;           /*!< accumulated time used for communication, and for collecting partial results, in tree-gravity */
  double CPU_Imbalance;         /*!< cumulative time lost accross all processors as work-load imbalance in gravitational tree */
  double CPU_HydCompWalk;       /*!< time used for actual SPH computations, including neighbour search */
  double CPU_HydCommSumm;       /*!< cumulative time used for communication in SPH, and for collecting partial results */
  double CPU_HydImbalance;      /*!< cumulative time lost due to work-load imbalance in SPH */
  double CPU_Hydro;             /*!< cumulative time spent for SPH related computations */
  double CPU_EnsureNgb;         /*!< time needed to iterate on correct neighbour numbers */
  double CPU_Predict;           /*!< cumulative time to drift the system forward in time, including dynamic tree updates */
  double CPU_TimeLine;          /*!< time used for determining new timesteps, and for organizing the timestepping, including kicks of active particles */
  double CPU_PM;                /*!< time used for long-range gravitational force */
  double CPU_Peano;             /*!< time required to establish Peano-Hilbert order */
  double CPU_Sinks;             /*!< time required for sink calculations */
#ifdef CHEMCOOL
  double CPU_Chemcool;          /*!< time used in chemistry & cooling subroutines */
#ifdef RAYTRACE
  double CPU_Raytrace;          /*!< time used in raytrace code */
#endif /* RAYTRACE */
#endif /* CHEMCOOL */

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                                     timestep is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep;       /*!< minimum allowed timestep. Normally, the simulation terminates if the
                                     timestep determined by the timestep criteria falls below this limit. */ 
  double MaxSizeTimestep;       /*!< maximum allowed timestep */

  double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
                                     in comoving coordinates.  To this end, the code computes the rms velocity
                                     of all particles, and limits the timestep such that the rms displacement
                                     is a fraction of the mean particle separation (determined from the
                                     particle mass and the cosmological parameters). This parameter specifies
                                     this fraction. */

  double CourantFac;		/*!< SPH-Courant factor */



  /* frequency of tree reconstruction/domain decomposition */

  double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */


  /* Gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening length).
   * Five groups of particles are supported 0="gas", 1="halo", 2="disk", 3="bulge", 4="stars", 5="bndry"
   */

  double MinGasHsmlFractional;  /*!< minimum allowed SPH smoothing length in units of SPH gravitational softening length */
  double MinGasHsml;            /*!< minimum allowed SPH smoothing length */


  double SofteningGas;          /*!< comoving gravitational softening lengths for type 0 */ 
  double SofteningHalo;         /*!< comoving gravitational softening lengths for type 1 */ 
  double SofteningDisk;         /*!< comoving gravitational softening lengths for type 2 */ 
  double SofteningBulge;        /*!< comoving gravitational softening lengths for type 3 */ 
  double SofteningStars;        /*!< comoving gravitational softening lengths for type 4 */ 
  double SofteningBndry;        /*!< comoving gravitational softening lengths for type 5 */ 

  double SofteningGasMaxPhys;   /*!< maximum physical softening length for type 0 */ 
  double SofteningHaloMaxPhys;  /*!< maximum physical softening length for type 1 */ 
  double SofteningDiskMaxPhys;  /*!< maximum physical softening length for type 2 */ 
  double SofteningBulgeMaxPhys; /*!< maximum physical softening length for type 3 */ 
  double SofteningStarsMaxPhys; /*!< maximum physical softening length for type 4 */ 
  double SofteningBndryMaxPhys; /*!< maximum physical softening length for type 5 */ 

  double SofteningTable[numtype];     /*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[numtype];     /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  double MassTable[numtype];          /*!< Table with particle masses for particle types with equal mass.
                                     If particle masses are all equal for one type, the corresponding entry in MassTable 
                                     is set to this value, allowing the size of the snapshot files to be reduced. */
  


  /* some filenames */

  char InitCondFile[MAXLEN_FILENAME];          /*!< filename of initial conditions */
  char OutputDir[MAXLEN_FILENAME];             /*!< output directory of the code */
  char SnapshotFileBase[MAXLEN_FILENAME];      /*!< basename to construct the names of snapshot files */
  char EnergyFile[MAXLEN_FILENAME];            /*!< name of file with energy statistics */
  char CpuFile[MAXLEN_FILENAME];               /*!< name of file with cpu-time statistics */
  char InfoFile[MAXLEN_FILENAME];              /*!< name of log-file with a list of the timesteps taken */
  char TimingsFile[MAXLEN_FILENAME];           /*!< name of file with performance metrics of gravitational tree algorithm */
  char RestartFile[MAXLEN_FILENAME];           /*!< basename of restart-files */
  /*SINK*/
  char SinkFile[MAXLEN_FILENAME];               /*!< name of file with sink statistics */
  char ResubmitCommand[MAXLEN_FILENAME];       /*!< name of script-file that will be executed for automatic restart */
  char OutputListFilename[MAXLEN_FILENAME];    /*!< name of file with list of desired output times */

  double OutputListTimes[MAXLEN_OUTPUTLIST];   /*!< table with desired output times */
  int OutputListLength;                        /*!< number of output times stored in the table of desired output times */

#ifdef POLYTROPE
#define MAX_SIZE_EOS_TABLE 1024

  /* These are input parameters */
  int WhichEOS;
  int EOSFullTableSize;  /* Size of our final table, after interpolation */  

  /* These are set by our choice of EOS */
  int   EOSInitTableSize;  /* Size of our initial table, before interpolation */
  FLOAT MinTabulatedDensity;
  FLOAT MaxTabulatedDensity;
  FLOAT PolyIndexLowDensity;
  FLOAT PolyIndexHighDensity;
  FLOAT EOSDensDel;
  FLOAT EOSDensity[MAX_SIZE_EOS_TABLE];
  FLOAT EOSPressure[MAX_SIZE_EOS_TABLE];
  FLOAT EOSEnergy[MAX_SIZE_EOS_TABLE];
#endif

   /* SINK: add variables for sink particles */

  double RInner;
  double ROuter;
  double HSinkCreate;
  double RInner2;
  double ROuter2;
  double SinkCriticalDens;
  double SinkCriticalRedshift;
  int MaxNumSinks;

  int      accrete_historyTIMESTEPNUM[100000][100];
  double   accrete_historyTIME[100000][100];
  int      accrete_historyID[100000][100];
  FLOAT    accrete_historyMASS[100000][100];
  int      sink_number_global;

  /* The array structure is as follows: time of accretion event, accreted particle id #, accreting particle id #, mass of accreted particle, mass of accreting particle.  It has room for 10000 mergers. */
  double   merger_historyTIME[10000];
  FLOAT    merger_historyMASSGROW[10000];
  FLOAT    merger_historyMASSGONE[10000];
  int      merger_historyIDGROW[10000];
  int      merger_historyIDGONE[10000];
  int      merger_number;
  int      accID_arr[100];
  int      sink_num; 



#ifdef CHEMCOOL
#define GADGET2_MAX(a,b)  ((a > b) ? a : b)
#define GADGET2_MIN(a,b)  ((a < b) ? a : b)

  int NeedAbundancesForOutput;

  double H2RefDustEff;
  double OxyAbund;
  double CarbAbund;
  double SiAbund;
  double DeutAbund;
  double MgAbund;
  double UVField;
  double PhiPAH;
  double InitDustTemp;
  double DustToGasRatio;
  double AVConversionFactor;
  double CosmicRayIonRate;
  double InitRedshift;
  double ExternalDustExtinction;
  int PhotochemApprox;
  int ChemistryNetwork;
  int ADRateFlag;
  int MNRateFlag;
  int AtomicFlag;
  int ThreeBodyFlagA;
  int ThreeBodyFlagB;
  int H3PlusRateFlag;
  double InitMolHydroAbund;
  double InitHPlusAbund;
  double InitDIIAbund;
  double InitHDAbund;
  double InitHeIIAbund;
  double InitHeIIIAbund;
/* Following relevent only in metal-enriched gas */ 
  double InitCIIAbund;
  double InitSiIIAbund;
  double InitOIIAbund;
  double InitCOAbund;
  double InitC2Abund;
  double InitOHAbund;
  double InitH2OAbund;
  double InitO2Abund;
  double InitHCOPlusAbund;
  double InitCHAbund;
  double InitCH2Abund;
  double InitSiIIIAbund;
  double InitCH3PlusAbund;
  double InitMgPlusAbund;
#endif

#ifdef RAYTRACE
  double RTXmin[3];
  double RTdX[3];
#endif

  long long int Time_last_raytrace;

  int ray_flag_sun;
  double max_dens;
  double ray_crit_dens;
  double ray_r_max_sink;
  double RefinementMass;

  double Teff, Tint, Prad_avg, Pres_avg;
  double lum_tot, sinkmass_sum_tot;
  double heat_ion[7];
  double star_mass, star_rad, mdot;
  double x_s, n_s, n_sink, r_s, t_s, t_s0, t_s0_acc, t_s0_sink, alpha0;
  double tacc, t_KH, t_acc;
  double r0, m0, r1, m1, mdot1, r2, m2, e2, alpha;
  double ch2, ch, v_bulk[3];
  int numtot, stod, trans1, trans1a, trans2, set1,star_read, n_near_sum, n_near_max;
  int ray_center_ID, flag_sink;
}
 All;                                          /*!< a container variable for global variables that are equal on all processors */



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  FLOAT Pos[3];			/*!< particle position at its current time */
  FLOAT Mass;			/*!< particle mass */
  FLOAT Vel[3];			/*!< particle velocity at its current time */
  FLOAT GravAccel[3];		/*!< particle acceleration due to gravity */
#ifdef PMGRID
  FLOAT GravPM[3];		/*!< particle acceleration due to long-range PM gravity force*/
#endif
#ifdef FORCETEST
  FLOAT GravAccelDirect[3];	/*!< particle acceleration when computed with direct summation */
#endif
  FLOAT Potential;		/*!< gravitational potential */
  FLOAT OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening criterion */
#ifndef LONGIDS
  int ID;                       /*!< particle identifier - sink code requires that it be signed */
#else
  long long ID;        /*!< particle identifier */
#endif

  int Type;		        /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  long long int Ti_endstep;               /*!< marks start of current timestep of particle on integer timeline */ 
  long long int Ti_begstep;               /*!< marks end of current timestep of particle on integer timeline */
#ifdef FLEXSTEPS
  long long int FlexStepGrp;		/*!< a random 'offset' on the timeline to create a smooth groouping of particles */
#endif
  float GravCost;		/*!< weight factor used for balancing the work-load */
#ifdef PSEUDOSYMMETRIC
  float AphysOld;               /*!< magnitude of acceleration in last timestep. Used to make a first order
                                     prediction of the change of acceleration expected in the future, thereby
                                     allowing to guess whether a decrease/increase of the timestep should occur
                                     in the timestep that is started. */
#endif
}
 *P,              /*!< holds particle data on local processor */
 *DomainPartBuf;  /*!< buffer for particle data used in domain decomposition */


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
#ifndef POLYTROPE
  FLOAT Entropy;                /*!< current value of entropy (actually entropic function) of particle */
#endif
  FLOAT Density;		/*!< current baryonic mass density of particle */
  FLOAT Hsml;			/*!< current smoothing length */
  FLOAT Left;                   /*!< lower bound in iterative smoothing length search */  
  FLOAT Right;                  /*!< upper bound in iterative smoothing length search */ 
  FLOAT NumNgb;                 /*!< weighted number of neighbours found */
  FLOAT Pressure;		/*!< current pressure */
  FLOAT Prad;
  FLOAT Prad_dir[3];
  FLOAT ystep;
#ifndef POLYTROPE
  FLOAT DtEntropy;              /*!< rate of change of entropy */
#endif
  FLOAT HydroAccel[3];		/*!< acceleration due to hydrodynamical force */
  FLOAT VelPred[3];		/*!< predicted SPH particle velocity at the current time */
  FLOAT VelB[3];
  FLOAT DivVel;			/*!< local velocity divergence */
  FLOAT CurlVel;		/*!< local velocity curl */
  FLOAT Rot[3];		        /*!< local velocity curl */
  FLOAT DhsmlDensityFactor;     /*!< correction factor needed in the equation of motion of the conservative entropy formulation of SPH */
  FLOAT MaxSignalVel;           /*!< maximum "signal velocity" occuring for this particle */
  FLOAT Gamma;
#ifdef CHEMCOOL
  FLOAT DtEntropyVisc;          /*!< rate of change of entropy due to viscous dissipation only */
  //FLOAT Gamma;                  /*!< adiabatic index of simulated gas (current time) */
  FLOAT TracAbund[TRAC_NUM];    /*!< current abundances of tracer variables */
  FLOAT EntropyOut;             /*!< value of entropy that will be output in snapshot file */
  FLOAT TracAbundOut[TRAC_NUM]; /*!< value of abundances of tracer variables that will be output in snapshot file */
  FLOAT DustTemp;               /*!< Temperature of dust, in K */
  FLOAT HM;
  FLOAT H2II;
#ifdef RAYTRACE
  FLOAT TotalColumnDensity[6];  /*!< Total H nuclei column density in the +/- x,y,z directions */
  FLOAT H2ColumnDensity[6];     /*!< Molecular hydrogen column density in the +/- x,y,z directions */
  FLOAT COColumnDensity[6];     /*!< CO column density in the +/- x,y,z directions */
#endif /* RAYTRACE */
#endif /* CHEMCOOL */
#ifdef RAYTRACE_TG
  FLOAT Ray_H_coeff;
  FLOAT Ray_He_coeff;
  FLOAT Ray_LW_coeff;
  FLOAT Ray_NH2;
#endif
#ifdef METALS_TG
  FLOAT Sigma;
  FLOAT const_A;
  FLOAT const_B;
  FLOAT Metallicity;
#endif
#ifdef SINKVAL
  FLOAT sink;                   // ARS: needed to use VB's sink algorithm
#endif
#ifdef BFF
  FLOAT bfield[3], dBdt[3], dScidt;
#endif
  FLOAT jacob[3][3], jacob2[3][3], Bsmooth[3], ksum, DivB, Sci, GradSci[3], VelRel[3];
}
 *SphP,                        	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;                 /*!< buffer for SPH particle data in domain decomposition */





/*  Variables for Tree
 */

extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */

extern struct NODE
{
  FLOAT len;			/*!< sidelength of treenode */
  FLOAT center[3];		/*!< geometrical center of node */
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  FLOAT maxsoft;                /*!< hold the maximum gravitational softening of particles in the 
                                     node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
#endif
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      FLOAT s[3];               /*!< center of mass of node */
      FLOAT mass;               /*!< mass of node */
      int bitflags;             /*!< a bit-field with various information on the node */
      int sibling;              /*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;             /*!< this gives the next node in case the current node needs to be opened */
      int father;               /*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
}
 *Nodes_base,                   /*!< points to the actual memory allocted for the nodes */
 *Nodes;                        /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
 				     gives the first allocated node */


extern int *Nextnode;	        /*!< gives next node in tree walk */
extern int *Father;	        /*!< gives parent node in tree    */


extern struct extNODE           /*!< this structure holds additional tree-node information which is not needed in the actual gravity computation */
{
  FLOAT hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
  FLOAT vs[3];			/*!< center-of-mass velocity */
}
 *Extnodes_base,                /*!< points to the actual memory allocted for the extended node information */
 *Extnodes;                     /*!< provides shifted access to extended node information, parallel to Nodes/Nodes_base */





/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[numtype];                        /*!< number of particles of each type in this file */
  double mass[numtype];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[numtype];          /*!< total number of particles of each type in this snapshot. This can be
                                            different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  int num_files;                       /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];	               /*!< fills to 256 Bytes */
}
 header;                               /*!< holds header for snapshot files */

   /*!< total number of defined information blocks for snapshot files.
        Must be equal to the number of entries in "enum iofields" */


#ifdef METALS_TG
   #ifdef CHEMCOOL
      #ifdef SINKVAL
         #ifdef RAYTRACE
             #ifdef BFF
               #define IO_NBLOCKS 19
             #else
               #define IO_NBLOCKS 18
             #endif
         #else
             #ifdef BFF
               #define IO_NBLOCKS 18
             #else
               #define IO_NBLOCKS 17
             #endif
         #endif /* RAYTRACE */
      #else
         #ifdef RAYTRACE
            #define IO_NBLOCKS  17
         #else
            #define IO_NBLOCKS  16
         #endif /* RAYTRACE */
      #endif
   #else /* CHEMCOOL */
      #ifdef POLYTROPE
         #define IO_NBLOCKS 15
      #else
         #define IO_NBLOCKS 14
      #endif /* POLYTROPE */
   #endif /* CHEMCOOL */

#else /* METALS_TG */
   #ifdef CHEMCOOL
      #ifdef SINKVAL
         #ifdef RAYTRACE
             #ifdef BFF
               #define IO_NBLOCKS 18
             #else
               #define IO_NBLOCKS 17
             #endif
         #else
             #ifdef BFF
               #define IO_NBLOCKS 17
             #else
               #define IO_NBLOCKS 16
             #endif
            //#define IO_NBLOCKS  14
         #endif /* RAYTRACE */
      #else
         #ifdef RAYTRACE
            #define IO_NBLOCKS  16
            #else
            #define IO_NBLOCKS  15
         #endif /* RAYTRACE */
      #endif
   #else
      #ifdef POLYTROPE
         #define IO_NBLOCKS 14
      #else
         #define IO_NBLOCKS 13
      #endif /* POLYTROPE */
   #endif /* CHEMCOOL */
#endif /* METALS_TG */


enum iofields           /*!< this enumeration lists the defined output blocks in snapshot files. Not all of them need to be present. */
{ 
  IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_U,
  IO_RHO,
  IO_HSML,
  IO_POT,
  IO_ACCEL,
  IO_DTENTR,
  IO_TSTP,
#ifdef POLYTROPE
  IO_PRESSURE,
#endif

#ifdef CHEMCOOL
  IO_CHEM,
  IO_GAMMA,
  #ifdef RAYTRACE
    IO_COLUMN,
  #endif /* RAYTRACE */

  #ifdef SINKVAL
    IO_SINK,
  #endif

  #ifdef BFF
    IO_BFF,
  #endif

#endif /* CHEMCOOL */

#ifdef OUTPUTCURL
  IO_CURL,
#endif

#ifdef OUTPUTDIVB
  IO_DIVB,
#endif


#ifdef METALS_TG
  IO_METALLICITY,
#endif
};


extern char Tab_IO_Labels[IO_NBLOCKS][4];   /*<! This table holds four-byte character tags used for fileformat 2 */


/* global state of system, used for global statistics
 */
extern struct state_of_system
{
  double Mass;
#ifdef TURBULENCE
  double EnergyDrv;
#endif /* TURBULENCE */
  double EnergyKin;
  double EnergyPot;
  double EnergyInt;
  double EnergyTot;
  double Momentum[4];
  double AngMomentum[4];
  double CenterOfMass[4];
  double MassComp[6];
#ifdef TURBULENCE
  double EnergyDrvComp[6];
#endif /* TURBULENCE */
  double EnergyKinComp[6];
  double EnergyPotComp[6];
  double EnergyIntComp[6];
  double EnergyTotComp[6];
  double MomentumComp[6][4]; 
  double AngMomentumComp[6][4]; 
  double CenterOfMassComp[6][4];
#ifdef CHEMCOOL
  double MolAbund[TRAC_NUM];
#endif /* CHEMCOOL */
}
 SysState;                       /*<! Structure for storing some global statistics about the simulation. */
 


/* Various structures for communication
 */
extern struct gravdata_in
{
  union
  {
    FLOAT Pos[3];
    FLOAT Acc[3];
    FLOAT Potential;
  }
  u;
#ifdef UNEQUALSOFTENINGS
  int Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  FLOAT Soft;
#endif
#endif
  union
  {
    FLOAT OldAcc;
    int Ninteractions;
  }
  w;
}
 *GravDataIn,                   /*!< holds particle data to be exported to other processors */
 *GravDataGet,                  /*!< holds particle data imported from other processors */
 *GravDataResult,               /*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;                  /*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

extern struct gravdata_index
{
  int Task;
  int Index;
  int SortIndex;
}
 *GravDataIndexTable;           /*!< the particles to be exported are grouped by task-number. This table allows the results to be disentangled again and to be assigned to the correct particle */



extern struct densdata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT bfield[3];
  FLOAT Hsml;
  FLOAT Sci;
#ifdef METALS_TG
  FLOAT Sigma;
  FLOAT Density;
#endif
  int Index;
  int Task;
}
 *DensDataIn,                   /*!< holds particle data for SPH density computation to be exported to other processors */
 *DensDataGet;                  /*!< holds imported particle data for SPH density computation */

extern struct densdata_out
{
  FLOAT Rho;
  FLOAT Div, Div_B, Rot[3], GSci[3], VRel[3];
  FLOAT DhsmlDensity;
  FLOAT Ngb;
#ifdef METALS_TG
  FLOAT Sigma;
  FLOAT const_A;
  FLOAT const_B;
#endif
}
 *DensDataResult,               /*!< stores the locally computed SPH density results for imported particles */
 *DensDataPartialResult;        /*!< imported partial SPH density results from other processors */



extern struct hydrodata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  FLOAT Mass;
  FLOAT Density;
  FLOAT Pressure;
  //FLOAT Prad;
  //FLOAT Prad_dir[3];
  FLOAT F1;
  FLOAT DhsmlDensityFactor;
  long long int   Timestep;
  long long int time0, time1;
  int   Task;
  int   Index;
#ifdef CHEMCOOL
  FLOAT Gamma;
#endif /* CHEMCOOL */
}
 *HydroDataIn,                  /*!< holds particle data for SPH hydro-force computation to be exported to other processors */
 *HydroDataGet;                 /*!< holds imported particle data for SPH hydro-force computation */

extern struct hydrodata_out
{
  FLOAT Acc[3];
#ifndef POLYTROPE
  FLOAT DtEntropy;
#endif /* POLYTROPE */
  FLOAT MaxSignalVel;
  FLOAT jacob[3][3];
  FLOAT jacob2[3][3];
  FLOAT Bsmooth[3], ksum;
}
 *HydroDataResult,              /*!< stores the locally computed SPH hydro results for imported particles */
 *HydroDataPartialResult;       /*!< imported partial SPH hydro-force results from other processors */

/* Various structure for communication during the sink 
 * accretion
 */
extern struct sinkdata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT GravAccel[3];
#ifdef PMGRID
  FLOAT GravPM[3];
#endif
  FLOAT Mass;
  int   ID;
} 
*SinkDataIn;

extern struct sinkdata_out
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT GravAccel[3];
#ifdef PMGRID
  FLOAT GravPM[3];
#endif
  FLOAT Mass;
  int   ID;
} 
*SinkDataResult,
*SinkDataPartialResult;

extern struct ngbsdata_in
{
  FLOAT Pos[3];
  FLOAT Mass;
} 
*NgbsDataIn;

#ifdef CHEMCOOL
extern struct{
  double temptab[NMD]; 
  double cltab[NMD][NCLTAB]; 
  double chtab[NMD][NCHTAB]; 
  double dtcltab[NMD][NCLTAB]; 
  double dtchtab[NMD][NCHTAB]; 
  double crtab[NCRTAB];
  double crphot[NCRPHOT];
  double phtab[NPHTAB];
  double cst[NCONST];
  double dtlog;
  double tdust; 
  double tmax; 
  double tmin; 
  double deff; 
  double abundc; 
  double abundo;
  double abundsi;
  double abundD;
  double abundmg; 
  double G0;
  double f_rsc;
  double phi_pah;
  double dust_to_gas_ratio;
  double AV_conversion_factor;
  double cosmic_ray_ion_rate;
  double redshift;
  double AV_ext;
  double HM;
  double H2II;
#ifdef METALS_TG
  double abratio_c;
  double abratio_o;
  double abratio_si;
#endif
#ifdef RAYTRACE_TG
  double ray_H_coeff;
  double ray_He_coeff;
  double ray_LW_coeff;
  double ray_NH2;
#endif
  double heat_ion[7];
}COOLR;

extern struct{
  int iphoto;
  int iflag_mn;
  int iflag_ad;
  int iflag_atom;
  int iflag_3bh2a;
  int iflag_3bh2b;
  int iflag_h3pra;
  int iflag_h2opc;
  int iflag_highT;
  int iflag_fixed_ion;
  int ipart_id;
#ifdef RAYTRACE_TG
  int ray_flag_sun;
#endif
}COOLI;

#ifdef RAYTRACE
extern FLOAT density_H[NCOL][NCOL][NCOL];
extern FLOAT density_H2[NCOL][NCOL][NCOL];
extern FLOAT RTbuffer[NCOL][NCOL][NCOL];
extern FLOAT column_H[6][NCOL][NCOL][NCOL];
extern FLOAT column_H2[6][NCOL][NCOL][NCOL];
#ifdef CO_SHIELDING
extern FLOAT density_CO[NCOL][NCOL][NCOL];
extern FLOAT column_CO[6][NCOL][NCOL][NCOL];
#endif
#endif /* RAYTRACE */

#endif /* CHEMCOOL */

#endif

#ifdef RAYTRACE_TG

#define N_theta 200						/* the resolution in theta direction */
#define N_shells 200						/* the number of radial shells */

extern struct RAY{
  int N_rays;								/* total number of rays */
  int flag_start;								/* indicates start of ray-tracing */
  int flag_continue;							/* indicates continuation of ray-tracing */
  int flag_continue_start;						/* 1 for continuation of run */
  int flag_continue_end;						/* 1 for writing position of I-front */
  int flag_switch;							/* a dummy flag */
  double time_start;							/* start of ray-tracing (as scale factor) */
  double time_end;							/* end of ray-tracing (as scale factor) */
  double center_x_local;						/* local x center of ray-tracing sphere */
  double center_y_local;						/* local y center of ray-tracing sphere */
  double center_z_local;						/* local z center of ray-tracing sphere */
  double center_x_local2;                                                /* local x center of ray-tracing sphere */
  double center_y_local2;                                                /* local y center of ray-tracing sphere */
  double center_z_local2;     
  double r_min_local;						/* local minimum size of ray-tracing sphere */
  double center_x;							/* global x center of ray-tracing sphere */
  double center_y;							/* global y center of ray-tracing sphere */
  double center_z;							/* global z center of ray-tracing sphere */
  double center_x2;                                                      /* global x center of ray-tracing sphere */
  double center_y2;                                                      /* global y center of ray-tracing sphere */
  double center_z2;   
  double r_min;								/* global minimum size of ray-tracing sphere */
  double r_max;							/* maximum size of ray-tracing sphere */
  double alpha_B_H;							/* HII and HeII recombination rate */
  double alpha_B_He;						/* HeIII recombination rate */
  double HI_ion_rate;						/* HI ionization rate */
  double HeI_ion_rate;						/* HeI ionization rate */
  double HeII_ion_rate;						/* HeII ionization rate */
  double LW_rate;							/* H2 dissociation rate */
  double temp_HII_estimate;					/* temperature estimate of HII region */
  double temp_HeIII_estimate;					/* temperature estimate of HeIII region */
  double lifetime;							/* lifetime of star */
  double N_H_photons;						/* total number of HI ionizing photons */
  double N_He_photons;						/* total number of HeII ionizing photons */
  double Q_H_ion;
  double Q_He_ion;
  double Qion_zams;
  int N_phi[N_theta];							/* number of phi directions for given theta */
  double theta[N_theta];						/* theta directions */
  double phi[N_theta][2*N_theta];				/* phi directions */
  double H_front[N_theta][2*N_theta];				/* position of HII front */
  double He_front[N_theta][2*N_theta];			/* position of HeIII front */
  double r[N_shells];							/* radial shells */
  double  nh_local[N_theta][2*N_theta][N_shells]; 	/* local n_h */
  double  nh_tot[N_theta][2*N_theta][N_shells]; 		/* global n_h */
  double  count_local[N_theta][2*N_theta][N_shells];	/* local normalization */
  double  count_tot[N_theta][2*N_theta][N_shells];	/* global normalization */
  //double  ne_local[N_theta][2*N_theta][N_shells];        /* local n_e  (electron density)*/
  //double  ne_tot[N_theta][2*N_theta][N_shells];          /* global n_e (electron density)*/
  double  NH2_local[N_theta][2*N_theta][N_shells];        /*local H2 column density*/
  double  NH2_tot[N_theta][2*N_theta][N_shells];        /*global H2 column density*/
}
ray;
#endif

