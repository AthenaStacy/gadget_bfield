#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <errno.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "proto.h"



/*! \file io.c
 *  \brief Routines for producing a snapshot file on disk.
 */

static int n_type[numtype];
static long long ntot_type_all[numtype];




/*! This function writes a snapshot of the particle distribution to one or
 *  several files using the selected file format.  If NumFilesPerSnapshot>1,
 *  the snapshot is distributed onto several files, several of them can be
 *  written simultaneously (up to NumFilesWrittenInParallel). Each file
 *  contains data from a group of processors.
 */
void savepositions(int num)
{
  double t0, t1;
  char buf[500];
  int i, j, *temp, n, filenr, gr, ngroups, masterTask, lastTask;

  t0 = second();

  if(ThisTask == 0)
    printf("\nwriting snapshot file %d... \n", num);

#if defined(SFR) || defined(BLACK_HOLES)
  rearrange_particle_sequence();
  /* ensures that new tree will be constructed */
  All.NumForcesSinceLastDomainDecomp = 1 + All.TreeDomainUpdateFrequency * All.TotNumPart;
#endif

  if(NTask < All.NumFilesPerSnapshot)
    {
      if(ThisTask == 0)
	printf("Fatal error.\nNumber of processors must be larger or equal than All.NumFilesPerSnapshot.\n");
      endrun(0);
    }
  if(All.SnapFormat < 1 || All.SnapFormat > 3)
    {
      if(ThisTask == 0)
	printf("Unsupported File-Format\n");
      endrun(0);
    }
#ifndef  HAVE_HDF5
  if(All.SnapFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun(0);
    }
#endif


  /* determine global and local particle numbers */
  for(n = 0; n < numtype; n++)
    n_type[n] = 0;

  for(n = 0; n < NumPart; n++)
    n_type[P[n].Type]++;

  /* because ntot_type_all[] is of type `long long', we cannot do a simple
   * MPI_Allreduce() to sum the total particle numbers 
   */
  temp = malloc(NTask * numtype * sizeof(int));
  MPI_Allgather(n_type, numtype, MPI_INT, temp, numtype, MPI_INT, MPI_COMM_WORLD);
  for(i = 0; i < numtype; i++)
    {
      ntot_type_all[i] = 0;
      for(j = 0; j < NTask; j++)
	ntot_type_all[i] += temp[j * numtype + i];
    }
  free(temp);


  /* assign processors to output files */
  distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

  fill_Tab_IO_Labels();


//ARS: Go here to change number of digits in snapshot number if we have over 1000 files!

  if(All.NumFilesPerSnapshot > 1)
    sprintf(buf, "%s%s_%03d.%d", All.OutputDir, All.SnapshotFileBase, num, filenr);
  else
//    sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);

//  if(num > 999)
    sprintf(buf, "%s%s_%04d", All.OutputDir, All.SnapshotFileBase, num);
//      sprintf(buf, "%s%s_%05d", All.OutputDir, All.SnapshotFileBase, num);

  ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
  if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
    ngroups++;

  for(gr = 0; gr < ngroups; gr++)
    {
      if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	write_file(buf, masterTask, lastTask);
      MPI_Barrier(MPI_COMM_WORLD);
    }


  if(ThisTask == 0)
    printf("done with snapshot.\n");

  t1 = second();

  All.CPU_Snapshot += timediff(t0, t1);

}



/*! This function fills the write buffer with particle data. New output blocks
 *  can in principle be added here.
 */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
  int n, k, pindex;
  //float *fp;
  double *fp;

#ifdef POLYTROPE
  FLOAT pressure, energy;
#endif /* POLYTROPE */

#ifdef LONGIDS
  long long *ip;
#else
  int *ip;
#endif

#ifdef PERIODIC
  FLOAT boxSize;
#endif
#ifdef PMGRID
  double dt_gravkick_pm = 0;
#endif
  double dt_gravkick, dt_hydrokick, a3inv = 1, fac1, fac2;

#ifdef CHEMCOOL
  double gamm1;
#endif

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      fac1 = 1 / (All.Time * All.Time);
#ifdef CHEMCOOL
      fac2 = 1;
#else
       fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
#endif
    }
  else
    a3inv = fac1 = fac2 = 1;

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkick_pm =
      get_gravkick_factor(All.PM_Ti_begstep,
			  All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkick_pm = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;
#endif



  fp = CommBuffer;
  ip = CommBuffer;

  pindex = *startindex;

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Pos[k];
#ifdef PERIODIC
		boxSize = All.BoxSize;
#ifdef LONG_X
		if(k == 0)
		  boxSize = All.BoxSize * LONG_X;
#endif
#ifdef LONG_Y
		if(k == 1)
		  boxSize = All.BoxSize * LONG_Y;
#endif
#ifdef LONG_Z
		if(k == 2)
		  boxSize = All.BoxSize * LONG_Z;
#endif
		while(fp[k] < 0)
		  fp[k] += boxSize;
		while(fp[k] >= boxSize)
		  fp[k] -= boxSize;
#endif
	      }
	    n++;
	    fp += 3;
	  }
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    if(All.ComovingIntegrationOn)
	      {
		dt_gravkick =
		  get_gravkick_factor(P[pindex].Ti_begstep,
				      All.Ti_Current) -
		  get_gravkick_factor(P[pindex].Ti_begstep,
				      (P[pindex].Ti_begstep + P[pindex].Ti_endstep) / 2);
		dt_hydrokick =
		  get_hydrokick_factor(P[pindex].Ti_begstep,
				       All.Ti_Current) -
		  get_hydrokick_factor(P[pindex].Ti_begstep,
				       (P[pindex].Ti_begstep + P[pindex].Ti_endstep) / 2);
	      }
	    else
	      dt_gravkick = dt_hydrokick =
		(All.Ti_Current - (P[pindex].Ti_begstep + P[pindex].Ti_endstep) / 2) * All.Timebase_interval;

	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Vel[k] + P[pindex].GravAccel[k] * dt_gravkick;
		if(P[pindex].Type == 0)
		  fp[k] += SphP[pindex].HydroAccel[k] * dt_hydrokick;
	      }
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += P[pindex].GravPM[k] * dt_gravkick_pm;
#endif
	    for(k = 0; k < 3; k++)
	      fp[k] *= sqrt(a3inv);

	    n++;
	    fp += 3;
	  }
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *ip++ = P[pindex].ID;
	    n++;
	  }
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Mass;
	    n++;
	  }
      break;

    case IO_U:			/* internal energy */
#ifdef POLYTROPE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    energy = get_energy(SphP[pindex].Density);
	    *fp++ = energy;
	    n++;
	  }
#else /* POLYTROPE */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#ifdef ISOTHERM_EQS
	    *fp++ = SphP[pindex].Entropy;
#else /* ISOTHERM_EQS */
#ifdef CHEMCOOL
	    gamm1 = SphP[pindex].Gamma - 1.0;
	    *fp++ =
	      dmax(All.MinEgySpec,
		   SphP[pindex].EntropyOut / gamm1 * pow(SphP[pindex].Density * a3inv, gamm1));
#else /* CHEMCOOL */
	    *fp++ =
	      dmax(All.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].Density * a3inv, GAMMA_MINUS1));
#endif /* CHEMCOOL */
#endif /* ISOTHERM_EQS */
	    n++;
	  }
#endif /* POLYTROPE */
      break;
    case IO_RHO:		/* density */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Density;
	    n++;
	  }
      break;

    case IO_HSML:		/* SPH smoothing length */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Hsml;
	    n++;
	  }
      break;


    case IO_POT:		/* gravitational potential */
#ifdef OUTPUTPOTENTIAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Potential;
	    n++;
	  }
#endif
      break;

    case IO_ACCEL:		/* acceleration */
#ifdef OUTPUTACCELERATION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = fac1 * P[pindex].GravAccel[k];
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += fac1 * P[pindex].GravPM[k];
#endif

#ifdef CHEMCOOL
	    if (All.ComovingIntegrationOn) {
              fac2 = 1 / pow(All.Time, 3 * SphP[pindex].Gamma - 2);
	    }
#endif
	    if(P[pindex].Type == 0)
	      for(k = 0; k < 3; k++)
		fp[k] += fac2 * SphP[pindex].HydroAccel[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_DTENTR:		/* rate of change of entropy */
#ifdef OUTPUTCHANGEOFENTROPY
#ifndef POLYTROPE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DtEntropy;
	    n++;
	  }
#endif
#endif
      break;

    case IO_TSTP:		/* timestep  */
#ifdef OUTPUTTIMESTEP

      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (P[pindex].Ti_endstep - P[pindex].Ti_begstep) * All.Timebase_interval;
	    n++;
	  }
#endif
      break;

#ifdef POLYTROPE
    case IO_PRESSURE:   /* Gas pressure */
#ifdef OUTPUTPRESSURE
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    pressure = get_pressure(SphP[pindex].Density);
	    *fp++ = pressure;
	    n++;
	  }
#endif /* OUTPUTPRESSURE */
      break;
#endif /* POLYTROPE */

#ifdef CHEMCOOL
    case IO_CHEM:
      for(n=0; n<pc; pindex++) {
	if (P[pindex].Type == type) {
	  for(k=0; k<TRAC_NUM; k++) {
  	    *fp++ = SphP[pindex].TracAbundOut[k];
	  }
	  n++;
	}
      }
      break;

    case IO_GAMMA:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Gamma;
	    n++;
	  }
      break;
#endif

#ifdef RAYTRACE 
    case IO_COLUMN:
#ifdef OUTPUTCOLUMN
      /* Convert to cgs units for output */
      for(n=0;n<pc;pindex++) {
	if (P[pindex].Type == type) {
	  for (k=0; k<6; k++) {
	    *fp++ = SphP[pindex].TotalColumnDensity[k] / (All.UnitLength_in_cm * All.UnitLength_in_cm);
	  }
	  for (k=0; k<6; k++) {
	    *fp++ = SphP[pindex].H2ColumnDensity[k] / (All.UnitLength_in_cm * All.UnitLength_in_cm);
	  }
#ifdef CO_SHIELDING
	  for (k=0; k<6; k++) {
	    *fp++ = SphP[pindex].COColumnDensity[k] / (All.UnitLength_in_cm * All.UnitLength_in_cm);
	  }
#endif
	  n++;
        }
      }
#endif
      break;
#endif

#ifdef METALS_TG
    case IO_METALLICITY:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Metallicity;
	    n++;
	  }
      break;
#endif

#ifdef SINKVAL
    case IO_SINK:
      //printf("writing sink values! \n");
      for(n = 0; n < pc; pindex++)
        if(P[pindex].Type == type)
          {
            *fp++ = SphP[pindex].sink;
            n++;
          }
      break;
#endif

#ifdef BFF
    case IO_BFF:
      //printf("writing bfield values! \n");
      for(n = 0; n < pc; pindex++)
        if(P[pindex].Type == type)
          {
           for (k=0; k<3; k++) {
            *fp++ = SphP[pindex].bfield[k];
            }
            n++;
          }
      break;
#endif

    case IO_CURL:
#ifdef OUTPUTCURL
      for(n = 0; n < pc; pindex++)
        {
        //if(n < 10) printf("writing curl values!\n");
        if(P[pindex].Type == type)
          {
            *fp++ = SphP[pindex].CurlVel;
            n++;
          }
         }
#endif
      break;

    case IO_DIVB:
#ifdef OUTPUTDIVB
      for(n = 0; n < pc; pindex++)
        {
        //if(n < 10) printf("writing curl values!\n");
        if(P[pindex].Type == type)
          {
            *fp++ = SphP[pindex].DivB;
            n++;
          }
         }
#endif
      break;

    }

}




/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file. If one wants to add a new output-block, this
 *  function should be augmented accordingly.
 */
int get_bytes_per_blockelement(enum iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
#ifdef BFF
    case IO_BFF:
#endif
      bytes_per_blockelement = 3 * sizeof(double);
      break;

    case IO_ID:
#ifdef LONGIDS
      bytes_per_blockelement = sizeof(long long);
#else
      bytes_per_blockelement = sizeof(int);
#endif
      break;

    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
#ifdef POLYTROPE
    case IO_PRESSURE:
#endif
#ifdef CHEMCOOL
    case IO_GAMMA:
#endif
#ifdef METALS_TG
    case IO_METALLICITY:
#endif
#ifdef SINKVAL
    case IO_SINK:
#endif
    case IO_CURL:
    case IO_DIVB:
      bytes_per_blockelement = sizeof(double);
      break;
#ifdef CHEMCOOL
    case IO_CHEM:
      bytes_per_blockelement = TRAC_NUM * sizeof(double);
      break;
#endif
#ifdef RAYTRACE
    case IO_COLUMN:
#ifdef CO_SHIELDING
      bytes_per_blockelement = 18 * sizeof(float);
#else
      bytes_per_blockelement = 12 * sizeof(float);
#endif
      break;
#endif
    }

  return bytes_per_blockelement;
}


/*! This function returns the type of the data contained in a given block of
 *  the output file. If one wants to add a new output-block, this function
 *  should be augmented accordingly.
 */
int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case IO_ID:
#ifdef LONGIDS
      typekey = 2;		/* native long long */
#else
      typekey = 0;		/* native int */
#endif
      break;

    default:
      typekey = 1;		/* native float */
      break;
    }

  return typekey;
}


/*! This function informs about the number of elements stored per particle for
 *  the given block of the output file. If one wants to add a new
 *  output-block, this function should be augmented accordingly.
 */
int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
#ifdef BFF
    case IO_BFF:
#endif
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
#ifdef POLYTROPE
    case IO_PRESSURE:
#endif
#ifdef CHEMCOOL
    case IO_GAMMA:
#endif
#ifdef METALS_TG
    case IO_METALLICITY:
#endif
#ifdef SINKVAL
    case IO_SINK:
#endif
    case IO_CURL:
    case IO_DIVB:
      values = 1;
      break;
#ifdef CHEMCOOL
    case IO_CHEM:
      values = TRAC_NUM;
      break;
#endif

#ifdef RAYTRACE
    case IO_COLUMN:
#ifdef CO_SHIELDING
      values = 18;
#else
      values = 12;
#endif
      break;
#endif
    }

  return values;
}


/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array. If one wants to
 *  add a new output-block, this function should be augmented accordingly.
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses, ngas, nstars;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < numtype; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  nstars = header.npart[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:

      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < numtype; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_DTENTR:
#ifdef POLYTROPE
    case IO_PRESSURE:
#endif
#ifdef CHEMCOOL
    case IO_CHEM:
    case IO_GAMMA:
#endif
#ifdef METALS_TG
    case IO_METALLICITY:
#endif
#ifdef SINKVAL
     case IO_SINK:
#endif
#ifdef BFF
     case IO_BFF:
#endif
#ifdef RAYTRACE
    case IO_COLUMN:
#endif
    case IO_CURL:
    case IO_DIVB:
      for(i = 1; i < numtype; i++)
	typelist[i] = 0;
      return ngas;
      break;
    }

  endrun(212);
  return 0;
}



/*! This function tells whether or not a given block in the output file is
 *  present, depending on the type of simulation run and the compile-time
 *  options. If one wants to add a new output-block, this function should be
 *  augmented accordingly.
 */
int blockpresent(enum iofields blocknr)
{

#ifndef OUTPUTPOTENTIAL
  if(blocknr == IO_POT)
    return 0;
#endif

#ifndef OUTPUTACCELERATION
  if(blocknr == IO_ACCEL)
    return 0;
#endif

#ifndef OUTPUTCHANGEOFENTROPY
  if(blocknr == IO_DTENTR)
    return 0;
#endif

#ifndef OUTPUTTIMESTEP
  if(blocknr == IO_TSTP)
    return 0;
#endif

#ifdef POLYTROPE
#ifndef OUTPUTPRESSURE
  if(blocknr == IO_PRESSURE)
    return 0;
#endif
#endif

#ifdef RAYTRACE
#ifndef OUTPUTCOLUMN
  if (blocknr == IO_COLUMN)
    return 0;
#endif
#endif

#ifndef OUTPUTCURL
  if(blocknr == IO_CURL)
    return 0;
#endif

#ifndef OUTPUTDIVB
  if(blocknr == IO_DIVB)
    return 0;
#endif

  return 1;			/* default: present */
}




/*! This function associates a short 4-character block name with each block
 *  number.  This is stored in front of each block for snapshot
 *  FileFormat=2. If one wants to add a new output-block, this function should
 *  be augmented accordingly.
 */
void fill_Tab_IO_Labels(void)
{
  enum iofields i;

  for(i = 0; i < IO_NBLOCKS; i++)
    switch (i)
      {
      case IO_POS:
	strncpy(Tab_IO_Labels[IO_POS], "POS ", 4);
	break;
      case IO_VEL:
	strncpy(Tab_IO_Labels[IO_VEL], "VEL ", 4);
	break;
      case IO_ID:
	strncpy(Tab_IO_Labels[IO_ID], "ID  ", 4);
	break;
      case IO_MASS:
	strncpy(Tab_IO_Labels[IO_MASS], "MASS", 4);
	break;
      case IO_U:
	strncpy(Tab_IO_Labels[IO_U], "U   ", 4);
	break;
      case IO_RHO:
	strncpy(Tab_IO_Labels[IO_RHO], "RHO ", 4);
	break;
      case IO_HSML:
	strncpy(Tab_IO_Labels[IO_HSML], "HSML", 4);
	break;
      case IO_POT:
	strncpy(Tab_IO_Labels[IO_POT], "POT ", 4);
	break;
      case IO_ACCEL:
	strncpy(Tab_IO_Labels[IO_ACCEL], "ACCE", 4);
	break;
      case IO_DTENTR:
	strncpy(Tab_IO_Labels[IO_DTENTR], "ENDT", 4);
	break;
      case IO_TSTP:
	strncpy(Tab_IO_Labels[IO_TSTP], "TSTP", 4);
	break;
#ifdef POLYTROPE
      case IO_PRESSURE:
	strncpy(Tab_IO_Labels[IO_PRESSURE], "PRES", 4);
	break;
#endif
#ifdef CHEMCOOL
      case IO_CHEM:
	strncpy(Tab_IO_Labels[IO_CHEM], "CHEM", 4);
	break;
      case IO_GAMMA:
	strncpy(Tab_IO_Labels[IO_GAMMA], "GAMM", 4);
	break;
#endif
#ifdef RAYTRACE
      case IO_COLUMN:
	strncpy(Tab_IO_Labels[IO_COLUMN], "COLN", 4);
	break;
#endif
#ifdef METALS_TG
      case IO_METALLICITY:
	strncpy(Tab_IO_Labels[IO_METALLICITY], "METL", 4);
	break;
#endif
#ifdef SINKVAL
      case IO_SINK:
        strncpy(Tab_IO_Labels[IO_SINK], "SINK", 4);
        break;
#endif
#ifdef BFF
      case IO_BFF:
        strncpy(Tab_IO_Labels[IO_BFF], "BFFF", 4);
        break;
#endif
#ifdef OUTPUTCURL
      case IO_CURL:
        strncpy(Tab_IO_Labels[IO_CURL], "CURL", 4);
        break;
#endif
#ifdef OUTPUTDIVB
      case IO_DIVB:
        strncpy(Tab_IO_Labels[IO_DIVB], "DIVB", 4);
        break;
#endif
      }
}

/*! This function returns a descriptive character string that describes the
 *  name of the block when the HDF5 file format is used.  If one wants to add
 *  a new output-block, this function should be augmented accordingly.
 *
 *  N.B. Names can be no more than 500 characters long! 
 */
void get_dataset_name(enum iofields blocknr, char *buf)
{

  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;
    case IO_HSML:
      strcpy(buf, "SmoothingLength");
      break;
    case IO_POT:
      strcpy(buf, "Potential");
      break;
    case IO_ACCEL:
      strcpy(buf, "Acceleration");
      break;
    case IO_DTENTR:
      strcpy(buf, "RateOfChangeOfEntropy");
      break;
    case IO_TSTP:
      strcpy(buf, "TimeStep");
      break;
#ifdef POLYTROPE
    case IO_PRESSURE:
      strcpy(buf, "Pressure");
      break;
#endif
#ifdef CHEMCOOL
    case IO_CHEM:
      strcpy(buf, "ChemicalAbundances");
      break;
    case IO_GAMMA:
      strcpy(buf, "Adiabatic index");
      break;
#endif
#ifdef RAYTRACE
    case IO_COLUMN:
      strcpy(buf, "Column densities");
      break;
#endif
#ifdef METALS_TG
    case IO_METALLICITY:
      strcpy(buf, "Metallicity");
      break;
#endif
#ifdef SINKVAL
    case IO_SINK:
      strcpy(buf, "SinkValue");
      break;
#endif
#ifdef BFF
    case IO_BFF:
      strcpy(buf, "Bfield");
      break;
#endif
#ifdef OUTPUTCURL
    case IO_CURL:
      strcpy(buf, "Curl");
      break;
#endif
#ifdef OUTPUTDIVB
    case IO_DIVB:
      strcpy(buf, "DivB");
      break;
#endif
    }
}



/*! This function writes an actual snapshot file containing the data from
 *  processors 'writeTask' to 'lastTask'. 'writeTask' is the one that actually
 *  writes.  Each snapshot file contains a header first, then particle
 *  positions, velocities and ID's.  Particle masses are written only for
 *  those particle types with zero entry in MassTable.  After that, first the
 *  internal energies u, and then the density is written for the SPH
 *  particles.  If cooling is enabled, mean molecular weight and neutral
 *  hydrogen abundance are written for the gas particles. This is followed by
 *  the SPH smoothing length and further blocks of information, depending on
 *  included physics and compile-time flags.  If HDF5 is used, the header is
 *  stored in a group called "/Header", and the particle data is stored
 *  separately for each particle type in groups called "/PartType0",
 *  "/PartType1", etc. The sequence of the blocks is unimportant in this case.
 */
void write_file(char *fname, int writeTask, int lastTask)
{
  int type, bytes_per_blockelement, npart, nextblock, typelist[numtype];
  int n_for_this_task, ntask, n, p, pc, offset = 0, task;
  int blockmaxlen, ntot_type[numtype], nn[numtype];
  enum iofields blocknr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank, pcsum = 0;
  char buf[500];
#endif

#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

  /* determine particle numbers of each type in file */

  if(ThisTask == writeTask)
    {
      for(n = 0; n < numtype; n++)
	ntot_type[n] = n_type[n];

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&nn[0], numtype, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
	  for(n = 0; n < numtype; n++)
	    ntot_type[n] += nn[n];
	}

      for(task = writeTask + 1; task <= lastTask; task++)
	MPI_Send(&ntot_type[0], numtype, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], numtype, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&ntot_type[0], numtype, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }



  /* fill file header */

  for(n = 0; n < numtype; n++)
    {
      header.npart[n] = ntot_type[n];
      header.npartTotal[n] = (unsigned int) ntot_type_all[n];
      header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }

  for(n = 0; n < numtype; n++)
    header.mass[n] = All.MassTable[n];

  header.time = All.Time;

  if(All.ComovingIntegrationOn)
    header.redshift = 1.0 / All.Time - 1;
  else
    header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;

#ifdef COOLING
  header.flag_cooling = 1;
#endif
#ifdef SFR
  header.flag_sfr = 1;
  header.flag_feedback = 1;
#ifdef STELLARAGE
  header.flag_stellarage = 1;
#endif
#ifdef METALS
  header.flag_metals = 1;
#endif
#endif

  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;


  /* open file and write header */

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  sprintf(buf, "%s.hdf5", fname);
	  hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

	  for(type = 0; type < numtype; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
		}
	    }

	  write_header_attributes_in_hdf5(hdf5_headergrp);
#endif
	}
      else
	{
	  if(!(fd = fopen(fname, "w")))
	    {
	      printf("can't open file `%s' for writing snapshot.\n", fname);
	      endrun(123);
	    }

	  if(All.SnapFormat == 2)
	    {
	      blksize = sizeof(int) + 4 * sizeof(char);
	      SKIP;
	      my_fwrite("HEAD", sizeof(char), 4, fd);
	      nextblock = sizeof(header) + 2 * sizeof(int);
	      my_fwrite(&nextblock, sizeof(int), 1, fd);
	      SKIP;
	    }

	  blksize = sizeof(header);
	  SKIP;
	  my_fwrite(&header, sizeof(header), 1, fd);
	  SKIP;
	}
    }

  ntask = lastTask - writeTask + 1;

  for(blocknr = 0; blocknr < IO_NBLOCKS; blocknr++)
    {
      if(blockpresent(blocknr))
	{
	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr);

	  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == writeTask)
		{

		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    {
		      if(All.SnapFormat == 2)
			{
			  blksize = sizeof(int) + 4 * sizeof(char);
			  SKIP;
			  my_fwrite(Tab_IO_Labels[blocknr], sizeof(char), 4, fd);
			  nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
			  my_fwrite(&nextblock, sizeof(int), 1, fd);
			  SKIP;
			}

		      blksize = npart * bytes_per_blockelement;
		      SKIP;

		    }
		}

	      for(type = 0; type < numtype; type++)
		{
		  if(typelist[type])
		    {
#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  switch (get_datatype_in_block(blocknr))
			    {
			    case 0:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
			      break;
			    case 1:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
			      break;
			    case 2:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
			      break;
			    }

			  dims[0] = header.npart[type];
			  dims[1] = get_values_per_blockelement(blocknr);
			  if(dims[1] == 1)
			    rank = 1;
			  else
			    rank = 2;

			  get_dataset_name(blocknr, buf);

			  hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
			  hdf5_dataset =
			    H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file,
				      H5P_DEFAULT);
			  pcsum = 0;
			}
#endif

		      for(task = writeTask, offset = 0; task <= lastTask; task++)
			{
			  if(task == ThisTask)
			    {
			      n_for_this_task = n_type[type];

			      for(p = writeTask; p <= lastTask; p++)
				if(p != ThisTask)
				  MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
			    }
			  else
			    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD,
				     &status);

			  while(n_for_this_task > 0)
			    {
			      pc = n_for_this_task;

			      if(pc > blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == task)
				fill_write_buffer(blocknr, &offset, pc, type);

			      if(ThisTask == writeTask && task != writeTask)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
					 TAG_PDATA, MPI_COMM_WORLD, &status);

			      if(ThisTask != writeTask && task == ThisTask)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask,
					  TAG_PDATA, MPI_COMM_WORLD);

			      if(ThisTask == writeTask)
				{
				  if(All.SnapFormat == 3)
				    {
#ifdef HAVE_HDF5
				      start[0] = pcsum;
				      start[1] = 0;

				      count[0] = pc;
				      count[1] = get_values_per_blockelement(blocknr);
				      pcsum += pc;

				      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							  start, NULL, count, NULL);

				      dims[0] = pc;
				      dims[1] = get_values_per_blockelement(blocknr);
				      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

				      hdf5_status =
					H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
						 hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

				      H5Sclose(hdf5_dataspace_memory);
#endif
				    }
				  else
				    my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
				}

			      n_for_this_task -= pc;
			    }
			}

#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && All.SnapFormat == 3 && header.npart[type] > 0)
			{
			  if(All.SnapFormat == 3)
			    {
			      H5Dclose(hdf5_dataset);
			      H5Sclose(hdf5_dataspace_in_file);
			      H5Tclose(hdf5_datatype);
			    }
			}
#endif
		    }
		}

	      if(ThisTask == writeTask)
		{
		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    SKIP;
		}
	    }
	}
    }

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      H5Gclose(hdf5_grp[type]);
	  H5Gclose(hdf5_headergrp);
	  H5Fclose(hdf5_file);
#endif
	}
      else
	fclose(fd);
    }
}




/*! This function writes the header information in case HDF5 is selected as
 *  file format.
 */
#ifdef HAVE_HDF5
void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total_HW", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);


  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  header.flag_entropy_instead_u = 0;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "Flag_Entropy_ICs", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_entropy_instead_u);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}
#endif





/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun(777);
    }
  return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun(778);
    }
  return nread;
}