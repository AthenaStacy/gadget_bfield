#ifdef TURBULENCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

void rsk_turbdriving()
{

  int i, j, k, ix, iy, iz, ixmin, ixmax, iymin, iymax, izmin, izmax;
  int nactive, nadvanced, ncritdrv, nactivesum, nadvancedsum;
  int index;
  double factor, aa, bb, cc, ampl, velrms, soundspeed;
  double *massf, *xvelf, *yvelf, *zvelf, *rubbish1, *rubbish2;
  double totalbef, tkinbef, tgravbef, ttermbef, angtobef;
  double totalaft, tkinaft, tgravaft, ttermaft, angtoaft;
  double totaldif, tkindif, tgravdif, ttermdif, angtodif;
  double tstart, tend, bbsum, velrmssum, aasum;
  double xmomsum,ymomsum,zmomsum,xmom, ymom, zmom,masstot, masssum;  
  double velrmsx, velrmssumx, velrmsy, velrmssumy, velrmsz, velrmssumz;
#ifdef CHEMCOOL
  double abe;
#endif 

  if(!(massf = malloc(sizeof(double)*IFIELDSIZE*IFIELDSIZE*IFIELDSIZE)))
  {
    printf("failed to allocate memory for massf.\n");
    endrun(2);
  }  
  
  if(!(xvelf = malloc(sizeof(double)*IFIELDSIZE*IFIELDSIZE*IFIELDSIZE)))
  {
    printf("failed to allocate memory for xvelf.\n");
    endrun(2);
  } 
   
  if(!(yvelf = malloc(sizeof(double)*IFIELDSIZE*IFIELDSIZE*IFIELDSIZE)))
  {
    printf("failed to allocate memory for xvelf.\n");
    endrun(2);
  } 

  if(!(zvelf = malloc(sizeof(double)*IFIELDSIZE*IFIELDSIZE*IFIELDSIZE)))
  {
    printf("failed to allocate memory for xvelf.\n");
    endrun(2);
  } 
  
  if(!(rubbish1 = malloc(sizeof(double)*IFIELDSIZE*IFIELDSIZE*IFIELDSIZE)))
  {
    printf("failed to allocate memory for rubbish1.\n");
    endrun(2);
  } 
   
  if(!(rubbish2 = malloc(sizeof(double)*IFIELDSIZE*IFIELDSIZE*IFIELDSIZE)))
  {
    printf("failed to allocate memory for rubbish2.\n");
    endrun(2);
  } 
  
  nactive = 0;
  /*  nadvanced = 0; */
  /*  ncritdrv = 0; */
  nactivesum = 0;
  /*  nadvancedsum = 0; */
  ixmin = iymin = izmin = IFIELDSIZE-1;
  ixmax = iymax = izmax = 0;

  for(i=0;i<NumPart;i++)
    {
      if((P[i].Type == 0) && (P[i].ID > 0))
	{
	  nactive++;
	}
	  /*	  if(P[i].ForceFlag != 0)
	   *        nadvanced++;
	   *      }
	   */
    }

  MPI_Allreduce(&nactive,&nactivesum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
/*  MPI_Allreduce(&nadvanced,&nadvancedsum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);*/
  /*  ncritdrv = (int)(nactivesum/40.0); */
  All.DeltaTimeDrv += All.TimeStep;

 /*  if(ThisTask==0) */
/*     { */
  /*    printf("\n%d :Driving:%g %g %d %d %g %g\n",ThisTask,All.Time,All.StartDriving,nadvanced,ncritdrv,All.DeltaTimeDrv,TDRVSTEPFACTOR*All.MaxSizeTimestep); */
/*      fflush(stdout); */
   /*  } */
   /*   if(All.Time > All.StartDriving && ((nadvancedsum > ncritdrv) || ((All.DeltaTimeDrv-All.StartDriving) > (TDRVSTEPFACTOR*All.MaxSizeTimestep)))) */
if(All.Time > All.StartDriving && ((All.DeltaTimeDrv-All.StartDriving) > (All.DrvTimestep)))
    {
#ifdef TURBULENCE_RANDOM
          icomp = 1;
	 /* iiseed = All.Seed0;*/
	  rsk_turbdriving_field();
	  icomp = 2;
	  /* iiseed = All.Seed1;*/
	  rsk_turbdriving_field();
	  icomp = 3;
	  /* iiseed = All.Seed2;*/
	  rsk_turbdriving_field();
	  
	  if(ThisTask==0)
	    {
	      printf("initialization of turbulence field finished.\n");
	      fflush(stdout);
	    }
#endif 
      All.DeltaTimeDrv -= All.StartDriving;
      compute_global_quantities_of_system(); /*maybe change to 2?*/
       
      totalbef = SysState.EnergyTot;
      tkinbef  = SysState.EnergyKin;
      tgravbef = SysState.EnergyPot;
      ttermbef = SysState.EnergyInt;
      angtobef = SysState.AngMomentum[3];
       
      factor = 1.0/All.BoxSize*IFIELDSIZE;

      for(k=0;k<IFIELDSIZE;k++)
	{
	  for(j=0;j<IFIELDSIZE;j++)
	    {
	      for(i=0;i<IFIELDSIZE;i++)
		{
                  index = i*IFIELDSIZE*IFIELDSIZE+j*IFIELDSIZE+k;
		  massf[index] = 0;
		  xvelf[index] = 0;
		  yvelf[index] = 0;
		  zvelf[index] = 0;
		}
	    }
	}
         
      for(i=0;i<NumPart;i++)
	{
	   if((P[i].Type == 0)&& (P[i].ID > 0)) 
 	    { 
	      ix = (int)(P[i].Pos[0]*factor);
	      iy = (int)(P[i].Pos[1]*factor);
	      iz = (int)(P[i].Pos[2]*factor);
    
  
	      if (ix < 0 || iy < 0 || iz < 0 || ix >= IFIELDSIZE || iy >= IFIELDSIZE || iz >= IFIELDSIZE)
		{
		  if(ThisTask==0)
		    {
		      printf("margins exceeded for 1. assignment: %d --- %d %d %d %g %g %g\n", i,ix,iy,iz,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
		      fflush(stdout);
		    }
		  if(ix < 0)
		    ix = 0;
		  if(iy < 0)
		    iy = 0;
		  if(iz < 0)
		    iz = 0;
		  if(ix >= IFIELDSIZE)
		    ix = IFIELDSIZE-1;
		  if(iy >= IFIELDSIZE)
		    iy = IFIELDSIZE-1;
		  if(iz >= IFIELDSIZE)
		    iz = IFIELDSIZE-1;
		}

	      if(ix < ixmin)
		ixmin = ix;
	      if(iy < iymin)
		iymin = iy;
	      if(iz < izmin)
		izmin = iz;
	      if(ix > ixmax)
		ixmax = ix;
	      if(iy > iymax)
		iymax = iy;
	      if(iz > izmax)
		izmax = iz;

              index = ix*IFIELDSIZE*IFIELDSIZE+iy*IFIELDSIZE+iz;

	      massf[index] += P[i].Mass;
	      xvelf[index] += P[i].Mass *P[i].Vel[0];
	      yvelf[index] += P[i].Mass *P[i].Vel[1];
	      zvelf[index] += P[i].Mass *P[i].Vel[2];

	    } 
	}

/* substract mean momentum*/
      
      xmomsum = 0.0;
      ymomsum = 0.0;
      zmomsum = 0.0;
      xmom = 0.0;
      ymom = 0.0;
      zmom = 0.0;
      masstot = 0.0;
      masssum = 0.0;

      for(k=izmin;k<=izmax;k++)
	{
	  for(j=iymin;j<=iymax;j++)
	    {
	      for(i=ixmin;i<=ixmax;i++)
		{
                  index = i*IFIELDSIZE*IFIELDSIZE+j*IFIELDSIZE+k;
		  xmom += massf[index] * xmatrix[index];
		  ymom += massf[index] * ymatrix[index];
		  zmom += massf[index] * zmatrix[index];
		  masstot += massf[index];		  
		}
	    }
	}

      
      MPI_Allreduce(&xmom, &xmomsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&ymom, &ymomsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&zmom, &zmomsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&masstot, &masssum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
      for(k=izmin;k<=izmax;k++)
	{
	  for(j=iymin;j<=iymax;j++)
	    {
	      for(i=ixmin;i<=ixmax;i++)
		{
                  index = i*IFIELDSIZE*IFIELDSIZE+j*IFIELDSIZE+k;
		  xmatrix[index] -= xmomsum/masssum;
		  ymatrix[index] -= ymomsum/masssum;
		  zmatrix[index] -= zmomsum/masssum;
		}
	    }
	}

      /*end of substracted momentum */


      for(k=izmin;k<=izmax;k++)
	{
	  for(j=iymin;j<=iymax;j++)
	    {
	      for(i=ixmin;i<=ixmax;i++)
		{
                  index = i*IFIELDSIZE*IFIELDSIZE+j*IFIELDSIZE+k;
		  rubbish1[index] = massf[index] * 
                                   (xmatrix[index]*xmatrix[index] + 
                                    ymatrix[index]*ymatrix[index] + 
                                    zmatrix[index]*zmatrix[index]);

		  rubbish2[index] = xmatrix[index]*xvelf[index] + 
                                    ymatrix[index]*yvelf[index] + 
                                    zmatrix[index]*zvelf[index]; 
		}
	    }
	}

      aa = 0.0;
      bb = 0.0;
      velrms = 0.0;
      bbsum = 0.0;
      aasum = 0.0;
      velrmssum = 0.0;

      for(k=izmin;k<=izmax;k++)
	{
	  for(j=iymin;j<=iymax;j++)
	    {
	      for(i=ixmin;i<=ixmax;i++)
		{
                  index = i*IFIELDSIZE*IFIELDSIZE+j*IFIELDSIZE+k;
		  aa += rubbish1[index];
		  bb += rubbish2[index];
		  if(isnan(aa) || isnan(bb) || isinf(aa) || isinf(bb))
		    printf("NANs:%d %d %d\n",i,j,k);
		 }
	    }
	}
          
      MPI_Allreduce(&bb, &bbsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&aa, &aasum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      bb=bbsum;
      aa=aasum;

      if (All.VelDispConst == 1)
	{
	  aa *= 0.5;
	  bb *= 1.0;
	  cc = -All.LDrv*(((double)nactivesum)/(int)All.TotN_gas)*All.DeltaTimeDrv ;
	  printf("new driving: %d %d\n",nactivesum,(int)All.TotN_gas);
	} 
      else
	{
	  aa *= 0.5;
	  bb *= 1.0;
	  cc = -All.LDrv*All.DeltaTimeDrv ;  
	}

      ampl = (-bb + sqrt(bb*bb - 4.0*aa*cc))/ (2.0*aa);
    /*   if(ThisTask==0) */
/* 	{ */
	 /*   printf("Amplitude:%d %g %g %g %g %d \n",ThisTask,ampl,bb,aa,cc,N_gas); */
/*  	  fflush(stdout); */
/* 	} */
      /* ------------------------------------------------------------------*/
      velrms = 0;
      velrmssum = 0;
      for(i=0;i<NumPart;i++)
	{
	  if((P[i].Type == 0) && (P[i].ID > 0))
	    {  
	      ix = (int)(P[i].Pos[0]*factor);
	      iy = (int)(P[i].Pos[1]*factor);
	      iz = (int)(P[i].Pos[2]*factor);
    
  
	      if (ix < 0 || iy < 0 || iz < 0 || ix >= IFIELDSIZE || iy >= IFIELDSIZE || iz >= IFIELDSIZE)
		{
		  if(ThisTask==0)
		    {
		      printf("margins exceeded for back assignment: %d --- %d %d %d %g %g %g\n", i,ix,iy,iz,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(ix < 0)
		    ix = 0;
		  if(iy < 0)
		    iy = 0;
		  if(iz < 0)
		    iz = 0;
		  if(ix >= IFIELDSIZE)
		    ix = IFIELDSIZE-1;
		  if(iy >= IFIELDSIZE)
		    iy = IFIELDSIZE-1;
		  if(iz >= IFIELDSIZE)
		    iz = IFIELDSIZE-1;
		}

              index = ix*IFIELDSIZE*IFIELDSIZE+iy*IFIELDSIZE+iz;

	      P[i].Vel[0] += ampl * xmatrix[index];
	      P[i].Vel[1] += ampl * ymatrix[index];
	      P[i].Vel[2] += ampl * zmatrix[index];
	      velrms += P[i].Vel[0]*P[i].Vel[0] 
                      + P[i].Vel[1]*P[i].Vel[1] 
                      + P[i].Vel[2]*P[i].Vel[2];
	    }
	}
       MPI_Allreduce(&velrms, &velrmssum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
       velrms = velrmssum;
       velrms = sqrt(velrms/nactivesum);
       if(ThisTask==0)
	{
#ifdef CHEMCOOL
	  /* XXX: This is WRONG -- need to use the current sound-speed, but
                  this varies throughtout the box; the right thing to do may
                  be to use the mean soundspeed, but I'm not sure it's 
                  really worth it */
          abe = initial_electron_fraction();
          soundspeed = (sqrt(GAS_CONST * All.InitGasTemp * 
                       (1.0 + ABHE - All.InitMolHydroAbund + abe) / 
                       (1.0 + 4.0 * ABHE))) / All.UnitVelocity_in_cm_per_s;
#else
          soundspeed = (sqrt(GAS_CONST*All.InitGasTemp/All.MWeight)) / 
                       All.UnitVelocity_in_cm_per_s;
#endif
	  printf("\n%d:%g Mach:%g\n",ThisTask,velrms,velrms/soundspeed);
	  fflush(stdout);
 	} 

      compute_global_quantities_of_system();

      totalaft = SysState.EnergyTot;
      tkinaft  = SysState.EnergyKin;
      tgravaft = SysState.EnergyPot;
      ttermaft = SysState.EnergyInt;
      angtoaft = SysState.AngMomentum[3];

      totaldif = totalaft - totalbef;
      tkindif  = tkinaft - tkinbef;
      tgravdif = tgravaft - tgravbef;
      ttermdif = ttermaft - ttermbef;
      angtodif = angtoaft - angtobef;


      SysState.EnergyDrvComp[0] += totaldif;
      for(i=0;i<6;i++)
	SysState.EnergyDrv += SysState.EnergyDrvComp[i];
      if(ThisTask==0)
	{
	  printf("Change due to driving:\ntotal energy:\t %g\nkinetic energy:\t %g\npotential energy:\t %g\ninternal energy:\t %g\ntotal angular momentum:\t %g\ndriving energy:\t %g\n\n",totaldif,tkindif,tgravdif,ttermdif,angtodif,SysState.EnergyDrvComp[0]);
	  fflush(stdout);
	}
  
      All.DeltaTimeDrv = All.StartDriving;
    }

  free(massf);
  free(xvelf);
  free(yvelf);
  free(zvelf);
  free(rubbish1);
  free(rubbish2);
}





#endif
