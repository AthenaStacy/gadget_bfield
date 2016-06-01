#ifdef TURBULENCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(double data[], int nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software 0(9p#3D0>4^c#>Z41.5.. */


/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */
#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define EPS (1.2E-07)
#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a<b)?a:b

double ran1(int *idum)
{
        int j,k;
        static int iv[NTAB],iy=0;
        void nrerror();
        static double NDIV = 1.0/(1.0+(IM-1.0)/NTAB);
        static double RNMX = (1.0-EPS);
        static double AM = (1.0/IM);

        if ((*idum <= 0) || (iy == 0)) {
                *idum = MAX(-*idum,*idum);
                for(j=NTAB+7;j>=0;j--) {
                        k = *idum/IQ;
                        *idum = IA*(*idum-k*IQ)-IR*k;
                        if(*idum < 0) *idum += IM;
                        if(j < NTAB) iv[j] = *idum;
                }
                iy = iv[0];
        }
        k = *idum/IQ;
        *idum = IA*(*idum-k*IQ)-IR*k;
        if(*idum<0) *idum += IM;
        j = iy*NDIV;
        iy = iv[j];
        iv[j] = *idum;
        return MIN(AM*iy,RNMX);
}
#undef IA 
#undef IM 
#undef IQ
#undef IR
#undef NTAB
#undef EPS 
#undef MAX
#undef MIN


void gen_Ak(double *Ak, int kx, int ky, int kz, double Pk, int dots)
{
    int nkx, nky, nkz;
    double phi, rnd1, rnd2, r;
    
    
    rnd1 = ran1(&iiseed);
    phi = 2 * PI * rnd1;                /* phase */

    rnd2 = ran1(&iiseed);
    r = sqrt(-2 * log(rnd2)*Pk);         /* amplitude */


    Ak[2*(kx*dots*dots + ky*dots + kz)] = r*cos(phi); /* real part */
    Ak[2*(kx*dots*dots + ky*dots + kz)+1] = r*sin(phi); /* imag. part */
    

    if((kx != 0) && (ky != 0) && (kz != 0))
      {
	  nkx = dots - kx;
	  nky = dots - ky;
	  nkz = dots - kz;

	  Ak[2*(nkx*dots*dots + nky*dots + nkz)] = r*cos(phi); /* real part */
	  Ak[2*(nkx*dots*dots + nky*dots + nkz)+1] = -r*sin(phi); /* imag. part */
      }
}



void rsk_turbdriving_field()
{
    int dots, kx, ky, kz, i,l, kmin, kmax, bytes; 
    int nn[4];
    double rcube, rcubehalf, xk, yk, zk, Pk, abs_k;
    /* double Ak[2*DOTS*DOTS*DOTS+1], matrix[DOTS][DOTS][DOTS]; */
    double *Ak, *matrix;
   
    if(!(Ak = malloc(bytes=sizeof(double)*(2*DOTS*DOTS*DOTS+1))))
    {
      printf("failed to allocate memory for Ak: %d.\n",bytes);
      endrun(2);
    }
    if(!(matrix = malloc(sizeof(double)*DOTS*DOTS*DOTS)))
    {
      printf("failed to allocate memory for matrix.\n");
      endrun(2);
    }
    dots = DOTS;
    rcube = 2.0;
    rcubehalf = rcube/2.0;

    kmin = All.kMin;
    kmax = All.kMax;

    All.DeltaTimeDrv = 0;
/* /\* ---------------------------------------------------------------------*\/ */
    
    ran1(&iiseed);

    nn[1] = dots; /* dimensions of array */
    nn[2] = dots;
    nn[3] = dots;

    if(kmin < 1)
	kmin = 1;
    for (kx = 0; kx < dots; kx++)
      {
	  xk = kx;
	  if(kx > dots/2)
	      xk = dots - xk;
	  for (ky = 0; ky < dots; ky++)
	    {
		yk = ky;
		if(ky > dots/2)
		    yk = dots - yk;
		for (kz = 0; kz < dots/2+1; kz++)
		  {
		      zk = kz;
		      abs_k = sqrt(xk*xk + yk*yk + zk*zk);
		      if((abs_k < kmin) || (abs_k > kmax))
			  Pk = 0;
		      else
			  Pk = 1/pow(abs_k,All.DrvIndx);
		      gen_Ak(Ak,kx,ky,kz,Pk,dots);   
		  }
	    }

      }

    Ak[2*(     0*dots*dots +      0*dots +      0)+1] = 0;
    Ak[2*(     0*dots*dots +      0*dots + dots/2)+1] = 0;
    Ak[2*(     0*dots*dots + dots/2*dots +      0)+1] = 0;
    Ak[2*(     0*dots*dots + dots/2*dots + dots/2)+1] = 0;
    Ak[2*(dots/2*dots*dots +      0*dots +      0)+1] = 0;
    Ak[2*(dots/2*dots*dots +      0*dots + dots/2)+1] = 0;
    Ak[2*(dots/2*dots*dots + dots/2*dots +      0)+1] = 0;
    Ak[2*(dots/2*dots*dots + dots/2*dots + dots/2)+1] = 0;

    for(i=2*dots*dots*dots;i>0;i--) /* array must start with Ak[1] for fourn */
      Ak[i]=Ak[i-1];

    fourn(Ak,nn,3,-1);

    for(i=1,l=0;i<2*dots*dots*dots;i+=2,l++)
	Ak[l] = Ak[i];

    l = 0;                           /*   assign the real part to matrix */
    for(kz = 0; kz<dots; kz++)
      {
	  for(ky = 0; ky<dots; ky++)
	    {
		for(kx = 0; kx<dots; kx++)
		  {
		      /* matrix[kx][ky][kz] = Ak[l]; */
		    matrix[kx*DOTS*DOTS+ky*DOTS+kz] = Ak[l];
		      l++;
		  }
	    }
      }

    for(l = 0; l<2*dots*dots*dots+1; l++) /*empty Ak*/
	 Ak[l] = 0;

    if(dots != IFIELDSIZE)
      {
	if(ThisTask==0)
	  {
	    printf("STOP: array sizes do nor agree\n");
	    printf("dots: %d\n",dots);
	    printf("ifieldsize: %d\n",IFIELDSIZE);
	    fflush(stdout);
	  }
        exit(1);
      }

    if(icomp == 1)             /*   copy matrix into xmatrix */
      {
	  for(kz = 0; kz<IFIELDSIZE; kz++)
	    {
		for(ky = 0; ky<IFIELDSIZE; ky++)
		  {
		      for(kx = 0; kx<IFIELDSIZE; kx++)
			{
			  xmatrix[kx*IFIELDSIZE*IFIELDSIZE+ky*IFIELDSIZE+kz] = matrix[kx*IFIELDSIZE*IFIELDSIZE+ky*IFIELDSIZE+kz];
			}
		  }
	    }
      }
    if(icomp == 2)              /*  copy matrix into ymatrix */
      {
	  for(kz = 0; kz<IFIELDSIZE; kz++)
	    {
		for(ky = 0; ky<IFIELDSIZE; ky++)
		  {
		      for(kx = 0; kx<IFIELDSIZE; kx++)
			{
			  ymatrix[kx*IFIELDSIZE*IFIELDSIZE+ky*IFIELDSIZE+kz] = matrix[kx*IFIELDSIZE*IFIELDSIZE+ky*IFIELDSIZE+kz];
			}
		  }
	    }
      }
    if(icomp == 3)             /*   copy matrix into zmatrix */
      {
	  for(kz = 0; kz<IFIELDSIZE; kz++)
	    {
		for(ky = 0; ky<IFIELDSIZE; ky++)
		  {
		      for(kx = 0; kx<IFIELDSIZE; kx++)
			{
			  zmatrix[kx*IFIELDSIZE*IFIELDSIZE+ky*IFIELDSIZE+kz] = matrix[kx*IFIELDSIZE*IFIELDSIZE+ky*IFIELDSIZE+kz];
			}
		  }
	    }
      }

    SysState.EnergyDrv = 0.0;
    for(i=0;i<6;i++)
      {
	SysState.EnergyDrvComp[i] = 0.0;
	SysState.EnergyDrv += SysState.EnergyDrvComp[i];
      }
    
/* ende */
    free(Ak);
    free(matrix);
    }  
#endif
