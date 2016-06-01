#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif


void rsk_turbdriving_field(void);
void gen_Ak(double *Ak, int kx, int ky, int kz, double Pk, int dots);
void fourn(double data[], int nn[], int ndim, int isign);
double ran1(int *idum);
void rsk_turbdriving(void);

