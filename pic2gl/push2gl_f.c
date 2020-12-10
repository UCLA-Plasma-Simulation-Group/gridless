/* C Library for Skeleton 2D Electrostatic Gridless PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr2_(double *part, double *vtx, double *vty, double *vdx,
             double *vdy, int *npx, int *npy, int *idimp, int *nop,
             int *nx, int *ny,
             int *ipbc);

void push2gl_(double *part, double complex *fxy, double complex *sctx,
              double *qbm, double *dt, double *ek, int *idimp, int *nop,
              int *nx, int *ny, int *nxhm, int *nyhm, int *nxvh,
              int *nyv, int *ipbc);

void dpost2gl_(double *part, double complex *q, double complex *sctx,
               double *qm, int *nop, int *idimp, int *nx, int *ny,
               int *nxhm, int *nyhm, int *nxvh, int *nyv);

void ffc2initgl_(double complex *ffc, double *ax, double *ay,
                 double *affp, int *nx, int *ny, int *nxhm, int *nyhm);

void pois22gl_(double complex *q, double complex *fxy, int *isign,
               double complex *ffc, double *ax, double *ay,
               double *affp, double *we, int *nx, int *ny, int *nxhm,
               int *nyhm, int *nxvh, int *nyv);

void evfield22gl_(double complex *fxy, double *gxy,
                  double complex *sctx, double *xp, double *yp, int *nx,
                  int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

/*--------------------------------------------------------------------*/
void cdistr2(double part[], double vtx, double vty, double vdx,
             double vdy, int npx, int npy, int idimp, int nop, int nx,
             int ny, int ipbc) {
   distr2_(part,&vtx,&vty,&vdx,&vdy,&npx,&npy,&idimp,&nop,&nx,&ny,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cpush2gl(double part[], double complex fxy[],
              double complex sctx[], double qbm, double dt, double *ek,
              int idimp, int nop, int nx, int ny, int nxhm, int nyhm, 
              int nxvh, int nyv, int ipbc) {
   push2gl_(part,fxy,sctx,&qbm,&dt,ek,&idimp,&nop,&nx,&ny,&nxhm,&nyhm,
            &nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdpost2gl(double  part[], double  complex q[],
               double  complex sctx[], double  qm, int nop, int idimp,
               int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   dpost2gl_(part,q,sctx,&qm,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,&nxvh,
             &nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cffc2initgl(double complex ffc[], double ax, double ay,
                 double affp, int nx, int ny, int nxhm, int nyhm) {
   ffc2initgl_(ffc,&ax,&ay,&affp,&nx,&ny,&nxhm,&nyhm);
   return;
}

/*--------------------------------------------------------------------*/
void cpois22gl(double complex q[], double complex fxy[], int isign,
               double complex ffc[], double ax, double ay, double affp, 
               double *we, int nx, int ny, int nxhm, int nyhm, int nxvh,
               int nyv) {
   pois22gl_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxhm,&nyhm,
             &nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cevfield22gl(double complex fxy[], double gxy[],
                  double complex sctx[], double xp, double yp, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   evfield22gl_(fxy,gxy,sctx,&xp,&yp,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}
