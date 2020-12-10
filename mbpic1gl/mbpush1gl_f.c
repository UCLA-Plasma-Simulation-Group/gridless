/* C Library for Skeleton 1-2/2D Electromagnetic Gridless OpenMP PIC */
/* Code Wrappers for calling the Fortran routines from a C main      */
/* program */

#include <complex.h>

double ranorm_();

double randum_();

void distr1h_(double *part, double *vtx, double*vty, double *vtz,
              double *vdx, double *vdy, double *vdz, int *npx,
              int *idimp, int *nop, int *nx, int *ipbc);

void mbpush13gl_(double *part, double complex *fxyz,
                 double complex *byz, double *omx, double *qbm,
                 double *dt, double *dtc, double *ek, int *idimp,
                 int *nop, int *nx, int *nxhm, int *nxvh);

void mrbpush13gl_(double *part, double complex *fxyz,
                  double complex *byz, double *omx, double *qbm,
                  double *dt, double *dtc, double *ci, double *ek,
                  int *idimp, int *nop, int *nx, int *nxhm, int *nxvh);

void mabpush13gl_(double *part, double complex *fxyz,
                  double complex *byz, double *omx, double *qbm,
                  double *dt, double *dtc, double *ek, int *idimp,
                  int *nop, int *nx, int *nxhm, int *nxvh);

void marbpush13gl_(double *part, double complex *fxyz,
                   double complex *byz, double *omx, double *qbm,
                   double *dt, double *dtc, double *ci, double *ek,
                   int *idimp, int *nop, int *nx, int *nxhm, int *nxvh);

void mearbpush13gl_(double *part, double complex *fxyz,
                    double complex *byz, double *omx, double *qbm,
                    double *dt, double *dtc, double *ci, double *ek,
                    int *idimp, int *nop, int *nx, int *nxhm,
                    int *nxvh);

void mdpost1gl_(double *part, double complex *q, double *qm, int *nop,
                int *idimp, int *nx, int *nxhm, int *nxvh);

void mdjpost1gl_(double *part, double complex *cu, double *cux0,
                 double *qm, double *dt, int *nop, int *idimp, int *nx,
                 int *nxhm, int *nxvh);

void mrdjpost1gl_(double *part, double complex *cu, double *cux0,
                  double *qm, double *dt, double *ci, int *nop,
                  int *idimp, int *nx, int *nxhm, int *nxvh);

void md2jpost1gl_(double *part, double complex *dcu, double *qm,
                  int *nop, int *idimp, int *nx, int *nxhm, int *nxvh);

void mrd2jpost1gl_(double *part, double complex *dcu, double *qm,
                   double *ci, int *nop, int *idimp, int *nx, int *nxhm,
                   int *nxvh);

void mdsjpost1gl_(double *part, double complex *cu, double complex *dcu,
                  double *qm, double *ci, int *nop, int *idimp, int *nx,
                  int *nxhm, int *nxvh);

void mrdsjpost1gl_(double *part, double complex *cu,
                   double complex *dcu, double *qm, double *ci,
                   int *nop, int *idimp, int *nx, int *nxhm, int *nxvh);

void ffc1initgl_(double complex *ffc, double *ax, double *affp, int *nx,
                 int *nxhm);

void pois1gl_(double complex *q, double complex *fx, int *isign,
              double complex *ffc, double *ax, double *affp, double *we,
              int *nx, int *nxhm, int *nxvh);

void ibpois13gl_(double complex *cu, double complex *byz,
                 double complex *ffc, double *ci, double *wm, int *nx,
                 int *nxhm, int *nxvh);

void amaxwel1gl_(double complex *eyz, double complex *byz,
                 double complex *cu, double complex *ffc, double *ci,
                 double *dt, double *wf, double *wm, int *nx, int *nxhm,
                 int *nxvh);

void emfield1gl_(double complex *fxyz, double complex *fx, 
                 double complex *eyz, double complex *ffc, int *nxhm,
                 int *nxvh);

void bmfield1gl_(double complex *fyz, double complex *eyz,
                 double complex *ffc, int *nxhm, int *nxvh);

void epois13gl_(double complex *dcu, double complex *eyz, int *isign,
                double complex *ffe, double *ax, double *affp,
                double *wp0, double *ci, double *wf, int *nx, int *nxhm,
                int *nxvh);

void maxwel1gl_(double complex *eyz, double complex *byz,
                double complex *cu, double complex *ffc, double *ci,
                double *dt, double *wf, double *wm, int *nx, int *nxhm,
                int *nxvh);
void mevfield13gl_(double complex *fxy, double *gxy, double *xp,
                   int *nx, int *nxhm, int *nxvh);

void elfield1gl_(double complex *q, double complex *fx,
                 double complex *ffc, double *we, int *nx, int *nxhm,
                 int *nxvh);

void poynt1gl_(double complex *q, double complex *eyz,
               double complex *byz, double complex *ffc, double *ex0,
               double *sx, double *sy, double *sz, int *nx, int *nxhm, 
               int *nxvh);

void ma0rbpush13gl_(double *part, double complex *fxyz,
                    double complex *byz, double *omx, double *qbm,
                    double *dt, double *dtc, double *ci, double *ek,
                    int *idimp, int *nop, int *nx, int *nxhm,
                    int *nxvh);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

/*--------------------------------------------------------------------*/
void cdistr1h(double part[], double vtx, double vty, double vtz,
              double vdx, double vdy, double vdz, int npx, int idimp,
              int nop, int nx, int ipbc) {
   distr1h_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&idimp,&nop,&nx,
            &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmbpush13gl(double part[], double complex fxyz[],
                 double complex byz[], double omx, double qbm,
                 double dt, double dtc, double *ek, int idimp, int nop,
                 int nx, int nxhm, int nxvh) {
   mbpush13gl_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,
               &nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmrbpush13gl(double part[], double complex fxyz[],
                  double complex byz[], double omx, double qbm,
                  double dt, double dtc, double ci, double *ek,
                  int idimp, int nop, int nx, int nxhm, int nxvh) {
   mrbpush13gl_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,
                &nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmabpush13gl(double part[], double complex fxyz[],
                  double complex byz[], double omx, double qbm,
                  double dt, double dtc, double *ek, int idimp, int nop,
                  int nx, int nxhm, int nxvh) {
   mabpush13gl_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,
                &nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmarbpush13gl(double part[], double complex fxyz[],
                   double complex byz[], double omx, double qbm,
                   double dt, double dtc, double ci, double *ek,
                   int idimp, int nop, int nx, int nxhm, int nxvh) {
   marbpush13gl_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,
                 &nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmearbpush13gl(double part[], double complex fxyz[],
                    double complex byz[], double omx, double qbm,
                    double dt, double dtc, double ci, double *ek,
                    int idimp, int nop, int nx, int nxhm, int nxvh) {
   mearbpush13gl_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,
                  &nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmdpost1gl(double part[], double complex q[], double qm, int nop,
                int idimp, int nx, int nxhm, int nxvh) {
   mdpost1gl_(part,q,&qm,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmdjpost1gl(double part[], double complex cu[], double *cux0,
                 double qm, double dt, int nop, int idimp, int nx,
                 int nxhm, int nxvh) {
   mdjpost1gl_(part,cu,cux0,&qm,&dt,&nop, &idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmrdjpost1gl(double part[], double complex cu[], double *cux0,
                  double qm, double dt, double ci, int nop, int idimp,
                  int nx, int nxhm, int nxvh) {
   mrdjpost1gl_(part,cu,cux0,&qm,&dt,&ci,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmd2jpost1gl(double part[], double complex dcu[], double qm,
                  int nop, int idimp, int nx, int nxhm, int nxvh) {
   md2jpost1gl_(part,dcu,&qm,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmrd2jpost1gl(double part[], double complex dcu[], double qm,
                   double ci, int nop, int idimp, int nx, int nxhm,
                   int nxvh) {
   mrd2jpost1gl_(part,dcu,&qm,&ci,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmdsjpost1gl(double part[], double complex cu[],
                  double complex dcu[], double qm, double ci, int nop,
                  int idimp, int nx, int nxhm, int nxvh) {
   mdsjpost1gl_(part,cu,dcu,&qm,&ci,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmrdsjpost1gl(double part[], double complex cu[], 
                   double complex dcu[], double qm, double ci, int nop,
                   int idimp, int nx, int nxhm, int nxvh) {
   mrdsjpost1gl_(part,cu,dcu,&qm,&ci,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cffc1initgl(double complex ffc[], double ax, double affp, int nx,
                 int nxhm) {
   ffc1initgl_(ffc,&ax,&affp,&nx,&nxhm);
   return;
}

/*--------------------------------------------------------------------*/
void cpois1gl(double complex q[], double complex fx[], int isign,
              double complex ffc[], double ax, double affp, double *we,
              int nx, int nxhm, int nxvh) {
   pois1gl_(q,fx,&isign,ffc,&ax,&affp,we,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cibpois13gl(double complex cu[], double complex byz[],
                 double complex ffc[], double ci, double *wm, int nx,
                 int nxhm, int nxvh) {
   ibpois13gl_(cu,byz,ffc,&ci,wm,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void camaxwel1gl(double complex eyz[], double complex byz[],
                 double complex cu[], double complex ffc[], double ci,
                 double dt, double *wf, double *wm, int nx, int nxhm,
                 int nxvh) {
   amaxwel1gl_(eyz,byz,cu,ffc,&ci,&dt,wf,wm,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cemfield1gl(double complex fxyz[], double complex fx[], 
                 double complex eyz[], double complex ffc[], int nxhm,
                 int nxvh) {
   emfield1gl_(fxyz,fx,eyz,ffc,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cbmfield1gl(double complex fyz[], double complex eyz[],
                 double complex ffc[], int nxhm, int nxvh) {
   bmfield1gl_(fyz,eyz,ffc,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cepois13gl(double complex dcu[], double complex eyz[], int isign,
                 double complex ffe[], double ax, double affp,
                 double wp0, double ci, double *wf, int nx, int nxhm,
                 int nxvh) {
   epois13gl_(dcu,eyz,&isign,ffe,&ax,&affp,&wp0,&ci,wf,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmaxwel1gl(double complex eyz[], double complex byz[],
                double complex cu[], double complex ffc[], double ci,
                double dt, double *wf, double *wm, int nx, int nxhm,
                int nxvh) {
   maxwel1gl_(eyz,byz,cu,ffc,&ci,&dt,wf,wm,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmevfield13gl(double complex fxy[], double gxy[], double xp,
                   int nx, int nxhm, int nxvh) {
   mevfield13gl_(fxy,gxy,&xp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void celfield1gl(double complex q[], double complex fx[],
                 double complex ffc[], double *we, int nx, int nxhm,
                 int nxvh) {
   elfield1gl_(q,fx,ffc,we,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cpoynt1gl(double complex q[], double complex eyz[],
               double complex byz[], double complex ffc[], double *ex0,
               double *sx, double *sy, double *sz, int nx, int nxhm, 
               int nxvh) {
   poynt1gl_(q,eyz,byz,ffc,ex0,sx,sy,sz,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cma0rbpush13gl(double part[], double complex fxyz[],
                    double complex byz[], double omx, double qbm,
                    double dt, double dtc, double ci, double *ek,
                    int idimp, int nop, int nx, int nxhm, int nxvh) {
   ma0rbpush13gl_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,
                  &nx,&nxhm,&nxvh);
   return;
}
