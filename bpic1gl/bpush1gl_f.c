/* C Library for Skeleton 1-2/2D Electromagnetic Gridless PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr1h_(double *part, double *vtx, double*vty, double *vtz,
              double *vdx, double *vdy, double *vdz, int *npx,
              int *idimp, int *nop, int *nx, int *ipbc);

void bpush13gl_(double *part, double complex *fxyz, double complex *byz,
                double complex *sctx, double *omx, double *qbm,
                double *dt, double *dtc, double *ek, int *idimp, 
                int *nop, int *nx, int *nxhm, int *nxvh);

void rbpush13gl_(double *part, double complex *fxyz,
                 double complex *byz, double complex *sctx,
                 double *omx, double *qbm, double *dt, double *dtc,
                 double *ci, double *ek, int *idimp, int *nop, int *nx,
                 int *nxhm, int *nxvh);

void abpush13gl_(double *part, double complex *fxyz,
                 double complex *byz, double complex *sctx, double *omx, 
                 double *qbm, double *dt, double *dtc, double *ek, 
                 int *idimp, int *nop, int *nx, int *nxhm, int *nxvh);

void arbpush13gl_(double *part, double complex *fxyz,
                  double complex *byz, double complex *sctx,
                  double *omx, double *qbm, double *dt, double *dtc,
                  double *ci, double *ek, int *idimp, int *nop, int *nx,
                  int *nxhm, int *nxvh);

void earbpush13gl_(double *part, double complex *fxyz,
                   double complex *byz, double complex *sctx,
                   double *omx, double *qbm, double *dt, double *dtc,
                   double *ci, double *ek, int *idimp, int *nop,
                   int *nx, int *nxhm, int *nxvh);

void dpost1gl_(double *part, double complex *q, double complex *sctx,
               double *qm, int *nop, int *idimp, int *nx, int *nxhm,
               int *nxvh);

void djpost1gl_(double *part, double complex *cu, double complex *sctx, 
                double *cux0, double *qm, double *dt, int *nop,
                int *idimp, int *nx, int *nxhm, int *nxvh);

void rdjpost1gl_(double *part, double complex *cu, double complex *sctx,
                 double *cux0, double *qm, double *dt, double *ci,
                 int *nop, int *idimp, int *nx, int *nxhm, int *nxvh);

void d2jpost1gl_(double *part, double complex *dcu,
                 double complex *sctx, double *qm, int *nop, int *idimp,
                 int *nx, int *nxhm, int *nxvh);


void rd2jpost1gl_(double *part, double complex *dcu,
                  double complex *sctx, double *qm, double *ci,
                  int *nop, int *idimp, int *nx, int *nxhm, int *nxvh);

void dsjpost1gl_(double *part, double complex *cu, double complex *dcu,
                 double complex *sctx, double *qm, double *ci, int *nop,
                 int *idimp, int *nx, int *nxhm, int *nxvh);

void rdsjpost1gl_(double *part, double complex *cu, double complex *dcu,
                  double complex *sctx, double *qm, double *ci,
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

void evfield13gl_(double complex *fxy, double *gxy,
                  double complex *sctx, double *xp, int *nx, int *nxhm,
                  int *nxvh);

void mjpost1gl_(double *part, double complex *amu, double complex *sctx,
                double *qm, int *nop, int *idimp, int *nx, int *nxhm,
                int *nxvh);

void rmjpost1gl_(double *part, double complex *amu,
                 double complex *sctx, double *qm, double *ci, int *nop,
                 int *idimp, int *nx, int *nxhm, int *nxvh);

void dcuperp13gl_(double complex *dcu, double complex *amu, int *nx,
                  int *nxhm, int *nxvh);

void maxwel1gl_(double complex *eyz, double complex *byz,
                double complex *cu, double complex *ffc, double *ci,
                double *dt, double *wf, double *wm, int *nx, int *nxhm,
                int *nxvh);

void elfield1gl_(double complex *q, double complex *fx,
                 double complex *ffc, double *we, int *nx, int *nxhm,
                 int *nxvh);

void psmooth1gl_(double complex *q, double complex *qs,
                 double complex *ffc, int *nxhm, int *nxvh);

void poynt1gl_(double complex *q, double complex *eyz,
               double complex *byz, double complex *ffc, double *ex0,
               double *sx, double *sy, double *sz, int *nx, int *nxhm, 
               int *nxvh);

void a0rbpush13gl_(double *part, double complex *fxyz,
                   double complex *byz, double complex *sctx,
                   double *omx, double *qbm, double *dt, double *dtc,
                   double *ci, double *ek, int *idimp, int *nop, 
                   int *nx, int *nxhm, int *nxvh);

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
void cbpush13gl(double part[], double complex fxyz[],
                double complex byz[], double complex sctx[], double omx,
                double qbm, double dt, double dtc, double *ek,
                int idimp, int nop, int nx, int nxhm, int nxvh) {
   bpush13gl_(part,fxyz,byz,sctx,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,
              &nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void crbpush13gl(double part[], double complex fxyz[],
                 double complex byz[], double complex sctx[],
                 double omx, double qbm, double dt, double dtc,
                 double ci, double *ek, int idimp, int nop, int nx,
                 int nxhm, int nxvh) {
   rbpush13gl_(part,fxyz,byz,sctx,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,
               &nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cabpush13gl(double part[], double complex fxyz[],
                 double complex byz[], double complex sctx[], 
                 double omx, double qbm, double dt, double dtc,
                 double *ek, int idimp, int nop, int nx, int nxhm,
                 int nxvh) {
   abpush13gl_(part,fxyz,byz,sctx,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,
              &nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void carbpush13gl(double part[], double complex fxyz[],
                  double complex byz[], double complex sctx[],
                  double omx, double qbm, double dt, double dtc,
                  double ci, double *ek, int idimp, int nop, int nx,
                  int nxhm, int nxvh) {
   arbpush13gl_(part,fxyz,byz,sctx,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,
                &nop,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cearbpush13gl(double part[], double complex fxyz[],
                   double complex byz[], double complex sctx[],
                   double omx, double qbm, double dt, double dtc,
                   double ci, double *ek, int idimp, int nop, int nx,
                   int nxhm, int nxvh) {
   earbpush13gl_(part,fxyz,byz,sctx,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,
                 &nop,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cdpost1gl(double part[], double complex q[], double complex sctx[],
               double qm, int nop, int idimp, int nx, int nxhm,
               int nxvh) {
   dpost1gl_(part,q,sctx,&qm,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cdjpost1gl(double part[], double complex cu[],
                double complex sctx[], double *cux0, double qm,
                double dt, int nop, int idimp, int nx, int nxhm,
                int nxvh) {
   djpost1gl_(part,cu,sctx,cux0,&qm,&dt,&nop, &idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void crdjpost1gl(double part[], double complex cu[],
                 double complex sctx[], double *cux0, double qm,
                 double dt, double ci, int nop, int idimp, int nx,
                 int nxhm, int nxvh) {
   rdjpost1gl_(part,cu,sctx,cux0,&qm,&dt,&ci,&nop,&idimp,&nx,&nxhm,
               &nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cd2jpost1gl(double part[], double complex dcu[],
                 double complex sctx[], double qm, int nop, int idimp,
                 int nx, int nxhm, int nxvh) {
   d2jpost1gl_(part,dcu,sctx,&qm,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void crd2jpost1gl(double part[], double complex dcu[],
                  double complex sctx[], double qm, double ci, int nop,
                  int idimp, int nx, int nxhm, int nxvh) {
   rd2jpost1gl_(part,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cdsjpost1gl(double part[], double complex cu[],
                 double complex dcu[], double complex sctx[], double qm,
                 double ci, int nop, int idimp, int nx, int nxhm, 
                 int nxvh) {
   dsjpost1gl_(part,cu,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void crdsjpost1gl(double part[], double complex cu[], 
                  double complex dcu[], double complex sctx[],
                  double qm, double ci, int nop, int idimp, int nx, 
                  int nxhm, int nxvh) {
   rdsjpost1gl_(part,cu,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&nxhm,&nxvh);
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
void cevfield13gl(double complex fxy[], double gxy[],
                  double complex sctx[], double xp, int nx, int nxhm,
                  int nxvh) {
   evfield13gl_(fxy,gxy,sctx,&xp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmjpost1gl(double part[], double complex amu[],
                double complex sctx[], double qm, int nop, int idimp,
                int nx, int nxhm, int nxvh) {
   mjpost1gl_(part,amu,sctx,&qm,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void crmjpost1gl(double part[], double complex amu[],
                 double complex sctx[], double qm, double ci, int nop,
                 int idimp, int nx, int nxhm, int nxvh) {
   rmjpost1gl_(part,amu,sctx,&qm,&ci,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cdcuperp13gl(double complex dcu[], double complex amu[], int nx,
                  int nxhm, int nxvh) {
   dcuperp13gl_(dcu,amu,&nx,&nxhm,&nxvh);
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
void celfield1gl(double complex q[], double complex fx[],
                 double complex ffc[], double *we, int nx, int nxhm,
                 int nxvh) {
   elfield1gl_(q,fx,ffc,we,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cpsmooth1gl(double complex q[], double complex qs[],
                 double complex ffc[], int nxhm, int nxvh) {
   psmooth1gl_(q,qs,ffc,&nxhm,&nxvh);
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
void ca0rbpush13gl(double part[], double complex fxyz[],
                   double complex byz[], double complex sctx[],
                   double omx, double qbm, double dt, double dtc,
                   double ci, double *ek, int idimp, int nop, int nx,
                   int nxhm, int nxvh) {
   a0rbpush13gl_(part,fxyz,byz,sctx,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,
                &nop,&nx,&nxhm,&nxvh);
   return;
}

