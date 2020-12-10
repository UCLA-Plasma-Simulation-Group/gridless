/* C Library for Skeleton 2-1/2D Electromagnetic Gridless PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr2h_(double *part, double *vtx, double *vty, double *vtz,
              double *vdx, double *vdy, double *vdz, int *npx, int *npy,
              int *idimp, int *nop, int *nx, int *ny, int *ipbc);

void bpush23gl_(double *part, double complex *fxy, double complex *bxy,
                double complex *sctx, double *qbm, double *dt,
                double *dtc, double *ek, int *idimp, int *nop, int *nx,
                int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv,
                int *ipbc);

void rbpush23gl_(double *part, double complex *fxy, double complex *bxy,
                 double complex *sctx, double *qbm, double *dt,
                 double *dtc, double *ci, double *ek, int *idimp, 
                 int *nop, int *nx, int *ny, int *nxhm, int *nyhm,
                 int *nxvh, int *nyv, int *ipbc);

void abpush23gl_(double *part, double complex *fxy, double complex *bxy,
                 double complex *sctx, double *qbm, double *dt,
                 double *dtc, double *ek, int *idimp, int *nop, int *nx,
                 int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv,
                 int *ipbc);
void arbpush23gl_(double *part, double complex *fxy,
                  double complex *bxy, double complex *sctx, 
                  double *qbm, double *dt, double *dtc, double *ci,
                  double *ek, int *idimp, int *nop, int *nx, int *ny,
                  int *nxhm, int *nyhm, int *nxvh, int *nyv, int *ipbc);

void earbpush23gl_(double *part, double complex *fxy,
                   double complex *bxy, double complex *sctx, 
                   double *qbm, double *dt, double *dtc, double *ci,
                   double *ek, int *idimp, int *nop, int *nx, int *ny,
                   int *nxhm, int *nyhm, int *nxvh, int *nyv,
                   int *ipbc);

void dpost2gl_(double *part, double complex *q, double complex *sctx,
               double *qm, int *nop, int *idimp, int *nx, int *ny,
               int *nxhm, int *nyhm, int *nxvh, int *nyv);

void djpost2gl_(double *part, double complex *cu, double complex *sctx,
                double *qm, double *dt, int *nop, int *idimp, int *nx,
                int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv,
                int *ipbc);

void rdjpost2gl_(double *part, double complex *cu, double complex *sctx,
                 double *qm, double *dt, double *ci, int *nop,
                 int *idimp, int *nx, int *ny, int *nxhm, int *nyhm,
                 int *nxvh, int *nyv, int *ipbc);

void d2jpost2gl_(double *part, double complex *dcu,
                 double complex *sctx, double *qm, int *nop, int *idimp,
                 int *nx, int *ny, int *nxhm, int *nyhm, int *nxvh,
                 int *nyv);

void rd2jpost2gl_(double *part, double complex *dcu,
                  double complex *sctx, double *qm, double *ci,
                  int *nop, int *idimp, int *nx, int *ny, int *nxhm, 
                  int *nyhm, int *nxvh, int *nyv);

void dsjpost2gl_(double *part, double complex *cu, double complex *dcu, 
                 double complex *sctx, double *qm, double *ci, int *nop,
                 int *idimp, int *nx, int *ny, int *nxhm, int *nyhm,
                 int *nxvh, int *nyv);

void rdsjpost2gl_(double *part, double complex *cu, double complex *dcu,
                  double complex *sctx, double *qm, double *ci,
                  int *nop, int *idimp, int *nx, int *ny, int *nxhm,
                  int *nyhm, int *nxvh, int *nyv);

void ffc2initgl_(double complex *ffc, double *ax, double *ay,
                 double *affp, int *nx, int *ny, int *nxhm, int *nyhm);

void pois23gl_(double complex *q, double complex *fxy, int *isign,
               double complex *ffc, double *ax, double *ay,
               double *affp, double *we, int *nx, int *ny, int *nxhm,
               int *nyhm, int *nxvh, int *nyv);

void cuperp2gl_(double complex *cu, int *nx, int *ny, int *nxhm,
                int *nyhm, int *nxvh, int *nyv);

void ibpois23gl_(double complex *cu, double complex *bxy,
                 double complex *ffc, double *ci, double *wm, int *nx,
                 int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv);

void amaxwel2gl_(double complex *exy, double complex *bxy,
                 double complex *cu, double complex *ffc, double *ci,
                 double *dt, double *wf, double *wm, int *nx, int *ny,
                 int *nxhm, int *nyhm, int *nxvh, int *nyv);

void emfield2gl_(double complex *fxy, double complex *exy,
                 double complex *ffc, int *isign, int *nxhm, int *nyhm,
                 int *nxvh, int *nyv);

void epois23gl_(double complex *dcu, double complex *exy, int *isign,
                double complex *ffe, double *ax, double *ay,
                double *affp, double *wp0, double *ci, double *wf,
                int *nx, int *ny, int *nxhm, int *nyhm, int *nxvh,
                int *nyv);

void evfield23gl_(double complex *fxy, double *gxy,
                  double complex *sctx, double *xp, double *yp, int *nx,
                  int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv);

void mjpost2gl_(double *part, double complex *amu, double complex *sctx,
                double *qm, int *nop, int *idimp, int *nx, int *ny,
                int *nxhm, int *nyhm, int *nxvh, int *nyv);

void rmjpost2gl_(double *part, double complex *amu,
                 double complex *sctx, double *qm, double *ci, int *nop,
                 int *idimp, int *nx, int *ny, int *nxhm, int *nyhm,
                 int *nxvh, int *nyv);

void dcuperp23gl_(double complex *dcu, double complex *amu, int *nx,
                  int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv);

void adcuperp23gl_(double complex *dcu, double complex *amu, int *nx,
                   int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv);

void maxwel2gl_(double complex *exy, double complex *bxy,
                double complex *cu, double complex *ffc, double *ci,
                double *dt, double *wf, double *wm, int *nx, int *ny,
                int *nxhm, int *nyhm, int *nxvh, int *nyv);

void elfield23gl_(double complex *q, double complex *fxy, 
                  double complex *ffc, double *we, int *nx, int *ny,
                  int *nxhm, int *nyhm, int *nxvh, int *nyv);

void psmooth2gl_(double complex *q, double complex *qs,
                 double complex *ffc, int *ny, int *nxhm, int *nyhm,
                 int *nxvh, int *nyv);

void poynt2gl_(double complex *q, double complex *exy,
               double complex *bxy, double complex *ffc, double *sx,
               double *sy, double *sz, int *nx, int *ny, int *nxhm,
               int *nyhm, int *nxvh, int *nyv);

void a0rbpush23gl_(double *part, double complex *fxy,
                   double complex *bxy, double complex *sctx, 
                   double *qbm, double *dt, double *dtc, double *ci,
                   double *ek, int *idimp, int *nop, int *nx, int *ny,
                   int *nxhm, int *nyhm, int *nxvh, int *nyv,
                   int *ipbc);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

/*--------------------------------------------------------------------*/
void cdistr2h(double part[], double vtx, double vty, double vtz,
              double vdx, double vdy, double vdz, int npx, int npy,
              int idimp, int nop, int nx, int ny, int ipbc) {
   distr2h_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&npy,&idimp,&nop,
            &nx,&ny,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cbpush23gl(double part[], double complex fxy[],
                double complex bxy[], double complex sctx[], double qbm,
                double dt, double dtc, double *ek, int idimp, int nop,
                int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv,
                int ipbc) {
   bpush23gl_(part,fxy,bxy,sctx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&ny,
              &nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void crbpush23gl(double part[], double complex fxy[],
                 double complex bxy[], double complex sctx[],
                 double qbm, double dt, double dtc, double ci,
                 double *ek, int idimp, int nop, int nx, int ny,
                 int nxhm, int nyhm, int nxvh, int nyv, int ipbc) {
   rbpush23gl_(part,fxy,bxy,sctx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,
               &ny,&nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cabpush23gl(double part[], double complex fxy[],
                 double complex bxy[], double complex sctx[],
                 double qbm, double dt, double dtc, double *ek,
                 int idimp, int nop, int nx, int ny, int nxhm, int nyhm,
                 int nxvh, int nyv, int ipbc) {
   abpush23gl_(part,fxy,bxy,sctx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&ny,
               &nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void carbpush23gl(double part[], double complex fxy[],
                  double complex bxy[], double complex sctx[],
                  double qbm, double dt, double dtc, double ci,
                  double *ek, int idimp, int nop, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv, int ipbc) {
   arbpush23gl_(part,fxy,bxy,sctx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,
                &ny,&nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cearbpush23gl(double part[], double complex fxy[],
                   double complex bxy[], double complex sctx[],
                   double qbm, double dt, double dtc, double ci,
                   double *ek, int idimp, int nop, int nx, int ny,
                   int nxhm, int nyhm, int nxvh, int nyv, int ipbc) {
   earbpush23gl_(part,fxy,bxy,sctx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,
                 &ny,&nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdpost2gl(double part[], double complex q[], double complex sctx[],
               double qm, int nop, int idimp, int nx, int ny, int nxhm,
               int nyhm, int nxvh, int nyv) {
   dpost2gl_(part,q,sctx,&qm,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,&nxvh,
             &nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cdjpost2gl(double part[], double complex cu[],
                double complex sctx[], double qm, double dt, int nop,
                int idimp, int nx, int ny, int nxhm, int nyhm, int nxvh,
                int nyv, int ipbc) {
   djpost2gl_(part,cu,sctx,&qm,&dt,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,
              &nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void crdjpost2gl(double part[], double complex cu[],
                 double complex sctx[], double qm, double dt, double ci,
                 int nop, int idimp, int nx, int ny, int nxhm, int nyhm,
                 int nxvh, int nyv, int ipbc) {
   rdjpost2gl_(part,cu,sctx,&qm,&dt,&ci,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,
               &nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cd2jpost2gl(double part[], double complex dcu[],
                 double complex sctx[], double qm, int nop, int idimp,
                 int nx, int ny, int nxhm, int nyhm, int nxvh, 
                 int nyv) {
   d2jpost2gl_(part,dcu,sctx,&qm,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,&nxvh,
               &nyv);
   return;
}

/*--------------------------------------------------------------------*/
void crd2jpost2gl(double part[], double complex dcu[],
                  double complex sctx[], double qm, double ci, int nop,
                  int idimp, int nx, int ny, int nxhm, int nyhm,
                  int nxvh, int nyv) {
   rd2jpost2gl_(part,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,
                &nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cdsjpost2gl(double part[], double complex cu[],
                 double complex dcu[], double complex sctx[], double qm,
                 double ci, int nop, int idimp, int nx, int ny,
                 int nxhm, int nyhm, int nxvh, int nyv) {
   dsjpost2gl_(part,cu,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,
               &nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void crdsjpost2gl(double part[], double complex cu[],
                  double complex dcu[], double complex sctx[],
                  double qm, double ci, int nop, int idimp, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   rdsjpost2gl_(part,cu,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&ny,&nxhm,
                &nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cffc2initgl(double complex ffc[], double ax, double ay,
                 double affp, int nx, int ny, int nxhm, int nyhm) {
   ffc2initgl_(ffc,&ax,&ay,&affp,&nx,&ny,&nxhm,&nyhm);
   return;
}

/*--------------------------------------------------------------------*/
void cpois23gl(double complex q[], double complex fxy[], int isign,
               double complex ffc[], double ax, double ay, double affp, 
               double *we, int nx, int ny, int nxhm, int nyhm, int nxvh,
               int nyv) {
   pois23gl_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxhm,&nyhm,
             &nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void ccuperp2gl(double complex cu[], int nx, int ny, int nxhm, int nyhm,
                int nxvh, int nyv) {
   cuperp2gl_(cu,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cibpois23gl(double complex cu[], double complex bxy[],
                 double complex ffc[], double ci, double *wm, int nx,
                 int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   ibpois23gl_(cu,bxy,ffc,&ci,wm,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void camaxwel2gl(double complex exy[], double complex bxy[],
                 double complex cu[], double complex ffc[], double ci,
                 double dt, double *wf, double *wm, int nx, int ny,
                 int nxhm, int nyhm, int nxvh, int nyv) {
   amaxwel2gl_(exy,bxy,cu,ffc,&ci,&dt,wf,wm,&nx,&ny,&nxhm,&nyhm,&nxvh,
               &nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cemfield2gl(double complex fxy[], double complex exy[],
                 double complex ffc[], int isign, int nxhm, int nyhm,
                 int nxvh, int nyv) {
   emfield2gl_(fxy,exy,ffc,&isign,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cepois23gl(double complex dcu[], double complex exy[], int isign,
                double complex ffe[], double ax, double ay,
                double affp, double wp0, double ci, double *wf, int nx,
                int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   epois23gl_(dcu,exy,&isign,ffe,&ax,&ay,&affp,&wp0,&ci,wf,&nx,&ny,
              &nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cevfield23gl(double complex fxy[], double gxy[],
                  double complex sctx[], double xp, double yp, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   evfield23gl_(fxy,gxy,sctx,&xp,&yp,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmjpost2gl(double part[], double complex amu[],
                double complex sctx[], double qm, int nop, int idimp,
                int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   mjpost2gl_(part,amu,sctx,&qm,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,&nxvh,
              &nyv);
   return;
}

/*--------------------------------------------------------------------*/
void crmjpost2gl(double part[], double complex amu[],
                 double complex sctx[], double qm, double ci, int nop,
                 int idimp, int nx, int ny, int nxhm, int nyhm,
                 int nxvh, int nyv) {
   rmjpost2gl_(part,amu,sctx,&qm,&ci,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,
               &nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cdcuperp23gl(double complex dcu[], double complex amu[], int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   dcuperp23gl_(dcu,amu,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cadcuperp23gl(double complex dcu[], double complex amu[], int nx,
                   int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   adcuperp23gl_(dcu,amu,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmaxwel2gl(double complex exy[], double complex bxy[],
                double complex cu[], double complex ffc[], double ci,
                double dt, double *wf, double *wm, int nx, int ny,
                int nxhm, int nyhm, int nxvh, int nyv) {
   maxwel2gl_(exy,bxy,cu,ffc,&ci,&dt,wf,wm,&nx,&ny,&nxhm,&nyhm,&nxvh,
              &nyv);
   return;
}

/*--------------------------------------------------------------------*/
void celfield23gl(double complex q[], double complex fxy[], 
                  double complex ffc[], double *we, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv) {
   elfield23gl_(q,fxy,ffc,we,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cpsmooth2gl(double complex q[], double complex qs[],
                 double complex ffc[], int ny, int nxhm, int nyhm,
                 int nxvh, int nyv) {
   psmooth2gl_(q,qs,ffc,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cpoynt2gl(double complex q[], double complex exy[],
               double complex bxy[], double complex ffc[], double *sx,
               double *sy, double *sz, int nx, int ny, int nxhm,
               int nyhm, int nxvh, int nyv) {
   poynt2gl_(q,exy,bxy,ffc,sx,sy,sz,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void ca0rbpush23gl(double part[], double complex fxy[],
                   double complex bxy[], double complex sctx[],
                   double qbm, double dt, double dtc, double ci,
                   double *ek, int idimp, int nop, int nx, int ny,
                   int nxhm, int nyhm, int nxvh, int nyv, int ipbc) {
   a0rbpush23gl_(part,fxy,bxy,sctx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,
                 &ny,&nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}
