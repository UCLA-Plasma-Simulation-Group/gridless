/* C Library for Skeleton 2-1/2D Electromagnetic Gridless OpenMP PIC */
/* Code                                                              */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr2h_(double *part, double *vtx, double *vty, double *vtz,
              double *vdx, double *vdy, double *vdz, int *npx, int *npy,
              int *idimp, int *nop, int *nx, int *ny, int *ipbc);

void mbpush23gl_(double *part, double complex *fxy, double complex *bxy,
                 double *qbm, double *dt, double *dtc, double *ek, 
                 int *idimp, int *nop, int *nx, int *ny, int *nxhm,
                 int *nyhm, int *nxvh, int *nyv, int *ipbc);

void mrbpush23gl_(double *part, double complex *fxy,
                  double complex *bxy, double *qbm, double *dt,
                  double *dtc, double *ci, double *ek, int *idimp,
                  int *nop, int *nx, int *ny, int *nxhm, int *nyhm,
                  int *nxvh, int *nyv, int *ipbc);

void mabpush23gl_(double *part, double complex *fxy,
                  double complex *bxy,  double *qbm, double *dt,
                  double *dtc, double *ek, int *idimp, int *nop,
                  int *nx, int *ny, int *nxhm, int *nyhm, int *nxvh, 
                  int *nyv, int *ipbc);

void marbpush23gl_(double *part, double complex *fxy,
                   double complex *bxy, double *qbm, double *dt,
                   double *dtc, double *ci, double *ek, int *idimp,
                   int *nop, int *nx, int *ny, int *nxhm, int *nyhm,
                   int *nxvh, int *nyv, int *ipbc);

void mearbpush23gl_(double *part, double complex *fxy,
                    double complex *bxy, double *qbm, double *dt,
                    double *dtc, double *ci, double *ek, int *idimp,
                    int *nop, int *nx, int *ny, int *nxhm, int *nyhm,
                    int *nxvh, int *nyv, int *ipbc);

void mdpost2gl_(double *part, double complex *q, double complex *sctx,
                double *qm, int *nop, int *idimp, int *nx, int *ny,
                int *nxhm, int *nyhm, int *nxvh, int *nyv);

void mdjpost2gl_(double *part, double complex *cu, double complex *sctx,
                 double *qm, double *dt, int *nop, int *idimp, int *nx,
                 int *ny, int *nxhm, int *nyhm, int *nxvh, int *nyv,
                 int *ipbc);

void mrdjpost2gl_(double *part, double complex *cu,
                  double complex *sctx, double *qm, double *dt,
                  double *ci, int *nop, int *idimp, int *nx, int *ny,
                  int *nxhm, int *nyhm, int *nxvh, int *nyv, int *ipbc);

void md2jpost2gl_(double *part, double complex *dcu,
                  double complex *sctx, double *qm, int *nop,
                  int *idimp, int *nx, int *ny, int *nxhm, int *nyhm,
                  int *nxvh, int *nyv);

void mrd2jpost2gl_(double *part, double complex *dcu,
                   double complex *sctx, double *qm, double *ci,
                   int *nop, int *idimp, int *nx, int *ny, int *nxhm, 
                   int *nyhm, int *nxvh, int *nyv);

void mdsjpost2gl_(double *part, double complex *cu, double complex *dcu, 
                  double complex *sctx, double *qm, double *ci,
                  int *nop, int *idimp, int *nx, int *ny, int *nxhm,
                  int *nyhm, int *nxvh, int *nyv);

void mrdsjpost2gl_(double *part, double complex *cu,
                   double complex *dcu, double complex *sctx, 
                   double *qm, double *ci, int *nop, int *idimp, 
                   int *nx, int *ny, int *nxhm, int *nyhm, int *nxvh,
                   int *nyv);

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

void m1bpush23gl_(double *part, double complex *fxy,
                  double complex *bxy, double complex *sctx,
                  double *qbm, double *dt, double *dtc, double *ek,
                  int *idimp, int *nop, int *nx, int *ny, int *nxhm,
                  int *nyhm, int *nxvh, int *nyv, int *ipbc);

void m1rbpush23gl_(double *part, double complex *fxy,
                   double complex *bxy, double complex *sctx,
                   double *qbm, double *dt, double *dtc, double *ci,
                   double *ek, int *idimp, int *nop, int *nx, int *ny,
                   int *nxhm, int *nyhm, int *nxvh, int *nyv,
                   int *ipbc);

void md4jpost2gl_(double *part, double complex *dcu,
                  double complex *dcu2, double complex *sctx,
                  double *qm, int *nop, int *idimp, int *nx, int *ny,
                  int *nxhm, int *nyhm, int *nxvh, int *nyv);

void mrd4jpost2gl_(double *part, double complex *dcu,
                   double complex *dcu2, double complex *sctx,
                   double *qm, double *ci, int *nop, int *idimp,
                   int *nx, int *ny, int *nxhm, int *nyhm, int *nxvh,
                   int *nyv);

void a4maxwel2gl_(double complex *exy, double complex *bxy,
                  double complex *cu, double complex *dcu,
                  double complex *dcu2, double complex *ffc,
                  double *ci, double *dt, double *wf, double *wm,
                  int *nx, int *ny, int *nxhm, int *nyhm, int *nxvh,
                  int *nyv);

void maxwel2gl_(double complex *exy, double complex *bxy,
                double complex *cu, double complex *ffc, double *ci,
                double *dt, double *wf, double *wm, int *nx, int *ny,
                int *nxhm, int *nyhm, int *nxvh, int *nyv);

void elfield23gl_(double complex *q, double complex *fxy,
                  double complex *ffc, double *we, int *nx, int *ny,
                  int *nxhm, int *nyhm, int *nxvh, int *nyv);

void mevfield23gl_(double complex *fxy, double *gxy,
                   double complex *sctx, double *xp, double *yp,
                   int *nx, int *ny, int *nxhm, int *nyhm, int *nxvh,
                   int *nyv);

void ma0rbpush23gl_(double *part, double complex *fxy,
                    double complex *bxy, double *qbm, double *dt,
                    double *dtc, double *ci, double *ek, int *idimp,
                    int *nop, int *nx, int *ny, int *nxhm, int *nyhm,
                    int *nxvh, int *nyv, int *ipbc);

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
void cmbpush23gl(double part[], double complex fxy[],
                 double complex bxy[], double qbm, double dt,
                 double dtc, double *ek, int idimp, int nop, int nx,
                 int ny, int nxhm, int nyhm, int nxvh, int nyv,
                 int ipbc) {
   mbpush23gl_(part,fxy,bxy,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&ny,&nxhm,
               &nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmrbpush23gl(double part[], double complex fxy[],
                  double complex bxy[], double qbm, double dt,
                  double dtc, double ci, double *ek, int idimp, int nop,
                  int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv,
                  int ipbc) {
   mrbpush23gl_(part,fxy,bxy,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,&ny,
                &nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmabpush23gl(double part[], double complex fxy[],
                  double complex bxy[],  double qbm, double dt,
                  double dtc, double *ek, int idimp, int nop, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv,
                  int ipbc) {
   mabpush23gl_(part,fxy,bxy,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&ny,&nxhm,
                &nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmarbpush23gl(double part[], double complex fxy[],
                   double complex bxy[], double qbm, double dt,
                   double dtc, double ci, double *ek, int idimp,
                   int nop, int nx, int ny, int nxhm, int nyhm,
                   int nxvh, int nyv, int ipbc) {
   marbpush23gl_(part,fxy,bxy,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,&ny,
                 &nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmearbpush23gl(double part[], double complex fxy[],
                    double complex bxy[], double qbm, double dt,
                    double dtc, double ci, double *ek, int idimp,
                    int nop, int nx, int ny, int nxhm, int nyhm,
                    int nxvh, int nyv, int ipbc) {
   mearbpush23gl_(part,fxy,bxy,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,&ny,
                  &nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmdpost2gl(double part[], double complex q[],
                double complex sctx[], double qm, int nop, int idimp,
                int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   mdpost2gl_(part,q,sctx,&qm,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,&nxvh,
              &nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmdjpost2gl(double part[], double complex cu[],
                 double complex sctx[], double qm, double dt, int nop,
                 int idimp, int nx, int ny, int nxhm, int nyhm, int nxvh,
                 int nyv, int ipbc) {
   mdjpost2gl_(part,cu,sctx,&qm,&dt,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,
               &nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmrdjpost2gl(double part[], double complex cu[],
                  double complex sctx[], double qm, double dt,
                  double ci, int nop, int idimp, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv, int ipbc) {
   mrdjpost2gl_(part,cu,sctx,&qm,&dt,&ci,&nop,&idimp,&nx,&ny,&nxhm,
                &nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmd2jpost2gl(double part[], double complex dcu[],
                  double complex sctx[], double qm, int nop, int idimp,
                  int nx, int ny, int nxhm, int nyhm, int nxvh, 
                  int nyv) {
   md2jpost2gl_(part,dcu,sctx,&qm,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,&nxvh,
                &nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmrd2jpost2gl(double part[], double complex dcu[],
                   double complex sctx[], double qm, double ci, int nop,
                   int idimp, int nx, int ny, int nxhm, int nyhm,
                   int nxvh, int nyv) {
   mrd2jpost2gl_(part,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,
                 &nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmdsjpost2gl(double part[], double complex cu[],
                  double complex dcu[], double complex sctx[],
                  double qm, double ci, int nop, int idimp, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   mdsjpost2gl_(part,cu,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&ny,&nxhm,
                &nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmrdsjpost2gl(double part[], double complex cu[],
                   double complex dcu[], double complex sctx[],
                   double qm, double ci, int nop, int idimp, int nx,
                   int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   mrdsjpost2gl_(part,cu,dcu,sctx,&qm,&ci,&nop,&idimp,&nx,&ny,&nxhm,
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
void cm1bpush23gl(double part[], double complex fxy[],
                  double complex bxy[], double complex sctx[],
                  double qbm, double dt, double dtc, double *ek,
                  int idimp, int nop, int nx, int ny, int nxhm,
                  int nyhm, int nxvh, int nyv, int ipbc) {
   m1bpush23gl_(part,fxy,bxy,sctx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&ny,
               &nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cm1rbpush23gl(double part[], double complex fxy[],
                   double complex bxy[], double complex sctx[],
                   double qbm, double dt, double dtc, double ci,
                   double *ek, int idimp, int nop, int nx, int ny,
                   int nxhm, int nyhm, int nxvh, int nyv, int ipbc) {
   m1rbpush23gl_(part,fxy,bxy,sctx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,
                 &ny,&nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmd4jpost2gl(double part[], double complex dcu[], 
                  double complex dcu2[], double complex sctx[],
                  double qm, int nop, int idimp, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv) {
   md4jpost2gl_(part,dcu,dcu2,sctx,&qm,&nop,&idimp,&nx,&ny,&nxhm,&nyhm,
                &nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmrd4jpost2gl(double part[], double complex dcu[],
                   double complex dcu2[], double complex sctx[],
                   double qm, double ci, int nop, int idimp, int nx,
                   int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   mrd4jpost2gl_(part,dcu,dcu2,sctx,&qm,&ci,&nop,&idimp,&nx,&ny,&nxhm,
                 &nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void ca4maxwel2gl(double complex exy[], double complex bxy[],
                  double complex cu[], double complex dcu[],
                  double complex dcu2[], double complex ffc[],
                  double ci, double dt, double *wf, double *wm, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   a4maxwel2gl_(exy,bxy,cu,dcu,dcu2,ffc,&ci,&dt,wf,wm,&nx,&ny,&nxhm,
                &nyhm,&nxvh,& nyv);
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
void cmevfield23gl(double complex fxy[], double gxy[],
                   double complex sctx[], double xp, double yp, int nx,
                   int ny, int nxhm, int nyhm, int nxvh, int nyv) {
   mevfield23gl_(fxy,gxy,sctx,&xp,&yp,&nx,&ny,&nxhm,&nyhm,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cma0rbpush23gl(double part[], double complex fxy[],
                    double complex bxy[], double qbm, double dt,
                    double dtc, double ci, double *ek, int idimp,
                    int nop, int nx, int ny, int nxhm, int nyhm,
                    int nxvh, int nyv, int ipbc) {
   ma0rbpush23gl_(part,fxy,bxy,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,&ny,
                  &nxhm,&nyhm,&nxvh,&nyv,&ipbc);
   return;
}
