/* header file for bpush2gl.c */

double ranorm_();

double randum_();

void cdistr2h(double part[], double vtx, double vty, double vtz,
              double vdx, double vdy, double vdz, int npx, int npy,
              int idimp, int nop, int nx, int ny, int ipbc);

void cbpush23gl(double part[], double complex fxy[],
                double complex bxy[], double complex sctx[], double qbm,
                double dt, double dtc, double *ek, int idimp, int nop,
                int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv,
                int ipbc);

void crbpush23gl(double part[], double complex fxy[],
                 double complex bxy[], double complex sctx[],
                 double qbm, double dt, double dtc, double ci,
                 double *ek, int idimp, int nop, int nx, int ny,
                 int nxhm, int nyhm, int nxvh, int nyv, int ipbc);

void cabpush23gl(double part[], double complex fxy[],
                 double complex bxy[], double complex sctx[],
                 double qbm, double dt, double dtc, double *ek, 
                 int idimp, int nop, int nx, int ny, int nxhm, int nyhm,
                 int nxvh, int nyv, int ipbc);

void carbpush23gl(double part[], double complex fxy[],
                  double complex bxy[], double complex sctx[],
                  double qbm, double dt, double dtc, double ci,
                  double *ek, int idimp, int nop, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv, int ipbc);

void cearbpush23gl(double part[], double complex fxy[],
                   double complex bxy[], double complex sctx[],
                   double qbm, double dt, double dtc, double ci,
                   double *ek, int idimp, int nop, int nx, int ny,
                   int nxhm, int nyhm, int nxvh, int nyv, int ipbc);

void cdpost2gl(double part[], double complex q[], double complex sctx[],
               double qm, int nop, int idimp, int nx, int ny, int nxhm,
               int nyhm, int nxvh, int nyv);

void cdjpost2gl(double part[], double complex cu[],
                double complex sctx[], double qm, double dt, int nop,
                int idimp, int nx, int ny, int nxhm, int nyhm, int nxvh,
                int nyv, int ipbc);

void crdjpost2gl(double part[], double complex cu[],
                 double complex sctx[], double qm, double dt, double ci,
                 int nop, int idimp, int nx, int ny, int nxhm, int nyhm,
                 int nxvh, int nyv, int ipbc);

void cd2jpost2gl(double part[], double complex dcu[],
                 double complex sctx[], double qm, int nop, int idimp,
                 int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv);

void crd2jpost2gl(double part[], double complex dcu[],
                  double complex sctx[], double qm, double ci, int nop,
                  int idimp, int nx, int ny, int nxhm, int nyhm,
                  int nxvh, int nyv);

void cdsjpost2gl(double part[], double complex cu[],
                 double complex dcu[], double complex sctx[], double qm,
                 double ci, int nop, int idimp, int nx, int ny,
                 int nxhm, int nyhm, int nxvh, int nyv);

void crdsjpost2gl(double part[], double complex cu[],
                  double complex dcu[], double complex sctx[],
                  double qm, double ci, int nop, int idimp, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cffc2initgl(double complex ffc[], double ax, double ay,
                 double affp, int nx, int ny, int nxhm, int nyhm);

void cpois23gl(double complex q[], double complex fxy[], int isign,
               double complex ffc[], double ax, double ay, double affp, 
               double *we, int nx, int ny, int nxhm, int nyhm, int nxvh,
               int nyv);

void ccuperp2gl(double complex cu[], int nx, int ny, int nxhm, int nyhm,
                int nxvh, int nyv);

void cibpois23gl(double complex cu[], double complex bxy[],
                 double complex ffc[], double ci, double *wm, int nx,
                 int ny, int nxhm, int nyhm, int nxvh, int nyv);

void camaxwel2gl(double complex exy[], double complex bxy[],
                 double complex cu[], double complex ffc[], double ci,
                 double dt, double *wf, double *wm, int nx, int ny,
                 int nxhm, int nyhm, int nxvh, int nyv);

void cemfield2gl(double complex fxy[], double complex exy[],
                 double complex ffc[], int isign, int nxhm, int nyhm,
                 int nxvh, int nyv);

void cepois23gl(double complex dcu[], double complex exy[], int isign,
                double complex ffe[], double ax, double ay,
                double affp, double wp0, double ci, double *wf, int nx,
                int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cevfield23gl(double complex fxy[], double gxy[],
                  double complex sctx[], double xp, double yp, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cmjpost2gl(double part[], double complex amu[],
                double complex sctx[], double qm, int nop, int idimp,
                int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv);

void crmjpost2gl(double part[], double complex amu[],
                 double complex sctx[], double qm, double ci, int nop,
                 int idimp, int nx, int ny, int nxhm, int nyhm,
                 int nxvh, int nyv);

void cdcuperp23gl(double complex dcu[], double complex amu[], int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cadcuperp23gl(double complex dcu[], double complex amu[], int nx,
                   int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cmaxwel2gl(double complex exy[], double complex bxy[],
                double complex cu[], double complex ffc[], double ci,
                double dt, double *wf, double *wm, int nx, int ny,
                int nxhm, int nyhm, int nxvh, int nyv);

void celfield23gl(double complex q[], double complex fxy[], 
                  double complex ffc[], double *we, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv);

void cpsmooth2gl(double complex q[], double complex qs[],
                 double complex ffc[], int ny, int nxhm, int nyhm,
                 int nxvh, int nyv);

void cpoynt2gl(double complex q[], double complex exy[],
               double complex bxy[], double complex ffc[], double *sx,
               double *sy, double *sz, int nx, int ny, int nxhm,
               int nyhm, int nxvh, int nyv);

void ca0rbpush23gl(double part[], double complex fxy[],
                   double complex bxy[], double complex sctx[],
                   double qbm, double dt, double dtc, double ci,
                   double *ek, int idimp, int nop, int nx, int ny,
                   int nxhm, int nyhm, int nxvh, int nyv, int ipbc);
