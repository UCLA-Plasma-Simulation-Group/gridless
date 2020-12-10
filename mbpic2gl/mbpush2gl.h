/* header file for mbpush2gl.c */

double ranorm_();

double randum_();

void cdistr2h(double part[], double vtx, double vty, double vtz,
              double vdx, double vdy, double vdz, int npx, int npy,
              int idimp, int nop, int nx, int ny, int ipbc);

void cmbpush23gl(double part[], double complex fxy[],
                 double complex bxy[],  double qbm, double dt,
                 double dtc, double *ek, int idimp, int nop, int nx,
                 int ny, int nxhm, int nyhm, int nxvh, int nyv,
                 int ipbc);

void cmrbpush23gl(double part[], double complex fxy[],
                  double complex bxy[], double qbm, double dt,
                  double dtc, double ci, double *ek, int idimp, int nop,
                  int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv,
                  int ipbc);

void cmabpush23gl(double part[], double complex fxy[],
                  double complex bxy[],  double qbm, double dt,
                  double dtc, double *ek, int idimp, int nop, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv,
                  int ipbc);

void cmarbpush23gl(double part[], double complex fxy[],
                   double complex bxy[], double qbm, double dt,
                   double dtc, double ci, double *ek, int idimp,
                   int nop, int nx, int ny, int nxhm, int nyhm,
                   int nxvh, int nyv, int ipbc);

void cmearbpush23gl(double part[], double complex fxy[],
                    double complex bxy[], double qbm, double dt,
                    double dtc, double ci, double *ek, int idimp,
                    int nop, int nx, int ny, int nxhm, int nyhm,
                    int nxvh, int nyv, int ipbc);

void cmdpost2gl(double part[], double complex q[], 
                double complex sctx[], double qm, int nop, int idimp,
                int nx, int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cmdjpost2gl(double part[], double complex cu[],
                 double complex sctx[], double qm, double dt, int nop,
                 int idimp, int nx, int ny, int nxhm, int nyhm,
                 int nxvh, int nyv, int ipbc);

void cmrdjpost2gl(double part[], double complex cu[],
                  double complex sctx[], double qm, double dt,
                  double ci, int nop, int idimp, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv, int ipbc);

void cmd2jpost2gl(double part[], double complex dcu[],
                  double complex sctx[], double qm, int nop, int idimp,
                  int nx, int ny, int nxhm, int nyhm, int nxvh,
                  int nyv);

void cmrd2jpost2gl(double part[], double complex dcu[],
                   double complex sctx[], double qm, double ci, int nop,
                   int idimp, int nx, int ny, int nxhm, int nyhm,
                   int nxvh, int nyv);

void cmdsjpost2gl(double part[], double complex cu[],
                  double complex dcu[], double complex sctx[],
                  double qm, double ci, int nop, int idimp, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cmrdsjpost2gl(double part[], double complex cu[],
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

void cm1bpush23gl(double part[], double complex fxy[],
                  double complex bxy[], double complex sctx[],
                  double qbm, double dt, double dtc, double *ek,
                  int idimp, int nop, int nx, int ny, int nxhm,
                  int nyhm, int nxvh, int nyv, int ipbc);

void cm1rbpush23gl(double part[], double complex fxy[],
                   double complex bxy[], double complex sctx[],
                   double qbm, double dt, double dtc, double ci,
                   double *ek, int idimp, int nop, int nx, int ny,
                   int nxhm, int nyhm, int nxvh, int nyv, int ipbc);

void cmd4jpost2gl(double part[], double complex dcu[], 
                  double complex dcu2[], double complex sctx[],
                  double qm, int nop, int idimp, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv);

void cmrd4jpost2gl(double part[], double complex dcu[],
                   double complex dcu2[], double complex sctx[],
                   double qm, double ci, int nop, int idimp, int nx,
                   int ny, int nxhm, int nyhm, int nxvh, int nyv);

void ca4maxwel2gl(double complex exy[], double complex bxy[],
                  double complex cu[], double complex dcu[],
                  double complex dcu2[], double complex ffc[],
                  double ci, double dt, double *wf, double *wm, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cmaxwel2gl(double complex exy[], double complex bxy[],
                double complex cu[], double complex ffc[], double ci,
                double dt, double *wf, double *wm, int nx, int ny,
                int nxhm, int nyhm, int nxvh, int nyv);

void celfield23gl(double complex q[], double complex fxy[],
                  double complex ffc[], double *we, int nx, int ny,
                  int nxhm, int nyhm, int nxvh, int nyv);

void cmevfield23gl(double complex fxy[], double gxy[],
                   double complex sctx[], double xp, double yp, int nx,
                   int ny, int nxhm, int nyhm, int nxvh, int nyv);

void cma0rbpush23gl(double part[], double complex fxy[],
                    double complex bxy[], double qbm, double dt,
                    double dtc, double ci, double *ek, int idimp,
                    int nop, int nx, int ny, int nxhm, int nyhm,
                    int nxvh, int nyv, int ipbc);
