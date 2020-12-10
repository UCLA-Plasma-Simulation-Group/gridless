/* header file for mbpush1gl.c */

double ranorm();

double randum();

void cdistr1h(double part[], double vtx, double vty, double vtz,
              double vdx, double vdy, double vdz, int npx, int idimp,
              int nop, int nx, int ipbc);

void cmbpush13gl(double part[], double complex fxyz[],
                 double complex byz[], double omx, double qbm,
                 double dt, double dtc, double *ek, int idimp, int nop,
                 int nx, int nxhm, int nxvh);

void cmrbpush13gl(double part[], double complex fxyz[],
                  double complex byz[], double omx, double qbm,
                  double dt, double dtc, double ci, double *ek,
                  int idimp, int nop, int nx, int nxhm, int nxvh);

void cmabpush13gl(double part[], double complex fxyz[],
                  double complex byz[], double omx, double qbm,
                  double dt, double dtc, double *ek, int idimp, int nop,
                  int nx, int nxhm, int nxvh);

void cmarbpush13gl(double part[], double complex fxyz[],
                   double complex byz[], double omx, double qbm,
                   double dt, double dtc, double ci, double *ek,
                   int idimp, int nop, int nx, int nxhm, int nxvh);

void cmearbpush13gl(double part[], double complex fxyz[],
                    double complex byz[], double omx, double qbm,
                    double dt, double dtc, double ci, double *ek,
                    int idimp, int nop, int nx, int nxhm, int nxvh);

void cmdpost1gl(double part[], double complex q[], double qm, int nop,
                int idimp, int nx, int nxhm, int nxvh);

void cmdjpost1gl(double part[], double complex cu[], double *cux0,
                 double qm, double dt, int nop, int idimp, int nx,
                 int nxhm, int nxvh);

void cmrdjpost1gl(double part[], double complex cu[], double *cux0,
                  double qm, double dt, double ci, int nop, int idimp,
                  int nx, int nxhm, int nxvh);

void cmd2jpost1gl(double part[], double complex dcu[], double qm, 
                  int nop, int idimp, int nx, int nxhm, int nxvh);

void cmrd2jpost1gl(double part[], double complex dcu[], double qm,
                   double ci, int nop, int idimp, int nx, int nxhm,
                   int nxvh);

void cmdsjpost1gl(double part[], double complex cu[],
                  double complex dcu[], double qm, double ci, int nop,
                  int idimp, int nx, int nxhm, int nxvh);

void cmrdsjpost1gl(double part[], double complex cu[], 
                   double complex dcu[], double qm, double ci, int nop,
                   int idimp, int nx, int nxhm, int nxvh);

void cffc1initgl(double complex ffc[], double ax, double affp, int nx,
                 int nxhm);

void cpois1gl(double complex q[], double complex fx[], int isign,
              double complex ffc[], double ax, double affp, double *we,
              int nx, int nxhm, int nxvh);

void cibpois13gl(double complex cu[], double complex byz[],
                 double complex ffc[], double ci, double *wm, int nx,
                 int nxhm, int nxvh);

void camaxwel1gl(double complex eyz[], double complex byz[],
                 double complex cu[], double complex ffc[], double ci,
                 double dt, double *wf, double *wm, int nx, int nxhm,
                 int nxvh);

void cemfield1gl(double complex fxyz[], double complex fx[], 
                 double complex eyz[], double complex ffc[], int nxhm,
                 int nxvh);

void cbmfield1gl(double complex fyz[], double complex eyz[],
                 double complex ffc[], int nxhm, int nxvh);

void cepois13gl(double complex dcu[], double complex eyz[], int isign,
                 double complex ffe[], double ax, double affp,
                 double wp0, double ci, double *wf, int nx, int nxhm,
                 int nxvh);

void cmaxwel1gl(double complex eyz[], double complex byz[],
                double complex cu[], double complex ffc[], double ci,
                double dt, double *wf, double *wm, int nx, int nxhm,
                int nxvh);

void cmevfield13gl(double complex fxy[], double gxy[], double xp,
                   int nx, int nxhm, int nxvh);

void celfield1gl(double complex q[], double complex fx[],
                 double complex ffc[], double *we, int nx, int nxhm,
                 int nxvh);

void cpoynt1gl(double complex q[], double complex eyz[],
               double complex byz[], double complex ffc[], double *ex0,
               double *sx, double *sy, double *sz, int nx, int nxhm, 
               int nxvh);

void cma0rbpush13gl(double part[], double complex fxyz[],
                    double complex byz[], double omx, double qbm,
                    double dt, double dtc, double ci, double *ek,
                    int idimp, int nop, int nx, int nxhm, int nxvh);
