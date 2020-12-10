/* header file for push1gl.c */

double ranorm();

double randum();

void cdistr1(double part[], double vtx, double vdx, int npx, int idimp,
             int nop, int nx, int ipbc);

void cpush1gl(double part[], double complex fx[], double complex sctx[],
              double qbm, double dt, double *ek, int idimp, int nop,
              int nx, int nxhm, int nxvh);

void cdpost1gl(double part[], double complex q[], double complex sctx[],
               double qm, int nop, int idimp, int nx, int nxhm,
               int nxvh);

void cpois1gl(double complex q[], double complex fx[], int isign,
              double complex ffc[], double ax, double affp, double *we,
              int nx, int nxhm, int nxvh);

void cffc1initgl(double complex ffc[], double ax, double affp, int nx,
                 int nxhm);
