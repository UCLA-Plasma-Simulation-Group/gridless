/* header file for mpush1gl.c */

double ranorm();

double randum();

void cdistr1(double part[], double vtx, double vdx, int npx, int idimp,
             int nop, int nx, int ipbc);

void cmpush1gl(double part[], double complex fx[], double qbm,
               double dt, double *ek, int idimp, int nop, int nx,
               int nxhm, int nxvh);

void cmdpost1gl(double part[], double complex q[], double qm, int nop,
                int idimp, int nx, int nxhm, int nxvh);

void cffc1initgl(double complex ffc[], double ax, double affp, int nx,
                 int nxhm);

void cpois1gl(double complex q[], double complex fx[], int isign,
              double complex ffc[], double ax, double affp, double *we,
              int nx, int nxhm, int nxvh);

void cmefield1gl(double complex fx[], double *gx, double xp, int nx,
                 int nxhm, int nxvh);
