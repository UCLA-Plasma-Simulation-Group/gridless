/* header file for push2gl.c */

double ranorm();

double randum();

void cdistr2(double part[], double vtx, double vty, double vdx,
             double vdy, int npx, int npy, int idimp, int nop, int nx,
             int ny, int ipbc);

void cpush2gl(double part[], double complex fxy[],
              double complex sctx[], double qbm, double dt, double *ek,
              int idimp, int nop, int nx, int ny, int nxhm, int nyhm,
              int nxvh, int nyv, int ipbc);

void cdpost2gl(double part[], double complex q[], double complex sctx[],
               double qm, int nop, int idimp, int nx, int ny, int nxhm,
               int nyhm, int nxvh, int nyv);

void cffc2initgl(double complex ffc[], double ax, double ay,
                 double affp, int nx, int ny, int nxhm, int nyhm);

void cpois22gl(double complex q[], double complex fxy[], int isign,
               double complex ffc[], double ax, double ay, double affp, 
               double *we, int nx, int ny, int nxhm, int nyhm, int nxvh,
               int nyv);

void cevfield22gl(double complex fxy[], double gxy[],
                  double complex sctx[], double xp, double yp, int nx,
                  int ny, int nxhm, int nyhm, int nxvh, int nyv);
