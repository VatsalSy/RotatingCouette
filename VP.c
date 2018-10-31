/* Title: Generalized Newtonian Fluid
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
char filename[80];
double tauy,mu_0,mumax;
double n;
#define imax (1e6)
int LEVEL;
#define error (1e-10)

int main()
{
  stokes = true;
  TOLERANCE = 1e-5;
  // The regularisation value of viscosity
  mumax=1e5;
  for (int counter = 0; counter < 4; counter++){
    switch (counter) {
    case 0: /* Newtonian */
      mu_0 = 1.; tauy = 0.; n = 1.; DT = (1e-2); LEVEL = 6;
      sprintf (filename, "%s", "lastNewt");
      break;
    case 1: /* Power-law (shear-thinning) n = 0.5 */
      mu_0 = 0.0625; tauy = 0.; n = 0.5; DT = (1e-5); LEVEL = 6;
      sprintf (filename, "%s", "lastShThn");
      break;
    case 2: /* Power-law (shear-thinning) n = 0.4 */
      mu_0 = 0.0625; tauy = 0.; n = 0.38; DT = (1e-5); LEVEL = 6;
      sprintf (filename, "%s", "lastShThn2");
      break;
    case 3: /* Bingham */
      mu_0 = 1.; tauy = 10.; n = 1.0; DT = (1e-6); LEVEL = 7;
      sprintf (filename, "%s", "lastBing");
      break;
    }
    fprintf(ferr, "mu0 = %g, tauy = %g, n = %g\n",mu_0, tauy, n);
    init_grid (1<<LEVEL);
    L0 = (1.0 + 8.0*1.0/(pow(2,LEVEL)));
    origin (-L0/2., -L0/2.);
    CFL = 0.05;
    run();
  }
}

/* fplot and D2plot save the embeded solid and the second principal invariant
  of the shear strain rate tensor for post-processing purposes.
*/
scalar un[], utheta[];
face vector muv[];

// initialization event
event init (t = 0) {
  // preparing viscosity to be used as Non-Newtonian fluid
  mu = muv;
  /*
  The geometry of the embedded boundary is defined as two cylinders of
  radii 0.5 and 0.25. */
  vertex scalar phi[];
  foreach_vertex()
    phi[] = difference (sq(0.5) - sq(x) - sq(y),
      sq(0.25) - sq(x) - sq(y));
  fractions (phi, cs, fs);
  /**
  The outer cylinder is fixed and the inner cylinder is rotating with
  an angular velocity unity. */
  u.n[embed] = dirichlet_embed (x*x + y*y > 0.14 ? 0. : - y);
  u.t[embed] = dirichlet_embed (x*x + y*y > 0.14 ? 0. :   x);
  /*
  We initialize the reference velocity field. This is used for checking for stationary solutions
  */
  foreach() {
    double theta = atan2(y, x); //, r = sqrt(x*x + y*y);
    utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
    un[] = utheta[];
  }
  boundary ((scalar *){utheta, un});
  dump (file = "start");
  boundary(all);
}

/**
We look for a stationary solution. */
event convergence (i++; i<= imax){
  static FILE * fp;
  foreach() {
    double theta = atan2(y, x); //, r = sqrt(x*x + y*y);
    utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
  }
  boundary ((scalar *){utheta});
  double du = change (utheta, un);
  fprintf(ferr, "i = %d: dt = %g, err = %g\n", i, dt, du);
  if (i > 100 && du < error){
    fp = fopen ("log", "a");
    fprintf (fp, "i = %d: dt = %g, err = %g\n", i, dt, du);
    fclose(fp);
    dump (file = filename);
    return 1; /* stop */
  }
  if (i == imax){
    fp = fopen ("log", "a");
    fprintf (fp, "i = %d: dt = %g, err = %g\n", i, dt, du);
    fclose(fp);
    dump (file = filename);
  }
}

/**
## Implementation of generalized Newtonian viscosity
*/

event properties(i++) {
  /*
   $$D_{11} = \frac{\partial u}{\partial x}$$
   $$D_{12} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{21} = \frac{1}{2}\left( \frac{\partial u}{\partial y}+ \frac{\partial v}{\partial x}\right)$$
   $$D_{22} = \frac{\partial v}{\partial y}$$
   The second invariant is $D_2=\sqrt{D_{ij}D_{ij}}$
   $$D_2^2= D_{ij}D_{ij}= D_{11}D_{11} + D_{12}D_{21} + D_{21}D_{12} + D_{22}D_{22}$$
   the equivalent viscosity is
   $$\mu_{eq}= \mu_0\left(\frac{D_2}{\sqrt{2}}\right)^{N-1} + \frac{\tau_y}{\sqrt{2} D_2 }$$
   **Note:** $\|D\| = D_2/\sqrt{2}$
   Finally, mu is the min of of $\mu_{eq}$ and a large $\mu_{max}$, then the fluid flows always, it is not a solid, but a very viscous fluid.
   $$ \mu = min\left(\mu_{eq}, \mu_{max}\right) $$
  */
  double muTemp = mu_0;
  foreach_face() {
    double D11 = (u.x[] - u.x[-1,0]);
    double D22 = ((u.y[0,1]-u.y[0,-1])+(u.y[-1,1]-u.y[-1,-1]))/4.0;
    double D12 = 0.5*(((u.x[0,1]-u.x[0,-1])+(u.x[-1,1]-u.x[-1,-1]))/4.0 + (u.y[] - u.y[-1,0]));
    double D2 = sqrt(sq(D11)+sq(D22)+2.0*sq(D12))/(Delta);

    if (D2 > 0.0) {
      double temp = tauy/(sqrt(2.0)*D2) + mu_0*exp((n-1.0)*log(D2/sqrt(2.0)));
      muTemp = min(temp, mumax);
    } else {
      if (tauy > 0.0 || n < 1.0){
        muTemp = mumax;
      } else {
        muTemp = (n == 1.0 ? mu_0 : 0.0);
      }
    }
    muv.x[] = fm.x[]*(muTemp);
  }
  boundary ((scalar *){muv});
}
