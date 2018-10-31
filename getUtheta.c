/* Title: Generalized Newtonian Fluid
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "embed.h"
#include "navier-stokes/centered.h"

char filename[80];
scalar fplot[];

event init(t = 0)
{
  restore (file = filename);
  boundary(all);
  scalar utheta[];
  foreach() {
    double theta = atan2(y, x); //, r = sqrt(x*x + y*y);
    utheta[] = - sin(theta)*u.x[] + cos(theta)*u.y[];
  }
  boundary ((scalar *){utheta});

  FILE * fp = ferr;
  int LEVEL = 6;
  for (double x = 0.25+1./pow(2.,LEVEL); x < 0.5; x += 1./pow(2.,LEVEL)){
    fprintf(ferr, "%g %g\n", x, interpolate(utheta, x, 0.0));
  }
  fflush (fp);
  fclose (fp);
}

int main(int a, char const *arguments[])
{
  u.n[embed] = dirichlet_embed (x*x + y*y > 0.14 ? 0. : - y);
  u.t[embed] = dirichlet_embed (x*x + y*y > 0.14 ? 0. :   x);
  sprintf (filename, "%s", arguments[1]);
  run();
}
