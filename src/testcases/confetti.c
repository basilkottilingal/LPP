/**
Copied from /sandbox/Antoonvh/
http://basilisk.fr/sandbox/Antoonvh/confetti.c
*/

#include "navier-stokes/centered.h"
#include "../lpp.h"

#define POSX (R1*sin(t/6))
#define POSY (R1*cos(t/6))

double R1 = 1, R2 = 0.2;

Particles * initial;

int main() {
  const face vector muc[] = {1./500., 1./500.};
  mu = muc;
  L0 = 10.;
  X0 = Y0 = -L0/2.;
  DT = 0.1;
  run();
}

event init (t = 0){
	Particles * P = particles_new(100);
	for(int p=0; p<100; ++p){
		double th = pi*noise();
		double r = 0.49*L0*noise();
		if(locate( (r*cos(th)) , (r*sin(th)) ).level<0) continue;
		Particle * p = add_particle(P);
		p->x[0] = r*cos(th);
		p->x[1] = r*sin(th);
	}
	initial = P;
}

event advance(i++){

	double dt = DT;
	Particles * P = initial;

	particles_advance(P, u);
}

event forcing (i++, t+=dt) {
  double angle = t/2.;
  double U = 1;
  foreach()
    if (sq(x - POSX) + sq(y - POSY) < sq(R2)) {
      u.x[] = U*sin(angle);
      u.y[] = U*cos(angle);
    }
}

event add (t += 2){
  //init_tp_square (xm = POSX, ym = POSY, l = R2);
	Particles * P = initial;
	int nadd = 30;
	particles_append(P, nadd);
	for(int p=0; p<1000; ++p){
		double xp = 0.5*R2*noise() + POSX;
		double yp = 0.5*R2*noise() + POSY;
		if(locate(xp, yp).level<0) continue;
		Particle * p = add_particle(P);
		p->x[0] = xp;
		p->x[1] = yp;
		if(nadd--) break;
	}
}

event adapt (i++)
  adapt_wavelet_particles ({u.x, u.y}, (double[]){0.1, 0.1}, 8, 6);
/*
event mov (t += 0.3) {
  foreach_P_in_list(tracer_particles) {
    scatter (P, pc = {sin(P), cos(P), sin(P*2.4)});
  }
  box();
  save ("movp.mp4");
  
  pstats ps = statsp (tracer_particles[1]);
  if (pid() == 0) 
    printf ("%g %g %g %g\n",t,
	    sqrt(sq(ps.stddev.x) + sq(ps.stddev.y)), //length of stddev vector
	    ps.avg.x, ps.avg.y);
}
*/
event dumpf (t += 1) {
	int tint = (int) t;
	char name[100];
	sprintf(name, "dumpdirectory/particle-%d", tint);
	particles_output(file = name, pmpi = 1);
}
event stop (t = 200);
