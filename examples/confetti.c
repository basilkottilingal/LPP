/**
Compile:
CC99='mpicc -std=c99' qcc -D_PTAG -D_MPI=1 -events -o a.out confetti.c -lm
Run
mpirun -np 6 ./a.out

Copied from /sandbox/Antoonvh/
http://basilisk.fr/sandbox/Antoonvh/confetti.c
*/

#include "navier-stokes/centered.h"
#include "../src/lpp.h"

#define POSX (R1*sin(t/6))
#define POSY (R1*cos(t/6))

double R1 = 1, R2 = 0.2;


int main() {
  const face vector muc[] = {1./500., 1./500.};
  mu = muc;
  L0 = 10.;
  X0 = Y0 = -L0/2.;
  DT = 0.1;
  run();
}

event init (t = 0){
	Particles * P = get_particles();  /*To get the default particles list*/
	/*Alternative
	 *  Particles * P = _particles;
	 */
	particles_append(P, 100); /*Memory for 100 particles*/
	for(int p=0; p<100; ++p){
		double th = pi*noise();
		double r = 0.49*L0*noise();
		if(locate( (r*cos(th)) , (r*sin(th)) ).level<0) continue;
		Particle * p = add_particle(P);
		p->x[0] = r*cos(th);
		p->x[1] = r*sin(th);
#if _PTAG
		p->tag = 0.;
#endif
	}
}

event advance(i++){

	double dt = DT;
	Particles * P = get_particles();

	particles_advance_eulerian(P, u);
}

event forcing (i++) {
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
	Particles * P = get_particles();
	int nadd = 30;
	particles_append(P, nadd); /* Creates memory to accomodate additional 30 particles*/
	for(int p=0; p<1000; ++p){
		double xp = 0.5*R2*noise() + POSX;
		double yp = 0.5*R2*noise() + POSY;
		if(locate(xp, yp).level<0) continue; /*Doesn't add the particle if the coordinate doesn't belong to pid()*/
		Particle * p = add_particle(P);
		p->x[0] = xp;
		p->x[1] = yp;
#if _PTAG
		p->tag = (long) t;
#endif
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
	sprintf(name, "particle-%d", tint);
#if _PTAG
	particles_output(file = name, ptag = 1);
#else
	particles_output(file = name, pmpi = 1);
#endif
}
event stop (t = 30);





/**
#Results
~~~gnuplot 
list=system('ls -1B particle-*')
do for [fname in list] {
        cmd=sprintf("echo %s | sed 's/particle-//g'",fname)
        number=int(system(cmd))
        tname=sprintf("confetti, %d", number);
        set title tname font ",14" textcolor rgbcolor "royalblue"
        set term png size 1024,1024
        set output sprintf('img-%d.png',number)
        set xlabel 'x'
        set ylabel 'y'
        set xrange [-5.:5.]
        set yrange [-5.:5.]
        set palette defined ( 0 "#B0B0B0", 0.333 "#FF0000", 0.666 "#0000FF", 1.0 "#000000" )
	set pointsize 1
        plot fname u 1:2:3 pt 7 lt 3 lc palette z
}
~~~
*/
