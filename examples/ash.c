/**
Compile:
CC99='mpicc -std=c99' qcc -D_PINERTIAL -D_MPI=1 -events -o a.out ash.c -lm
Run:
mpirun -np 6 ./a.out


http://basilisk.fr/sandbox/Antoonvh/ash.c
*/


/**
![Turbudity currents are driven by particles. Image via [Octavio E. Sequeiros et al. (2010)](https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0000200)](https://ascelibrary.org/cms/attachment/a5c1fe40-ea5f-4d00-a9d1-3a8694cd7137/7.gif)

# Vulcanic ash settling in water

We use 64000 spherical sand particles to model the descend of ganules
in water.

![You can "see" the two-way coupling in action](ash/locs.mp4)

The movie displays a surprisingly pleasing bview bug as is renders only part of the grid slice...
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define TWO_WAY 1
#include "../src/lpp.h"
//#include "view.h"

u.t[bottom] = dirichlet (0.);
u.r[bottom] = dirichlet (0.);


#define MUZ 1e-3

int main() {
  const face vector muc[] = {MUZ, MUZ, MUZ}; //mu Water
  periodic (left);
  mu = muc;
  G.y = -9.81;
  L0 = 0.1;
  X0 = Y0 = Z0 = -L0/2; 
  N = 64;
  run();
}
/**
Initially, the flow is steady and seeded with particles near the top
of the domain.
 */

event init (t = 0) {
  
	const scalar rhof[] = 1000;
  	rho = rhof;
  	int pni = 40;
  
	Particles * P = get_particles();
	particles_append(P, pni*pni*pni); /* create enoigh memory to add 40x40x40 particles*/
	int n = 0;
	for (int i = 0; i < pni; i++) 
		for (int j = 0; j < pni; j++)  
			for (int k = 0; k < pni; k++) {
	coord _p = { (double) i/(100*pni),  (double) j/(100*pni) + L0/4.,  (double) k/(100*pni)};
	if(locate(_p.x, _p.y, _p.z).level >= 0){
		Particle * p = add_particle(P);
		p->x[0] = _p.x;
		p->x[1] = _p.y;
		p->x[2] = _p.z;
		p->tag = i*pni*pni + j*pni + k;
		p->rho = 2000; //approx denisty of silicon
		p->r = 0.05e-3 + 0.01e-3*noise();
	}
			}
}
/**
## Bottom boundary

Particles bounce at the bottom boundary. They should be able to rest
 on there.
*/
double damp = 2;
int  bottom_boundary (Particle * particle) {
	particle_var()
	if (yp < Y0) {
		particle->x[1] = Y0 + Y0 - yp;
		//foreach_dimension()
		//p().u.y /= damp;
      		particle->u[0] /= damp;
		particle->u[1] *= (-1/damp);  
      		particle->u[2] /= damp;
	}
	return 1;
}

event defaults(i=0){
	particle_boundary = bottom_boundary;
}

event adapt (i++) 
	adapt_wavelet_particles ((scalar*){u}, (double[]){0.01, 0.01, 0.01}, 7, 4);

int nout = 0;
event mov (t += .02) {
	char name[100];
	sprintf(name, "particle-%d", nout);
#if _PTAG
	particles_output(file = name, ptag = 1);
#else
	particles_output(file = name, pmpi = 1);
#endif
	++nout;
}
/*
event mov (t += .02) {
  view (theta = t/10. - 0.5);
  scatter (sand, pc = {1, 1, 1});
  box();
  translate (z = -L0/2)
    cells();
  save ("locs.mp4");
}
*/

event end (t = 10);

/**
#Results
~~~gnuplot 
list=system('ls -1B particle-*')
do for [fname in list] {
        cmd=sprintf("echo %s | sed 's/particle-//g'",fname)
        number=int(system(cmd))
        tname=sprintf("randomly moving particles, %d", number);
        set title tname font ",14" textcolor rgbcolor "royalblue"
        set term png size 1024,1024
        set output sprintf('img-%d.png',number)
	set view 60,75
	set xlabel 'x'
	set ylabel 'z'
	set zlabel 'y'
	set xrange [-0.05:0.05]
	set yrange [-0.05:0.05]
	set zrange [-0.05:0.05]
	set palette defined ( 0 "#B0B0B0", 0.333 "#FF0000", 0.666 "#0000FF", 1.0 "#000000" )
	set pointsize 1
	splot fname u 1:3:2:4 pt 7 lt 2 lc palette z
}
~~~
*/
