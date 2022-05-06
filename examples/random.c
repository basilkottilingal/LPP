

/**
Compile:
CC99='mpicc -std=c99' qcc -D_PTAG -D_MPI=1 -events -o a.out random.c -lm
Run
mpirun -np 6 ./a.out


Particles let move freely inside a closed box;
*/

#include "run.h"
#include "../src/lpp.h"

#define LEVEL 3
#define CFL 0.1

/*
#define reflection(_u) {\
	_u[0] = (xp + _u[0]) < X0 ? (2*X0 - (2*xp + _u[0])) :  (xp + _u[0]) > (X0+L0) ? (2*(X0+L0) - (2*xp + _u[0])) : _u[0]; \
	_u[1] = (yp + _u[1]) < Y0 ? (2*Y0 - (2*yp + _u[1])) :  (yp + _u[1]) > (Y0+L0) ? (2*(Y0+L0) - (2*yp + _u[1])) : _u[1]; \
}
*/
#define reflection(_x) {\
	_x[0] = _x[0] < X0 ? (2*X0 - _x[0]) : _x[0] > (X0+L0) ? (2*(X0+L0) - _x[0]) : _x[0]; \
	_x[1] = _x[1] < Y0 ? (2*Y0 - _x[1]) : _x[1] > (Y0+L0) ? (2*(Y0+L0) - _x[1]) : _x[1]; \
}


int main() {
 	L0 = 1.;
 	DT = 1;
	init_grid(1<<LEVEL);
 	run();
}

event init (t = 0){
	Particles * _P = get_particles();  /*To get the default particles list*/
	int np = tree->leaves.n; /*number of leavf cells */
	particles_append(_P, np); /* create enoigh memory to add  one particle in each leaf cell*/
       	
	foreach(){
		Particle * p = add_particle(_P);
		p->x[0] = x+0.5*Delta*noise();
		p->x[1] = y+0.5*Delta*noise();
#if _PTAG
		p->tag = (1<<LEVEL)*(point.j-GHOSTS) + point.i-GHOSTS;
#endif
	}
}

int reflection_boundary(Particle * particle){
	particle_var();
	reflection(particle->x);
	return 1;
}

event defaults(i=0){
	particle_boundary = reflection_boundary;
}

event advance(i++){

	dt = DT;
	Particles * P = get_particles();

	double d = CFL*L0/(1<<LEVEL);
	foreach_particle(P){
		double dx[2] = {d*noise(), d*noise()};
		particle->x[0] += dx[0];
		particle->x[1] += dx[1];
	}

	particles_boundary(P);
}

event dumpf (t += dt) {
	int tint = (int) t;
	char name[100];
	sprintf(name, "particle-%d", tint);
#if _PTAG
	particles_output(file = name, ptag = 1);
#else
	particles_output(file = name, pmpi = 1);
#endif
}
event end (i = 1000);

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
        set xlabel 'x'
        set ylabel 'y'
        set xrange [0.:1.]
        set yrange [0.:1.]
        set palette defined ( 0 "#B0B0B0", 0.333 "#FF0000", 0.666 "#0000FF", 1.0 "#000000" )
	set pointsize 3
        plot fname u 1:2:3 pt 7 lt 3 lc palette z
}
~~~
*/
