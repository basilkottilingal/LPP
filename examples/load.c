/**
Compile:
CC99='mpicc -std=c99' qcc -D_PTAG -D_MPI=1 -events -o a.out load.c -lm
Run:
mpirun -np 2 ./a.out
mpirun -np 4 ./a.out
mpirun -np 8 ./a.out
mpirun -np 16 ./a.out


Load test for npe() = 2,4,8,16
Particles let move freely inside a closed box;

*/

#include "run.h"
#include "../src/lpp.h"

#define LEVEL 8
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

double _time = 0.;

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

	MPI_Barrier(MPI_COMM_WORLD);
	_time = MPI_Wtime() - _time;
	
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

event end (i = 1000) {
	MPI_Barrier(MPI_COMM_WORLD);
	_time = MPI_Wtime() - _time;
	if(!pid()){
		FILE * fp = fopen("time.dat", "a");
		fprintf(fp, "%d %f\n", npe(), _time);
		fclose(fp);
	}
}

/**
#Results
~~~gnuplot 
#comment
tname=sprintf("time vs npe()");
set title tname font ",14" textcolor rgbcolor "royalblue"
set term png size 1024,1024
set output sprintf('load.png')
set xlabel 'num of procs'
set ylabel 'time'
set logscale xy
#set xrange [-5.:5.]
#set yrange [-5.:5.]
#set palette defined ( 0 "#B0B0B0", 0.333 "#FF0000", 0.666 "#0000FF", 1.0 "#000000" )
#set pointsize 1
#plot 'time.dat' u 1:2 w l lt 3 pt 7 lc rgb "red"
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2 pointtype 7 pointsize 1.5
plot 'time.dat' with linespoints linestyle 1
~~~
*/
