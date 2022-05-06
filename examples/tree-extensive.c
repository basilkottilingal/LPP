/**
Compile:
CC99='mpicc -std=c99' qcc -D_PDEBUG -D_MPI=1 -events -o a.out tree-extensive.c -lm
Run
mpirun -np 6 ./a.out

checking lpp-tree
(adapt_wavelet_particles, coarsening/refine etc.)
*/
#include "../src/lpp-common.h"
#include "../src/lpp-mpi.h"
#include "../src/lpp-tree.h"
#include "run.h"

#include "fractions.h"

#define ParticlesMax 1000
Particles * random_particles = NULL;

int main(){
	size(1.0);
  init_grid (1 << 7);
	
	run();		
	return 0;
}

event init(i=0){
	Particles * P = particles_new(ParticlesMax);
	random_particles = P;
	for(int i=0; i<ParticlesMax; ++i){
		double xi[2] = {0.5 *(1 + noise()), 0.5 *(1 + noise())};
		if(locate_coord(xi).level < 0) continue;
		Particle * p = add_particle(P);
		memcpy(p->x, xi, sizeof(double)*dimension);
	}

	foreach_particle(P)
		assert(locate_coord(particle->x).level>= 0);

	particles_output(file = "dumpdirectory/particle-0", pmpi =1); 

}

event adapt(i++; i<9){

	scalar f[];
	double xc = 0.5 *(1 + noise()), yc = 0.5 *(1 + noise());
	double rc = 0.2;
	fraction(f, sq(x-xc)+sq(y-yc) - rc*rc);
		
	adapt_wavelet_particles({f}, (double []) {0.05}, maxlevel = 7, minlevel = 4);
	char name[100];
	sprintf(name, "dumpdirectory/particle-%d", i);
	particles_output(file = name, pmpi =1); 
}

event end(i=10){
	particles_free(random_particles);
	random_particles = NULL;
}
