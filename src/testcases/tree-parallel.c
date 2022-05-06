#include "../lpp-common.h"
#include "../lpp-mpi.h"
#include "../lpp-tree.h"
#include "run.h"

#define ParticlesMax 40000
Particles * random_particles = NULL;

int main(){
	size(1.0);
  init_grid (1 << 7);
	
	run();		
	return 0;
}

scalar f[];
event init(i=0){
	Particles * P = particles_new(ParticlesMax);
	for(int i=0; i<ParticlesMax; ++i){
		double xi[2] = {0.5 *(1 + noise()), 0.5 *(1 + noise())};
		if(locate_coord(xi).level < 0) continue;
		Particle * p = add_particle(P);
		memcpy(p->x, xi, sizeof(double)*dimension);
	}

	foreach_particle(P)
		assert(locate_coord(particle->x).level>= 0);

particles_output(file = "dumpdirectory/particle-0", pmpi =1); 

	foreach()
		f[] = (x<0.75) + (x>0.25);
	adapt_wavelet_particles({f}, (double []) {0.05}, maxlevel = 7, minlevel = 6);

particles_output(file = "dumpdirectory/particle-1", pmpi =1); 

	adapt_wavelet_particles({f}, (double []) {0.05}, maxlevel = 7, minlevel = 7);

particles_output(file = "dumpdirectory/particle-2", pmpi =1); 

	random_particles = P;

	particles_free(random_particles);
	random_particles = NULL;
}
