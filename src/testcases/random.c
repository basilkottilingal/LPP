#include "../lpp-common.h"
#include "run.h"

#define ParticlesMax 100
Particles * random_particles = NULL;

int main(){
	size(1.0);
  init_grid (1 << 7);
	
	run();		
	return 0;
}

event init(i=0){
	Particles * P = particles_new(ParticlesMax);
	for(int i=0; i<ParticlesMax; ++i){
		Particle * p = add_particle(P);
		for(int d=0; d<dimension; ++d)
			p->x[d] = 0.5*L0 *(1 + noise());
	}
	random_particles = P;
}

event advance(i++){
	double Dt = 0.5*L0/(1<<7);
	Particles * P = random_particles;

	foreach_particle(P)
		for(int d=0; d<dimension; ++d)
			particle->x[d] += Dt*noise();

	particles_remove_nonlocal(P);
}

event poutput(i++){
	char name[100]; sprintf(name, "dumpdirectory/particle-%d", i);
	particles_output(name);
}

event freemem(i=500){
	particles_free(random_particles);
}

