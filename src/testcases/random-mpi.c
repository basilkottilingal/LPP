#include "../lpp-common.h"
#include "../lpp-mpi.h"
#include "run.h"

#define ParticlesMax 200
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
		double xi[2] = {0.5 *(1 + noise()), 0.5 *(1 + noise())};
		if(locate_coord(xi).level < 0) continue;
		Particle * p = add_particle(P);
		memcpy(p->x, xi, sizeof(double)*dimension);
	}

	foreach_particle(P)
		assert(locate_coord(particle->x).level>= 0);

	random_particles = P;
}

//int onlyonce = 0;
event advance(i++; i<=49){

	double Dt = 0.5*L0/(1<<7);
	Particles * P = random_particles;
			
	foreach_particle(P)
		for(int d=0; d<dimension; ++d){
			double uj = Dt*noise();
			particle->x[d] += uj;
		}
	
	particles_redistribute(P);

}

event output(i++){
	char file[100];
	sprintf(file, "dumpdirectory/particle-%d", i);
	if(random_particles)
		particles_output(file);
}

event freemem(i=50){
	particles_free(random_particles);
	random_particles = NULL;
}

