#include "../lpp-common.h"
#include "../lpp-mpi.h"
#include "../lpp-tree.h"
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
		double xi[3] = {0.5 *(1 + noise()), 0.5 *(1 + noise()), 0.5 *(1+noise())};
		if(locate_coord(xi).level < 0) continue;
		Particle * p = add_particle(P);
		memcpy(p->x, xi, sizeof(double)*dimension);
	}

	foreach_particle(P)
		assert(locate_coord(particle->x).level>= 0);

	particles_output(file = "dumpdirectory/particle-0", pmpi =1); 

}

event advection(i++, i<100){
	double Dt = 0.5*L0/(1<<7);
	Particles * P = random_particles;
			
	foreach_particle(P)
		for(int d=0; d<dimension; ++d){
			double uj = Dt*noise();
			particle->x[d] += uj;
		}
	
	particles_redistribute(P);

	char name[100];
	sprintf(name, "dumpdirectory/particle-%d", i);
	particles_output(file = name, pmpi =1); 
} 

event adapt(i+=10){
	scalar f[];
	double xc = 0.5 *(1 + noise()), yc = 0.5 *(1 + noise()), zc = 0.5 *(1+noise());
	double rc = 0.2;
	fraction(f, sq(x-xc)+sq(y-yc)+sq(z-zc)- rc*rc);
		
	adapt_wavelet_particles({f}, (double []) {0.05}, maxlevel = 7, minlevel = 4);
}

event end(i=100){
	particles_free(random_particles);
	random_particles = NULL;
}
