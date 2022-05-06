#include "../lpp-common.h"
#include "../lpp-mpi.h"
#include "../lpp-cache.h"
#include "run.h"

#define ParticlesMax 200

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


	//order particles in z-order
	particles_cache_update(P);

	//checking nparticles[], zparticles
	double np=0;
	foreach(){
		assert(np == zparticles[]);
		np += nparticles[];
	}
	//check lenghth of cache
	assert(P->len == particles_cache()->len);

	//Output
	particles_output("dumpdirectory/particle-0");

	//testing foreach_particle_cache()
	for(int p=0; p<npe(); ++p){
		MPI_Barrier(MPI_COMM_WORLD);
		if(p == pid()) {
			FILE * fp = fopen("dumpdirectory/particle-1",  pid() ? "a" : "w");
			foreach_particle_cache()
				fprintf(fp, "%f %f\n", xp, yp);
			fclose(fp);
		}
	}

	particles_free(P);
}
