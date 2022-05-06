#include "../lpp-common.h"
#include "../lpp-cache.h"
#include "../lpp-tree.h"
#include "run.h"

#define LEVEL 7

#define ParticlesMax 10000
Particles * particles = NULL;

int main(){
	size(1.0);
  init_grid (1 << LEVEL);
	
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
	particles_cache_update(P);
#if _PDEBUG
	foreach()
		is_particle_in_cell(point);
#endif
	particles = P;
}

scalar f[];

event adapt(i=1){

			FILE * fp = fopen("dumpdirectory/particle-1",  "w");
			//foreach()
				foreach_particle_cache()
					fprintf(fp, "%f %f\n", xp, yp);
				//if(nparticles[])
					//fprintf(fp, "%d %d\n", (int) nparticles[], (int) zparticles[]);
			fclose(fp); fp = NULL;

	foreach()
		f[] = x>0.5;
		//f[] = sign(noise());
	adapt_wavelet({f}, (double []) {0.05}, maxlevel = LEVEL, minlevel = LEVEL-1);

			fp = fopen("dumpdirectory/particle-2",  "w");
			//foreach() {
				foreach_particle_cache()
					fprintf(fp, "%f %f\n", xp, yp);
				//if(nparticles[])
					//fprintf(fp, "%d %d\n", (int) nparticles[], (int) zparticles[]);
			fclose(fp);

	adapt_wavelet({f}, (double []) {0.05}, maxlevel = LEVEL, minlevel = LEVEL);

			fp = fopen("dumpdirectory/particle-3",  "w");
			//foreach() {
				foreach_particle_cache()
					fprintf(fp, "%f %f\n", xp, yp);
				//if(nparticles[])
					//fprintf(fp, "%d %d\n", (int) nparticles[], (int) zparticles[]);
			fclose(fp);
}

event end(i=2){
	particles_free(particles);
}

