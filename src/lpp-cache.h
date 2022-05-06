#include "lpp-common.h"

typedef struct{
	Particle ** cache;
	unsigned int len, max;
}Particles_cache;

Particles_cache * Pcache = NULL;

Particles_cache * particles_cache(){
	if(!Pcache) {
		Pcache = (Particles_cache *)malloc(sizeof(Particles_cache));
		Pcache->cache = (Particle **)malloc(10*sizeof(Particle *));
		Pcache->len = 0;	
		Pcache->max = 10;//default	
	}		
	return Pcache;
}

void particles_cache_realloc(int l){
	Particles_cache * pcache = particles_cache();
	pcache->len = 0;
	if (l < pcache->max)
		return;
	pcache->cache = (Particle **)realloc((void *)(pcache->cache), (l+1)*sizeof(Particle *));
	pcache->max = l+1;
}

void particles_cache_free(){
	if(Pcache){
		free(Pcache->cache);
		free(Pcache);
		Pcache =  NULL;
	}
}

@def foreach_particle_cache() 
	{
		Particles_cache * pcache = particles_cache();	
		Particle ** particles = pcache->cache, * particle;
		NOT_UNUSED(particle);		NOT_UNUSED(particles);
		for (int ipart=0; ipart<pcache->len; ++ipart, ++particles) {
			particle = *particles;
			particle_var()
@
@def end_foreach_particle_cache()
		}
	}
@

/**
Like all other foreach_***() macros, this macro also ..
.. cannot be used inside foreach(){}
*/
scalar nparticles[], zparticles[];
@def foreach_particle_in_cell(zp, np) 
	{
		Particles_cache * pcache = particles_cache();	
		Particle ** particles = pcache->cache + zp, * particle;
		NOT_UNUSED(particle);		NOT_UNUSED(particles);
		int npart = np;
		for (int ipart=0; ipart<npart; ++ipart, ++particles) {
			particle = *particles;
			particle_var()
@
@def end_foreach_particle_in_cell()
		}
	}
@

#if _PDEBUG
int is_particle_in_cell(Point point) {
	foreach_particle_in_cell((int) zparticles[], (int) nparticles[])
		foreach_dimension()
			if(!(xp >= (x-0.5*Delta) && xp <  (x+0.5*Delta)))
				return 0;
	return 1;
}
#endif

void particles_cache_update(Particles * P){
	
	Particles_cache * pcache = particles_cache();	
	assert(pcache);

	particles_cache_realloc(P->len);

	foreach()
		nparticles[] = 0.;

	int n= 0;
	foreach_particle(P){
assert(locate_particle().level>=0);
		n++;
		Point point = locate_particle();
		nparticles[] += 1.0;
	}

	pcache->len = n;

	double cp = 0.;
	foreach() {
		cp += nparticles[];
		zparticles[] = cp;
	}

	foreach_particle(P){
		Point point = locate_particle();
		zparticles[] -= 1.0;
		pcache->cache[(int) zparticles[]] = particle;
	}

#if _PDEBUG
	foreach()
		assert(is_particle_in_cell(point));
#endif

}
