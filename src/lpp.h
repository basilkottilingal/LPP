#include "lpp-common.h"
#include "lpp-mpi.h"
#include "lpp-tree.h"

int (* particle_boundary)(Particle * p) = NULL;

void particles_boundary(Particles * P) {
	if(particle_boundary != NULL)
		foreach_particle(P)
			particle_boundary(particle);
#if _MPI	
	particles_redistribute(P);
#else
	particles_remove_nonlocal(P);
#endif
}

#if _PINERTIAL
#include "lpp-inertial.h"

/**
For advancing particle with lagrangian velocity vector, particle->u.
Explicit advection.
*/
static inline void particle_advance_lagrangian(Particle * particle){
	/**
	Assume particle velocity is updated in
	particle->u[];
	*/
	for(int i=0; i<dimension; ++i)
		particle->x[i] += dt*(particle->u)[i];
}

void particles_advance_lagrangian(Particles * P){
	foreach_particle(P)
		particle_advance_lagrangian(particle);
	particles_boundary(P);
}

#endif

//include stability condition here
/**
For advancing particle with eulerian velocity vector, u.
Explicit advection.
*/
static inline void particle_advance_eulerian(Particle * particle, vector u){
	particle_var()
	Point point = locate_particle();
	assert(point.level >= 0);
	double * q = particle->x;
	q[0] += dt*interpolate_linear(point, (struct _interpolate){u.x, xp, yp, zp});
	q[1] += dt*interpolate_linear(point, (struct _interpolate){u.y, xp, yp, zp});
#if dimension == 3
	q[2] += dt*interpolate_linear(point, (struct _interpolate){u.z, xp, yp, zp});
#endif
}

void particles_advance_eulerian(Particles * P, vector u){
	foreach_particle(P)
		particle_advance_eulerian(particle, u);
	particles_boundary(P);
}

