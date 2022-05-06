/**
 http://basilisk.fr/sandbox/Antoonvh/stokes-particles.h
 */

extern vector u;        //Fluid medium flow
extern face vector mu;  //Fluid medium dynamic viscosity
extern scalar rho;      //Fluid medium density
//extern coord G;                //Gravity acceleration vector (should it be external) ?
coord G;                //Gravity acceleration vector (should it be external) ?
scalar * _automatics_ = NULL;

event defaults (i= 0) {
	foreach_dimension()
		_automatics_ = list_add (_automatics_ , u.x);
	if (!is_constant (rho))
		_automatics_ = list_add (_automatics_, rho);
	if (!is_constant (mu.x)) {
		foreach_dimension()
			_automatics_ = list_add (_automatics_, mu.x);
	}
}

/*acceleration of particle*/
void particle_acceleration(Particle *  particle, double Dt){
	particle_var();
  	double muc = is_constant (mu.x) ? constant(mu.x) : interpolate (mu.x, xp, yp, zp);
  	double tau  = rp ? (muc > 0 ? rhop * sq(2*rp)/(18*muc) : HUGE) : HUGE;
  	//double tau  = rp ? (muc > 0 ? rhop * sq(2*rp)/(18*muc) : HUGE) : pa.u2.z ? pa.u2.z : HUGE;
  	double itau = tau > 0 ? 1/tau : HUGE;
  	double idt  = dt > 0 ? 1/Dt: HUGE;
  	//fixme: replace with foreach_dimension()
	up[0] = (idt*up[0] + itau*interpolate(u.x, xp, yp ,zp))/(idt + itau);
	up[1] = (idt*up[1] + itau*interpolate(u.y, xp, yp ,zp))/(idt + itau);
#if dimension ==3 
	up[2] = (idt*up[2] + itau*interpolate(u.z, xp, yp ,zp))/(idt + itau);
#endif
  	// Gravity
  	double rhof = is_constant(rho) ? constant(rho) : interpolate (rho, xp, yp, zp);
    	up[0] += rhop ? (G.x*(rhop - rhof)/rhop)/(idt + itau) : 0;
    	up[1] += rhop ? (G.y*(rhop - rhof)/rhop)/(idt + itau) : 0;
#if dimension ==3 
    	up[2] += rhop ? (G.z*(rhop - rhof)/rhop)/(idt + itau) : 0;
#endif
}

double dtprev = 0.;

event intertial_particles_step1 (i += 2, last) {
  
	dtprev = dt;
	double Dt = i > 0 ? 2*dt : dt;
	boundary (_automatics_);
	Particles * P = get_particles();
	foreach_particle(P) {
		particle_acceleration (particle, Dt);
	}
}

event intertial_particles_step2 (i = 1; i += 2) {
	double Dt = dt + dtprev;
	Particles * P = get_particles();
	foreach_particle(P)
		for(int i=0; i<dimension; ++i)
			particle->x[i] += Dt*(particle->u)[i];

	particles_boundary(P);
}

#if TWO_WAY
event defaults (i = 0) {  
	if (is_constant(a.x)) {
		a = new face vector;
		foreach_face()
			a.x[] = 0.;
		boundary ((scalar *){a});
  	}
}

void add_force (Particles * P, vector F) {
	//particle_boundary (P);
	//double rhof = constant(rho);
	boundary ((scalar*){u});
	foreach_particle(P) {
		double muc = is_constant (mu.x) ? constant(mu.x) : interpolate (mu.x, xp, yp, zp);
		double pref = 6*pi*muc*rp;
		Point point = locate_particle ();
#if dimension == 2
		coord _up = {up[0], up[1]};
#else
		coord _up = {up[0], up[1], up[2]};
#endif
		foreach_dimension()
			F.x[] += pref*(_up.x - interpolate(u.x, xp, yp, zp));
	}
}

event acceleration (i++) {
	face vector av = a;
	vector Fp[];
	foreach() 
		foreach_dimension()
			Fp.x[] = 0;
	Particles * P = get_particles();
   	add_force (P, Fp);  
  	
	double rhof = constant(rho);
  
	boundary ((scalar*){Fp});
  	foreach_face() 
    		av.x[] += face_value(Fp.x, 0)/(dv()*rhof);
}
#endif

