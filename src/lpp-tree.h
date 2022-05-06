#if  TREE
#include "lpp-mpi.h"
#include "lpp-cache.h"

Particles_cache * Pcache_temp = NULL;

Particle ** particles_cache_temp(unsigned int l){
#if _PDEBUG
	assert(l>0);
#endif
	if(!Pcache_temp) {
		Particles_cache * p = (Particles_cache *)malloc(sizeof(Particles_cache));
		p->cache = (Particle **)malloc((l+1)*sizeof(Particle *));
		Pcache_temp = p;
	}
	else{
		Particles_cache * p = Pcache_temp;
		if(l >= p->max)
			p->cache = (Particle **)realloc((void *)(p->cache), (l+1)*sizeof(Particle *));
	}
	Particles_cache * p = Pcache_temp;
	p->len = l;	
	p->max = l+1;	
	return (p->cache);
}

/**
WARNING: Using adapt_wavelet function outside 
	event adapt() will crash the program.
	And you can use max ONCE in each timestep.
	Because resitribution of particles(between processors)
	during adapt() has to be done after each adapt_wavelet()
*/
#if _MPI
scalar rparticles[]; //to store rank of pid()

void update_rparticles(){
	foreach()
		rparticles[] = pid();
}
/**
	refine, coarsen fns for rparticles[]
*/

static inline void refine_rparticles (Point point, scalar s)
{
	double pid = s[];
	foreach_child()
		s[] = pid;
}

static inline void coarsen_rparticles (Point point, scalar s)
{
	int pid = 0;
	foreach_child()
		pid += cell.pid;
#if dimension ==  2
	s[] = pid/4.;
#elif dimension == 3
	s[] = pid/8.;
#endif
}
#endif

/**
	Refine, coarsen fns for nparticles
*/
static inline void refine_nparticles (Point point, scalar s)
{
	foreach_child() s[] = 0.;

	if(!s[]) return;

	double X = x, Y = y;
#if dimension == 3
	double Z = z;
#endif

	int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS,
#if dimension == 3
	_k = 2*point.k - GHOSTS,
#endif
	_l = point.level + 1;

	foreach_particle_in_cell((int) zparticles[], (int) nparticles[]) {
		point.i = _i + (xp >= X);  point.j = _j + (yp >= Y);
#if dimension == 3
	  point.k = _k + (zp >= Z);
#endif
		point.level = _l;
		s[] += 1.;
	}

  point.level = _l - 1;
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
#if dimension == 3
  point.k = (_k + GHOSTS)/2;
#endif

}

static inline void coarsen_nparticles (Point point, scalar s)
{
	double val = 0;
	foreach_child() {
		if(cell.pid == pid())
    	val += s[];
	}
	s[] = val;
}

/**
	Refine, coarsen fns for zparticles
*/
static inline void refine_zparticles (Point point, scalar s)
{
#if _PDEBUG
	assert(is_particle_in_cell(point));
#endif

	refine_nparticles(point, nparticles);

	double zchild = s[];
	foreach_child() {
		zchild += nparticles[];
    s[] = zchild;
	}

	if(!nparticles[]) return;
	double X = x, Y = y;
#if dimension == 3
	double Z = z;
#endif

  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS,
#if dimension == 3
	    _k = 2*point.k - GHOSTS,
#endif
      _l = point.level + 1;

	int nps = (int) nparticles[], zparent = (int) s[];
	Particle ** temp = particles_cache_temp(nps);

	foreach_particle_in_cell(zparent, nps) {
		point.i = _i + (xp >= X);  point.j = _j + (yp >= Y);
#if dimension == 3
	  point.k = _k + (zp >= Z);
#endif
		point.level = _l;
		POINT_VARIABLES;
		s[] -= 1.;
		temp[(int) s[] - zparent] = particle;
	}

  point.level = --_l;
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
#if dimension == 3
  point.k = (_k + GHOSTS)/2;
#endif

	Particle ** start = particles_cache()->cache + zparent;
	memcpy(start, temp, nps*sizeof(Particle *));

#if _PDEBUG
	int check = (int) s[];
	foreach_child(){
		assert(check == ((int) s[])); 
		check += (int) nparticles[];
	}
	assert( check == ((int) (s[] + nparticles[])) );

	foreach_child()
		assert(is_particle_in_cell(point));
#endif

}

static inline void coarsen_zparticles (Point point, scalar s)
{
	coarsen_nparticles(point, nparticles);

  bool f=false;
	double val = 0;
	foreach_child()
		if (!f) 
			if(cell.pid == pid()){
    		val = s[];
				f = true;
			}
	assert(f);
	s[] = val;
#if _PDEBUG
	//foreach_child()
	assert(is_particle_in_cell(point));
#endif
}

event defaults(i=0) {
#if _MPI
	rparticles.prolongation = refine_rparticles;
	rparticles.refine       = refine_rparticles;
	rparticles.restriction  = coarsen_rparticles;
#endif

	zparticles.prolongation = no_restriction;
	zparticles.refine       = refine_zparticles;
	zparticles.restriction  = coarsen_zparticles;

	nparticles.prolongation = no_restriction;
	nparticles.refine       = no_restriction;
	nparticles.restriction  = no_restriction;
}

#if _MPI

#define foreach_proc_serial_output(fname)\
	{\
		FILE  * fp;\
		if(pid() == 0) {\
			fp = fopen(fname, "w");\
			fclose(fp); fp =  NULL;\
		}\
		for(int p=0; p<npe(); ++p){\
			if(npe() > 20) break;\
			if(p == pid()) {\
				fp = fopen(fname, "a");

#define end_foreach_proc_serial_output()\
				fclose(fp); fp = NULL;\
			}\
			MPI_Barrier(MPI_COMM_WORLD);\
		}\
	}

/**
Functions that returns respectively first and last leaf cells of this pid()
*/ 

Point first_leafpoint() {
	Index * p = tree->leaves.p;
#if dimension == 2
	Point point = { p->i, p->j, p->level};
#elif dimension == 3
	Point point = { p->i, p->j, p->k, p->level};
#endif
	return point;
}

Point last_leafpoint() {
	Index * p = tree->leaves.p + (tree->leaves.n-1);
#if dimension == 2
	Point point = { p->i, p->j, p->level};
#elif dimension == 3
	Point point = { p->i, p->j, p->k, p->level};
#endif
	return point;
}

Array * particles_to_left(){

	Array * l = array_new();
	int in, out = 0;
	Particles_cache * pcache = particles_cache();

	if(pid() < npe()-1){
		Point point = last_leafpoint();
		in = ((double) pid())< rparticles[];
		MPI_Request Req;	
		MPI_Isend (&in, 1, MPI_INT, pid()+1, 12, MPI_COMM_WORLD, &Req);
	}
	if(pid() > 0)
  	MPI_Recv (&out, 1, MPI_INT, pid()-1, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if(out) {
		Point point = first_leafpoint();
		double pid = (double) pid(), rk = rparticles[];
		assert (rk >= pid);
		int lim = 0;
		if(rk == pid)
			lim = (int) zparticles[];
		else if(rk >= (pid+1.))
			lim = pcache->len;
		else {
			Particle ** particles = pcache->cache + pcache->len, * particle;
			for(int i=pcache->len-1; i>=0; --i){
				--particles; particle = *particles;
				if(locate_coord(particle->x).level < 0) {
					lim = i + 1; break;
				}
			}
		}
		Particle ** particles = pcache->cache, * particle;
		for(int i=0; i<lim; ++i){
			particle = *particles;
#if _PDEBUG
			assert(locate_coord(particle->x).level<0);
#endif
			array_append(l, particle, sizeof(Particle));
			particle->pid = -1;
			++particles;
		}
	}

	return l;	
}

Array * particles_to_right(){
	Array * r = array_new();
	int in, out=0;
	Particles_cache * pcache = particles_cache();

	if(pid() >0) {
		Point point = first_leafpoint();
		in = (double) pid() > rparticles[];	
		MPI_Request Req;	
		MPI_Isend (&in, 1, MPI_INT, pid()-1, 21, MPI_COMM_WORLD, &Req);
	}
	if(pid() < npe()-1)
  	MPI_Recv (&out, 1, MPI_INT, pid()+1, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if(out) {
		Point point = last_leafpoint();
		double pid = (double) pid(), rk = rparticles[];
		assert (rk <= pid);
		int lim = pcache->len;
		if(rk == pid)
			lim = (int) (zparticles[]+nparticles[]);
		else if(rk <= (pid-1.))
			lim = 0;
		else {
			Particle ** particles = pcache->cache, * particle;
			for(int i=0; i<pcache->len; ++i){
				particle = *particles;
				if(locate_coord(particle->x).level < 0) {
					lim = i; break;
				}
				particles++;
			}
		}
		Particle ** particles = pcache->cache + pcache->len, * particle;
		for(int i=pcache->len-1; i>=lim; --i){
			--particles; particle = *particles;
#if _PDEBUG
			assert(locate_coord(particle->x).level<0);
#endif
			array_append(r, particle, sizeof(Particle));
			particle->pid = -1;
		}
	}

	return r;	
}

void array_accomodate(Array * a, size_t len){
  if (len >= a->max) {
    a->max = max (len, a->max + 4096);
    a->p = realloc (a->p, a->max);
  }
  a->len = len;
}

void particles_add_local(Particles * P, Array * a){
	int np = (a->len)/sizeof(Particle);
	particles_append(P, np);
	Particle * p = (Particle *) (a->p);
	p += np;
	for(int i=np-1; i>=0; --i){
		--p;
		if(locate_coord(p->x).level < 0) break;
		Particle * pnew = add_particle(P);
		memcpy(pnew, p, sizeof(Particle));
		pnew->pid = pid();
		a->len -= sizeof(Particle);	
	}
}


astats adapt_wavelet_particles (struct Adapt p)
{
	//fixme: SO UGLY. But written in a style which I can optimie in future.
	// refer to parallel/fmpi/fmpi-rebalance.h (master branch)
	/**
	update rparticles, nparticles, zparticles before adapt
	*/
	Particles * P = get_particles(); 
	update_rparticles();
	particles_cache_update(P);

#if _PDEBUG
	unsigned int nparts = P->len;
#endif

	astats st = adapt_wavelet (p); 

	/**
	Redistribute particles to (pid()-1, pid()+1) during balance()*/
	update_cache();
	Array * r = particles_to_right(), * l = particles_to_left();


	int lpid = pid()-1, rpid = pid()+1, done = false;
	MPI_Request reqr[2], reql[2];
	for(int p=0; p<npe() -1; ++p){
		//send
		if(lpid >= 0) {
  		MPI_Isend (&(l->len), 1, MPI_LONG, lpid, 123, MPI_COMM_WORLD, reql);
  		if (l->len > 0) 
    		MPI_Isend (l->p, l->len, MPI_BYTE, lpid, 123, MPI_COMM_WORLD, reql+1); 
		}

		if(rpid < npe()){
  		MPI_Isend (&(r->len), 1, MPI_LONG, rpid, 123, MPI_COMM_WORLD, reqr);
  		if (r->len > 0) 
    		MPI_Isend (r->p, r->len, MPI_BYTE, rpid, 123, MPI_COMM_WORLD, reqr+1); 
		}
		//receive	

		particles_remove_nonlocal(P);

		//fixme: MPI_Wait just after MPI_Isend makes thenon-blocking communicn ineffective. Improve
		if(lpid >= 0) {
			MPI_Wait(reql, MPI_STATUS_IGNORE);
  		if (l->len > 0) 
				MPI_Wait(reql+1, MPI_STATUS_IGNORE);
		}
		if(rpid < npe()){
			MPI_Wait(reqr, MPI_STATUS_IGNORE);
  		if (r->len > 0) 
				MPI_Wait(reqr+1, MPI_STATUS_IGNORE);
		}

		l->len = r->len = 0;
  	long len;
		if(rpid < npe()){
  		mpi_recv_check (&len, 1, MPI_LONG, rpid, 123,
		  	MPI_COMM_WORLD, MPI_STATUS_IGNORE, "return_pmpi (len)");
  		if (len > 0) {
				array_accomodate(l, len);
  			mpi_recv_check (l->p, len, MPI_BYTE, rpid, 123,
		  	     MPI_COMM_WORLD, MPI_STATUS_IGNORE, "return_pmpi (data)");
				particles_add_local(P, l);	
			}
		}
		if(lpid >= 0) {
  		mpi_recv_check (&len, 1, MPI_LONG, lpid, 123,
		  	MPI_COMM_WORLD, MPI_STATUS_IGNORE, "return_pmpi (len)");
  		if (len > 0) {
				array_accomodate(r, len);
  			mpi_recv_check (r->p, len, MPI_BYTE, lpid, 123,
		  	     MPI_COMM_WORLD, MPI_STATUS_IGNORE, "return_pmpi (data)");
				particles_add_local(P, r);	
			}
		}
		int notyet = (r->len) || (l->len);
		mpi_all_reduce(notyet, MPI_INT, MPI_LOR);
		if(!notyet){
			done = true; break;
		}
	}

	assert(done);

	array_free(r);
	array_free(l);

#if _PDEBUG
	foreach_particle(P)
		assert(locate_particle().level>=0);

	nparts -= P->len;
	mpi_all_reduce(nparts, MPI_INT, MPI_SUM);
	assert(!nparts);
#endif

	return st;
}

#else
astats adapt_wavelet_particles (struct Adapt p)
{
	astats st = adapt_wavelet (p); 
	return st;
}
#endif

#endif
