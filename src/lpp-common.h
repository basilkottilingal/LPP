#if _PINERTIAL
#ifndef _PTAG
#define _PTAG 1
#endif
#endif
/**
	Particle or the lagrangian particle
*/ 
typedef struct{
	double x[dimension];
#if _MPI
	int pid;
#endif
#if _PINERTIAL
	double u[dimension];
	double rho, r;
#endif
#if _PTAG
	long tag;
#endif
}Particle;

typedef struct {
	Particle * p;
	unsigned int len, max;
}Particles;

@def particle_var() 
			double xp = particle->x[0]; NOT_UNUSED(xp);
			double yp = particle->x[1]; NOT_UNUSED(yp);
#if dimension == 2	
			double zp = 0;              NOT_UNUSED(zp);
#else
			double zp = particle->x[2]; NOT_UNUSED(zp);
#endif
#if _PINERTIAL
			double rhop = particle->rho; NOT_UNUSED(rhop);
			double rp = particle->r;     NOT_UNUSED(rp);
			double * up = particle->u;   NOT_UNUSED(up); 
#endif
@

/**
foreach_particle(fr) is macro that iterates through each marker particles of the front*/
@def foreach_particle_all(P) 
	{
		Particle * particle = P->p;
		int len = P->len;
		NOT_UNUSED(particle);
		for (int ipart=0; ipart<len; ++ipart, ++particle) {
			particle_var()
@
@def end_foreach_particle_all()
		}
	}
@

/**
To see if the particle is local */
#if _MPI
#define local_particle() (pid()==particle->pid)
#endif

@def foreach_particle(P) 
	foreach_particle_all(P)
#if _MPI
		if(local_particle()) {
#endif
@
@def end_foreach_particle()
#if _MPI
		}
#endif
	end_foreach_particle_all()
@


#if dimension == 2
		//fixme: fixfperiodic
#define distance(x, y) sqrt((x[0]-y[0])*(x[0]-y[0]) \
                       + (x[1]-y[1])*(x[1]-y[1]))
#endif


/**
Default
*/
Particles * _particles;

Particles * get_particles(){
	assert(_particles);
	return _particles;
}

Particles * particles_new(int l){
		
	Particles * p = (Particles *)malloc(sizeof(Particles));
	p->p   = (Particle *)malloc((l+1)*sizeof(Particle));
	p->len = 0;
	p->max = l+1;

_particles = p;	

	return(p);
}

event defaults(i=0){
	particles_new(1000);
}

Particle * add_particle(Particles * P){
	Particle * p = &(P->p[P->len++]);
#if _MPI
	p->pid = pid();
#endif
	assert(P->len < P->max);
	return p;
}

//copy from Array
void particles_shrink(Particles * p){
/*
	p->particle = realloc(p->particle, p->len*sizeof(Particle));
	p->max = p-> 
*/
}

void particles_append(Particles * P, unsigned int r){
	//r : extra requirement
	unsigned int l = r + P->len;
	if(l < P->max)
		return;
	P->p = (Particle *)realloc((void *)P->p, (l+1)*sizeof(Particle));
	P->max = l+1;
}

/*
Particles * get_particles(){
#if dimension==2
	if (fr_circle == NULL)
		fr_circle = front_new();
	return fr_circle;
#elif dimension==3		
	if (fr_sphere == NULL)
		fr_sphere = front_new();
	return fr_sphere;
#endif
}
*/
#if dimension == 2
#define locate_coord(xp)       locate(xp[0], xp[1])
#define locate_particle()      locate(xp, yp)
#elif dimension == 3
#define locate_coord(xp)       locate(xp[0], xp[1], xp[2])
#define locate_particle()      locate(xp, yp, zp)
#endif

void particles_remove_nonlocal(Particles * P){
	//Particles * particles = particles_get();
	Particle * pa = P->p, * pb;
	int len = P->len, i = 0;
	bool f = true;

	for(; i<len; ++i, ++pa)
#if _MPI
		if(pa->pid != pid())
#else
		if(locate_coord(pa->x).level<0)
#endif
		{
			f = false;
			break;
		}

	if(f) return; //All particles belong to this processor domain

	pb = pa+1;

	for(++i; i<len; ++i, pb++)
#if _MPI
		if(pb->pid == pid())
#else
		if(locate_coord(pb->x).level>=0)
#endif
		{
			memcpy(pa, pb, sizeof(Particle)); 
			++pa;
		}
	P->len = (int) (pa - P->p);
}

/**
Free memory*/
void particles_free(Particles * P){
	assert(P != NULL);
	free(P->p);
	free(P);
	P = NULL;	
}

int locate_rank (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = { .level = l };
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + GHOSTS;
#if dimension >= 2
    point.j = (p.y - Y0)/L0*n + GHOSTS;
#endif
#if dimension >= 3
    point.k = (p.z - Z0)/L0*n + GHOSTS;
#endif
    if (point.i >= GHOSTS && point.i < n + GHOSTS
#if dimension >= 2
	&& point.j >= GHOSTS && point.j < n + GHOSTS
#endif
#if dimension >= 3
	&& point.k >= GHOSTS && point.k < n + GHOSTS
#endif
	) {
      if (allocated(0) && is_leaf(cell))
	return cell.pid;
    }
    else
			return -1;
			/*This arises when Point is out of Computational box*/
  }
	fprintf(stdout, "\nFAILED %f %f", p.x, p.y); fflush(stdout);
	assert(false);
	/*This assert arises when CFL of moving points > 1.0*/

	return 1;
}

#if dimension == 2
#define locate_particle_rank() locate_rank(xp, yp)
#else
#define locate_particle_rank() locate_rank(xp, yp, zp)
#endif

/**
For advancing particle.
Explicit advection.
*/

struct PDump{
	char * file;
	int pmpi, ptag;
};

#if _PTAG
#define particle_color() ( ptag ? (double) particle->tag : 0. )
#elif _MPI
#define particle_color() ( pmpi ? (double) pid() : 0. )
#else
#define particle_color()  0.
#endif


void particles_output(struct PDump p){

//if(pid() == 1) assert(0);
	//fixme : textfile
	Particles * P = get_particles();
	char dfile[] = "dumpdirectory/particle-default",
	   * file = p.file ? p.file : dfile;
	bool pmpi = p.pmpi ? 1 : 0, ptag = p.ptag ? 1 : 0;

	int nlocal = 0;
	foreach_particle(P) {
#if _MPI
		particle->pid = locate_particle_rank();
		if(local_particle())
#else
		if(locate_particle().level>=0)
#endif
			nlocal++;
	}

	//fixme : WHy yu need to delete non-local particles?
	particles_remove_nonlocal(P);

	FILE * fp = fopen(file, "w");
	int pstart = 0;
	int nspace  = (dimension+1)*10;
#if _MPI
	int * nps = (int *)malloc(npe()*sizeof(int));
	MPI_Allgather(&nlocal, 1,  MPI_INT, nps, 1, MPI_INT, MPI_COMM_WORLD);
	for(int p=0; p<pid(); ++p)
		pstart += nps[p];
	free(nps);

#endif
	long offset = (nspace+1)*pstart;

	//fixme : separate files	
	foreach_particle(P) 
#if _MPI
		if(local_particle())
#else
		if(locate_particle().level>=0)
#endif
		{
			fseek (fp, offset, SEEK_SET);
			for(int j=0; j<nspace; ++j)
				fprintf(fp, " ");
			fprintf(fp, "\n");
			fseek (fp, offset, SEEK_SET);
			for(int j=0; j<dimension; ++j)
				fprintf(fp, "%f ", particle->x[j]);
			double col = particle_color();
			fprintf(fp, "%.3f ", col);
			offset += nspace+1;
		}

	fclose(fp);
} 
