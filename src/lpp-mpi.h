#if _MPI

typedef struct {
	MPI_Request s[2];
	int pid;
	Array * q;
}Pmpi_sr;

typedef struct{
	int npid;
	Pmpi_sr * sr;
	int * mappid;
}Pmpi;

@def foreach_pmpi_sr(f)
{
	Pmpi_sr * sr = f->sr; 	
	for(int p=0; p < f->npid; ++p, ++sr) {
@
@def end_foreach_pmpi_sr()
	}
}
@

Pmpi * pmpi_new() {
		
	MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
	RcvPid * rcvpid = (mpi->mpi_level).snd;

	int npid = rcvpid->npid;
	//assert(npid >0);
	Pmpi_sr * sr = (Pmpi_sr *)malloc(npid*sizeof(Pmpi_sr)), * s;
	s = sr;
	for(int p=0; p<npid; ++p, ++s) {
		s->pid = (rcvpid->rcv[p]).pid;
		s->q = array_new();
	}

	//fixme : mappid slows down
	int * mappid = (int *)malloc(npe()*sizeof(int));
	for(int p=0; p<npe(); p++) 
		mappid[p] = -1;
	for(int p=0; p<npid; ++p)
		mappid[(rcvpid->rcv[p]).pid] = p;	

	Pmpi * f = (Pmpi *)malloc(sizeof(Pmpi));
	f->mappid = mappid;
	f->npid = npid;
	f->sr = sr;
	return f;
}

void pmpi_free(Pmpi * f){
	foreach_pmpi_sr(f)
		array_free(sr->q);
	free(f->sr);
	free(f->mappid);
	free(f);
}

//To reassign particles moving out of domain bdry
void particles_redistribute(Particles * P){
	/**
Make sure CFL < 1.0
*/

	Pmpi * f = pmpi_new();

	foreach_particle(P) {
		int rpid = locate_particle_rank();
		particle->pid = rpid;
		if(rpid == pid() || rpid == -1) continue;
		int map  = f->mappid[rpid]; 
		array_append(f->sr[map].q, particle, sizeof(Particle));
	} 

	foreach_pmpi_sr(f){
		Array * s = sr->q;
  	MPI_Isend (&(s->len), 1, MPI_LONG, sr->pid, 123, MPI_COMM_WORLD, sr->s);
  	if (s->len > 0) 
    	MPI_Isend (s->p, s->len, MPI_BYTE, sr->pid, 123, MPI_COMM_WORLD, sr->s + 1); 
	}
	
	particles_remove_nonlocal(P);

	foreach_pmpi_sr(f){
  	long len;
  	mpi_recv_check (&len, 1, MPI_LONG, sr->pid, 123,
		  MPI_COMM_WORLD, MPI_STATUS_IGNORE, "return_pmpi (len)");
  	if (!len) continue;
		unsigned int np = len/sizeof(Particle);
		particles_append(P, np);
		Particle * pa = P->p + P->len; 
    mpi_recv_check (pa, len, MPI_BYTE, sr->pid, 123,
		    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "return_pmpi (p)");
		P->len += np;
	}

#if _PDEBUG
	foreach_particle(P)
		assert(local_particle());
	foreach_particle(P)
		assert(locate_particle().level>=0);
#endif

	pmpi_free(f);
}

#endif

