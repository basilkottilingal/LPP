/**

C99='mpicc -std=c99' qcc  -D_MPI=1 -events -o a.out attempt1.c -lm
*/
#include "grid/octree.h"
//#include "axi.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "contact.h"
#include "vof.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "curvature.h"
#define _PTAG 1
#include "../../../src/lpp.h"

//#include "view.h"



//Dimensionless quantities

#define Oh 2.185e-03  // MU1/sqrt(Rho1*SIGMA*Radius)  // Ohnesorge number for water/air
//#define Bo 3.171    // Rho1*(A-g)*Req*Req/SIGMA  // Acceleration Bond number (gamma = 4.2)
#define bo 2.
#define Rratio 1000.0  // Rho1/Rho2
#define MUratio 100.0  //MU1/MU2

//System properties (SI)
#define Rho1 998.0  // liquid
#define MU1 0.00102

#define Rho2 Rho1/Rratio  // gas
#define MU2 MU1/MUratio

#define Radius 3e-03 // Radius of the hemispherical droplet (3mm) V=56.5µL // Req = 0.79*Radius

#define SIGMA (sq(MU1)/(Rho1*Radius*sq(Oh))) // 72.8 mN/m for water/air
#define A (bo*SIGMA/(Rho1*0.79*Radius*0.79*Radius)) // acceleration of the plate in the frame of reference
#define Tc sqrt(Rho1*0.5*Radius*Radius*Radius/SIGMA) // capillary time 13.6ms


//System characteristics (dimensionless)
//We rescale everything according to Radius=1, Rho2_nd=Rho2 and Tc constant

#define Radius_nd (Radius/Radius)
#define Rho2_nd Rho2
#define Rho1_nd (Rratio*Rho2_nd)

#define SIGMA_nd (Rho1_nd*0.5*pow(Radius_nd,3)/pow(Tc,2))

#define A_nd (bo*SIGMA_nd/(Rho1_nd*0.79*Radius_nd*0.79*Radius_nd))

#define MU1_nd (Oh*sqrt(Rho1_nd*SIGMA_nd*Radius_nd))
#define MU2_nd (MU1_nd/MUratio)

double theta_d = 90.*(pi/180);
#define tend 120*Tc //time upto which simulation runs or movie is generated


vector h[];
h.t[back] = contact_angle (theta_d); // impose a constant contact angle
h.r[back] = contact_angle (theta_d);
u.t[back] = dirichlet(0); // no slip (comment these lines for free slip BC)
u.r[back] = dirichlet(0); // no slip (comment these lines for free slip BC)
u.n[back] = dirichlet(0);
u.n[front] = u.n[] > 0 ? neumann(0) : 0;
/*
Particles __new;
Particles Piclb;
Particles Piclc;
Particles Picld;
Particles Piil;
Particles Piol;
Particles Pbtt;
*/
int LEVEL = 6;

int main() {
  /*
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */
  size (10*Radius_nd);
  f.height = h;
  
  rho1 = Rho1_nd;
  mu1 = MU1_nd;
  rho2 = Rho2_nd;
  mu2 = MU2_nd;
  f.sigma = SIGMA_nd;              
  

  N = 1 << 6;
  init_grid(N);
  
    //for (bo=1.0; bo<= 7.0; bo += 1.0)
    run();
/*
  for (Bo=1; Bo<= 5; Bo += 1)
    run();
}
*/
}


event init (t = 0)
{
  if (!restore (file = "restart")){
		//fraction (f, - (sq(x) + sq(y) - sq(Radius_nd)));
 		fraction (f, - (sq(x) + sq(y) + sq(z) - sq(Radius_nd)));
 		//fprintf(ferr,"bo \t%e\t%e\t%e\n",bo,Rho1*A*Radius*Radius/SIGMA,Rho1_nd*A_nd*Radius_nd*Radius_nd/SIGMA_nd);
  	fprintf(ferr,"Oh \t%e\t%e\t%e\n",Oh,MU1/sqrt(Rho1*SIGMA*Radius),MU1_nd/sqrt(Rho1_nd*Radius_nd*SIGMA_nd));
  	fprintf(ferr,"Tc \t%e\t%e\t%e\n",Tc, sqrt(Rho1*pow(Radius,3)/SIGMA),sqrt(Rho1_nd*pow(Radius_nd,3)/SIGMA_nd));
  	fprintf(ferr,"Val \t%e\t%e\t%e\n",A_nd,SIGMA_nd,MU1_nd);
  	fprintf(ferr,"Volume in microliters \t%e\n",(1e9)*(4/3)*(pi)*pow(Radius*0.79 , 3));
  	fprintf(ferr,"Actual Volume in microliters \t %e",(1e9)*(4/3)*(pi)*pow(Radius , 3)*0.5);  
  }

	Particles * P = get_particles();
	for (int _k=0; _k<5; ++_k){
		int np = 100;
		particles_append(P, np); /*Memory for 100 particles*/
		for(int ip=0; ip<np; ++ip){
			double Rt = 0.9 + 0.1*noise();
			double th = (2*_k+1)*pi/16. + (pi/16.)*noise();
			double xp = Rt*cos(th), yp = Rt*sin(th), zp = 0.1*noise();
			if(locate(xp,yp,zp).level >= 0) {
		Particle * p = add_particle(P);
		p->x[0] = xp;		p->x[1] = yp; p->x[2] = zp;
#if _PTAG
		p->tag = (double) _k;
#endif
			}
		}
	}
 
	boundary({f});
	dump("dump-initial");
		
}

event advance(i++){

	double dt = DT;
	Particles * P = get_particles();
	particles_advance_eulerian(P, u);
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(z)
    av.z[] += A_nd;      
}

event adapt (i++){
  adapt_wavelet_particles ((scalar*){f,u}, (double[]){0.01, 0.1, 0.1,0.1}, LEVEL);
}

event snapshots (t += 0.1*Tc ; t<=tend)
{
  char name[80];
  sprintf (name, "s-%g-%d", t/Tc,LEVEL);
  scalar pid[];
  scalar omega[];
  vorticity (u, omega);
  scalar kappa[];
  curvature (f, kappa);
  foreach(){
      pid[] = fmod(pid()*(npe() + 37), npe());
      omega[] *= 0.1;
    }
  boundary ({pid,omega});
  
/*
h.x.nodump = true;
h.y.nodump = true;
h.z.nodump = true;    

g.x.nodump = true;
g.y.nodump = true;
g.z.nodump = true;


rhov.nodump = true;
//f.nodump = true;

u.x.nodump = true;
u.y.nodump = true;
u.z.nodump = true;
*/
dump (name);
    output_facets (f, stdout);

}
/*
event img (t = 0 ; t += 0.1*Tc ; t<=tend)
{
 char name1[80];
  sprintf (name1, "img-kappa-%g.png",t/Tc);
   scalar omega[];
  vorticity (u, omega);
  foreach()
    omega[] *= 0.1;
   scalar kappa[];
        curvature (f, kappa);
  clear();
 view (fov = 15., quat = {0.474458,0.144142,0.234923,0.836017},
	  tx = 0., ty = -0.3 , bg = {1,1,1},
	  width = 300, height = 1000);
  //box(notics = true , lw =2);
  //squares("omega");
    draw_vof ("f", color = "kappa");
 scatter (__new, s = 5, pc = {0., 4., 4.});
    scatter (Piclb, s = 5, pc = {4., 0., 4.});
    scatter (Piclc, s = 5, pc = {0., 4., 0.});
    scatter (Picld, s = 5, pc = {0., 0., 4.});
    scatter (Piil, s = 5, pc = {4., 0., 0.});
    scatter (Piol, s = 5, pc = {0., 0., 0.});
    scatter (Pbtt, s = 5, pc = {4., 4., 0.});
    cells (lc = {1,0,0});

    mirror (n = {1,0,0}) {
    draw_vof ("f", color = "kappa");
 scatter (__new, s = 5, pc = {0., 4., 4.});
    scatter (Piclb, s = 5, pc = {4., 0., 4.});
    scatter (Piclc, s = 5, pc = {0., 4., 0.});
    scatter (Picld, s = 5, pc = {0., 0., 4.});
    scatter (Piil, s = 5, pc = {4., 0., 0.});
    scatter (Piol, s = 5, pc = {0., 0., 0.});
    scatter (Pbtt, s = 5, pc = {4., 4., 0.});
    cells (lc = {1,0,0});
    }
    mirror (n = {0,1,0}) {
      draw_vof ("f", color = "kappa");
      //draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
      mirror (n = {1,0,0}) {
	draw_vof ("f", color = "kappa");
	//draw_vof ("f", edges = true);
	cells (lc = {1,0,0});
      }
    }

  save (name1); 

 char name2[80];
  sprintf (name2, "uz-%g.png",t/Tc);
  clear();
 view (fov = 15., quat = {0.474458,0.144142,0.234923,0.836017},
	  tx = 0., ty = -0.3 , bg = {1,1,1},
	  width = 300, height = 1000);
  //box(notics = true , lw =2);
  //squares("omega");
    draw_vof ("f", color = "u.z");
 scatter (__new, s = 5, pc = {0., 4., 4.});
    scatter (Piclb, s = 5, pc = {4., 0., 4.});
    scatter (Piclc, s = 5, pc = {0., 4., 0.});
    scatter (Picld, s = 5, pc = {0., 0., 4.});
    scatter (Piil, s = 5, pc = {4., 0., 0.});
    scatter (Piol, s = 5, pc = {0., 0., 0.});
    scatter (Pbtt, s = 5, pc = {4., 4., 0.});
    cells (lc = {1,0,0});

    mirror (n = {1,0,0}) {
    draw_vof ("f", color = "u.z");
 scatter (__new, s = 5, pc = {0., 4., 4.});
    scatter (Piclb, s = 5, pc = {4., 0., 4.});
    scatter (Piclc, s = 5, pc = {0., 4., 0.});
    scatter (Picld, s = 5, pc = {0., 0., 4.});
    scatter (Piil, s = 5, pc = {4., 0., 0.});
    scatter (Piol, s = 5, pc = {0., 0., 0.});
    scatter (Pbtt, s = 5, pc = {4., 4., 0.});
    cells (lc = {1,0,0});
    }
    mirror (n = {0,1,0}) {
      draw_vof ("f", color = "kappa");
      //draw_vof ("f", edges = true);
      cells (lc = {1,0,0});
      mirror (n = {1,0,0}) {
	draw_vof ("f", color = "kappa");
	//draw_vof ("f", edges = true);
	cells (lc = {1,0,0});
      }
    }

  save (name2); 

}
*/
