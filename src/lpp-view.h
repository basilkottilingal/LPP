/**
# Copied from  scatter2.h (sandbox/Antoonvh)
 */
#include "lpp.h"

#if dimension == 3 
void glPointParameterfv(GLenum pname, const GLfloat * params);
#endif

struct _lpp_view {
  float s, pc[3], coefs[3]; // point size, colour and distance attenuation coefs.
};

trace
void lpp_view (struct _lpp_view p){
  Particles * P = get_particles();    // particles
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
#else // Dimension == 3
  if (!p.coefs[0]){ // A guess:
    p.coefs[0] = 0.01;
    p.coefs[1] = 0.2;
    p.coefs[2] = 0.5;
  }
  glPointParameterfv(GL_POINT_DISTANCE_ATTENUATION, p.coefs);
#endif
  glEnable (GL_BLEND);
  glEnable (GL_POINT_SMOOTH);
#if _MPI	
	double cmap[127][3];
	int icmap = ((int)  (128*pid()/npe())) % 128;
	jet(cmap);
	//fixme:  Overrridng p.pc	
	glColor3f(cmap[icmap][0], cmap[icmap][1], cmap[icmap][2]);
#else
	glColor3f(p.pc[0], p.pc[1], p.pc[2]);
#endif		
  if (p.s <= 0.)
    p.s = 20;
  glPointSize(p.s);

  glBegin (GL_POINTS);
  foreach_particle(P) {
#if dimension == 2    
    glvertex2d (view, xp, yp);
#else
    glvertex3d (view, xp, yp, zp);
#endif
  }
  glEnd();

  view->ni++; 
}
