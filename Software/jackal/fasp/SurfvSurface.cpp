#include "source.h"
/*
#include <stdlib.h>
#include <assert.h>
#include "surfvsurface.h"
#include "icosahedron.h"
#include "chain.h"
#include "app.h"
*/

//using namespace Troll;


/*
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
*/

SurfvSurface::SurfvSurface() {
sp=0;
rp2=0;
gb=0;
probe_radius=1.4;
atom=0;
}

SurfvSurface::SurfvSurface(Atom**a,int spherelevel,float pro)
{
int i,n;

atom=a;

int tt=0;
while(atom[tt])tt++;
if(tt==0) return;

probe_radius=pro;

sp=new Icosahedron(spherelevel);

gb=new Grid(atom,2.0+probe_radius);

//n=atom.size();
n=tt;
rp2=new float[n]; 
Atom **ap=atom;//atom.begin();
for(i=0;i<n;i++) {
  float r=ap[i]->radius+probe_radius;
  rp2[i]=r*r*.99999;
  }
}

void SurfvSurface::ready(Atom**a,int spherelevel,float pro) {

int i,n;

atom=a;

int tt=0;
while(atom[tt])tt++;
if(tt==0) return;

probe_radius=pro;

if(sp==0)sp=new Icosahedron(spherelevel);

gb=new Grid(atom,2.0+probe_radius);

//n=atom.size();
n=tt;
rp2=new float[n];
Atom **ap=atom;//atom.begin();
for(i=0;i<n;i++) {
  float r=ap[i]->radius+probe_radius;
  rp2[i]=r*r*.99999;
}

}


/*
void SurfvSurface::Build(unsigned long flag)

{
int i;

int natoms=0;//=atom.size();

for(i=0;i<natoms;i++) Surface(atom[i]->pos,atom[i]->radius+probe_radius);

}
*/



SurfvSurface::~SurfvSurface()

{

if(gb) delete gb;
if(sp) delete sp;
if(rp2) delete rp2;

}



float SurfvSurface::Surface(Vector center,float rad)
{
int j,k,l,nacc;
int iv,ie,ifc,ivm;
register float r,cx,cy,cz,dx,dy,dz;
int satom[4][3];
int atom1, atom2, atom3;
int a1, a2, a3, a4, a5, a6;
int *vncl;
Vector *v,u;
ObjectBox box;
Atom **ap;
Grid *gridbox=0;
Icosahedron *sphere=0;
float *rplusp2=0;

if(gridbox==NULL) gridbox=gb;
if(sphere==NULL) sphere=sp;
if(rplusp2==NULL) rplusp2=rp2;

vncl=sphere->vncl;

//ap=atom.begin();
ap=atom;

cx=center.x;
cy=center.y;
cz=center.z;

vncl=sphere->vncl;

nacc=0;


// See if first twelve vertices are buried.

for(iv=0;iv<12;++iv) {
  
  v=&sphere->ver[iv];
  
  u.x=cx+v->x*rad;
  u.y=cy+v->y*rad;
  u.z=cz+v->z*rad;


// Get box vertex is in.

  box=(*gridbox)[u];
  

// Does one of atoms in grid box bury vertex?
  
  for(k=0;k<box.ne;k++) {

    j=box.element[k];
    v=&ap[j]->pos;
    r=rplusp2[j];
      
    dx=u.x-v->x;
    dy=u.y-v->y;
    dz=u.z-v->z;
    
    if((dx*=dx)>r) continue;
    if((dy*=dy)>r) continue;
    if((dz*=dz)>r) continue;
      
    if((dx+dy+dz)<r) { vncl[iv]=j; break; }

    }

  if(k==box.ne) { vncl[iv]=-1; nacc++; }

  }


// See edge vertices are buried.

for(ie=0,ivm=12;ie<sphere->ne;ie++,ivm++) {

  a1=vncl[sphere->edg2[ie][0]];
  a2=vncl[sphere->edg2[ie][1]];
  
  if(a1>0) if(a1==a2) {
    vncl[ivm]=a1;
    continue;
    }

  v=&sphere->ver[ivm];
  u.x=cx+v->x*rad;
  u.y=cy+v->y*rad;
  u.z=cz+v->z*rad;
  
  if(a1>=0) {
    v=&ap[a1]->pos;
    dx=u.x-v->x;
    dy=u.y-v->y;
    dz=u.z-v->z;

    if((dx*dx+dy*dy+dz*dz)<rplusp2[a1]) {
      vncl[ivm]=a1;
      continue;
      }
    }

  if(a2>=0) {
    v=&ap[a2]->pos;
    dx=u.x-v->x;
    dy=u.y-v->y;
    dz=u.z-v->z;

    if((dx*dx+dy*dy+dz*dz)<rplusp2[a2]) {
      vncl[ivm]=a2;
      continue;
      }
    }

  // See if nearby atom buries vertex.
  
  box=(*gridbox)[u];

  for(k=0;k<box.ne;k++) {

    j=box.element[k];
    v=&ap[j]->pos;
    r=rplusp2[j];
      
    dx=u.x-v->x;
    dy=u.y-v->y;
    dz=u.z-v->z;
    
    if((dx*=dx)>r) continue;
    if((dy*=dy)>r) continue;
    if((dz*=dz)>r) continue;
      
    if((dx+dy+dz)<r) { vncl[ivm]=j; break; }

    }

  if(k==box.ne) { vncl[ivm]=-1; nacc++; }

  }
  

// Check midpoints of faces.

for(ifc=0;ifc<sphere->nf;ifc++) {

  // Get index of nearest neighbor to three outer vertices.
  
  a1=vncl[sphere->fc[ifc][0]];
  a2=vncl[sphere->fc[ifc][1]];
  a3=vncl[sphere->fc[ifc][2]];

  
  // If three outer vertices of face buried by same atom, done.
  
  if(a1>=0) if(a1 == a2 && a1 == a3) continue;

  a4=vncl[sphere->fcedg[ifc][0]];
  a5=vncl[sphere->fcedg[ifc][1]];
  a6=vncl[sphere->fcedg[ifc][2]];

  satom[0][0]=a1;
  satom[0][1]=a4;
  satom[0][2]=a6;
  satom[1][0]=a2;
  satom[1][1]=a4;
  satom[1][2]=a5;
  satom[2][0]=a6;
  satom[2][1]=a3;
  satom[2][2]=a5;
  satom[3][0]=a4;
  satom[3][1]=a5;
  satom[3][2]=a6;
  

  // Check if four midpoints are buried.
  
  for(j=0;j<4;j++) {
    atom1=satom[j][0];
    atom2=satom[j][1];
    atom3=satom[j][2];
    
    // If surrounding vertices buried by same atom, don't check middle vertex.
    
    v=&sphere->fcemd[ifc][j];
    u.x=cx+v->x*rad;
    u.y=cy+v->y*rad;
    u.z=cz+v->z*rad;

    if(atom1>=0) { if(atom1==atom2) if(atom2==atom3) continue; }
    else if(atom2<0) if(atom3<0)  { nacc++; continue; } 
    
    // If the vertex is buried by any of its neighbors, go to next vertex.
    
    if(atom1>=0) {
      v=&ap[atom1]->pos;
      dx=u.x-v->x;
      dy=u.y-v->y;
      dz=u.z-v->z;

      if((dx*dx+dy*dy+dz*dz)<rplusp2[atom1]) continue;
      }

    if(atom2>=0) {
      v=&ap[atom2]->pos;
      dx=u.x-v->x;
      dy=u.y-v->y;
      dz=u.z-v->z;

      if((dx*dx+dy*dy+dz*dz)<rplusp2[atom2]) continue;
      }

    if(atom3>=0) {
      v=&ap[atom3]->pos;
      dx=u.x-v->x;
      dy=u.y-v->y;
      dz=u.z-v->z;

      if((dx*dx+dy*dy+dz*dz)<rplusp2[atom3]) continue;
      }

    

    // If not buried by nearest neighbors, use grid to look at nearby atoms.
    
    box=(*gridbox)[u];

    for(k=0;k<box.ne;k++) {

      l=box.element[k];
      v=&ap[l]->pos;
      r=rplusp2[l];
      
      dx=u.x-v->x;
      dy=u.y-v->y;
      dz=u.z-v->z;
    
      if((dx*=dx)>r) continue;
      if((dy*=dy)>r) continue;
      if((dz*=dz)>r) continue;
      
      if((dx+dy+dz)<r) break;
      }
    
    if(k==box.ne) nacc++;

    }
  }

return rad*rad*(float)nacc*12.56637/(float)(sphere->nv+4*sphere->nf);

} 


