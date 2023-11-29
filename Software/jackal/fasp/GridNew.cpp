#include "source.h"

void GridNew::initial() {
	gridvalue=0;
	gridsize=3.4;
	resolution=0;
	leg=5.;
	npoints=0;
	level=1;
	nx=ny=nz=0;
	slicesize=0;
	data=0;
}

GridNew::GridNew() {
initial();
}

void GridNew::ready(Pdb *pd){
Atom **data0=pd->getAtomAll();
setGridNew(data0);
}

void GridNew::ready(Atom **data0) {
setGridNew(data0);
}

void GridNew::setGridNew(Atom **data0) 
{
int i;
int nelements;

/*
int j,k,l,o;
  int x_ll, x_ul, y_ll, y_ul, z_ll, z_ul;
  int *temp;
*/
data=data0;
nelements=0;
while(data[nelements])nelements++;
if(nelements==0) return;

defaultvalue.ne=0;

float sz=gridsize;

origin=data[0]->pos;
extent=data[0]->pos;

origin=origin-Vector(sz,sz,sz);
extent=extent+Vector(sz,sz,sz);


for(i=0;i<nelements;i++) {
  if(data[i]->pos.x<origin.x) origin.x=data[i]->pos.x;
  if(data[i]->pos.y<origin.y) origin.y=data[i]->pos.y;
  if(data[i]->pos.z<origin.z) origin.z=data[i]->pos.z;
  if(data[i]->pos.x>extent.x) extent.x=data[i]->pos.x;
  if(data[i]->pos.y>extent.y) extent.y=data[i]->pos.y;
  if(data[i]->pos.z>extent.z) extent.z=data[i]->pos.z;
}

resolution=1.0/gridsize;

// Determine number of grid boxes.

origin.x-=gridsize+leg;
origin.y-=gridsize+leg;
origin.z-=gridsize+leg;
extent.x+=gridsize+leg;
extent.y+=gridsize+leg;
extent.z+=gridsize+leg;

nx=(int)(1+(extent.x-origin.x)/gridsize);
ny=(int)(1+(extent.y-origin.y)/gridsize);
nz=(int)(1+(extent.z-origin.z)/gridsize);


// Allocate memory for grid boxes;

npoints=nx*ny*nz+1;
gridvalue=new ObjectBox[npoints];
if(gridvalue==NULL) 
{
cerr<<"GridNew::GridNew(), too many grid boxes"<<endl;
}

for(i=0;i<npoints;i++) {
  gridvalue[i].ne=0;
  gridvalue[i].max=30;
  gridvalue[i].element=NULL;
  gridvalue[i].element=0;
}


// Fill the boxes with objects.

for(i=0;i<nelements;i++) {
  add(i);
}
  
}

void GridNew::add(int i) {
	
  int j,k,l,o;
  int x_ll, x_ul, y_ll, y_ul, z_ll, z_ul;
  int *temp;
  Atom *aa=data[i];
  Vector v=aa->pos;
  j=(int)((v.x-origin.x)*resolution);
  k=(int)((v.y-origin.y)*resolution);
  l=(int)((v.z-origin.z)*resolution);
  x_ll=max(0,j-level);
  x_ul=min(nx-1,j+level);
  y_ll=max(0,k-level);
  y_ul=min(ny-1,k+level);
  z_ll=max(0,l-level);
  z_ul=min(nz-1,l+level);

  for(j=x_ll;j<=x_ul; j++) for(k=y_ll;k<=y_ul;k++) for(l=z_ll;l<=z_ul;l++) {
    o=l*nx*ny + k*nx + j;
    if(gridvalue[o].ne==0) gridvalue[o].element=new int[30];
    if(gridvalue[o].ne==gridvalue[o].max) {
      gridvalue[o].max+=10;
      temp=new int[gridvalue[o].max];
      assert(temp!=NULL);
      memmove(temp,gridvalue[o].element,(gridvalue[o].max-10)*sizeof(int));
      delete gridvalue[o].element;
      gridvalue[o].element=temp;
    }
    gridvalue[o].element[gridvalue[o].ne++]=aa->id;
  }
}

void GridNew::putoff(Res *a) {
for(Atm *b=a->atm;b;b=b->next) putoff(b);
}

void GridNew::putoff(Atm *a) {
  int i,j,k,l;
  int x_ll, x_ul, y_ll, y_ul, z_ll, z_ul;
  int o;
  i=a->id0;
  Vector v=data[i]->pos;
  j=(int)((v.x-origin.x)*resolution);
  k=(int)((v.y-origin.y)*resolution);
  l=(int)((v.z-origin.z)*resolution);
  
  x_ll=max(0,j-level);
  x_ul=min(nx-1,j+level);
  y_ll=max(0,k-level);
  y_ul=min(ny-1,k+level);
  z_ll=max(0,l-level);
  z_ul=min(nz-1,l+level);
  
  for(j=x_ll;j<=x_ul;j++) for(k=y_ll;k<=y_ul;k++) for(l=z_ll;l<=z_ul;l++) {
    
    o=l*nx*ny + k*nx + j;
    if(gridvalue[o].ne==0) continue;
    gridvalue[o].remove(i);
  }
}

void GridNew::puton(Atm *a) {
  int i,j,k,l;
  int x_ll, x_ul, y_ll, y_ul, z_ll, z_ul;
  int *temp;
  int o;
  i=a->id0;
  Vector v=data[i]->pos;

  /*
  data[i]->pos.x=a->xyz[0];
  data[i]->pos.y=a->xyz[1];
  data[i]->pos.z=a->xyz[2];
  data[i]->radius=a->tatm->eng->radius;
  */

  j=(int)((v.x-origin.x)*resolution);
  k=(int)((v.y-origin.y)*resolution);
  l=(int)((v.z-origin.z)*resolution);
  x_ll=max(0,j-level);
  x_ul=min(nx-1,j+level);
  y_ll=max(0,k-level);
  y_ul=min(ny-1,k+level);
  z_ll=max(0,l-level);
  z_ul=min(nz-1,l+level);

  for(j=x_ll;j<=x_ul; j++) for(k=y_ll;k<=y_ul;k++) for(l=z_ll;l<=z_ul;l++) {
    o=l*nx*ny + k*nx + j;
    if(gridvalue[o].exist(i)) continue;
    if(gridvalue[o].ne==0) gridvalue[o].element=new int[30];
    if(gridvalue[o].ne==gridvalue[o].max) {
      gridvalue[o].max+=10;
      temp=new int[gridvalue[o].max];
      assert(temp!=NULL);
      memmove(temp,gridvalue[o].element,(gridvalue[o].max-10)*sizeof(int));
      delete gridvalue[o].element;
      gridvalue[o].element=temp;
    }
    //gridvalue[o].add(i);
    gridvalue[o].element[gridvalue[o].ne++]=i;
  }
}


void GridNew::puton(Res *a) {
for(Atm *b=a->atm;b;b=b->next) puton(b);
}




ObjectBox& GridNew::operator[](Vector v)
{
int i, j, k;

i=(int)((v.x-origin.x)*resolution);
j=(int)((v.y-origin.y)*resolution);
k=(int)((v.z-origin.z)*resolution);
if(i<0 || i>=nx) return defaultvalue;
if(j<0 || j>=ny) return defaultvalue;
if(k<0 || k>=nz) return defaultvalue;

i=k*nx*ny + j*nx + i;

return gridvalue[i];

}


GridNew::~GridNew()
{
int i;

if(npoints) for(i=0;i<npoints;i++) if(gridvalue[i].ne) if(gridvalue[i].element!=NULL) 
  delete gridvalue[i].element;
if(gridvalue!=NULL) delete gridvalue;

}

int GridNew::gettotal() {

        int n=0;
        for(int i=0;i<npoints;i++) {
                n+=gridvalue[i].ne;
        }
        return n;
}

void GridNew::puton(Res *a,int n) {
        for(Res *s=a;s&&s->id0<n;s=s->next) puton(s);
}

void GridNew::puton(Res *a,int n,int m) {
        for(Res *s=a;s;s=s->next) {
                if(s->id0>=n&&s->id0<m) puton(s);
        }
}

void GridNew::putoff(Res *a,int n,int m) {
        for(Res *s=a;s;s=s->next) {
                if(s->id0>=n&&s->id0<m) putoff(s);
        }
}

void GridNew::putoff(Res *a,int n) {
        for(Res *s=a;s&&s->id0<n;s=s->next) putoff(s);
}

void GridNew::puton(Chn *a) {
	for(Res *s=a->res;s;s=s->next) puton(s);
}
void GridNew::puton(Pdb *a) {

	for(Chn *s=a->chn;s;s=s->next) puton(s);
}
void GridNew::putoff(Pdb *a) {

        for(Chn *s=a->chn;s;s=s->next) putoff(s);
}

void GridNew::putoff(Chn *a){
	for(Res *s=a->res;s;s=s->next) putoff(s);
}


void GridNew::putoff() {

	for(int i=0;i<npoints;i++) {
		if(gridvalue[i].ne) {
			if(gridvalue[i].element) delete gridvalue[i].element;
			gridvalue[i].element=0;
			gridvalue[i].max=30;
			gridvalue[i].ne=0;
		}
		else {
			gridvalue[i].max=30;
			gridvalue[i].element=0;
			gridvalue[i].ne=0;
		}
	}
}
