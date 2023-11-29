#include "source.h"

 
Icosahedron::Icosahedron(int nlvl)

{
int i, j, k;
double alpha;

ver=new Vector[NVER];

fcemd=new Vector*[NFCE];
for(i=0;i<NFCE;i++) fcemd[i]=new Vector[4];

edg2=new int*[NEDGE];
for(i=0;i<NEDGE;i++) edg2[i]=new int[2];

fc=new int*[NFCE];
for(i=0;i<NFCE;i++) fc[i]=new int[3];

fcedg=new int*[NFCE];
for(i=0;i<NFCE;i++) fcedg[i]=new int[6];

edg=new int*[NFCE];
for(i=0;i<NFCE;i++) edg[i]=new int[6];

nedg=new int[NVER];

ver[0].x=0.0;
ver[0].y=0.0;
ver[0].z=1.0;
ver[11].x=0.0;
ver[11].y=0.0;
ver[11].z=-1.0;

for(i=1;i<6;i++) ver[i].z=(float)0.447214;
for(i=6;i<11;i++) ver[i].z=(float)-0.447214;
  

nv=12;
ne=0;
nf=0;

for(i=1;i<6;i++) {
  alpha = ((i+1)*3.141593/2.5)-3.141593;

  ver[i].x=cos(alpha)*(float).894427;
  ver[i].y=sin(alpha)*(float).894427;

  j=i+1;
  if(j>=6) j=1;

  fc[nf+i-1][0]=0;
  fc[nf+i-1][1]=i;
  fc[nf+i-1][2]=j;

  }

nf+=5;

for (i=1;i<6;i++) {
  j=i+1;
  if(j>=6) j=1;
  k=i+5;
  fc[nf+i-1][0]=i;
  fc[nf+i-1][1]=j;
  fc[nf+i-1][2]=k;
  }

nf+=5;

for(i=6;i<11;i++) {
  j=i+1;
  if(j>=11) j=6;
  k=i-4;
  if(k>=6) k=1;
  fc[nf+i-6][0]=i;
  fc[nf+i-6][1]=j;
  fc[nf+i-6][2]=k;
  }

nf += 5;

for(i=6;i<11;i++) {
  alpha=(i-6)*3.141593/2.5;
  ver[i].x=cos(alpha)*(float).894427;
  ver[i].y=sin(alpha)*(float).894427;
  j=i+1;
  if(j>=11) j=6;
  fc[nf+i-6][0]=i;
  fc[nf+i-6][1]=j;
  fc[nf+i-6][2]=11;
  }

nf+=5;

face(nlvl);

vncl=new int[nv];

delete nedg;

} 



void Icosahedron::face(int nlvl)

{
int ilvl, i, k, iv;
int fct[NFCE][3];
int ifc,ntf;
int iv1,iv2,iv3,imd1,imd2,imd3;

for(ilvl=1;ilvl<=nlvl;ilvl++) {
  ntf=0;

  for(iv=0;iv<nv;iv++) 
    nedg[iv]=0;

  for(ifc=0;ifc<nf;ifc++)
    {
    iv1=fc[ifc][0];
    iv2=fc[ifc][1];
    iv3=fc[ifc][2];
    imd1=divide(iv1,iv2);
    imd2=divide(iv2,iv3);
    imd3=divide(iv3,iv1);
    fct[ntf][0]=iv1;
    fct[ntf][1]=imd1;
    fct[ntf][2]=imd3;
    ++ntf;
    fct[ntf][0]=imd1;
    fct[ntf][1]=iv2;
    fct[ntf][2]=imd2;
    ++ntf;
    fct[ntf][0]=imd3;
    fct[ntf][1]=imd2;
    fct[ntf][2]=iv3;
    ++ntf;
    fct[ntf][0]=imd2;
    fct[ntf][1]=imd3;
    fct[ntf][2]=imd1;
    ++ntf;
    fcedg[ifc][0]=imd1;
    fcedg[ifc][1]=imd2;
    fcedg[ifc][2]=imd3;
    } 

  if(ilvl==nlvl) break;

  for(i=0;i<ntf;i++) for(k=0;k<3;k++)
    fc[i][k]=fct[i][k];

  nf = ntf;
  }

trimid();
return;

}


int Icosahedron::divide(int a,int b)

{
int j;
static int midver[NVER][6];
float xm, ym, zm, dis;


for(j=0;j<nedg[a];j++) if(edg[a][j]==b) return midver[a][j]; 

xm=(ver[a].x+ver[b].x)/2.0;
ym=(ver[a].y+ver[b].y)/2.0;
zm=(ver[a].z+ver[b].z)/2.0;

dis=sqrt(xm*xm+ym*ym+zm*zm);

ver[nv].x=xm/dis;
ver[nv].y=ym/dis;
ver[nv].z=zm/dis;
edg[a][nedg[a]]=b;
edg[b][nedg[b]]=a;
midver[a][nedg[a]]=nv;
nedg[a]++;
midver[b][nedg[b]]=nv;
nedg[b]++;


edg2[ne][0]=(a<b ? a : b);
edg2[ne][1]=(a>b ? a : b);

ne++;

assert(nedg[a]<7);
assert(nedg[b]<7);

return nv++;

}


void Icosahedron::trimid()

{
Vector v1,v2,v3,v4;
int ie1, ie2, ie3, iv1, iv2, iv3, ifc;

for(ifc=0;ifc<nf;ifc++)
  {
  iv1=fc[ifc][0];
  iv2=fc[ifc][1];
  iv3=fc[ifc][2];

  ie1=fcedg[ifc][0];
  ie2=fcedg[ifc][1];
  ie3=fcedg[ifc][2];

  v1=(ver[iv1]+ver[ie1]+ver[ie3])/3.0;
  v2=(ver[iv2]+ver[ie1]+ver[ie2])/3.0;
  v3=(ver[iv3]+ver[ie2]+ver[ie3])/3.0;
  v4=(ver[ie1]+ver[ie2]+ver[ie3])/3.0;

  v1=v1/v1.norm();
  v2=v2/v2.norm();
  v3=v3/v3.norm();
  v4=v4/v4.norm();

  fcemd[ifc][0]=v1;
  fcemd[ifc][1]=v2;
  fcemd[ifc][2]=v3;
  fcemd[ifc][3]=v4;
  }

return;

}



Icosahedron::~Icosahedron()

{
int i;

delete ver;

for(i=0;i<NFCE;i++) delete fcemd[i];
delete fcemd;

for(i=0;i<NEDGE;i++) delete edg2[i];
delete edg2;

for(i=0;i<NFCE;i++) delete fc[i];
delete fc;

for(i=0;i<NFCE;i++) delete fcedg[i];
delete fcedg;

for(i=0;i<NFCE;i++) delete edg[i];
delete edg;

delete vncl;

}
