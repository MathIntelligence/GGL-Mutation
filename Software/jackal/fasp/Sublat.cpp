#include"source.h"

Sublat::Sublat(Lattice *s)
{
 int i;
 lattice=s;
 atm=new Atm;
 obtain=0;
 next=0;
 grdsiz=0;
 for(i=0;i<3;i++) mlen[i]=0;
 oncell=0;
 total=0;
 flag=0;
 nget=0;
}

Sublat::~Sublat()
{
 int i;
 if(atm)    delete atm;
 if(obtain) delete [] obtain;
 if(oncell) delete [] oncell;
 lattice=0;
 atm=0;
 obtain=0;
 next=0;
 grdsiz=0;
 for(i=0;i<3;i++) mlen[i]=0;
 oncell=0;
 total=0;
 flag=0;
 nget=0;
}

void Sublat::ready()
{
 Tres *t;
 float d;
 d=-100;
 for(t=&TRES;t;t=t->next)if(t->mlen>d) d=t->mlen;
 ready(d); 
}

void Sublat::ready(float d)
{
int i;
if(grdsiz==0) {cerr<<"no grid step for second lattice"<<endl;exit(0);}
if(d==0) {cerr<<"second lattice maxlen is zero"<<endl;exit(0);}
if(oncell) {delete [] oncell;oncell=0;}
if(obtain) {delete [] obtain;obtain=0;}
for(i=0;i<3;i++) mlen[i]=2*d;
for(i=0;i<3;i++) mgrd[i]=mlen[i]/grdsiz+0.5;
total=mgrd[0]*mgrd[1]*mgrd[2];
oncell=new int[total];
memset(oncell,0,total*sizeof(int));
obtain=new Atm*[lattice->hash/2];
}

void Sublat::ready(Chn *chn)
{
float co1[3],co2[3],len;
Res *r;
Atm *a;
int i;
Tres *t;

if(grdsiz==0) {cerr<<"no grid step for second lattice"<<endl;exit(0);}
if(oncell) {delete [] oncell;oncell=0;}
if(obtain) {delete [] obtain;obtain=0;}

len=-100;
for(t=&TRES;t;t=t->next)if(t->mlen>len) len=t->mlen;
len+=3;

for(i=0;i<3;i++) {co1[i]=100000;co2[i]=-100000;}

for(r=chn->res;r;r=r->next)
for(a=r->atm;a;a=a->next)
{
  if(a->tatm->id>4) continue;
  for(i=0;i<3;i++) 
  {
     if(co1[i]>a->xyz[i]) co1[i]=a->xyz[i];
     if(co2[i]<a->xyz[i]) co2[i]=a->xyz[i];
  }
}

for(i=0;i<3;i++)atm->xyz[i]=co1[i]-len;
for(i=0;i<3;i++)mlen[i]=co2[i]-co1[i]+2*len;
for(i=0;i<3;i++)mgrd[i]=mlen[i]/grdsiz+0.5;
total=mgrd[0]*mgrd[1]*mgrd[2];
oncell=new int[total];
memset(oncell,0,total*sizeof(int));
obtain=new Atm*[lattice->hash/2];
}

int Sublat::getcell(Atm *s,float d,int n)
{
  int i,k,m;
  m=lattice->flag;
  lattice->flag=n;
  k=lattice->getcell(s,d);
  lattice->flag=m;
  nget=lattice->nget;
  for(i=0;i<nget;i++) obtain[i]=lattice->obtain[i]->atm;
  return k;
}

int Sublat::getcell(Chn *s)
{
  int i;
  Res *r;
  Atm *a;
  i=0;
  for(r=s->res;r;r=r->next)
  for(a=r->atm;a;a=a->next)
  {
    obtain[i++]=a;
  }
  return i;
}
void Sublat::addmap(Atm *a,float cut,int ff)
{
 int i,yz;  
 float xyz[3],d,r,co[3];
 int kx[3],ky[3],kz[3];
 yz=mgrd[1]*mgrd[2];

 for(i=0;i<3;i++) xyz[i]=a->xyz[i]-atm->xyz[i]; 
  
 for(i=0;i<3;i++)
 {
   kx[i]=(xyz[i]-cut)/grdsiz-0.5;
   kx[i]=max(0,kx[i]);
   ky[i]=(xyz[i]+cut)/grdsiz+0.5;
   ky[i]=min(ky[i],mgrd[i]);
 }

 r=cut*cut;

 for(kz[0]=kx[0];kz[0]<=ky[0];kz[0]++)
 for(kz[1]=kx[1];kz[1]<=ky[1];kz[1]++)
 for(kz[2]=kx[2];kz[2]<=ky[2];kz[2]++)
 {
   for(i=0;i<3;i++)co[i]=kz[i]*grdsiz+atm->xyz[i]-a->xyz[i];
   d=co[0]*co[0]+co[1]*co[1]+co[2]*co[2];
   if(d>r) continue;
   i=kz[0]*yz+kz[1]*mgrd[2]+kz[2];
   oncell[i]+=ff;
 }
}

void Sublat::getmap(Chn *s,int ff)
{
  Res *r;
  for(r=s->res;r;r=r->next)
  getmap(r,ff);
}

void Sublat::getmap(Res *r,int ff)
{
 Atm *a;
 float cut;
 for(a=r->atm;a;a=a->next)
 { 
   cut=(1.8+a->tatm->eng->radius)*0.7;
   addmap(a,cut,ff);
 }
}

void Sublat::getmap(int ff)
{
  Atm *a;
  int j;
  float cut;
  memset(oncell,0,total*sizeof(int));

  for(j=0;j<nget;j++)
  {
    a=obtain[j];
    cut=0.7*(a->tatm->eng->radius+1.8);
    addmap(a,cut,ff);
  }
}

int Sublat::clash(Atm *a)
{
int e;
int yz,kk[3],i,nn[3],j;
yz=mgrd[2]*mgrd[1];
e=0;
for(i=0;i<3;i++) kk[i]=(a->xyz[i]-atm->xyz[i])/grdsiz;
for(nn[0]=kk[0];nn[0]<kk[0]+2;nn[0]++)
for(nn[1]=kk[1];nn[1]<kk[1]+2;nn[1]++)
for(nn[2]=kk[2];nn[2]<kk[2]+2;nn[2]++)
{
  j=0;
  for(i=0;i<3;i++) if(nn[i]<0||nn[i]>mgrd[i]-1) j=1;
  if(j) continue;
  i=nn[0]*yz+nn[1]*mgrd[2]+nn[2];
  if(oncell[i]>0) break;
  e++;
}
if(e==8) return 1;
return 0;
}

int Sublat::clash(Res *r,int begin,int end)
{
int e;
Atm *a;
for(a=r->atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') break;
  if(a->tatm->id<begin||a->tatm->id>end)continue;
  e=clash(a);
  if(e==1) return 1;
}
return 0;
}

int Sublat::clash(Res *r,int ntat)
{
float e;
Atm *a;
for(a=r->atm;a;a=a->next)
{
if(a->tatm->name[1]=='H') break;
if(a->tatm->rotm!=ntat)continue;
e=clash(a);
if(e==1) return 1;
}
return 0;
}

