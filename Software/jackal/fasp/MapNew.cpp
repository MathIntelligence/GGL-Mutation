#include"source.h"

MapNew::MapNew()
{
int i;
pdb=0;
node=0;
for(i=0;i<3;i++) {offmin[i]=0;offmax[i]=0;edge[i]=0;}
grdsiz=2;
probe=0;
leg=0;
}

void MapNew::initial() {
int i;
pdb=0;
node=0;
for(i=0;i<3;i++) {offmin[i]=0;offmax[i]=0;edge[i]=0;}
grdsiz=2;
probe=0;
leg=0;
}

MapNew::~MapNew()
{
int i,j;
pdb=0;
grdsiz=0;
if(node)
{
   for(i=0;i<edge[0];i++)
   {
     for(j=0;j<edge[1];j++)
     {
       delete [] node[i][j];
     }
     delete [] node[i];
   }
   delete [] node;
   node=0;
 }
}

void MapNew::deletenode() {
int i,j;
pdb=0;
grdsiz=0;
if(node)
{
   for(i=0;i<edge[0];i++)
   {
     for(j=0;j<edge[1];j++)
     {
       delete [] node[i][j];
     }
     delete [] node[i];
   }
   delete [] node;
   node=0;
 }
}

void MapNew::clear()
{
 int i,j,k;
 //set initial
 for(i=0;i<edge[0];i++)
   for(j=0;j<edge[1];j++)
     for(k=0;k<edge[2];k++)
       node[i][j][k]=0;
}

void MapNew::restart() {
deletenode();
initial();
}

void MapNew::ready(Pdb *s)
{
 int i,j,k;
 Chn *chn;
 Res *r;
 Atm *a;
 pdb=s;

 if(grdsiz==0) 
 {
   cerr<<"setting grdsiz please..."<<endl;grdsiz=2;
 }
 if(node) 
 {
   for(i=0;i<edge[0];i++) 
   {
     for(j=0;j<edge[1];j++)
     {
       delete [] node[i][j];
     }
     delete [] node[i];
   }
   delete [] node;
   node=0;
 }

 for(i=0;i<3;i++) {offmin[i]=1000000;offmax[i]=-100000;}

 float x=probe+leg;
 for(chn=pdb->chn;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next)
 {
   for(a=r->atm;a;a=a->next) 
   {
      for(i=0;i<3;i++)
      {
        offmin[i]=min(offmin[i],a->xyz[i]-a->tatm->eng->radius-x);
	offmax[i]=max(offmax[i],a->xyz[i]+a->tatm->eng->radius+x);
	/*
        if(a->xyz[i]-a->tatm->eng->radius-probe-leg<offmin[i]) 
           offmin[i]=a->xyz[i]-a->tatm->eng->radius-probe-leg;
        if(a->xyz[i]+a->tatm->eng->radius+probe+leg>offmax[i]) 
           offmax[i]=a->xyz[i]+a->tatm->eng->radius+probe+leg;
	*/
      }
   }
 }
 for(i=0;i<3;i++)
 {
   edge[i]=(int)((offmax[i]-offmin[i])/grdsiz+2);
 }
 node=new int**[edge[0]];
 for(i=0;i<edge[0];i++) 
 { 
    node[i]=new int*[edge[1]];
    for(j=0;j<edge[1];j++) node[i][j]=new int[edge[2]];
 }

 //set initial
 for(i=0;i<edge[0];i++) 
   for(j=0;j<edge[1];j++)
     for(k=0;k<edge[2];k++)
       node[i][j][k]=0;
}

void MapNew::puton()
{
 Chn *chn;
 Res *r;
 Atm *a;
 for(chn=pdb->chn;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next)
 {
   for(a=r->atm;a;a=a->next)
   puton(a);
 }
}

void MapNew::putoff()
{
 Chn *chn;
 Res *r;
 Atm *a;
 for(chn=pdb->chn;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next)
 {
   for(a=r->atm;a;a=a->next)
   putoff(a);
 }
}

void MapNew::puton(Chn *chn)
{
 Res *r;
 Atm *a;
 if(chn->pdb!=pdb) {cerr<<" different pdb"<<endl; return ;}
 for(r=chn->res;r;r=r->next)
 {
   for(a=r->atm;a;a=a->next)
   puton(a);
 }
}

void MapNew::putoff(Chn *chn)
{
 Res *r;
 Atm *a;
 if(chn->pdb!=pdb) {cerr<<" different pdb"<<endl; return ;}
 for(r=chn->res;r;r=r->next)
 {
   for(a=r->atm;a;a=a->next)
   putoff(a);
 }
}


void MapNew::puton(Res *s,int to)
{
 Res *r;
 Atm *a;
 if(s->chn->pdb!=pdb){cerr<<"wrong chain in mapping"<<endl;exit(0);}
 for(r=s;r;r=r->next)
 {
   if(r->id0>to)break;
   for(a=r->atm;a;a=a->next)
   puton(a);
 } 
}

void MapNew::putoff(Res *s,int to)
{
 Res *r;
 Atm *a;
 if(s->chn->pdb!=pdb){cerr<<"wrong chain in mapping"<<endl;exit(0);}
 for(r=s;r;r=r->next)
 {
   if(r->id0>to)break;
   for(a=r->atm;a;a=a->next)
   putoff(a);
 }
}

void MapNew::puton(Atm *a) 
{
int i,j,k;
float d,xyz[3];

int rr0[3],rr[3];

for(i=0;i<3;i++)
{
  rr0[i]=(int)((a->xyz[i]-(probe+a->tatm->eng->radius)-offmin[i])/grdsiz-1);
  if(rr0[i]<0) rr0[i]=0;
  rr[i] =(int)((a->xyz[i]+(probe+a->tatm->eng->radius)-offmin[i])/grdsiz+1);
  if(rr[i]>edge[i]-1) rr[i]=edge[i]-1;
}

float x=0;
for(i=rr0[0];i<=rr[0];i++)
for(j=rr0[1];j<=rr[1];j++)
for(k=rr0[2];k<=rr[2];k++)
{
  xyz[0]=i*grdsiz+offmin[0]-a->xyz[0];
  xyz[1]=j*grdsiz+offmin[1]-a->xyz[1];
  xyz[2]=k*grdsiz+offmin[2]-a->xyz[2];
  d=xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2];
  x=a->tatm->eng->radius+probe;
  if(d<=x*x) node[i][j][k]++;
}
}

void MapNew::putoff(Atm *a)
{
int i,j,k;
float d,xyz[3];

int rr0[3],rr[3];

for(i=0;i<3;i++)
{
  rr0[i]=(int)((a->xyz[i]-(probe+a->tatm->eng->radius)-offmin[i])/grdsiz-1);
  if(rr0[i]<0) rr0[i]=0;
  rr[i] =(int)((a->xyz[i]+(probe+a->tatm->eng->radius)-offmin[i])/grdsiz+1);
  if(rr[i]>edge[i]-1) rr[i]=edge[i]-1;
}

float x=0;
for(i=rr0[0];i<=rr[0];i++)
for(j=rr0[1];j<=rr[1];j++)
for(k=rr0[2];k<=rr[2];k++)
{
  xyz[0]=i*grdsiz+offmin[0]-a->xyz[0];
  xyz[1]=j*grdsiz+offmin[1]-a->xyz[1];
  xyz[2]=k*grdsiz+offmin[2]-a->xyz[2];
  d=xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2];
  x=a->tatm->eng->radius+probe;
  if(d<=x*x&&node[i][j][k]>0) node[i][j][k]--;
}
}



void MapNew::reform(float *xyz,int x,int y,int z)
{
  xyz[0]=x*grdsiz+offmin[0];
  xyz[1]=x*grdsiz+offmin[1];
  xyz[2]=x*grdsiz+offmin[2];
}
 
int MapNew::cover(float *xyz)
{
return cover(xyz[0],xyz[1],xyz[2]);
}
int MapNew::cover(float x,float y,float z)
{
int i,j,k;
//int ii0[3],ii[3];
int a,b,c;
/*
for(i=0;i<3;i++)
{
  ii0[i]=(xyz[i]-offmin[i])/grdsiz;
  if(ii0[i]<0) ii0[i]=0;
  //ii[i] =(xyz[i]-offmin[i])/grdsiz+1;
  //if(ii[i]>edge[i]-1) ii[i]=edge[i]-1;
}
*/

a=max(0,(int)((x-offmin[0])/grdsiz));
b=max(0,(int)((y-offmin[1])/grdsiz));
c=max(0,(int)((z-offmin[2])/grdsiz));


for(i=a;i<=a+1;i++)
for(j=b;j<=b+1;j++)
for(k=c;k<=c+1;k++)
{
 if(node[i][j][k]==0) return 0;
}
return 1;
}

int MapNew::cover(Atm *a)
{
int i,j,k;
int ii0[3],ii[3];

for(i=0;i<3;i++)
{
  ii0[i]=(int)((a->xyz[i]-a->tatm->eng->radius-offmin[i])/grdsiz);
  if(ii0[i]<0) ii0[i]=0;
  ii[i] =(int)((a->xyz[i]+a->tatm->eng->radius-offmin[i])/grdsiz+1);
  if(ii[i]>edge[i]-1) ii[i]=edge[i]-1;
}

for(i=ii0[0];i<=ii[0];i++)
for(j=ii0[1];j<=ii[1];j++)
for(k=ii0[2];k<=ii[2];k++)
{
  if(node[i][j][k]==0) return 0;
}

return 1;

}

