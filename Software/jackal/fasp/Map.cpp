#include"source.h"

Map::Map()
{
int i;
pdb=0;
node=0;
for(i=0;i<3;i++) {offmin[i]=0;offmax[i]=0;edge[i]=0;}
grdsiz=2;
}

Map::~Map()
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

void Map::clear()
{
 int i,j,k;
 //set initial
 for(i=0;i<edge[0];i++)
   for(j=0;j<edge[1];j++)
     for(k=0;k<edge[2];k++)
       node[i][j][k]=0;
}
void Map::ready(Pdb *s)
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
 for(chn=pdb->chn;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next)
 {
   for(a=r->atm;a;a=a->next) 
   {
      for(i=0;i<3;i++)
      {
        if(a->xyz[i]-a->tatm->eng->radius<offmin[i]) 
           offmin[i]=a->xyz[i]-a->tatm->eng->radius;
        if(a->xyz[i]+a->tatm->eng->radius>offmax[i]) 
           offmax[i]=a->xyz[i]+a->tatm->eng->radius;
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

void Map::puton()
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

void Map::puton(Chn *chn)
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


void Map::puton(Res *s,int to)
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

void Map::puton(Atm *a) 
{
int i,j,k;
int ii0[3],ii[3];
float d,xyz[3];

for(i=0;i<3;i++)
{
  ii0[i]=(int)((a->xyz[i]-a->tatm->eng->radius*0.70710678-offmin[i])/grdsiz+1);
  if(ii0[i]<0) ii0[i]=0;
  ii[i] =(int)((a->xyz[i]+a->tatm->eng->radius*0.70710678-offmin[i])/grdsiz);
  if(ii[i]>edge[i]-1) ii[i]=edge[i]-1;
}

for(i=ii0[0];i<=ii[0];i++)
for(j=ii0[1];j<=ii[1];j++)
for(k=ii0[2];k<=ii[2];k++)
{
  node[i][j][k]=1;
}  

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
  if(node[i][j][k]==1) continue;
  xyz[0]=i*grdsiz+offmin[0]-a->xyz[0];
  xyz[1]=j*grdsiz+offmin[1]-a->xyz[1];
  xyz[2]=k*grdsiz+offmin[2]-a->xyz[2];
  d=xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2];
  if(d<=a->tatm->eng->radius*a->tatm->eng->radius)
  node[i][j][k]=1;
}
}

void Map::reform(float *xyz,int x,int y,int z)
{
  xyz[0]=x*grdsiz+offmin[0];
  xyz[1]=x*grdsiz+offmin[1];
  xyz[2]=x*grdsiz+offmin[2];
}
 
int Map::cover(float *xyz)
{
int i,j,k;
int ii0[3],ii[3];

for(i=0;i<3;i++)
{
  ii0[i]=(int)((xyz[i]-offmin[i])/grdsiz);
  if(ii0[i]<0) ii0[i]=0;
  ii[i] =(int)((xyz[i]-offmin[i])/grdsiz+1);
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

int Map::cover(Atm *a)
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

