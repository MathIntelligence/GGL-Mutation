#include"source.h"
Tatm::Tatm()
{
sidecenter=0;
int i;
for(i=0;i<3;i++) xyz[i]=0;
for(i=0;i<4;i++) {bond[i]=0; bondlen[i]=0;}
strcpy(name," ???\0");
strcpy(type," ???\0");
next=0;
id=0;
balance=0;
hbond=0;
rotate=0;
polarkeep=1;
keep=1;
nonp=0;
eng=0;
nbond=0;
eng=new Eng(this);
ishn=0;
area=0;
rotm=0;
ntat=0;
ispolar=0;
ring=0;
sim=0;
}

Tatm::~Tatm() 
{
if(eng) delete eng;
if(next) delete next; 
eng=0;
next=0;
}

void Tatm::bondleng()
{
  int i;
  for(i=0;i<4;i++)
  {
    if(bond[i]==0)
    {
      bondlen[i]=0;
    }
    else
    {
      bondlen[i]=tres->distance(xyz,bond[i]->xyz);
    }
  }
}

int Tatm::isbackbone()
{
  if(name[1]=='H'&&bond[0]->id<=3) return 1; 
  if(name[1]!='H'&&id<=3) return 1;
  return 0;
}

int Tatm::getnumberhbonds() {

	int n=0;
	for(int i=0;i<4;i++) {
		if(bond[i]&&bond[i]->name[1]=='H') n++;	
	}
	return n;
}

int Tatm::isnear(Tatm *s,int dis)
{
//this program is used to calculate the bond distance between
//two atoms. the maximum distance is 4
int i,j,k,m;
Tatm *a0,*a1,*a2,*a3,*a4,*tatm;

tatm=this;
a0=this;
i=s->tres->id-a0->tres->id;

if(abs(i)>1) return -1;

if(a0==s) return 0;

 
for(i=0;i<tatm->nbond;i++)
{
  a1=bond[i];
  if(a1==0) continue;
  if(a1==a0)continue;
  if(a1==s) return 1;
}
if(dis==1) return -1;

for(i=0;i<tatm->nbond;i++)
{
  a1=bond[i];
  if(a1==0) continue;
  if(a1==a0)continue;
  for(j=0;j<a1->nbond;j++)
  {
    a2=a1->bond[j];
    if(a2==0) continue;
    if(a2==a0||a2==a1) continue;
    if(a2==s) return 2;
  }
}
if(dis==2) return -1;

for(i=0;i<tatm->nbond;i++)
{
  a1=bond[i];
  if(a1==0) continue;
  if(a1==a0)  continue;
  for(j=0;j<a1->nbond;j++)
  {
    a2=a1->bond[j];
    if(a2==0) continue;
    if(a2==a0||a2==a1) continue;

    for(k=0;k<a2->nbond;k++)
    { 
     a3=a2->bond[k];
     if(a3==0) continue;
     if(a3==a1||a3==a0||a3==a2) continue;
     if(a3==s) return 3;
    }
  }
}
if(dis==3) return -1;

for(i=0;i<tatm->nbond;i++)
{
  a1=bond[i];
  if(a1==0) continue;
  if(a1==a0) continue;
  for(j=0;j<a1->nbond;j++)
  {
    a2=a1->bond[j];
    if(a2==0) continue;
    if(a2==a0||a2==a1) continue;

    for(k=0;k<a2->nbond;k++)
    {
     a3=a2->bond[k];
     if(a3==0) continue;
     if(a3==a1||a3==a2||a3==a0) continue;

     for(m=0;m<a3->nbond;m++)
     {
       a4=a3->bond[k];
       if(a4==0) continue;
       if(a4==s) return 4;
     }
    }
  }
}
return -1;;
}
