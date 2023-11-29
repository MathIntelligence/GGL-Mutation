#include"source.h"

void Lattice::initial() {
ishet=0;
pdb=0;
busket=0;
atom=0;
obtain=0;
hash=0;
nget=0;
grdsiz=0;
radall=100;
radmax=0;
flag=0;
hydrogen=1;
offset=0;
//atm=0;
}
Lattice::Lattice()
{
initial();
}

Lattice::~Lattice()
{
int i;
if(atom)
{
for(i=0;i<hash/2;i++) if(atom[i])delete atom[i];
if(atom)delete [] atom;atom=0;
}
if(busket) delete [] busket;
if(obtain) delete [] obtain;
hash=0;
busket=0;
pdb=0;
atom=0;
obtain=0;
nget=0;
radall=0;
grdsiz=0;
flag=0;
//atm=0;
}
 
void Lattice::clean()
{
int i;
if(atom)
{
for(i=0;i<hash/2;i++) if(atom[i])delete atom[i];
delete [] atom;
}
if(busket) delete [] busket;
if(obtain) delete [] obtain;
hash=0;
obtain=0;
busket=0;
pdb=0;
atom=0;
nget=0;
//radall=0;
//grdsiz=0;
//flag=0;
}


void Lattice::ready(Pdb *s)
{
Chn *chn;
Res *res;
Atm *atm;
int i,j;

//new bugs found and fixed at 1/1/2002
clean();
//

if(grdsiz<=0) {cerr<<"set the grdsiz for lattice!\n";exit(0);}
i=1;
pdb=s; 
radmax=0; 
offset=s->chn->res->atm->id0;

for(chn=s->chn;chn;chn=chn->next)
for(res=chn->res;res;res=res->next)
for(atm=res->atm;atm;atm=atm->next)
{
if(atm->tatm&&atm->tatm->eng->radius>radmax) radmax=atm->tatm->eng->radius;
i=atm->id0;
}

for(chn=s->hetatm;chn&&ishet;chn=chn->next)
for(res=chn->res;res;res=res->next)
for(atm=res->atm;atm;atm=atm->next)
{
if(atm->tatm&&atm->tatm->eng->radius>radmax) radmax=atm->tatm->eng->radius;
i=atm->id0;
}

i=i-offset+1;
hash=i*2;
if(i==0) {
	cerr<<"offset =0"<<endl;exit(0);
}
if(hash>0) atom=new Cell*[hash/2];
for(j=0;j<hash/2;j++)atom[j]=new Cell;
if(hash>0) busket=new Cell*[hash];
if(i>0) obtain=new Cell*[i];
//if(TRES.logg>3) s->writeold("tt.pdb");
for(chn=s->chn;chn;chn=chn->next)
for(res=chn->res;res;res=res->next)
for(atm=res->atm;atm;atm=atm->next)
{
 i=atm->id0-offset;
 //if(TRES.logg>3) cerr<<"lattice .."<<atm->id0<<endl;
 if(i<0) {cerr<<"offset <0"<<endl;exit(0);}
 if(atom) atom[i]->atm=atm;
}

for(chn=s->hetatm;chn&&ishet;chn=chn->next)
for(res=chn->res;res;res=res->next)
for(atm=res->atm;atm;atm=atm->next)
{
 i=atm->id0-offset;
 //if(TRES.logg>3) cerr<<"lattice .."<<atm->id0<<endl;
 if(i<0) {cerr<<"offset <0"<<endl;exit(0);}
 if(atom) atom[i]->atm=atm;
}


for(i=0;i<hash;i++) if(busket) busket[i]=0;
}

int Lattice::indx(float *s)
{
int i,j,k,n;
float x;

x=grdsiz/2;

if(s[0]>=0) i=(int)((s[0]+x)/grdsiz);
else i=(int)((s[0]-x)/grdsiz);

if(s[1]>=0) j=(int)((s[1]+x)/grdsiz);
else j=(int)((s[1]-x)/grdsiz);

if(s[2]>=0) k=(int)((s[2]+x)/grdsiz);
else k=(int)((s[2]-x)/grdsiz);

n=i*250000+j*500+k;

return n;
}

int Lattice::indx(float s)
{
float x;
int i;

x=grdsiz/2;

if(s>=0) i=(int)((s+x)/grdsiz);
else i=(int)((s-x)/grdsiz);

return i;
}

void Lattice::puton(Res *s,int a)
{
  Res *r;
  for(r=s;r;r=r->next)
  {
    if(r->id0>a) break;
    puton(r);
  }
}

void Lattice::puton(Res *s)
{
  puton(s,0,100);
}
void Lattice::puton(Chn *s)
{
 Res *r;
 if(s->pdb!=pdb)return;
 for(r=s->res;r;r=r->next) puton(r);
}

void Lattice::puton(Pdb *s)
{
 Chn *chn;
 Res *r;
 if(s!=pdb) return;
 for(chn=s->chn;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next) puton(r);
}

void Lattice::putonhetatm(Pdb *s)
{
 if(ishet==0) return;
 Chn *chn;
 Res *r;
 if(s!=pdb) return;
 Atm *a;
 for(chn=s->hetatm;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next) 
 for(a=r->atm;a;a=a->next) puton(a); 
}
 
void Lattice::puton()
{
 Chn *chn;
 Res *r;
 Pdb *s=pdb;
 if(s!=pdb) return;
 for(chn=s->chn;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next) puton(r);
 //
 Atm *a;
 for(chn=s->hetatm;chn&&ishet;chn=chn->next)
 for(r=chn->res;r;r=r->next) 
 for(a=r->atm;a;a=a->next) puton(a);
}

void Lattice::putonhetatm()
{
 if(ishet==0) return;
 Chn *chn;
 Res *r;
 Pdb *s=pdb;
 if(s!=pdb) return;
 //for(chn=s->chn;chn;chn=chn->next)
 //for(r=chn->res;r;r=r->next) puton(r);
 //
 Atm *a;
 for(chn=s->hetatm;chn&&ishet;chn=chn->next)
 for(r=chn->res;r;r=r->next) 
 for(a=r->atm;a;a=a->next) puton(a);
}

void Lattice::puton(Chn *s,int a,int b,int c)
{
Res *res;
if(s->pdb!=pdb)return;
for(res=s->res;res;res=res->next)
puton(res,a,b,c);
}

void Lattice::puton(Pdb *s,int a,int b,int c)
{
Res *res;
if(s!=pdb)return;
Chn *chn;
for(chn=s->chn;chn;chn=chn->next)
for(res=chn->res;res;res=res->next)
puton(res,a,b,c);
}

void Lattice::puton(Chn *s,int a,int b)
{
Res *res;
if(s->pdb!=pdb) return;
for(res=s->res;res;res=res->next)
puton(res,a,b);
}

void Lattice::puton(Pdb *s,int a,int b)
{
Res *res;
if(s!=pdb) return;
Chn *chn;
for(chn=s->chn;chn;chn=chn->next)
for(res=chn->res;res;res=res->next)
puton(res,a,b);
}

void Lattice::puton(Res *s,int a,int b,int c)
{
Atm *atm,*atm1;
int i;
if(s->chn->pdb!=pdb) return;
for(atm=s->atm;atm;atm=atm->next)
{
  if(atm->tatm->name[1]=='H') continue;
  if(atm->tatm->rotm<a||atm->tatm->rotm>b) continue;
  puton(atm);
  for(i=0;i<atm->tatm->nbond;i++)
  {
    atm1=atm->bond[i];
    if(atm1==0) continue;
    if(atm1->tatm->name[1]!='H') continue;
    if(hydrogen==0) continue;
    puton(atm1);
  }
}
}

void Lattice::puton(Res *s,int a,int b)
{
Atm *atm,*atm1;
int i;
if(s->chn->pdb!=pdb) return;
for(atm=s->atm;atm;atm=atm->next)
{
  if(atm->tatm->name[1]=='H') continue;
  if(atm->tatm->id<a||atm->tatm->id>b) continue;
  puton(atm);
  for(i=0;i<atm->tatm->nbond;i++)
  {
    atm1=atm->bond[i];
    if(atm1==0) continue;
    if(atm1->tatm->name[1]!='H') continue;
    if(hydrogen==0) continue;
    puton(atm1); 
  }
}
}

void Lattice::puton(Atm *s)
{
int i,j,n;
Cell *cell_temp;
if(s->res->chn->pdb!=pdb) return;
j=s->id0-offset;
if(atom[j]->used!=-1) return; //already on
if(s->good==0) return;  //this atom should not be used!
n=indx(s->xyz);
atom[j]->id=n;
i=(n%hash+hash)%(hash-1);
atom[j]->used=i;
if(busket[i]==0) busket[i]=atom[j];
else
{
cell_temp=busket[i];
while(cell_temp->next)cell_temp=cell_temp->next;
cell_temp->next=atom[j];
}
}

int Lattice::putoff()
{
int i,j,k;
j=hash/2;
k=0;
for(i=0;i<j;i++)
{
 atom[i]->next=0;
 if(atom[i]->used!=-1)k++;
 atom[i]->id=0;
 atom[i]->used=-1;
}
for(i=0;i<hash;i++) busket[i]=0;
return k;
}

int Lattice::putoff(Res *s,int n)
{
  Res *r;
  int i;
  i=0;
  for(r=s;r;r=r->next)
  {
    if(r->id0>n) break;
    i+=putoff(r);
  }
  return i;
}
 
int Lattice::putoff(Res *s)
{
  return putoff(s,0,1000000);
}

int Lattice::putoff(Chn *s,int a, int b,int c)
{
Res *res;
int j,n;
n=0;
if(s->pdb!=pdb) return 0;
if(a==0&&b>50)
{
 j=putoff(s,a,b);
 return j;
}
else
{
for(res=s->res;res;res=res->next)n+=putoff(res,a,b,c);
}
return n;
}


int Lattice::putoff(Pdb *s,int a, int b,int c)
{
Chn *chn;
Res *res;
int j,n;
n=0;
if(s!=pdb) return 0;
if(a==0&&b>50)
{
 j=putoff(s,a,b);
 return j;
}
else
{
for(chn=s->chn;chn;chn=chn->next)
for(res=chn->res;res;res=res->next)n+=putoff(res,a,b,c);
}
return n;
}

int Lattice::putoff(Chn *s,int a,int b)
{
Res *res;
int i,j,n;
if(s->pdb!=pdb) return 0;
n=0;
if(a==0&&b>50)
{
j=hash/2;
for(i=0;i<j;i++)
{
 atom[i]->next=0;
 atom[i]->id=0;
 atom[i]->used=-1;
}
for(i=0;i<hash;i++) busket[i]=0;
return j;
}
else
{
for(res=s->res;res;res=res->next)n+=putoff(res,a,b);
}
return n;
}

int Lattice::putoff(Pdb *s,int a,int b)
{
Chn *chn;
Res *res;
int i,j,n;
if(s!=pdb) return 0;
n=0;
if(a==0&&b>50)
{
j=hash/2;
for(i=0;i<j;i++)
{
 atom[i]->next=0;
 atom[i]->id=0;
 atom[i]->used=-1;
}
for(i=0;i<hash;i++) busket[i]=0;
return j;
}
else 
{
for(chn=s->chn;chn;chn=chn->next)
for(res=chn->res;res;res=res->next)n+=putoff(res,a,b);
}
return n;
}

int Lattice::putoff(Res *s, int a,int b,int c)
{
Atm *atm;
int i;
if(s->chn->pdb!=pdb) return 0;
i=0;
for(atm=s->atm;atm;atm=atm->next)
{
if(atm->tatm->name[1]=='H'&&hydrogen==0) continue;
if(atm->tatm->rotm>=a&&atm->tatm->rotm<=b)
i+=putoff(atm);
}
return i;
}

int Lattice::putoff(Res *s,int a,int b)
{
Atm *f,*atm;
int i;
if(s->chn->pdb!=pdb) return 0;
i=0;
for(atm=s->atm;atm;atm=atm->next)
{
  if(atm->tatm->name[1]=='H') continue;
  if(atm->tatm->id<a||atm->tatm->id>b) continue;
  i+=putoff(atm);
  for(i=0;i<atm->tatm->nbond;i++)
  {
    f=atm->bond[i];
    if(f==0) continue;
    if(f->tatm->name[1]!='H') continue;
    if(hydrogen==0) continue;
    i+=putoff(f);
  }
}
return i;
}

int Lattice::putoff(Atm *s)
{
int i,j;
Cell *cell0,*cell1;
if(s->res->chn->pdb!=pdb) return 0;
j=s->id0-offset;
i=atom[j]->used;
atom[j]->id=0;
atom[j]->used=-1;
if(i==-1) return 0;
cell0=atom[j];
cell1=busket[i];
if(cell0==cell1) busket[i]=cell0->next;
else
{
 while(cell1->next!=cell0) cell1=cell1->next;
 cell1->next=cell0->next;
}
cell0->next=0;
return 1;
}

int Lattice::exist(Atm *s)
{
int i,j;
if(s->res->chn->pdb!=pdb) return 0;
j=s->id0-offset;
i=atom[j]->used;
if(i==-1) return 0;
return 1;
}

void Lattice::printoutLattice() {

	int i;
	for(i=0;i<nget;i++) {
		Atm *a=obtain[i]->atm;
		cerr<<a->res->name<<a->res->id0<<" "<<a->name<<" "<<a->id0<<" : "<<obtain[i]->id<<endl;
	}
	pdb->write("se");
} 
int Lattice::getcell(float x,float y,float z,float d) {
float xyz[3];

xyz[0]=x;
xyz[1]=y;
xyz[2]=z;
return getcell(xyz,d);

}
int Lattice::getcell(float *co,float d)
{
int i;
Atm s;
s.transfer(co,1);
i=getcell(&s,d);
return i;
}

int Lattice::getcell(Atm *s,float d)
{
Cell *temp;
int i0,h,f,i,j,k,n,m,close;
float d1;
int x[6],co[3],ff;
nget=0;
close=0;
i0=hash/2-1;
ff=0;
//atm=s; 
if(d>radall) { ff=1;d=3.; }
d1=d+0.001; 
for(i=0;i<3;i++)
{
j=2*i;
x[j]  =indx(s->xyz[i]-d1);
x[j+1]=indx(s->xyz[i]+d1);
co[i]=(x[j]+x[j+1])/2;
}

for(i=0;i<6;i++){
	if(fabs(x[i])>100000) {
		cerr<<"warning! atom position corrupted: "<<s->res->name<<s->res->id0<<":"<<s->name<<"-->"; 
		cerr<<s->xyz[0]<<" "<<s->xyz[1]<<" "<<s->xyz[2]<<endl;
		nget=0;
		obtain[nget]=0;
		return nget;
	}	
}
for(i=x[0];i<=x[1];i++)
for(j=x[2];j<=x[3];j++)
for(k=x[4];k<=x[5];k++)
{
  n=i*250000+j*500+k;
  m=(n%hash+hash)%(hash-1);
  h=co[0]-i;
  f=1;
  if(h<=1&&h>=-1)
  {
    h=co[1]-j;
    if(h<=1&&h>=-1) 
    {
      h=co[2]-k; 
      if(h<=1&&h>=-1) f=0;
    }
  }
  if(f==0)
  {
     for(temp=busket[m];temp;temp=temp->next)
     {
       if(temp->id!=n) continue;
       if(flag==0)if(s->res->id0==temp->atm->res->id0) continue;
       if(flag==1)if(s->id==temp->atm->id) continue;
       if(flag==2&&s->res->id0==temp->atm->res->id0)
       {
         if(temp->atm->res->chn->ishet==0&&s->res->chn->ishet==0&&s->isnear(temp->atm,3)!=-1) continue;
       }
       if(flag==3)if(abs(s->res->id-temp->atm->res->id)<=1) continue;

       obtain[i0-close]=temp;
       if(ff==1) temp->flag=1;
       close++;
     }
  }
  else
  {
     for(temp=busket[m];temp;temp=temp->next)
     {
       if(temp->id!=n) continue;
       if(flag==0) if(s->res->id0==temp->atm->res->id0) continue;
       if(flag==1) if(s->id==temp->atm->id) continue;
       if(flag==2&&s->res->id0==temp->atm->res->id0)
       {
         if(temp->atm->res->chn->ishet==0&&s->res->chn->ishet==0&&s->isnear(temp->atm,3)!=-1) continue;
       }
       if(flag==3)if(abs(s->res->id-temp->atm->res->id)<=1) continue;
       obtain[nget]=temp;
       if(ff==1) temp->flag=1;
       nget++;
     }
  }
}
 
for(i=i0-close+1;i<=i0;i++)
{
  obtain[nget++]=obtain[i];
}
obtain[nget]=0;
if(ff==0) return nget;
for(m=0;m<hash;m++)
{
  for(temp=busket[m];temp;temp=temp->next)
  {
    if(temp->flag==1) {temp->flag=0;continue;}
    if(flag==0) if(s->res->id0==temp->atm->res->id0) continue;
    if(flag==1&&s->id==temp->atm->id) continue;
    if(flag==2&&s->res->id0==temp->atm->res->id0)
    {
      if(temp->atm->res->chn->ishet==0&&s->res->chn->ishet==0&&s->isnear(temp->atm,3)!=-1) continue;
    }
    if(flag==3)if(abs(s->res->id-temp->atm->res->id)<=1) continue;
    obtain[nget++]=temp;
  }
}
obtain[nget]=0;
return nget;
}

int Lattice::cutoff(float *co,float d)
{
int i;
Atm s;
s.transfer(co,1);
i=cutoff(&s,d);
return i;
}

int Lattice::cutoff(float *co)
{
int i;
Atm s;
s.transfer(co,1);
i=cutoff(&s);
return i;
}

int Lattice::cutoff(Atm *s,float d)
{
int i,j;
float c;
j=0;
if(d>radall) return 0;
d=d*d;
for(i=0;i<nget;i++)
{
 c=TRES.distsqr(s->xyz,obtain[i]->atm->xyz);
 if(c>d) continue;
 obtain[j]=obtain[i];
 j++;
}
i=nget-j;
nget=j;
return i;
}

int Lattice::cutoff(Atm *s)
{
int i,j;
float c,d;
j=0;
for(i=0;i<nget;i++)
{
 d=s->tatm->eng->radius+obtain[i]->atm->tatm->eng->radius;
 d=d*d;
 c=TRES.distsqr(s->xyz,obtain[i]->atm->xyz);
 if(c>d) continue;
 obtain[j]=obtain[i];
 j++;
}
i=nget-j;
nget=j;
return i;
}
 
int Lattice::cover(float *s,float d0)
{
Cell *a;
int i;
float c,d;
for(i=0;i<nget;i++)
{
 d=d0+obtain[i]->atm->tatm->eng->radius;
 d=d*d;
 c=TRES.distsqr(s,obtain[i]->atm->xyz);
 if(c>d) continue;
 a=obtain[0];
 obtain[0]=obtain[i];
 obtain[i]=a;
 return 1;
}
return 0;
}

int Lattice::resonly()
{
int i,j,k,n;

if(nget==0)
{
 return 0; 
}
j=1;

for(i=1;i<nget;i++)
{
  n=0;
  for(k=0;k<j;k++)
  if(obtain[i]->atm->res->id0==obtain[k]->atm->res->id0) {n=1;break;}
  if(n==1) continue;
  else {obtain[j]=obtain[i];j++;} 
}
i=nget-j;
nget=j;
return i;
}

int  Lattice::getrescutcharge(Atm *s) {
	
	if(s==0) {
		nget=0;
		return 0;
	}
	int *reso=new int[2*hash+100];

	int i,m;
	for(i=0;i<2*hash+100;i++) reso[i]=0;
	
	for(i=nget-1;i>=0;i--)
	{
		Atm *b=obtain[i]->atm;
		if(b==0) continue;
		if(s->res->id0==b->res->id0) continue;
		reso[b->res->id0]=1;
	}				
	

	Chn *c;
	Res *r;
	
	for(c=pdb->chn;c;c=c->next) 
	for(r=c->res;r;r=r->next) {

		if(reso[r->id0]==1) continue;
		if(r->id0==s->res->id0) continue;
		Tatm *t;
		float e=0;
		for(t=r->tres->tatm;t;t=t->next) {
			e+=t->eng->charge;
		}
		if(fabs(e)>0.01) reso[r->id0]=1;
	}
	
		
	Cell *temp;


	nget=0;
	for(m=0;m<hash;m++)
	{
  		for(temp=busket[m];temp;temp=temp->next)
  		{
			Atm *b=temp->atm;
			if(b==0) continue;
			if(reso[b->res->id0]==0) continue;
			if(s->res->id0==b->res->id0) continue;
    			obtain[nget++]=temp;
			//float d=TRES.distance(b,s);
			//if(d<3.0) cerr<<b->res->name<<b->res->oid<<" "<<b->name<<":"<<s->res->name<<s->oid<<" "<<s->name<<" "<<d<<endl;
  		}
	}

	obtain[nget]=0;
	delete [] reso;
	return nget;
}

int  Lattice::getrescut(Atm *s) {
	
	if(s==0) {
		nget=0;
		return 0;
	}
	int *reso=new int[2*hash+100];

	int i,m;
	for(i=0;i<2*hash+100;i++) reso[i]=0;
	
	for(i=nget-1;i>=0;i--)
	{
		Atm *b=obtain[i]->atm;
		if(b==0) continue;
		if(s->res->id0==b->res->id0) continue;
		reso[b->res->id0]=1;
	}				
	
	/*
	Chn *c;
	Res *r;
	
	for(c=pdb->chn;c;c=c->next) 
	for(r=c->res;r;r=r->next) {

		if(reso[r->id0]==1) continue;
		if(r->id0==s->res->id0) continue;
		Tatm *t;
		float e=0;
		for(t=r->tres->tatm;t;t=t->next) {
			e+=t->eng->charge;
		}
		if(fabs(e)>0.01) reso[r->id0]=1;
	}
	*/
		
	Cell *temp;


	nget=0;
	for(m=0;m<hash;m++)
	{
  		for(temp=busket[m];temp;temp=temp->next)
  		{
			Atm *b=temp->atm;
			if(b==0) continue;
			if(reso[b->res->id0]==0) continue;
			if(s->res->id0==b->res->id0) continue;
    			obtain[nget++]=temp;
			//float d=TRES.distance(b,s);
			//if(d<3.0) cerr<<b->res->name<<b->res->oid<<" "<<b->name<<":"<<s->res->name<<s->oid<<" "<<s->name<<" "<<d<<endl;
  		}
	}

	obtain[nget]=0;
	delete [] reso;
	return nget;
}


void Lattice::addifnotexist(Atm *a) {

	int i;
	for(i=0;i<nget;i++) {
		if(obtain[i]==0) continue;
		Atm *aa1=obtain[i]->atm;
		if(a==aa1) return;
	}
	i=a->id0-offset;
	obtain[nget]=atom[i];
	nget++;
}
