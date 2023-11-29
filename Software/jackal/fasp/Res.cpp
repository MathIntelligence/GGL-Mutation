#include"source.h"
void Res::initial()
{
rotenergy=0;
rotdegree=0;
rotchance=0;
addhanyway=0;
chn=0;
next=0;
more=0;
atm=0;
tres=0;
name='?';
id=0;
id0=0;
number=0;
nhydr=0;
nummore=0;
regular=1;
area=0;
sec='-';
temp=0;
nemp=0;
flag=0;
hbond=0;
oid=0;
last=0;
ambgt=0;
energy=0;
rescard[0]='\0';
}

Res::Res()
{
initial(); 
} 

Res::Res(Res *s) // create the identical residue
{

//set initial
 
initial();

//copy from residue s
rotdegree=s->rotdegree;
rotenergy=s->rotenergy;
tres=s->tres;
name=s->name;
id=s->id;
id0=s->id0;
oid=s->oid;
sec=s->sec;
strcpy(rescard,s->rescard);
number=s->number;
nhydr=s->nhydr;
nummore=s->nummore;
regular=s->regular;
rotchance=s->rotchance;
//create all Atoms
if(s->atm) atm=s->atm->cloneatm();

//set atom's residues
Atm *a;
for(a=atm;a;a=a->next) {
 a->res=this;
}

//create more residues
if(s->more) more=new Res(s->more);
}


Res *Res::cloneres() {
  Res *t=new Res(this);
  if(next) t->next=next->cloneres();
  else t->next=0;
  return t;
}

Res::~Res()
{
if(atm)  {delete atm;atm=0;}
if(next) {delete next;next=0;}
if(more) {delete more;more=0;}
if(temp) {delete [] temp;nemp=0;temp=0;}
if(hbond) {delete hbond;hbond=0;}
if(ambgt) {delete ambgt;ambgt=0;}
}


void Res::write(FILE *fp,int f)
{
char line[256];
Atm *a;

if(f==0) {sprintf(line,"%d%c",id,name);fprintf(fp,"%s",line);}
if(f==1) {sprintf(line,"%c %d%c",chn->id,id,name);fprintf(fp,"%s",line);}
if(f==2) {sprintf(line,"%s %c %d%c",chn->pdb->name,chn->id,id,name);fprintf(fp,"%s",line);}

if(f==4) 
{
for(a=atm;a;a=a->next) a->write(fp,4);
fflush(fp);
return; 
}
if(f==40) {
for(a=atm;a;a=a->next) a->write(fp,40);
fflush(fp);
return;
}
if(f==5)
{
for(a=atm;a;a=a->next) a->write(fp,5);
fflush(fp);
return;
}
fflush(fp);
return;
}


void Res::write(FILE *fp)
{
Atm *atm_temp;
for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
atm_temp->write(fp);
}
 
void Res::writeold(FILE *fp)
{
Atm *atm_temp;
for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
atm_temp->writeold(fp);
}

void Res::writerescard(FILE *fp)
{
Atm *atm_temp;
for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
atm_temp->writerescard(fp);
}

void Res::writeoldmore(FILE *fp)
{
Atm *atm_temp;
Res *r0;

int nmore=0;
for(r0=this;r0;r0=r0->more) nmore++;
if(TRES.logg==-1) cerr<<name<<oid<<":" <<nmore<<endl;
int nn=1;


for(r0=this;r0;r0=r0->more) {
	r0->oid=oid;
	r0->chn=chn;
	if(r0!=this) nn++;
	char ch=' ';
	if(nn==0) ch=' ';
	else if(nn>=1&&nn<=26) ch=64+nn;
	else if(nn>=27&&nn<=52) ch=nn-26+96;
	else ch='*';
	if(nmore==1) ch=' ';
	for(atm_temp=r0->atm;atm_temp;atm_temp=atm_temp->next)
	{
		Atm *s=this->isatm(atm_temp->name);
		if(s==0) continue;
		atm_temp->oid=s->oid;
		atm_temp->writeold(fp,ch);
	}
}
}



void Res::write(FILE *fp,int a,int b)
{
Atm *atm_temp;
for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
if(atm_temp->tatm->id>=a&&atm_temp->tatm->id<=b)
atm_temp->write(fp);
}


Atm *Res::operator[](int k)
{
Atm *temp;
for(temp=atm;temp;temp=temp->next) if(temp->tatm->id==k) return temp;
return 0;
}

Atm *Res::isatmid(int k)
{
Atm *temp;
for(temp=atm;temp;temp=temp->next) if(temp->tatm->id==k) return temp;
return 0;
}

void Res::setbackbonetorsion() {

Atm *temp;
for(temp=atm;temp;temp=temp->next) {
	if(temp->tatm->id>=4) break;
	if(temp->tatm->id==1||temp->tatm->id==2) {
		temp->chi=temp->gettorsionangle();
	}
}
}
Atm *Res::operator[](char *s)
{
Atm *temp;
for(temp=atm;temp;temp=temp->next) 
if(strncmp(temp->tatm->name,s,4)==0) return temp;
return 0;
}

Atm *Res::isatm(char *s)
{
Atm *temp;
for(temp=atm;temp;temp=temp->next)
if(strncmp(temp->name,s,4)==0) return temp;
return 0;
}

Atm *Res::isatm(int n)
{
Atm *temp;
for(temp=atm;temp;temp=temp->next)
if(temp->tatm->ntat==n) return temp;
return 0;
}

void Res::dihedral() {
  dihedral(0);
}

void Res::dihedral(FILE *fp)
{
Atm *atm_temp;
//if(regular==0) return;
if(fp&&chn)fprintf(fp,"%c%5i %c  ",name,id+chn->start,sec);
else if(fp) fprintf(fp,"%c%5i %c  ",name,-1000,sec);
for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
atm_temp->dihedral(fp);
if(fp)fprintf(fp,"\n");
}

void Res::dihedraloption(FILE *fp,int m)
{
Atm *atm_temp;
//if(regular==0) return;
if(fp)fprintf(fp,"%c%5i %c  ",name,oid,chn->id);
for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next){
	if((m==0||m==4)&&atm_temp==atm) {//omega
		float d=atm_temp->gettorsionangle();
		if(fp)fprintf(fp,"  %6.1f",d);
	}
	else if((m==2||m==0||m==1)&&atm_temp->tatm->id>=1&&atm_temp->tatm->id<=2) {//psi,phi
		float d=atm_temp->gettorsionangle();
		if(fp)fprintf(fp,"  %6.1f",d);
	}
	else if((m==3||m==0||m==1)&&atm_temp->tatm->id>3&&atm_temp->tatm->rotate==1) {//side-chain
		float d=atm_temp->dihedral(0);
		if(fp)fprintf(fp,"  %6.1f",d);
	}
}
if(fp)fprintf(fp,"\n");
}

void Res::dihedral(FILE *fp,int m)
{
Atm *atm_temp;
int i;
float u=0;
u=180;
if(m==-4&&fp)
{
  fprintf(fp,"%s    ",tres->name3);
  i=0;
  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
  {
    if(atm_temp->tatm->id<4||atm_temp->tatm->rotate!=1) continue;
    i++;
    atm_temp->dihedral(0);
    fprintf(fp,"%6.1f ",atm_temp->chi);
  }
  if(i==0)
  fprintf(fp,"%6.1f",u);
  //fprintf(fp,"\n");
  return;
}

if(m==-5&&fp)
{
  fprintf(fp,"%s    ",tres->name3);
  i=0;
  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
  {
    if(atm_temp->tatm->id>=4||atm_temp->tatm->rotate!=1) continue;
    i++;
    atm_temp->dihedral(0);
    fprintf(fp,"%6.1f ",atm_temp->chi);
  }
  if(i==0) fprintf(fp,"%6.1f",u);
  fprintf(fp,"\n");
  return;
}


//if(regular==0) return;
if(fp&&m==-3)fprintf(fp,"%c ",name);
else if(fp)fprintf(fp,"%c%5i  ",name,id);
for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
atm_temp->dihedral(fp,m);
if(fp)fprintf(fp,"\n");
}


void Res::transfer(Res *s)
{
Atm *atm_temp,*atm_temp0;
int i;
if(s->name!=name) return;
for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
{
  atm_temp0=(*s)[atm_temp->tatm->id];
  if(atm_temp0==0) continue;
  for(i=0;i<3;i++) atm_temp->xyz[i]=atm_temp0->xyz[i];
}
}


int Res::getid()
{
Res *a,*b;
int i;
i=0;
b=(*chn)[id0];
for(a=b;a;a=a->more)
{
if(a==this) return i;
i++;
}
return -1;
}

Res *Res::isres(int n)
{
Res *s,*a;
int i;
i=0;
a=(*this->chn)[id0];
for(s=a;s;s=s->more)
if(n==i) return s;
else i++;
return 0;
}

void Res::transfer(float s)
{
  Atm *a;
  for(a=atm;a;a=a->next)a->transfer(s);
}

void Res::transfer(float *co,int f)
{
int i,j;
Atm *a; 
for(a=atm;a;a=a->next)
{
i=a->tatm->id;
for(j=0;j<3;j++) 
if(f) a->xyz[j]=co[i*3+j];
else co[i*3+j]=a->xyz[j];
}
}

void Res::transfertemp(Res *co)
{
if(temp) delete [] temp;
temp=0;
if(temp==0) temp=new float[9];
int j;
for(j=0;j<9;j++) temp[j]=co->temp[j];
}

void Res::transfertemp(float *co,int f)
{
int j;
for(j=0;j<9;j++) 
if(f) temp[j]=co[j];
else co[j]=temp[j];
}



void Res::transfer(Atm *at, float *co,int f)
{
int i,j,k;
Atm *a,*c;
i=0;
for(a=at;a;a=a->next)
{
if(a->tatm->name[1]=='H') break;
i=a->tatm->id;
for(j=0;j<3;j++)
if(f) a->xyz[j]=co[i*3+j];
else co[i*3+j]=a->xyz[j];
for(k=1;k<a->tatm->nbond;k++)
{
 c=a->bond[k];
 if(c==0) continue;
 i=c->tatm->id;
 if(c->tatm->name[1]!='H') continue;
 for(j=0;j<3;j++)
 if(f) c->xyz[j]=co[i*3+j];
 else co[i*3+j]=c->xyz[j];
}
}
}

int Res::rotable(Atm **s,int n1,int n2)
{
int i;
Atm *a;
i=0;
for(a=atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') break;
  if(a->tatm->id>=n1&&a->tatm->id<=n2)
  if(a->tatm->rotate==1) s[i++]=a;
}
  s[i]=0;
  return i;
}


int Res::rotable(int n1,int n2)
{
int i;
Atm *a;
i=0;
for(a=atm;a;a=a->next)
{
  if(a->tatm->id>=n1&&a->tatm->id<=n2)
  if(a->tatm->rotate==1) i++;
}
  return i;
}


int Res::getuse(char *s,int m,int n)
{
Atm *a;
int i;
i=0;
for(a=atm;a;a=a->next)
if(strncmp(a->name+m,s+m,n)==0) i++;
return i;
}

int Res::getuseambgt(char *s,int m,int n)
{
Atm *a;
int i;
i=0;
for(a=ambgt;a;a=a->next)
if(strncmp(a->name+m,s+m,n)==0) i++;
return i;
}


int Res::good(int b,int e,int m)
{
Atm *a;
int i;
i=0;
for(a=atm;a;a=a->next)
{
if(a->tatm->id<b&&a->tatm->id>e) continue;
if(a->good!=m) i++;
a->good=m;
}
return i;
}

float Res::bury(int n1,int n2)
{
Atm *a;
float e,e1;
e=0;
e1=0;
for(a=atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') 
  {
    if(a->tatm->bond[0]->id<n1||a->tatm->bond[0]->id>n2)continue;
  }
  else
  {
    if(a->tatm->id<n1||a->tatm->id>n2) continue;
  }
  e+=a->area;
  e1+=a->tatm->area;
}
if(e>e1) e=e1;
if(e1==0) return -1;
return 1-e/e1;
}

float Res::clash(int f)//internal clashes
{
float e,x;
Atm *a,*b,*c;
int i;
int n1,n2,m1,m2,near;
int k1,k2; 

if(f%10==0) //all against all interaction
{
 n1=0;n2=100;
 m1=0;m2=100;
}
else if(f%10==1) //mainchain against sidechain
{
 n1=0;n2=3;
 m1=4;m2=100;
}
else if(f%10==2)//all against sidechain
{
n1=0;n2=100;
m1=4;m2=100;
}
else if(f%10==3)//neighboring backbone against sidechain
{
n1=0;n2=100;
m1=4;m2=100;
k1=0;k2=4;
}
else return 0;

e=0;

//calculate torsional energy;

for(a=atm;a;a=a->next)
{
 if(a->tatm->name[1]=='H') break;
 if(a->tatm->id<m1||a->tatm->id>m2) continue;
 if(f/10>0)
 {
  if(a->tatm->rotm>=f/10) continue; 
 }
 e+=a->torsion(-1);
}

//calculate internal vdw clashes

for(a=atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') continue;
  if(a->tatm->id<m1||a->tatm->id>m2) continue;
  if(f/10>0)
  {
    if(a->tatm->rotm>=f/10+2) continue;
  }
  for(b=atm;b;b=b->next)
  {
    if(b->tatm->name[1]=='H') continue;
    if(b->tatm->id<n1||b->tatm->id>n2) continue;
    if(f/10>0)
    {
       if(b->tatm->rotm>=f/10+2) continue;
    }
    near=a->isnear(b,3);
    if(near>=0&&near<3) continue;
    if(near==-1) 
    {
      if(b->tatm->id>=m1&&b->tatm->id<=m2) x=0.5;
      else x=1.;
      e+=x*a->clash(b);
    }
    
    for(i=0;i<b->tatm->nbond;i++) 
    {
     c=b->bond[i];
     if(c==0) continue;
     if(c->tatm->name[1]!='H') continue;
     e+=a->clash(c);
    }
    
  }
}

//interacting with neighboring residues
if(f%10!=3) return e;
if(chn==0) return e;
Res *r1=(*chn)[id0-1];
if(r1==0) r1=this;
Res *r2;
for(r2=r1;r2;r2=r2->next) {
if(r2==this) continue;
if(r2->id0-id0>1) break;
for(a=atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') continue;
  if(a->tatm->id<m1||a->tatm->id>m2) continue;
  if(f/10>0)
  {
    if(a->tatm->rotm>=f/10+2) continue;
  }
  for(b=r2->atm;b;b=b->next)
  {
    if(b->tatm->name[1]=='H') continue;
    if(b->tatm->id<k1||b->tatm->id>k2) continue;
    if(f/10>0)
    {
       if(b->tatm->rotm>=f/10+2) continue;
    }
    near=a->isnear(b,3);
    if(near>=0&&near<3) continue;
    if(near==-1)
    {
      x=1.;
      e+=x*a->clash(b);
    }

    for(i=0;i<b->tatm->nbond;i++)
    {
     c=b->bond[i];
     if(c==0) continue;
     if(c->tatm->name[1]!='H') continue;
     e+=a->clash(c);
    }
  }
}
}

return e;
}

float Res::clashnoring(int f,int gg)//internal clashes
{
//gg remove torsion in the ring
float e,x;
Atm *a,*b,*c;
int i;
int n1,n2,m1,m2,near;
int k1,k2; 

if(f%10==0) //all against all interaction
{
 n1=0;n2=100;
 m1=0;m2=100;
}
else if(f%10==1) //mainchain against sidechain
{
 n1=0;n2=3;
 m1=4;m2=100;
}
else if(f%10==2)//all against sidechain
{
n1=0;n2=100;
m1=4;m2=100;
}

else if(f%10==3)//neighboring backbone against sidechain
{
n1=0;n2=100;
m1=4;m2=100;
k1=0;k2=4;
}
else if(f%10==4)//main against main
{
n1=0;n2=3;
m1=0;m2=3;
}
else if(f%10==5)//side against side 
{
n1=4;n2=100;
m1=4;m2=100;
}
else if(f%10==6)//side against side
{
n1=100;n2=100;
m1=100;m2=100;
}
else return 0;

e=0;

//calculate torsional energy;

for(a=atm;a;a=a->next)
{
 if(a->tatm->name[1]=='H') break;
 if(a->tatm->id<m1||a->tatm->id>m2) continue;
 if(f/10>0)
 {
  if(a->tatm->rotm>=f/10) continue; 
 }
 if(gg&&a->tatm->ring) continue;
 e+=a->torsion(-1);
}

if(TRES.selfvdw==0) return e;

//calculate internal vdw clashes

for(a=atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') continue;
  if(a->tatm->id<m1||a->tatm->id>m2) continue;
  if(f/10>0)
  {
    if(a->tatm->rotm>=f/10+2) continue;
  }
  for(b=atm;b;b=b->next)
  {
    if(b->tatm->name[1]=='H') continue;
    if(b->tatm->id<n1||b->tatm->id>n2) continue;
    if(a->tatm->ring&&b->tatm->ring&&a->res==b->res) continue;
    if(f/10>0)
    {
       if(b->tatm->rotm>=f/10+2) continue;
    }
    near=a->isnear(b,3);
    if(near>=0&&near<3) continue;
    if(near==-1) 
    {
      if(b->tatm->id>=m1&&b->tatm->id<=m2) x=0.5;
      else x=1.;
      e+=x*a->clash(b);
      if(TRES.logg==-1) cerr<<a->res->name<<a->res->id0<<" "<<a->name<<" "<<b->res->name<<b->res->id0<<" "<<b->name<<endl;
    }
    
    for(i=0;i<b->tatm->nbond;i++) 
    {
     c=b->bond[i];
     if(c==0) continue;
     if(c->tatm->name[1]!='H') continue;
     e+=a->clash(c);
     if(TRES.logg==-1) cerr<<a->res->name<<a->res->id0<<" "<<a->name<<" "<<c->res->name<<c->res->id0<<" "<<c->name<<endl;
    }
    
  }
}

//interacting with neighboring residues
if(f%10!=3) return e;
if(chn==0) return e;
Res *r1=(*chn)[id0-1];
if(r1==0) r1=this;
Res *r2;
for(r2=r1;r2;r2=r2->next) {
if(r2==this) continue;
if(r2->id0-id0>1) break;
for(a=atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') continue;
  if(a->tatm->id<m1||a->tatm->id>m2) continue;
  if(f/10>0)
  {
    if(a->tatm->rotm>=f/10+2) continue;
  }
  for(b=r2->atm;b;b=b->next)
  {
    if(b->tatm->name[1]=='H') continue;
    if(b->tatm->id<k1||b->tatm->id>k2) continue;
    if(a->tatm->ring&&b->tatm->ring&&a->res==b->res) continue;
    if(f/10>0)
    {
       if(b->tatm->rotm>=f/10+2) continue;
    }
    near=a->isnear(b,3);
    if(near>=0&&near<3) continue;
    if(near==-1)
    {
      x=1.;
      e+=x*a->clash(b);
    }

    for(i=0;i<b->tatm->nbond;i++)
    {
     c=b->bond[i];
     if(c==0) continue;
     if(c->tatm->name[1]!='H') continue;
     e+=a->clash(c);
    }
  }
}
}

return e;
}



float Res::clash(Atm *s,int f)
{
Atm *a;
int m1,m2;
float e;

if(s->res->id==id&&s->res->chn==chn)return 0;

if(f==0) //all against all interaction
{
 m1=0;m2=100;
}
else if(f==1) //mainchain against sidechain
{
 m1=4;m2=100;
}
else if(f==2)//all against sidechain
{
m1=4;m2=100;
}
else return 0;

e=0;
for(a=atm;a;a=a->next)
{
 if(a->tatm->name[1]=='H') continue;
 if(a->tatm->id<m1||a->tatm->id>m2) continue; 
 e+=a->clash(s);
}
return e;
}

float Res::rmsd(Res *s,int n1,int n2)
{
Atm *a1, *a2;
float x;
int n;
x=0;n=0;
if(s->name!=name) return -1;
for(a1=atm;a1;a1=a1->next)
{
 if(a1->tatm->name[1]=='H') continue;
 if(a1->tatm->id<n1||a1->tatm->id>n2) continue;
 a2=(*s)[a1->tatm->id];
 if(a2==0) continue; 
 x+=TRES.distsqr(a1->xyz,a2->xyz);
 n++; 
}
if(n==0) n=1;
x=sqrt(x/n);
return x+n*1000000;
}

float Res::directrmsd(Res *s,int n1,int n2)
{
Atm *a1, *a2;
float x;
int n;
x=0;n=0;
if(s->name!=name) return -1;
for(a1=atm;a1;a1=a1->next)
{
 if(a1->tatm->name[1]=='H') continue;
 if(a1->tatm->id<n1||a1->tatm->id>n2) continue;
 a2=(*s)[a1->tatm->id];
 if(a2==0) continue;
 x+=TRES.distsqr(a1->xyz,a2->xyz);
 n++;
}
x=sqrt(x/n);
return x;
}

float Res::directrmsdanyway(Res *s,int n1,int n2)
{
Atm *a1, *a2;
float x;
int n;
x=0;n=0;
//if(s->name!=name) return -1;
for(a1=atm;a1;a1=a1->next)
{
 if(a1->tatm->name[1]=='H') continue;
 if(a1->tatm->id<n1||a1->tatm->id>n2) continue;
 a2=(*s)[a1->tatm->id];
 if(a2==0) continue;
 x+=TRES.distsqr(a1->xyz,a2->xyz);
 n++;
}
x=sqrt(x/n);
return x;
}


void Res::copytemp(Res *s)
{
int i,j;

if(nemp>s->nemp) j=s->nemp;
else j=nemp;
nemp=j;
for(i=0;i<j;i++) temp[i]=s->temp[i]; 
}

void Res::copytemp(float *s,int flg)
{
  if(nemp==0) {cerr<<"no copytemp acted"<<endl;return;}
  if(flg==0) TRES.copy(temp,s,nemp);
  else       TRES.copy(s,temp,nemp);
}

void Res::giveid()
{
Res *r;
int i;
i=0;
for(r=more;r;r=r->more)
{
  r->id0=i++;//r->id=id;
}
}


Res *Res::ismore(int n)
{
Res *r;
int m;
m=0; 
for(r=more;r;r=r->more)
{
  if(m==n) return r;
  m++;
}
return 0;
}

Res *Res::findbestmore(float a,float b)
{
Res *r,*r0;
float dd=10000;
r0=0;
for(r=more;r;r=r->more)
{
	float d=fabs(r->atm->next->chi-a)+fabs(r->atm->next->next->chi-b);
	if(d<dd) {
		r0=r;
		dd=d;
	}
}
return r0;
}

 
int Res::ssbond(Res *r1, Res *r2,float s_s,float s_cb)
{
  Atm *a3,*a1,*a2;
  float d1;
  float d2;
  d2=0.5;
  a1=(*r1)[" SG "];a2=(*r2)[" SG "];
  d1=TRES.distance(a1->xyz,a2->xyz);
  if(d1>s_s+d2||d1<s_s-d2) return 0;
  a3=(*a2->res)[" CB "];
  d1=TRES.distance(a1->xyz,a3->xyz);
  if(d1>s_cb+d2||d1<s_cb-d2) return 0;
  a3=(*a1->res)[" CB "];
  d1=TRES.distance(a2->xyz,a3->xyz);
  if(d1>s_cb+d2||d1<s_cb-d2) return 0;
  return 1;
}

int Res::many()
{
  Res *r;
  int k;
  k=0;
  for(r=more;r;r=r->more)k++;
  return k;
}

Res::Res(Tres *tr)
{
Tatm *a1;//,*a2;
Atm  *s1;
int i;
initial();

name=tr->name;
tres=tr;
for(a1=tr->tatm;a1;a1=a1->next)
{
  if(atm==0)
  {
      atm=new Atm;
      s1=atm;
  }
  else
  {
      s1->next=new Atm;
      s1=s1->next;
  }
  s1->res=this;
  strcpy(s1->name,a1->name);
  s1->tatm=a1;s1->id=a1->id;
  for(i=0;i<3;i++) s1->xyz[i]=a1->xyz[i];
}

}

void Res::bondleng()
{
  Atm *a;
  for(a=atm;a;a=a->next)
  {
    a->bondleng();
  }
}

void Res::transfer(Res *s,int n)
{
  Res *r,*r0;
  int i;
  r0=s;r=this;i=0;
  for(; ;)
  {
    if(r0->name!=r->name)
    {cerr<<"warning: not the same residue copying!"<<endl;}   
    r->transfer(r0);
    i++;
    if(i==n) break;
    r0=r0->next;
    r=r->next;
    if(r0==0||r==0) break;
  }
}

float *Res::center(int fa,int fb)
{
int i,m;
float *xyz;
Atm *a;

xyz=new float[3];
for(i=0;i<3;i++) xyz[i]=0;

m=0;
for(a=atm;a;a=a->next)
{
  if(fa==0&&a->tatm->id!=1) continue;
  if(fa==1&&a->tatm->id>=4) continue;
  if(fa==2&&a->tatm->name[1]=='H') continue;
  for(i=0;i<3;i++) xyz[i]+=a->xyz[i];
  m++;
}

for(i=0;i<3;i++) xyz[i]=xyz[i]/m; 

if(fb) return xyz;

for(a=atm;a;a=a->next)
{
  for(i=0;i<3;i++)a->xyz[i]-=xyz[i];
}
return xyz;
}

void Res::configure() {
  
  Tatm *tatm_temp;
  Atm *atm_temp,*atm0[30];
  int i,j;

  //re-order
  for(i=0;i<30;i++)atm0[i]=0;
  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
  {
    atm0[atm_temp->tatm->id]=atm_temp;
  }
  j=0;
  for(i=0;i<tres->number;i++)
  {
   if(atm0[i]==0) continue;
   atm0[j++]=atm0[i];
  }
  atm=atm0[0];
  for(i=0;i<j-1;i++)
  {
    atm0[i]->next=atm0[i+1];
  }
  atm0[j-1]->next=0;

  //atom id and number and nhydr
  number=0;nhydr=0;
  int  n=0,nn=0;
  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
  {
    number++;
    atm_temp->id0=n++;
    //n++;
    if(atm_temp->name[1]=='H') nhydr++;
  }
 
  //set residues of atoms
  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next) atm_temp->res=this;

  //atom id with tres complete supposed
  nn=0;
  for(tatm_temp=tres->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  {
    atm_temp=(*this)[tatm_temp->id];
    if(atm_temp)atm_temp->id=nn;
    nn++;
  }
 
  
// assign topology to the pdb read

  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
  {
     for(i=0;i<4;i++)
     {
       tatm_temp=atm_temp->tatm->bond[i];
       if(!tatm_temp) continue;
       atm_temp->bond[i]=(*this)[tatm_temp->id];
     }
     if(atm_temp->tatm->id==0)  atm_temp->bond[0]=0;
     if(atm_temp->tatm->id==2)  atm_temp->bond[1]=0;
  }
}

void Res::configureaddhatom() {
  
  Tatm *tatm_temp;
  Atm *atm_temp,*atm0[30];
  int i,j;

  //re-order
  for(i=0;i<30;i++)atm0[i]=0;
  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
  {
    atm0[atm_temp->tatm->id]=atm_temp;
  }
  j=0;
  for(i=0;i<tres->number;i++)
  {
   if(atm0[i]==0) continue;
   atm0[j++]=atm0[i];
  }
  atm=atm0[0];
  for(i=0;i<j-1;i++)
  {
    atm0[i]->next=atm0[i+1];
  }
  atm0[j-1]->next=0;

  //atom id and number and nhydr
  number=0;nhydr=0;
  int  n=0,nn=0;
  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
  {
    number++;
    atm_temp->id0=n++;
    //n++;
    if(atm_temp->name[1]=='H') nhydr++;
  }
 
  //set residues of atoms
  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next) atm_temp->res=this;

  //atom id with tres complete supposed
  nn=0;
  for(tatm_temp=tres->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  {
    atm_temp=(*this)[tatm_temp->id];
    if(atm_temp)atm_temp->id=nn;
    nn++;
  }
 
  
// assign topology to the pdb read

  for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next)
  {
     for(i=0;i<4;i++)
     {
       if(atm_temp->tatm->id==0&&i==0) continue;
       if(atm_temp->tatm->id==2&&i==1) continue;
       tatm_temp=atm_temp->tatm->bond[i];
       if(!tatm_temp) continue;
       atm_temp->bond[i]=(*this)[tatm_temp->id];
     }
     //if(atm_temp->tatm->id==0)  atm_temp->bond[0]=0;
     //if(atm_temp->tatm->id==2)  atm_temp->bond[1]=0;
  }
}
void Res::attachAtm(Atm *a) {

   Atm *b;
   for(b=atm;b->next;b=b->next);
   if(b==0) atm=a;
   else b->next=a;
   a->next=0;
} 

void Res::addsidechain() {

   Atm *a,*b;
   Tatm *t;
   int cb=1;

   if((*this)[4]==0) cb=0;

   int n=0;
   for(t=tres->tatm;t;t=t->next) {
	a=(*this)[t->name];
        if(a) continue;
	if(t->id<4) {
		cerr<<"residue:"<<name<<id+chn->start<<" missing mainchains"<<endl;
		continue;
        }
	if(t->name[1]=='H'&&t->bond[0]->id<4) continue;
	Atm *a0=new Atm(t);
	a0->tatm=t;
	attachAtm(a0);
	n++;
   }

   if(n) configure();
   else return;

   Rotate rot; 

   if(name!='P') {
	if(cb==0) rot.hook(this);
	else     rot.hook(this,4);
   }
   else {   
        Res *res0=0; 
	float dd=1000;
        for(Res *res_temp=TRES.rotamer->next->chn->isres(name,0)->more;res_temp;res_temp=res_temp->more) {
		res0=res_temp;
		rot.hook(res_temp,this);	
		a=(*this)[" CD "];
		b=(*this)[" N  "];
 		float dt=TRES.distance(a,b);
		if(dd>fabs(dt-1.5)) {
			dd=fabs(dt-1.5);
			res0=res_temp;
		}
        }
	if(res0) {
		 rot.hook(res0,this);
		 a=(*this)[" CD "];
                 b=(*this)[" N  "];
		 float dt=TRES.distance(a,b);
		 if(TRES.logg) cerr<<"proline residue:"<<id<<" checked with CD to N distance:"<<dt<<endl;
	}
   }
}
void Res::addmissingatms() {

   Atm *a;
   Tatm *t;
  
   int n=0;
   for(t=tres->tatm;t;t=t->next) {
	a=(*this)[t->name];
        if(a) continue;
	cerr<<"add missing atom: "<<name<<oid<<" "<<t->name<<endl;
	Atm *a0=new Atm(t);
	a0->tatm=t;
	attachAtm(a0);
	n++;
   }

   if(n) configure();
   return;  
}
void Res::hooksidechain() {

   Atm *a,*b;
   //Tatm *t;
   int cb=1;

   /*
   if((*this)[4]==0) cb=0;
	
   int n=0;
   for(t=tres->tatm;t;t=t->next) {
	a=(*this)[t->name];
        if(a) continue;
	if(t->id<4) {
		cerr<<"residue:"<<name<<id+chn->start<<" missing mainchains"<<endl;
		continue;
        }
	if(t->name[1]=='H'&&t->bond[0]->id<4) continue;
	attachAtm(new Atm(t));
	n++;
   }

   if(n) configure();
   else return;
   */

   a=isatm(" CA ");
   b=isatm(" CB ");
   
   if(TRES.distance(a,b)>2) cb=0;	

   a=isatm(" N  ");
   if(TRES.distance(a,b)>2) cb=0;	

   a=isatm(" C  ");
   if(TRES.distance(a,b)>2) cb=0;	

   Rotate rot; 

   if(name!='P') {
	if(cb==0) rot.hook(this);
	else     rot.hook(this,4);
   }
   else {   
        Res *res0=0; 
	float dd=1000;
        for(Res *res_temp=TRES.rotamer->next->chn->isres(name,0)->more;res_temp;res_temp=res_temp->more) {
		res0=res_temp;
		rot.hook(res_temp,this);	
		a=(*this)[" CD "];
		b=(*this)[" N  "];
 		float dt=TRES.distance(a,b);
		if(dd>fabs(dt-1.5)) {
			dd=fabs(dt-1.5);
			res0=res_temp;
		}
        }
	if(res0) {
		 rot.hook(res0,this);
		 a=(*this)[" CD "];
                 b=(*this)[" N  "];
		 float dt=TRES.distance(a,b);
		 cerr<<"proline residue:"<<id<<" checked with CD to N distance:"<<dt<<endl;
	}
   }
}

void Res::addhatoms(int f) {
   addbackbonehatoms(f);
   addsidechainhatoms(f);
}


void Res::addbackbonehatoms(int f) {

   Atm *a;
   Tatm *t;

   for(t=tres->tatm;t;t=t->next) {
        if(t->name[1]!='H') continue;
	if(t->bond[0]==0||t->bond[0]->id>=4) continue;
        a=(*this)[t->name];
	if(a) {
		if(f) continue;
	}
	else {
		a=new Atm(t);
		attachAtm(a);
		configureaddhatom();
	}
        int ff=setbackbonehatomxyz(a);
	if(ff==0) {
		takeoffatm(a);
		configureaddhatom();
	}
   }
}

void Res::addhnatoms(int f) {

	//configure ?
	Tatm *t=tres->isatm(" HN ");
	if(t==0) return; //not configured
	//exist?
	Atm *a=isatm(" HN ");
	
	if(a) {
                if(f) return;
        }
        else {
                a=new Atm(t);
                attachAtm(a);
                configureaddhatom();
        }

	//int ff=sethnatomxyz(a);
	int ff=setbackbonehatomxyz(a);

	if(ff==0) {
                takeoffatm(a);
                configureaddhatom();
        }
}

int Res::isbackbonecomplete() {
	Atm *a;
	int i;
	for(i=0;i<4;i++) {
		a=isatmid(i);
		if(a==0) return 0;
	}
	return 1;
}

int Res::issidechaincomplete() {
        Atm *a;
	Tatm *ta;
	for(ta=tres->tatm;ta;ta=ta->next) {
		if(ta->name[1]=='H'&&ta->bond[0]->id<=3) continue;
		if(ta->name[1]!='H'&&ta->id<=3) continue;	
		a=isatmid(ta->id);
		if(a==0) return 0;
	}
        return 1;
}



void Res::takeoffatm(Atm *a) {

   Atm *b,*b0;
   b0=0;
   for(b=atm;b;b=b->next) {
	if(b==a) {
		if(b0==0) {
			atm=b->next;
			b->next=0;
			delete b;
			return;
		}	
		else {
			b0->next=b->next;
			b->next=0;
			delete b;
			return;
		}
	}
	b0=b;
   }
   return;
} 

void Res::addsidechainhatoms(int f) {

   Atm *a;
   Tatm *t;

   for(t=tres->tatm;t;t=t->next) {
        if(t->name[1]!='H') continue;
        if(t->bond[0]==0||t->bond[0]->id<4) continue;
        a=(*this)[t->name];
        if(a) {
                if(f) continue;
        }
        else {
                a=new Atm(t);
                attachAtm(a);
                configureaddhatom();
        }
        int ff=setsidechainhatomxyz(a);
	if(ff==0) {
		takeoffatm(a);
		configureaddhatom();
	}
   }
}


int Res::setbackbonehatomxyz(Atm *aa) {

  AtmGeom g;
  //Tatm *t2;
  Atm *a;
  float d;
  //int yes=0;

  //t=aa->tatm;
  if(strcmp(aa->name," HN ")==0) {
	
	float *n,*ca,*o,*c;
	float *th,*tn,*to;

	//h
	//h=aa->xyz;
	th=aa->tatm->xyz; 

	//n
	a=(*this)[0];
	if(a==0) return 0;
	n=a->xyz;
	tn=a->tatm->xyz;

	//ca
	
	a=(*this)[1];
        if(a==0) return 0;
        ca=a->xyz;
        //tca=a->tatm->xyz;
	

	//o
	o=temp+6;
	to=tres->head+6;

	//c
	c=temp+3;
        //tc=tres->head+3;

	//add n
	d=TRES.distance(th,tn);
	g.addbounds(n,d,0.0001);
	float dd=d;

	//add ca
	float ange=TRES.angle(c,n,ca);	
	ange=(360-ange)/2;

	//float ange=TRES.angle(th,tn,tca);
	
	d=dd*dd+TRES.distsqr(n,ca)-2*dd*TRES.distance(n,ca)*cos(ange*3.14/180);
	d=sqrt(d);
	g.addbounds(ca,d,0.001);
		

	//add c
	//ange=TRES.angle(th,tn,to)-TRES.angle(c,n,o);
	//ange=TRES.angle(th,tn,tc);
        d=dd*dd+TRES.distsqr(n,c)-2*dd*TRES.distance(n,c)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(c,d,0.001);
	//float angee=ange;
	float de=d;

	//add o
	
	//ange=TRES.angle(th,tn,tc)+TRES.angle(c,n,o);
	//ange=TRES.angle(th,tn,to);//+angee;
        //d=dd*dd+TRES.distsqr(n,o)-2*dd*TRES.distance(n,o)*cos(ange*3.14/180);
        //d=sqrt(d);
        //g.addbounds(o,d,0.001);
	//g.findcoo();
  	g.findmincoo(2,0.0005);
	if(TRES.logg==-51) cerr<<"nh.."<<TRES.angle(th,tn,to)<<"  "<<TRES.angle(g.coo,n,o)<<endl;
	if(TRES.logg==-51) cerr<<"nhd..."<<de<<" "<<TRES.distance(g.coo,c)<<" "<<g.evaluate()<<"  "<<g.evaluaterange()<<endl;
  }
  else if(strcmp(aa->name," HA ")==0) {


	float *n,*ca,*cb,*c;
	float *th,*tn,*tca,*tcb,*tc;

	//h
	//h=aa->xyz;
	th=aa->tatm->xyz;

	//n
	a=(*this)[0];
	if(a==0) return 0;
	n=a->xyz;
	tn=a->tatm->xyz;

	//ca
	a=(*this)[1];
        if(a==0) return 0;
        ca=a->xyz;
        tca=a->tatm->xyz;
	
	//cb
	a=(*this)[4];
	if(a==0) return 0;
	cb=a->xyz;
	tcb=a->tatm->xyz;
 
	//c
	a=(*this)[2];
	if(a==0) return 0;
	c=a->xyz;
	tc=a->tatm->xyz;
 

	//add ca
	d=TRES.distance(th,tca);
	g.addbounds(ca,d,0.05,10);
	float dd=d;

	//add n
	float ange=TRES.angle(th,tca,tn);
	d=dd*dd+TRES.distsqr(n,ca)-2*dd*TRES.distance(n,ca)*cos(ange*3.14/180);
	d=sqrt(d);
	g.addbounds(n,d,0.3);

	//add c
	ange=TRES.angle(th,tca,tc);
        d=dd*dd+TRES.distsqr(ca,c)-2*dd*TRES.distance(ca,c)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(c,d,0.3);

	//add cb
	ange=TRES.angle(th,tca,tcb);
        d=dd*dd+TRES.distsqr(ca,cb)-2*dd*TRES.distance(ca,cb)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(cb,d,0.3);
 
  	g.findmincoo(0.06);
  }
  else if(strcmp(aa->name,"1HA ")==0) {
	
	
	float *n,*ca,*c;
	float *th,*tn,*tca,*tc;

	//h
	//h=aa->xyz;
	th=aa->tatm->xyz;

	//n
	a=(*this)[0];
	if(a==0) return 0;
	n=a->xyz;
	tn=a->tatm->xyz;

	//ca
	a=(*this)[1];
        if(a==0) return 0;
        ca=a->xyz;
        tca=a->tatm->xyz;
	
	/*
	//cb
	a=(*this)[4];
	if(a==0) return 0;
	cb=a->xyz;
	tcb=a->tatm->xyz;
 	*/

	//c
	a=(*this)[2];
	if(a==0) return 0;
	c=a->xyz;
	tc=a->tatm->xyz;
 

	//add ca
	d=TRES.distance(th,tca);
	g.addbounds(ca,d,0.05,20);
	float dd=d;

	//add n
	float ange=TRES.angle(th,tca,tn);
	d=dd*dd+TRES.distsqr(n,ca)-2*dd*TRES.distance(n,ca)*cos(ange*3.14/180);
	d=sqrt(d);
	g.addbounds(n,d,0.3);

	//add c
	ange=TRES.angle(th,tca,tc);
        d=dd*dd+TRES.distsqr(ca,c)-2*dd*TRES.distance(ca,c)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(c,d,0.3);

	/*//add cb
	ange=TRES.angle(th,tca,tcb);
        d=dd*dd+TRES.distsqr(ca,cb)-2*dd*TRES.distance(ca,cb)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(cb,d,0.3);
 	*/
  	g.findmincoo(0.06);
  }
  else if(strcmp(aa->name,"2HA ")==0) {

	float *n,*ca,*cb,*c;
	float *th,*tn,*tca,*tcb,*tc;

	//h
	//h=aa->xyz;
	th=aa->tatm->xyz;

	//n
	a=(*this)[0];
	if(a==0) return 0;
	n=a->xyz;
	tn=a->tatm->xyz;

	//ca
	a=(*this)[1];
        if(a==0) return 0;
        ca=a->xyz;
        tca=a->tatm->xyz;
	
	//cb
	a=(*this)["1HA "];
	if(a==0) return 0;
	cb=a->xyz;
	tcb=a->tatm->xyz;
 
	//c
	a=(*this)[2];
	if(a==0) return 0;
	c=a->xyz;
	tc=a->tatm->xyz;
 

	//add ca
	d=TRES.distance(th,tca);
	g.addbounds(ca,d,0.05,20);
	float dd=d;

	//add n
	float ange=TRES.angle(th,tca,tn);
	d=dd*dd+TRES.distsqr(n,ca)-2*dd*TRES.distance(n,ca)*cos(ange*3.14/180);
	d=sqrt(d);
	g.addbounds(n,d,0.3);

	//add c
	ange=TRES.angle(th,tca,tc);
        d=dd*dd+TRES.distsqr(ca,c)-2*dd*TRES.distance(ca,c)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(c,d,0.3);

	//add cb
	ange=TRES.angle(th,tca,tcb);
        d=dd*dd+TRES.distsqr(ca,cb)-2*dd*TRES.distance(ca,cb)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(cb,d,0.5);
  	g.findmincoo(0.06);
  }
  else {
	cerr<<"not backbone h atoms"<<endl;
	return 0;
  }
  if(TRES.logg) cerr<<aa->name<<" "<<tres->name<<id0<<" "<<g.num<<endl;
  d=g.evaluatemin();
  
  if(d>0.4&&addhanyway==0) {
	cerr<<"the atom position of "<<aa->name<<" is not accurate standard! "<< d<<endl;
	return 0;
  }
  if(d>0.4) {
	cerr<<"the atom position of "<<name<<oid<<"--"<<aa->name<<" is not accurate standard! "<< d<<endl;
  }
  for(int i=0;i<3;i++) aa->xyz[i]=g.coo[i];
  return 1;
}	

int Res::sethnatomxyz(Atm *aa) {

  AtmGeom g;
  //Tatm *t2;
  Atm *a;
  float d;
  //int yes=0;

  //t=aa->tatm;
  if(strcmp(aa->name," HN ")==0) {
	
	float *n,*ca,*o,*c;
	float *th,*tn,*to;

	//h
	//h=aa->xyz;
	th=aa->tatm->xyz; 

	//n
	a=(*this)[0];
	if(a==0) return 0;
	n=a->xyz;
	tn=a->tatm->xyz;

	//ca	
	a=(*this)[1];
        if(a==0) return 0;
        ca=a->xyz;
        //tca=a->tatm->xyz;
	
	Res *ro=chn->isres0(id0-1);
	if(ro==0||id-ro->id!=1) return 0;
	//o
	Atm *ao=ro->isatmid(3);
	if(ao==0) return 0;
	o=ao->xyz;//o=temp+6;
	to=tres->head+6;

	//c
	Atm *ac=ro->isatmid(2);
	if(ac==0) return 0;
	c=ac->xyz;//c=temp+3;
	if(TRES.distance(c,n)>2) return 0;
        //tc=tres->head+3;

	//add n
	d=TRES.distance(th,tn);
	g.addbounds(n,d,0.0001);
	float dd=d;

	//add ca
	float ange=TRES.angle(c,n,ca);	
	ange=(360-ange)/2;

	//float ange=TRES.angle(th,tn,tca);
	
	d=dd*dd+TRES.distsqr(n,ca)-2*dd*TRES.distance(n,ca)*cos(ange*3.14/180);
	d=sqrt(d);
	g.addbounds(ca,d,0.001);
		

	//add c
	//ange=TRES.angle(th,tn,to)-TRES.angle(c,n,o);
	//ange=TRES.angle(th,tn,tc);
        d=dd*dd+TRES.distsqr(n,c)-2*dd*TRES.distance(n,c)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(c,d,0.001);
	//float angee=ange;
	float de=d;

	//add o
	
	//ange=TRES.angle(th,tn,tc)+TRES.angle(c,n,o);
	//ange=TRES.angle(th,tn,to);//+angee;
        //d=dd*dd+TRES.distsqr(n,o)-2*dd*TRES.distance(n,o)*cos(ange*3.14/180);
        //d=sqrt(d);
        //g.addbounds(o,d,0.001);
	//g.findcoo();
  	g.findmincoo(2,0.0005);
	if(TRES.logg) cerr<<"nh.."<<TRES.angle(th,tn,to)<<"  "<<TRES.angle(g.coo,n,o)<<endl;
	if(TRES.logg) cerr<<"nhd..."<<de<<" "<<TRES.distance(g.coo,c)<<" "<<g.evaluate()<<"  "<<g.evaluaterange()<<endl;
  }
 
  else {
	cerr<<"not backbone h atoms"<<endl;
	return 0;
  }
  if(TRES.logg>3) cerr<<aa->name<<" "<<tres->name<<id0<<" "<<g.num<<endl;
  d=g.evaluatemin();
  
  if(d>0.4&&addhanyway==0) {
	cerr<<"the atom position of "<<aa->name<<" is not accurate standard! "<< d<<endl;
	return 0;
  }
  if(d>0.4) {
        cerr<<"the atom position of "<<name<<oid<<"--"<<aa->name<<" is not accurate standard! "<< d<<endl;
  }

  for(int i=0;i<3;i++) aa->xyz[i]=g.coo[i];
  return 1;
}	


int Res::setsidechainhatomxyz(Atm *a) {

   //find the parent;

   if(a==0) {
	//cerr<<"the atom:"<<a->name<<" does not bonded"<<endl;
	return 0;
   }

   if(a->tatm==0) {
	cerr<<"the atom:"<<a->res->name<<a->res->oid<<" is not complete standard residue"<<endl;
	cerr<<"protons not added:"<<a->name<<endl;
	cerr<<endl;
        return 0;
   }

   Tatm *t0=a->tatm->bond[0];
   Atm  *s0=a->bond[0];

   if(t0==0||s0==0) {
	cerr<<"the atom:"<<a->res->name<<a->res->oid<<" is not complete standard residue"<<endl;
	cerr<<"protons not added:"<<a->name<<endl;
	cerr<<endl;
        return 0;
   }

   AtmGeom g;


   float d=TRES.distance(t0->xyz,a->tatm->xyz);
   float dd=d;
   int i;
   //int n=0; 

   if(t0->ring==1&&strchr("YFW",tres->name)) {
	
	Tatm *tf=t0->tres->findfarringatm(t0);
	if(tf==0) {
		cerr<<"the atom:"<<a->name<<" is not complete standard residue"<<endl;
		cerr<<"protons not added:"<<a->name<<endl;
		cerr<<endl;
        	return 0;
	}
 	Atm *sf=s0->res->isatm(tf->name);
	if(sf==0) {
		cerr<<"the atom:"<<a->name<<" is not complete standard residue"<<endl;
		cerr<<"protons not added:"<<a->name<<endl;
		cerr<<endl;
        	return 0;
	}
	//cerr<<name<<" "<<s0->name<<" "<<sf->name<<endl;
	float de=TRES.distance(sf->xyz,s0->xyz);
	if(de<0.1) {
		cerr<<"the atom:"<<a->name<<" is not complete standard residue"<<endl;
		cerr<<"protons not added:"<<a->name<<endl;
		cerr<<endl;
        	return 0;
	}
	int ie;
	for(ie=0;ie<3;ie++) a->xyz[ie]=s0->xyz[ie]+(dd)/de*(s0->xyz[ie]-sf->xyz[ie]);
	return 1; 
   }
   else if(t0->ring&&strchr("HYFW",tres->name)) {
	
	Tatm **tf=t0->tres->findtwofarringatm(t0);
	if(tf==0||tf[0]==0||tf[1]==0) {
		cerr<<"the atom:"<<a->name<<" is not complete standard residue"<<endl;
		cerr<<"protons not added:"<<a->name<<endl;
		cerr<<endl;
        	return 0;
	}
 	Atm *sf[2];
	sf[0]=s0->res->isatm(tf[0]->name);
	sf[1]=s0->res->isatm(tf[1]->name);
	if(tf) delete [] tf;tf=0;
	if(sf[0]==0||sf[1]==0) {
		cerr<<"the atom:"<<a->name<<" is not complete standard residue"<<endl;
		cerr<<"protons not added:"<<a->name<<endl;
		cerr<<endl;
        	return 0;
	}
	//cerr<<name<<" "<<s0->name<<" "<<sf->name<<endl;
	int ie;
	float sxyz[3];
	for(ie=0;ie<3;ie++) sxyz[ie]=(sf[0]->xyz[ie]+sf[1]->xyz[ie])/2;
	float de=TRES.distance(sxyz,s0->xyz);
	if(de<0.1) {
		cerr<<"the atom:"<<a->name<<" is not complete standard residue"<<endl;
		cerr<<"protons not added:"<<a->name<<endl;
		cerr<<endl;
        	return 0;
	}
	
	for(ie=0;ie<3;ie++) a->xyz[ie]=s0->xyz[ie]+(dd)/de*(s0->xyz[ie]-sxyz[ie]);
	return 1; 
   }
   else if(t0->ring&&strchr("YFWH",tres->name)) {
	 
	if(s0->bond[0]==0||s0->bond[1]==0) return 0;
   	g.addbounds(s0->xyz,dd,0.001);
	float ange=TRES.angle(s0->bond[0]->xyz,s0->xyz,s0->bond[1]->xyz);
	ange=(360-ange)/2;

	//bond[0]
	d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.1415926/180);
	d=sqrt(d);
	g.addbounds(s0->bond[0]->xyz,d,0.0001);
	
	//bond[1]
	d=dd*dd+TRES.distsqr(s0->bond[1]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[1]->xyz,s0->xyz)*cos(ange*3.1415926/180);
        d=sqrt(d);
        g.addbounds(s0->bond[1]->xyz,d,0.0001);
	g.findmincoo(5,0.0005);
	 
   }
   else if(tres->name=='R'&&strstr(a->tatm->name," HE ")) {

	g.addbounds(s0->xyz,dd,0.001);
	float angee=TRES.angle(s0->bond[0]->xyz,s0->xyz,s0->bond[1]->xyz);
	float ange=(360-angee)/2;
        d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(s0->bond[0]->xyz,d,0.001);

	ange=120;
        d=dd*dd+TRES.distsqr(s0->bond[1]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[1]->xyz,s0->xyz)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(s0->bond[1]->xyz,d,0.001);

	g.findmincoo(2,0.0005);

   }
   else if(tres->name=='R'&&strstr(a->tatm->name,"HH")) {

	g.addbounds(s0->xyz,dd,0.001);
	if(strstr(a->tatm->name,"1HH1")) {
		float ange;
		ange=120;
		d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.14/180);
        	d=sqrt(d);
        	g.addbounds(s0->bond[0]->xyz,d,0.001);

		Atm *nh2;
		nh2=isatm(" NH2");
		if(nh2==0) return 0;
		ange=120+TRES.angle(s0->bond[0]->xyz,s0->xyz,nh2->xyz);
		d=dd*dd+TRES.distsqr(nh2->xyz,s0->xyz)-2*dd*TRES.distance(nh2->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(nh2->xyz,d,0.001);
	}
	else if(strstr(a->tatm->name,"2HH1")) {
		float ange;
                ange=120;
                d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(s0->bond[0]->xyz,d,0.001);

                Atm *nh2;
                nh2=isatm(" NH2");
                if(nh2==0) return 0;
                ange=120-TRES.angle(s0->bond[0]->xyz,s0->xyz,nh2->xyz);
                d=dd*dd+TRES.distsqr(nh2->xyz,s0->xyz)-2*dd*TRES.distance(nh2->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(nh2->xyz,d,0.001);
        }
	else if(strstr(a->tatm->name,"1HH2")) {
		float ange;
		ange=120;
                d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(s0->bond[0]->xyz,d,0.001);
                Atm *nh1;
                nh1=isatm(" NH1");
                if(nh1==0) return 0;
                ange=120+TRES.angle(s0->bond[0]->xyz,s0->xyz,nh1->xyz);
                d=dd*dd+TRES.distsqr(nh1->xyz,s0->xyz)-2*dd*TRES.distance(nh1->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(nh1->xyz,d,0.001);
	}
	else if(strstr(a->tatm->name,"2HH2")) {
		float ange;
                ange=120;
                d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(s0->bond[0]->xyz,d,0.001);

                Atm *nh1;
                nh1=isatm(" NH1");
                if(nh1==0) return 0;
                ange=120-TRES.angle(s0->bond[0]->xyz,s0->xyz,nh1->xyz);
                d=dd*dd+TRES.distsqr(nh1->xyz,s0->xyz)-2*dd*TRES.distance(nh1->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(nh1->xyz,d,0.001);
        }
	g.findmincoo(2,0.0005);
   }
   else if(s0->name[1]=='N'&&s0->bond[0]&&s0->bond[1]&&s0->bond[0]->name[1]!='H'&&s0->bond[1]->name[1]!='H') {

	g.addbounds(s0->xyz,dd,0.001);
        float angee=TRES.angle(s0->bond[0]->xyz,s0->xyz,s0->bond[1]->xyz);
        float ange=(360-angee)/2;
        d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(s0->bond[0]->xyz,d,0.001);

        ange=120;
        d=dd*dd+TRES.distsqr(s0->bond[1]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[1]->xyz,s0->xyz)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(s0->bond[1]->xyz,d,0.001);
        g.findmincoo(2,0.0005);
   }
   else if(s0->name[1]=='N'&&strchr("NQ",name)) {
	g.addbounds(s0->xyz,dd,0.001);
	
	Atm *o;
	if(name=='N')  o=isatm(" OD1");
	else	       o=isatm(" OE1");
	if(o==0) return 0;
	
	float ange=120;
	if(strstr(a->name,"1H")) {
		d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.14/180);
        	d=sqrt(d);
        	g.addbounds(s0->bond[0]->xyz,d,0.001);

		float angee=TRES.angle(s0->bond[0]->xyz,s0->xyz,o->xyz);
		ange=120-angee;
		d=dd*dd+TRES.distsqr(o->xyz,s0->xyz)-2*dd*TRES.distance(o->xyz,s0->xyz)*cos(ange*3.14/180);
		d=sqrt(d);
                g.addbounds(o->xyz,d,0.001);
	}
	else {

		d=dd*dd+TRES.distsqr(s0->bond[0]->xyz,s0->xyz)-2*dd*TRES.distance(s0->bond[0]->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(s0->bond[0]->xyz,d,0.001);

                float angee=TRES.angle(s0->bond[0]->xyz,s0->xyz,o->xyz);
                ange=120+angee;
                d=dd*dd+TRES.distsqr(o->xyz,s0->xyz)-2*dd*TRES.distance(o->xyz,s0->xyz)*cos(ange*3.14/180);
                d=sqrt(d);
                g.addbounds(o->xyz,d,0.001);
	}
	g.findmincoo(2,0.0005);
   }
   else if(t0->getnumberhbonds()==3&&s0->getnumberhbonds()==1&&strncmp(a->tatm->name,"1H",2)==0) {
		
      g.addbounds(s0->xyz,dd,0.001,100);
      float angee=0;
      for(i=0;i<1;i++) {  //find the parent's child
           Tatm *t=t0->bond[i];
           if(t==0) continue;
           if(t==a->tatm) continue;
           Atm *s=s0->bond[i];
           if(s==0) continue;
           float ange=TRES.angle(a->tatm->xyz,t0->xyz,t->xyz);
           d=dd*dd+TRES.distsqr(s->xyz,s0->xyz)-2*dd*TRES.distance(s->xyz,s0->xyz)*cos(ange*3.14/180);
           d=sqrt(d);
	   if(t==t0->bond[0]) angee=ange;
           //d=TRES.distance(a->tatm->xyz,t->xyz);
           if(t->name[1]!='H') g.addbounds(s->xyz,d,0.001,1);
           else  g.addbounds(s->xyz,d,0.001,1);
      }
      if(name=='L'||name=='I'||name=='K'||name=='T'||name=='V'||name=='A'||name=='M') {
	   
	   Atm *hg;
	   for(i=0;i<4;i++) {
		hg=0;
		if(s0->bond[0]==0) continue;
		hg=s0->bond[0]->bond[i];
		if(hg==0) continue;
		if(hg->tatm->name[1]=='H') break;
		hg=0;
	   }
	   if(name=='M') hg=isatm(" CG "); 
	   if(hg) {
		float ange=TRES.angle(hg->xyz,s0->xyz,s0->bond[0]->xyz)+angee; //109;
		d=dd*dd+TRES.distsqr(hg->xyz,s0->xyz)-2*dd*TRES.distance(hg->xyz,s0->xyz)*cos(ange*3.14/180);
		d=sqrt(d);
		g.addbounds(hg->xyz,d,0.001,1);
	   }
      }
      
      g.findmincoo(2,0.0005);

   }
   else {
      g.addbounds(s0->xyz,dd,0.001,100);
      for(i=0;i<4;i++) {  //find the parent's child
	   Tatm *t=t0->bond[i];
	   if(t==0) continue;
	   if(t==a->tatm) continue;
	   Atm *s=s0->bond[i];
	   if(s==0) continue;
	   float ange=TRES.angle(a->tatm->xyz,t0->xyz,t->xyz);
	   d=dd*dd+TRES.distsqr(s->xyz,s0->xyz)-2*dd*TRES.distance(s->xyz,s0->xyz)*cos(ange*3.14/180);	
	   d=sqrt(d);
	   //d=TRES.distance(a->tatm->xyz,t->xyz);
	   if(t->name[1]!='H') g.addbounds(s->xyz,d,0.3,1);
	   else  g.addbounds(s->xyz,d,0.5,1);
      }
      g.findmincoo(5,0.06);
   }
   //g.findcoo();
   d=g.evaluate();
   if(d>0.4&&addhanyway==0){
	cerr<<"the atom position of "<<a->name<<" is not accurate standard!"<<endl;
	return 0;
   }
   if(d>0.4) {
        cerr<<"the atom position of "<<name<<oid<<"--"<<a->name<<" is not accurate standard! "<< d<<endl;
  }

   for(i=0;i<3;i++) a->xyz[i]=g.coo[i];
   return 1;
}

void Res::mutateResidue(char c){
	
   if(TRES[c]==0) {
	cerr<<"residue: "<<c<<" is not standard!"<<endl;
   }

   if(c==name) return; 

   for(Res *r=more;r;r=r->more) {
	r->next=0;
   }

   if(more) delete more; more=0;
   
   Atm *aa[30];
   Atm *bb[30];

   int i;
   for(i=0;i<30;i++) {aa[i]=0;bb[i]=0;}

   int n=0;int j=0; 

   tres=TRES[c];
   name=tres->name;

   Atm *a;
   for(a=atm;a;a=a->next) {
	Tatm *t=tres->istatm(a->name);
	if(a->tatm->name[1]!='H'&&a->tatm->id<=4&&t) {
		aa[n++]=a;
		a->tatm=t;
	}
	else if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id<=4&&t) {
		aa[n++]=a;
		a->tatm=t;
	}
	else bb[j++]=a;
   }  
   
   for(i=0;i<30;i++) {
	if(aa[i]) aa[i]->next=0; 
        if(bb[i]) bb[i]->next=0;
   }
   
   for(i=0;i<j;i++) delete [] bb[j];

   atm=aa[0];a=atm;
   for(i=0;i<n-1;i++) {
	a->next=aa[i+1];
	a=a->next;
   }
   
   addsidechain();
   configure();
}

float Res::colonycoeff() {
   
   	int n=0;float d=0;
        for(Res *r2=more;r2;r2=r2->more) {
                for(Res *r3=more;r3;r3=r3->more) {
                        if(r3==r2) continue;
                        d+=r2->directrmsd(r3,4,100);
                        n++;
                }
        }

	if(n==0) return 0;
	d=d/n;
	d=log(2)/pow(d,3);
	return d;
}

void Res::buryangle(FILE *fp,int from,int to) {

    float d=bury(from,to);

    Atm *atm_temp;
    if(fp)fprintf(fp,"%c%5i  ",name,id);
    for(atm_temp=atm;atm_temp;atm_temp=atm_temp->next) atm_temp->dihedral(fp);
    fprintf(fp," %f",d);
    if(fp)fprintf(fp,"\n");
}

void Res::addhbondlist(Atm *a,Atm *b,float e,int t) {

    if(hbond==0) hbond=new HBondList();
	
    hbond->add(a,b,e,t);
}

int Res::ishbondedwithres(int n){

    HBondList *s=hbond;
	
    Res *r=chn->isres(n);
    if(r==0) return 0;

    for(s=hbond;s;s=s->next) {
	if(s->donor->res==r||s->acceptor->res==r) return 1;
    }

    return 0;
}

int Res::ishbondedwithres(Res *r){

    if(r==0) return 0;

    HBondList *s=hbond;

    for(s=hbond;s;s=s->next) {
        if(s->donor->res==r||s->acceptor->res==r) return 1;
    }

    return 0;
}


int Res::ishbondedwithatm(Atm *a,Atm *b){

    HBondList *s=hbond;

    for(s=hbond;s;s=s->next) {
        if(s->donor==a&&s->acceptor==b) return 1;
	if(s->donor==b&&s->acceptor==a) return 1;
    }

    return 0;


}

void Res::writehbondlist(FILE *fp) {

    HBondList *s=hbond;
    if(fp==0) fp=stdout;
    for(s=hbond;s;s=s->next) {
	if(s->donor->name[1]=='S'||s->acceptor->name[1]=='S') continue;
	if(s->donor->res!=this) continue;
	float d=TRES.distance(s->donor,s->acceptor);
	fprintf(fp,"%c%5i %s %c%5i %s %c %5.2f\n",s->donor->res->name,s->donor->res->oid,s->donor->name,s->acceptor->res->name,s->acceptor->res->oid,s->acceptor->name, s->donor->res->chn->id, d);
    }
}

void Res::writessbondlist(FILE *fp) {

    if(name!='C') return;
    HBondList *s=hbond;
    if(fp==0) fp=stdout;
    for(s=hbond;s;s=s->next) {
	
        float d=TRES.distance(s->donor,s->acceptor);
 	if(s->donor->name[1]!='S'||s->acceptor->name[1]!='S') continue;
	if(s->donor->res!=this) continue;
        fprintf(fp,"%c%5i %s %c%5i %s %c %5.2f\n",s->donor->res->name,s->donor->res->oid,s->donor->name,s->acceptor->res->name,s->acceptor->res->oid,s->acceptor->name, s->donor->res->chn->id,d);
    }
}


int Res::ishbondedsource(Res *s){

	Atm *n=isatmid(0); 
	if(s==0||n==0) return 0;
	Atm *o=s->isatmid(3);

	if(o==0) return 0;

	return ishbondedwithatm(n,o);
	

}

void Res::setflg(int f) {

	for(Atm *a=atm;a;a=a->next) {
		a->flag=f;
	}
}

void Res::setoid(int n) {
	for(Atm *a=atm;a;a=a->next) {
                a->oid=a->id0+n;
        }
}

float Res::getarea() {
float e=0;
for(Atm *a=atm;a;a=a->next) {
        e+=a->area; 
}
return e;
}

void Res::addambgt(char *line) {

	
	Tres *t=TRES[line+17];

	if(t!=tres) return;

	if(t==0) return;

	int i=getuseambgt(line+12,1,3);

	if(i>0) return;

	Tatm * tatm_temp=tres->isatm(line+12,1,3,i);

	Atm *a=new Atm();

	for(i=0;i<3;i++) a->xyz[i]=atof(line+30+8*i);

	
	strncpy(a->name,tatm_temp->name,4);//strncpy(atm_temp->name,line+12,4)

	a->name[4]='\0';

	if(ambgt==0) ambgt=a;
	else {
		Atm *b;
		for(b=ambgt;b->next;b=b->next);
		b->next=a;
	}
	a->res=this;
	a->tatm=tatm_temp;
	return;
}


float Res::ambgtrmsd(int n) {

	if(ambgt==0) return 0;

	Atm *a;
	
	int nn=0;
	float dt=0;
	for(a=ambgt;a;a=a->next) {
		Atm *b=isatm(a->name);
		if(b==0) continue;
		if(n==0) { //all
			float d=TRES.distance(a->xyz,b->xyz);
			dt+=d*d;
			nn++;
		}
		else if(n==1) { //only sidechain
			
			if(b->tatm->name[1]=='H'&&b->tatm->bond[0]->id<=4) continue;
			else if(b->tatm->id<=4) continue;
			float d=TRES.distance(a->xyz,b->xyz);
			dt+=d*d;
			nn++;
		}
		else if(n==2) { //only backbone
			if(b->tatm->name[1]=='H'&&b->tatm->bond[0]->id>4) continue;
			else if(b->tatm->id>4) continue;
			 
			float d=TRES.distance(a->xyz,b->xyz);
			dt+=d*d;
			nn++;
		}
		else if(n==3) { //only up to chi angle 2
			
			if(b->tatm->name[1]=='H'&&(b->tatm->bond[0]->rotm>4||b->tatm->bond[0]->rotm<=2)) continue;
			else if(b->tatm->rotm>4||b->tatm->rotm<=2) continue;

			if(b->tatm->name[1]=='H'&&b->tatm->bond[0]->id<=4) continue;
                        else if(b->tatm->id<=4) continue;
			 

			float d=TRES.distance(a->xyz,b->xyz);
			dt+=d*d;
			nn++;
		}
		else if(n==4) { //only beyond chi angle 2
			
			if(b->tatm->name[1]=='H'&&b->tatm->bond[0]->rotm<=4) continue;
			else if(b->tatm->rotm<=4) continue;
			 
			float d=TRES.distance(a->xyz,b->xyz);
			dt+=d*d;
			nn++;
		}
	}
	if(nn==0) return 0;
	else return sqrt(dt/nn);
}

void Res::setfreenextonmore(){
	for(Res *rr=more;rr;rr=rr->more) rr->next=0;
}

void Res::write(char *s){
	FILE *fp=fopen(s,"w");
	if(fp==0) return;
	write(fp);
	fclose(fp);
}

void Res::transferbackbone(Res *s){
	
	Atm *a1,*a2;
	a1=s->isatm(" N  ");
	a2=isatm(" N  ");
	if(a1&&a2)a2->transfer(a1);

	a1=s->isatm(" CA ");
        a2=isatm(" CA ");
        if(a1&&a2)a2->transfer(a1);

	a1=s->isatm(" C  ");
        a2=isatm(" C  ");
        if(a1&&a2)a2->transfer(a1);

	a1=s->isatm(" O  ");
        a2=isatm(" O  ");
        if(a1&&a2)a2->transfer(a1);

	a1=s->isatm(" CB ");
        a2=isatm(" CB ");
        if(a1&&a2) a2->transfer(a1);

	a1=s->isatm(" HN ");
        a2=isatm(" HN ");
        if(a1&&a2) a2->transfer(a1);

	a1=s->isatm(" HA ");
        a2=isatm(" HA ");
        if(a1&&a2) a2->transfer(a1);

	a1=s->isatm(" HA ");
        a2=isatm(" HA ");
        if(a1&&a2) a2->transfer(a1);

	a1=s->isatm("1HA ");
        a2=isatm("1HA ");
        if(a1&&a2) a2->transfer(a1);


	a1=s->isatm("2HA ");
        a2=isatm("2HA ");
        if(a1&&a2) a2->transfer(a1);

	a1=s->isatm("1HB ");
        a2=isatm("1HB ");
        if(a1&&a2) a2->transfer(a1);

	a1=s->isatm("2HB ");
        a2=isatm("2HB ");
        if(a1&&a2) a2->transfer(a1);

	a1=s->isatm("3HB ");
        a2=isatm("3HB ");
        if(a1&&a2) a2->transfer(a1);

}

void Res::deletecontact(){
	for(Atm *a=atm;a;a=a->next) {
		a->deletecontact();
	}
}
/*
int Res::iscomplete(Res *r) {

	Tatm *a;
	for(a=tres->tatm;a;a=a->next) {
		if(a->name[1]=='H') continue;
		Atm *b=isatmid(a->id);
		if(
	}
}
*/
