#include"source.h"

void Atm::initial()
{
int i;
res=0;
next=0;
for(i=0;i<4;i++) {bond[i]=0;bondlen[i]=0;}
for(i=0;i<2;i++) hbond[i]=0;
near=0;
nnear=0;
tatm=0;
strcpy(name," ???\0");
id=0;
id0=0;
for(i=0;i<3;i++) xyz[i]=0;
occup=1;
bfact=0;
chi=999;
area=0;
good=1;
flag=0;
temp=0;
nemp=0;
oid=0;
mask=0;
contact=0;
energy=0;
}

Atm::Atm()
{
initial();
}

Atm::Atm(Atm *s)
{
 int i;

 initial();

 tatm=s->tatm;
 strcpy(name,s->name);
 id=s->id;
 id0=s->id0;
 oid=s->oid;
 for(i=0;i<3;i++) xyz[i]=s->xyz[i];
 occup=s->occup;
 bfact=s->bfact;
 mask=s->mask;
}

Atm::Atm(Tatm *s)
{
 int i;

 initial();

 tatm=s;
 strcpy(name,s->name);
 id=s->id;
 for(i=0;i<3;i++) xyz[i]=s->xyz[i];
}

Atm::~Atm()
{
if(next) {delete next;next=0;}
if(near) {delete [] near;near=0;}
if(nnear){delete [] nnear;nnear=0;}
if(temp) {delete [] temp;temp=0;}
if(contact) {delete contact;contact=0;}
}
 
void Atm::write(FILE *fp,int f)
{
char line[256];
if(res->chn->pdb&&res->chn->pdb->name==0)
{
 res->chn->pdb->name=new char[10];
 strcpy(res->chn->pdb->name,"unknown");
}
if(f==0) sprintf(line,"%d%s",id0,name);
if(f==1) sprintf(line,"%d%c %d%s",res->id,res->name,id0,name);
if(f==2) sprintf(line,"%c %d%c %d%s",res->chn->id,res->id,res->name,id0,name); 
if(f==3) sprintf(line,"%s %c %d%c %d%s",res->chn->pdb->name,
                       res->chn->id, res->id,res->name,id0,name);
if(f==4) sprintf(line,"%s %c %8.3f %8.3f %8.3f\n",tatm->name,res->tres->name,
                       xyz[0],xyz[1],xyz[2]);
if(f==5) 
{
 if(!(res->tres->name=='P'||res->tres->name=='A'))
 {
   if(tatm->name[1]=='H')
   {
     if(tatm->bond[0]->id>=4) return;
   }
   else
   {
     if(tatm->id>4) return;
   }  
 }
 sprintf(line,"%s %c %8.3f %8.3f %8.3f\n",tatm->name,res->tres->name,
                       xyz[0],xyz[1],xyz[2]);
}
fprintf(fp,"%s",line);
return;
}

void Atm::write(FILE *fp)
{
char line[216],line0[20];
if(res->chn->pdb&&res->chn->pdb->name==0)
{
 res->chn->pdb->name=new char[10];
 strcpy(res->chn->pdb->name,"unknown");
}
if(res->chn->pdb&&res->chn->pdb->name) strncpy(line0,res->chn->pdb->name,20);
else strcpy(line0,"unknown");
line0[19]='\0';
line[0]='\0';
if(res->chn->ishet==1)  sprintf(line,"HETATM%5d",id0+1);
else sprintf(line,"ATOM  %5d",id0+1);

sprintf(line,"%s %s ",line,name);
if(res->chn->ishet==1) sprintf(line,"%s%s ",line,res->tres->name3);
else sprintf(line,"%s%s ",line,TRES.swap(res->name));
sprintf(line,"%s%c",line,res->chn->id);
sprintf(line,"%s%4d",line,res->chn->start+res->id);
sprintf(line,"%s    %8.3f%8.3f%8.3f%6.2f%6.2f    %s\n",line,xyz[0],
                    xyz[1],xyz[2], occup,bfact,line0);
fprintf(fp,line);
}

void Atm::writeold(FILE *fp)
{
char line[216],line0[20];
if(res->chn->pdb&&res->chn->pdb->name==0)
{
 res->chn->pdb->name=new char[10];
 strcpy(res->chn->pdb->name,"unknown");
}
if(res->chn->pdb&&res->chn->pdb->name) strncpy(line0,res->chn->pdb->name,20);
else strcpy(line0,"unknown");
line0[19]='\0';
line[0]='\0';
if(res->chn->ishet==1) sprintf(line,"HETATM%5d",oid);
else sprintf(line,"ATOM  %5d",oid);
sprintf(line,"%s %s ",line,name);
if(res->chn->ishet==1) sprintf(line,"%s%s ",line,res->tres->name3);
else sprintf(line,"%s%s ",line,TRES.swap(res->name));
sprintf(line,"%s%c",line,res->chn->id);
sprintf(line,"%s%4d",line,res->oid);
sprintf(line,"%s    %8.3f%8.3f%8.3f%6.2f%6.2f    %s\n",line,xyz[0],
                    xyz[1],xyz[2], occup,bfact,line0);
fprintf(fp,line);
}

void Atm::writerescard(FILE *fp)
{
char line[216],line0[20];
if(res->chn->pdb&&res->chn->pdb->name==0)
{
 res->chn->pdb->name=new char[10];
 strcpy(res->chn->pdb->name,"unknown");
}
if(res->chn->pdb&&res->chn->pdb->name) strncpy(line0,res->chn->pdb->name,20);
else strcpy(line0,"unknown");
line0[19]='\0';
line[0]='\0';
if(res->chn->ishet==1) sprintf(line,"HETATM%5d",oid);
else sprintf(line,"ATOM  %5d",oid);
sprintf(line,"%s %s ",line,name);
if(res->chn->ishet==1) sprintf(line,"%s%s ",line,res->tres->name3);
else sprintf(line,"%s%s ",line,TRES.swap(res->name));
sprintf(line,"%s%c",line,res->chn->id);
//sprintf(line,"%s%4d",line,res->oid);
sprintf(line,"%s%s",line,res->rescard);
sprintf(line,"%s   %8.3f%8.3f%8.3f%6.2f%6.2f    %s\n",line,xyz[0],
                    xyz[1],xyz[2], occup,bfact,line0);
fprintf(fp,line);
}

void Atm::writeold(FILE *fp,char opt)
{
char line[216],line0[20];
if(res->chn->pdb&&res->chn->pdb->name==0)
{
 res->chn->pdb->name=new char[10];
 strcpy(res->chn->pdb->name,"unknown");
}
if(res->chn->pdb&&res->chn->pdb->name) strncpy(line0,res->chn->pdb->name,20);
else strcpy(line0,"unknown");
line0[19]='\0';
line[0]='\0';
if(res->chn->ishet==1) sprintf(line,"HETATM%5d",oid);
else sprintf(line,"ATOM  %5d",oid);
sprintf(line,"%s %s%c",line,name,opt);
if(res->chn->ishet==1) sprintf(line,"%s%s ",line,res->tres->name3);
else sprintf(line,"%s%s ",line,TRES.swap(res->name));
sprintf(line,"%s%c",line,res->chn->id);
sprintf(line,"%s%4d",line,res->oid);
sprintf(line,"%s    %8.3f%8.3f%8.3f%6.2f%6.2f    %s\n",line,xyz[0],
                    xyz[1],xyz[2], occup,bfact,line0);
fprintf(fp,line);
}

float Atm::dihedral() {
  return dihedral(0);
}

float Atm::dihedral(FILE *fp)
{
Atm *atm0[4];
int i;
if(tatm->rotate==0) return 999;
atm0[0]=bond[1];
atm0[1]=this;
atm0[2]=bond[0];
if(atm0[2]&&atm0[2]->bond[0]!=this) atm0[3]=atm0[2]->bond[0];
else {chi=999;goto re200;}
for(i=0;i<4;i++) if(atm0[i]==0) {chi=999;goto re200;}
i=atm0[0]->good*atm0[1]->good*atm0[2]->good*atm0[3]->good;
if(i==0) {chi=999;goto re200;}
//make changes so that a1->a2->a3->a4 in a stream
//cout<<"old";
chi=TRES.dihedral(atm0[0]->xyz,atm0[1]->xyz,atm0[2]->xyz,atm0[3]->xyz);
//chi=TRES.dihedral(atm0[3]->xyz,atm0[2]->xyz,atm0[1]->xyz,atm0[0]->xyz);
re200:
if(fp)fprintf(fp,"  %6.1f",chi);
return chi; 
}

float Atm::gettorsionangle()
{
Atm *atm0[4];
//int i;
//if(tatm->rotate==0) return 999;
atm0[0]=bond[1];
if(atm0[0]==0) return 999;
atm0[1]=this;
atm0[2]=bond[0];
if(atm0[2]&&atm0[2]->bond[0]!=this) atm0[3]=atm0[2]->bond[0];
else {return chi=999;}
//for(i=0;i<4;i++) if(atm0[i]==0) {chi=999;goto re200;}
//i=atm0[0]->good*atm0[1]->good*atm0[2]->good*atm0[3]->good;
//if(i==0) {chi=999;goto re200;}
//make changes so that a1->a2->a3->a4 in a stream
//cout<<"old";
if(atm0[3]==0) return 999;
chi=TRES.dihedral(atm0[0]->xyz,atm0[1]->xyz,atm0[2]->xyz,atm0[3]->xyz);
//chi=TRES.dihedral(atm0[3]->xyz,atm0[2]->xyz,atm0[1]->xyz,atm0[0]->xyz);
//re200:
//if(fp)fprintf(fp,"  %6.1f",chi);
return chi;
}



float Atm::dihedral(FILE *fp,int m)
{
Atm *a1,*a2,*a3;
int i,j,k;
float d;
a1=bond[0];
if(m==-2||m==-3) {if(fp&&tatm->rotate)fprintf(fp,"  %6.1f",chi);return chi;}

if(tatm->eng->mchi==0)  {d=999;goto re200;}
if(a1==0)  {d=999;goto re200;}
if(good==0||a1->good==0)  {d=999;goto re200;}

k=-1; 
for(i=0;i<a1->tatm->nbond;i++)
{
  a2=a1->bond[i];
  if(a2==this) continue;
  for(j=1;j<tatm->nbond;j++)
  {
    k++;
    a3=bond[j];
    if(a2==0||a3==0) continue;
    if(a2->good==0||a3->good==0) continue;
    if(m==-1) 
    {
      if(fp) 
      {
        d=TRES.dihedral(a3->xyz,xyz,a1->xyz,a2->xyz);
        if(fp)fprintf(fp,"%s %s %s %s: %6.1f\n",a3->name,name,
                                          a1->name,a2->name,d);
      }
    }
    else if(k==m)   
    {
        d=TRES.dihedral(a3->xyz,xyz,a1->xyz,a2->xyz);
        if(fp) 
        fprintf(fp,"%s %s %s %s: %6.1f\n",a3->name,name,
                                          a1->name,a2->name,d);
        break; 
    }
  }
}
re200:
if(fp)fprintf(fp,"  %6.1f",d);
return d;
}

Atm *Atm::cloneatm() {
   
   Atm *t=new Atm(this);
   if(next) t->next=next->cloneatm();	
   else t->next=0;
   return t;
}

int Atm::isbond(Atm *s)
{
int i;
for(i=0;i<4;i++)
{
if(bond[i]) if(bond[i]->id0==s->id0) return 1;
}
return 0;
}

int Atm::isbond(char *s)
{
int i,j;
i=strlen(s);
for(j=0;j<4;j++)
if(bond[j]&&strncmp(s,bond[j]->name,i)==0) return 1;
return 0;
}

void Atm::transfer(float s)
{
int i;
for(i=0;i<3;i++)xyz[i]=s;
}

void Atm::transfer(Atm *s)
{
int i;
for(i=0;i<3;i++) xyz[i]=s->xyz[i];
}

void Atm::transfer(float *s,int f) 
{
int i;
for(i=0;i<3;i++)
{
if(f) xyz[i]=s[i];
else s[i]=xyz[i];
}
} 

int Atm::isnear(Atm *s) {

	int j;

	if(near==0) return 0;

	for(j=0;j<10000&&near[j];j++) {
		if(near[j]==s) return 1;
	}

	return 0;
}

int Atm::isbondnear(Atm *s) {

        int j;

	if(s==this) return 0;
        if(near==0) return 0;

        for(j=0;j<10000&&near[j];j++) {
                if(near[j]==s) return nnear[j];
        }

        return 0;
}

int Atm::isnewbondnear(Atm *s) {

        int j;

        if(s==this) return 0;
        if(near==0) return -1;

        for(j=0;j<10000&&near[j];j++) {
                if(near[j]==s) return nnear[j];
        }

        return -1;
}

int Atm::isnear(Atm *s,int dis)
{
//this program is used to calculate the bond distance between
//two atoms. the maximum distance is 4
int i,j,k,m;
Atm *a0,*a1,*a2,*a3,*a4;

a0=this;
i=s->res->id-a0->res->id;

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
  for(j=0;j<a1->tatm->nbond;j++)
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
  for(j=0;j<a1->tatm->nbond;j++)
  {
    a2=a1->bond[j];
    if(a2==0) continue;
    if(a2==a0||a2==a1) continue;

    for(k=0;k<a2->tatm->nbond;k++)
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
  for(j=0;j<a1->tatm->nbond;j++)
  {
    a2=a1->bond[j];
    if(a2==0) continue;
    if(a2==a0||a2==a1) continue;

    for(k=0;k<a2->tatm->nbond;k++)
    {
     a3=a2->bond[k];
     if(a3==0) continue;
     if(a3==a1||a3==a2||a3==a0) continue;

     for(m=0;m<a3->tatm->nbond;m++)
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

int Atm::dsspdefinedhbond(Atm *s) {

Atm *n,*o;
Atm *a=this;

if(a==0||s==0) return 0;
if(s->res->chn!=a->res->chn) return 0; //not same chain
if(fabs(s->res->id-a->res->id)<=1) return 0;
int nn=a->tatm->hbond*s->tatm->hbond;
if(nn!=2&&nn!=3&&nn!=9&&nn!=6) return 0;

if(a->tatm->hbond==2&&s->tatm->hbond==1)
{
	n=s; o=a;
}
else if(a->tatm->hbond==1&&s->tatm->hbond==2)
{
	n=a; o=s;
}
else if(a->tatm->hbond==1&&s->tatm->hbond==3)
{
	return a->ishbondexist(s,0,3.5,30);
}
else if(a->tatm->hbond==3&&s->tatm->hbond==1)
{
	return a->ishbondexist(s,0,3.5,30);
}
else if(a->tatm->hbond==2&&s->tatm->hbond==3)
{
	return a->ishbondexist(s,0,3.5,30);
}
else if(a->tatm->hbond==3&&s->tatm->hbond==2)
{
	return a->ishbondexist(s,0,3.5,30);
}
else
{
   return 0;
}

if(n->tatm->id!=0) return a->ishbondexist(s,0,3.5,30);

Atm *c=o->bond[0];
if(c==0) return 0;
float h[3];

Atm *u=n->res->isatm(" HN ");
if(u==0){
	float dx,dy,dz;
	if(n->res->temp==0) return 0;
	dx=n->res->temp[3]-n->res->temp[6];
	dy=n->res->temp[4]-n->res->temp[7];	
	dz=n->res->temp[5]-n->res->temp[8];	
	float dist=dx*dx + dy*dy + dz*dz;
	float dco=sqrt(dist);
	h[0]=n->xyz[0]+dx/dco;
	h[1]=n->xyz[1]+dy/dco;
	h[2]=n->xyz[2]+dz/dco;
}
else {
	for(int i=0;i<3;i++) h[i]=u->xyz[i];
}

float dnc=TRES.distance(n->xyz,c->xyz);
float dno=TRES.distance(n->xyz,o->xyz);
float dho=TRES.distance(h,o->xyz);
float dhc=TRES.distance(h,c->xyz);
float Q=-27888.;
float HBLOW=-9900;
float HBHIGH=-500;

float hbe;
if(dho<0.5||dhc<0.5||dnc<0.5||dno<0.5) hbe=HBLOW;
else {
	hbe = (long)floor(Q / dho - Q / dhc + Q / dnc - Q / dno + 0.5);		
}	
if(hbe <= HBLOW) hbe=HBLOW;
		
if(hbe>HBHIGH) return 0;
else return 1;
}

int Atm::ishbondexist(Atm *s) {
	return ishbondexist(s,0,3.5,30);
}

int Atm::ishbondexist(Atm *s,float dmin,float dmax,float ang) {

Atm *a,*a0,*s0;
float theta,phita;
float d,anga,angs;

if(s==0) return 0;

//dis=0.5;
//ang=30;
a=this;
a0=a->bond[0];
s0=s->bond[0];

if(a0==0||s0==0) return 0; //not complete atom ancetor

if(s->res->chn!=a->res->chn) return 0; //not same chain
if(fabs(s->res->id-a->res->id)<=1) return 0;
//if(TRES.distance(a->xyz,s->xyz)>dmax) return 0;
int nn=a->tatm->hbond*s->tatm->hbond;
if(nn!=2&&nn!=3&&nn!=9&&nn!=6) return 0;

if(a->tatm->hbond==0||s->tatm->hbond==0)
{
   return 0;
}
else if(a->tatm->hbond==5||s->tatm->hbond==5)
{
   return 0;
}
else if(a->tatm->hbond==2&&s->tatm->hbond==1&&s->tatm->id!=0){
   theta=109;phita=120;
}
else if(s->tatm->hbond==2&&a->tatm->hbond==1&&a->tatm->id!=0){
   theta=120;phita=110;
}
else if(a->tatm->hbond==2&&s->tatm->hbond==1) 
{ 
   theta=150;phita=120;
}
else if(a->tatm->hbond==1&&s->tatm->hbond==2)
{
   theta=120;phita=150;
}
else if(a->tatm->hbond==1&&s->tatm->hbond==3)
{
   theta=120;phita=108;
} 
else if(a->tatm->hbond==3&&s->tatm->hbond==1)
{
   theta=108;phita=120;
}
else if(a->tatm->hbond==2&&s->tatm->hbond==3)
{
   theta=110;phita=109;
}
else if(a->tatm->hbond==3&&s->tatm->hbond==2)
{
   theta=109;phita=110;
}
else
{
   return 0;
}

//calculate the energy
//x=3.14/180;
d=TRES.distance(a->xyz,s->xyz); 
if(d>dmax||d<dmin) return 0; //not meet distance criteria
anga=TRES.angle(a0->xyz,a->xyz,s->xyz);
if(anga<theta-ang||anga>theta+ang) return 0;
angs=TRES.angle(s0->xyz,s->xyz,a->xyz);
if(angs<phita-ang||angs>phita+ang) return 0;
//d=fabs(d-3);
//dd=cos(x*(anga-theta))*cos(x*(angs-phita))*exp(-d);
//if(dd<0) dd=0;
return 1;
}

float Atm::ishbond(Atm *s,float dis,float  ang,float eng)
{
Atm *a,*a0,*s0;
float theta,phita;
float d,anga,angs,dd,x;

a=this;
if(s==0) return 0;
a0=a->bond[0];
s0=s->bond[0]; 
if(a0==0||s0==0) return 0; //not complete atom ancetor

//if(s->res->chn!=a->res->chn) return 0; //not same chain
if(fabs(s->res->id-a->res->id)<=1&&s->res->chn==a->res->chn) return 0;
int nn=a->tatm->hbond*s->tatm->hbond;
if(nn!=2&&nn!=3&&nn!=9&&nn!=6) return 0;

//find the angle for hbond!

if(a->tatm->hbond==0||s->tatm->hbond==0)
{
   return 0;
}
else if(a->tatm->hbond==5||s->tatm->hbond==5)
{
   return 0;
}
else if(a->tatm->hbond==2&&s->tatm->hbond==1&&s->tatm->id!=0){
   theta=110;phita=120;
}
else if(a->tatm->hbond==1&&s->tatm->hbond==2&&a->tatm->id!=0)
{
   theta=120;phita=110;
}
else if(a->tatm->hbond==2&&s->tatm->hbond==1) 
{ 
   theta=150;phita=120;
}
else if(a->tatm->hbond==1&&s->tatm->hbond==2)
{
   theta=120;phita=150;
}
else if(a->tatm->hbond==1&&s->tatm->hbond==3)
{
   theta=120;phita=108;
} 
else if(a->tatm->hbond==3&&s->tatm->hbond==1)
{
   theta=108;phita=120;
}
else if(a->tatm->hbond==2&&s->tatm->hbond==3)
{
   theta=110;phita=109;
}
else if(a->tatm->hbond==3&&s->tatm->hbond==2)
{
   theta=109;phita=110;
}
else
{
   return 0;
}

//calculate the energy
x=3.14/180;
d=TRES.distance(a->xyz,s->xyz);
if(d>3.0+dis||d<3.0-dis) return 0; //not meet distance criteria
anga=TRES.angle(a0->xyz,a->xyz,s->xyz);
if(anga<theta-ang||anga>theta+ang) return 0;
angs=TRES.angle(s0->xyz,s->xyz,a->xyz);
if(angs<phita-ang||angs>phita+ang) return 0;
d=fabs(d-3);
dd=cos(x*(anga-theta))*cos(x*(angs-phita))*exp(-3*d);
if(dd<0) dd=0;
return dd*eng;
}

void Atm::writehbondparm(Atm *s)
{
Atm *a,*a0,*s0;
//float theta,phita;
float d,anga,angs;

a=this;
if(s==0) return;
a0=a->bond[0];
s0=s->bond[0]; 
if(a0==0||s0==0) return ; //not complete atom ancetor
 
//calculate the energy
//float x=3.14/180;
d=TRES.distance(a->xyz,s->xyz);
cerr<<res->name<<res->id<<" "<<name<<" "<<s->res->name<<s->res->id<<" "<<s->name<<":"<<d<<" ";
d=TRES.distance(a->xyz,s0->xyz);
cerr<<d<<" ";
d=TRES.distance(a0->xyz,s->xyz);
cerr<<d<<" ";
d=TRES.distance(a0->xyz,s0->xyz);
cerr<<d<<" ";
anga=TRES.angle(a0->xyz,a->xyz,s->xyz);
cerr<<anga<<" ";
 
angs=TRES.angle(s0->xyz,s->xyz,a->xyz);
cerr<<angs<<" "<<endl;
 
 
}
float Atm::clash(Atm *a)
{
Atm *s;
int i,j;
float e,d,r,c;
float ta,tb,tc,tn;

i=4*TRES.smt;
ta=TRES.smooth[i+0];
tb=TRES.smooth[i+1];
tc=TRES.smooth[i+2];
tn=TRES.smooth[i+3];

e=0;
s=this;
r=a->tatm->eng->radius+s->tatm->eng->radius;
j=s->tatm->hbond*a->tatm->hbond;
if(j==2||j==3||j==6||j==9) r=3.0;
//if(s->tatm->hbond==5) r-=s->tatm->eng->radius;
//if(a->tatm->hbond==5) r-=a->tatm->eng->radius;
if(a->tatm->name[1]=='H'&&s->tatm->name[1]=='H'&&r>1.)r=1.;
if(r<=0) return 0;
d=TRES.distsqr(a->xyz,s->xyz);
//float de=d;
d=d/r/r;
if(d<0.1) d=0.1;
else if(d>4)   return 0;
c=sqrt(s->tatm->eng->epslon*a->tatm->eng->epslon);
c=ta*c*exp(-d*tb);
d=1/d;
d=pow(d,tn/2);
e=c*d*(d-tc);
return e;
}


float Atm::torsion(int f)
{
Atm *a0,*a1,*a2,*a3,*s;
int i,j,k;
float d,e,e1,e2,e3;
s=this;
a0=s;
a1=a0->bond[0];
if(a1==0) return 0;
if(a1==0||s->tatm->eng->mangle==0) return 0;
k=-1;
e=0;
for(i=0;i<a1->tatm->nbond;i++)
{
  a2=a1->bond[i];
  if(a2==0) continue;
  if(a2==a0) continue;
  for(j=1;j<a0->tatm->nbond;j++)
  {
    k++;
    a3=a0->bond[j];
    if(a3==0) continue;
    if(f==-2)
    {
      if(a3->name[1]=='H'||a2->name[1]=='H') continue;
    }
    else if(f>=0)
    {
      if(k!=f) continue;
    }
    d=TRES.dihedral(a3->xyz,a0->xyz,a1->xyz,a2->xyz);
    e3=d*s->tatm->eng->nchi[k]-s->tatm->eng->chi0[k];
    e3=e3/180*3.1416;
    e1=s->tatm->eng->kchi[k]*(1+cos(e3));
    if(s->tatm->eng->nchi2[k]>9999) e2=10000;
    else e2=s->tatm->eng->kchi2[k]*(1+cos(d*s->tatm->eng->nchi2[k]-
            s->tatm->eng->chi02[k]));
    if(e1>e2) e1=e2; 
    e+=e1;
  }
}
return e;
}
  


void Atm::bondleng()
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
      bondlen[i]=TRES.distance(xyz,bond[i]->xyz);
    }
  }
}

int Atm::allnear(int dis,int flg)
{
//flg==0; itself plus neighboring residues
//flg==1; itself onlu
//flg==-1; neighboring only

  Atm *bb[100];
  int i,j,aa[100];
  Res *r;
  Atm *a;
  
  if(near) {delete [] near;near=0;}
  if(nnear) {delete [] nnear;nnear=0;}
  i=0;

  if(flg==0||flg==-1)
  {
    r=res->next;
  
    if(r&&r->id-res->id==1)
    {
      for(a=r->atm;a;a=a->next)
      {
        j=isnear(a,dis);
        if(j>0&&j<=dis) {bb[i]=a;aa[i]=j;i++;}
      }
    }
  }

  if(flg==0||flg==1)
  {
    for(a=res->atm;a;a=a->next)
    {
      j=isnear(a,dis);
      if(j>0&&j<=dis) {bb[i]=a;aa[i]=j;i++;}
    }
  }

  if(flg==0||flg==-1)
  {
    r=res->chn->isres(res->id-1);
    if(r)
    {
      for(a=r->atm;a;a=a->next)
      {
        j=isnear(a,dis);
        if(j>0&&j<=dis) {bb[i]=a;aa[i]=j;i++;}
      }
    }
  }

  near=new Atm*[i+1];
  nnear=new int[i+1];
  for(j=0;j<i+1;j++) near[j]=0;
  for(j=0;j<i+1;j++) nnear[j]=0;
  for(j=0;j<i;j++) {near[j]=bb[j]; nnear[j]=aa[j];}

  return i-1;
}

float Atm::surface(Lattice *lattice,int m,float probe)
{
Tres *tr;
Tatm *ta;
Atm *a;
float xyz[3];
float t1,t2,f1,f2,c1,c2,c3;
int i,j,k;
int flag0;

area=0;

//add the probe to atoms radius;

for(tr=&TRES;tr;tr=tr->next)
for(ta=tr->tatm;ta;ta=ta->next)
ta->eng->radius+=probe;

lattice->radmax=lattice->radmax+probe;
if(m==0) return 0;
t1=3.1415926/(m+1);
t2=2*3.1415926/(2*m+1);
c3=t1*t2;

//set up the lattice

flag0=lattice->flag;
lattice->flag=1;

a=this;
lattice->getcell(a,a->tatm->eng->radius+lattice->radmax);
lattice->cutoff(a);
for(i=0;i<m+1;i++)
for(j=0;j<2*m+1;j++)
{  
  f1=i*t1;f2=j*t2;
  c1=a->tatm->eng->radius*sin(f1);
  c2=a->tatm->eng->radius*a->tatm->eng->radius;
  xyz[0]=c1*cos(f2)+a->xyz[0];
  xyz[1]=c1*sin(f2)+a->xyz[1];
  xyz[2]=a->tatm->eng->radius*cos(f1)+a->xyz[2];
  k=lattice->cover(xyz,0);
  if(k) continue;
  a->area+=c3*c2*sin(f1);
}

lattice->flag=flag0;
lattice->radmax=lattice->radmax-probe;
for(tr=&TRES;tr;tr=tr->next)
for(ta=tr->tatm;ta;ta=ta->next)
ta->eng->radius-=probe;
return area;
}

int Atm::getnumberhbonds() {

        int n=0;
        for(int i=0;i<4;i++) {
                if(bond[i]&&bond[i]->name[1]=='H') n++;
        }
        return n;
}

void Atm::deletecontact(){
	if(contact) delete contact;contact=0;	
}
Atm *Atm::ischildatm(char *s){
	//Atm *temp;
	int i;
	for(i=0;i<4;i++) {
		if(bond[i]==0) continue;
		if(strncmp(bond[i]->name,s,4)==0) return bond[i];
	}
	return 0;
}
