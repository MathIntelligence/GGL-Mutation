#include"source.h"

Scap::Scap()
{
perfil=70;
prbrad=1.6;
fastmodeH=0;
fastmode=0;
layermode=0;
proring=0;
hydroeng=0.025;
hydroflg=0;
hbondeng=-1.0;
niterate=5;
allcharges=0;
mutationonly=0;
pdb=0;
rotnum=0;
rotamer=0;
flag=0;
cut=6.0;
ener=0;
lattice0=new Lattice;
lattice=new Lattice;
lattice1=new Lattice;
seqn=0;
mine=100000;
bmax=1;
strcpy(force,"14");
arbt=0;
nummore=30;
rott=0;
cbetaflg=0;
cbeta=0;
tormax=1;
nout=0;
cutbb=5;
colony=0;
armsd=0;
ncolony=2;
colonyline=0;
resultout=1;
bordertorsion=0;
singletorsion=0;
initside=0;
ring=0;
disulfide=0;
torsionring=0;
coleng=0;
nncolony=2;
coleffect=0.5;
solventeff=0;
solventorder=0;
solventbury=0;
//solventdirectbury=0;
solventcolony=0;
iterateweight=1.0;
iweight=0;
dielectric=1;
outinitial=0;
dsqrt=0;
outenforce=0;
outforcemag=0;
outcharge=0;
solventorderrandom=0;
printouteng=0;
//ermsd=0;
pdborg=0;
hydrophobic=0;
adaptvdw=0;
emilcharge=0;
latsurfv=0;
emilhydro=0;
emilnewcharge=0;
delphi=0;
delfile=0;
indi=4.0;
exdi=80;
delphiself=1;
delphilinit=30;
delphiscale=0.5;
delphigsize=45;
delphionself=1;
includeself=0;
}

Scap::~Scap() 
{
if(lattice0) {delete lattice0;lattice0=0;}
if(lattice)  {delete lattice; lattice=0;}
if(lattice1) {delete lattice1; lattice1=0;}
if(latsurfv) {delete latsurfv;latsurfv=0;}
if(solventbury) delete [] solventbury; solventbury=0;
if(solventorder) delete [] solventorder; solventorder=0;
if(delfile) delete [] delfile;delfile=0;
if(ener) delete [] ener;
if(armsd) delete [] armsd;
if(cbetaflg) delete [] cbetaflg;
for(int i=0;rotamer&&i<pdb->getresnum();i++)
{
 for(int j=0;rotnum&&j<rotnum[i];j++)
 {
   if(rotamer&&rotamer[i][j]) {delete [] rotamer[i][j];rotamer[i][j]=0;}
 }
 if(rotamer&&rotamer[i]) delete [] rotamer[i];rotamer[i]=0;
}
if(rotnum)delete [] rotnum;

}

void Scap::sidelib(char *filnam)
{
FILE *fp;
Tres *tres;
Res *res_temp,*res_temp0;
Atm *atm_temp;
char *resmer,name[200],line[200];
float x;
int i,k,n,j,m,dd[30],c;
Chn *chn;

//calculate how many rotamers for each kind of residues

if(rott==0) goto re200;
resmer=RCS["library"];
strcpy(name,filnam);
sprintf(name,"%s/%s",resmer,filnam);
i=strlen(resmer);
while(((fp=fopen(name,"r"))==NULL))
{
  resmer=resmer+i+1;
  i=strlen(resmer);
  if(*resmer=='\0')
  {cerr<<"warning:could not open: "<<filnam<<endl;exit(0);}
  sprintf(name,"%s/%s",resmer,filnam);
}

//decide the matrix should be opened

for(i=0;i<30;i++) 
{ 
 if(flag%10==0)dd[i]=2;
 else dd[i]=500; 
}

if(flag%10==0)
{
 while(fgets(line,256,fp)!=NULL)
 {
  c=TRES.swap(line);
  if(c=='?') continue;
  j=c-'A';
  dd[j]++;
 }
} 
 
//allocating space 

k=flag/10000;
i=0;
for(tres=&TRES;tres;tres=tres->next)
{
  n=tres->name-'A';
  j=tres->numrot-2;
  j=min(j,k);
  dd[n]*=(int) pow(3,j);
  i+=dd[n];
}


cout<<"the total rotamer is: "<<i<<endl;

re200:

if(rott==0)
{
  for(tres=&TRES;tres;tres=tres->next)
  dd[tres->name-'A']=1;
} 

rotnum=0;
rotamer=0;

n=pdb->getresnum();
rotnum=new int[n];
rotamer=new float**[n];

for(i=0;i<n;i++) {rotnum[i]=0;rotamer[i]=0;}
if(rott==0) goto re300;
for(chn=pdb->chn;chn;chn=chn->next)
for(res_temp=chn->res;res_temp;res_temp=res_temp->next)
{
 if(res_temp->flag<-9999)continue;
 i=res_temp->id0;
 if(flag%10==0)
 {
  n=res_temp->tres->numrot-1;
  if(n==0) n=1; 
 }
 else 
 {
  n=6;
 }
 j=res_temp->name-'A';
 m=dd[j];
 rotnum[i]=0;
 rotamer[i]=new float*[m];
 for(j=0;j<m;j++)
 {
  rotamer[i][j]=new float[n];
  for(k=0;k<n;k++)
  {
   rotamer[i][j][k]=0;
  }
 }
}


fclose(fp);
if(flag%10==0)
{
bbind(name);
}
else 
{
bbdep(name);
}

//including its own conformation

if(flag%100>=10)
{
  pdb->next->chn->dihedral();
  for(Chn *chn=pdb->chn;chn;chn=chn->next)
  for(res_temp=chn->res;res_temp;res_temp=res_temp->next)
  {
    if(res_temp->flag<-9999) continue;
    i=rotnum[res_temp->id0];
    j=0;
    res_temp0=chn->pdb->next->chn->isres(res_temp->id0);
    for(atm_temp=res_temp0->atm;atm_temp;atm_temp=atm_temp->next)
    {
      if(atm_temp->tatm->id<4) continue;
      if(atm_temp->tatm->rotate==1)
      {
        x=0;
        rotamer[res_temp->id0][i][j]=atm_temp->chi+x;
        j++;
      }
    }
    rotnum[res_temp->id0]++;
  }
}

re300:
lattice->putoff();
lattice->flag=0;
lattice->ishet=1; //modifed to account for hetatms
lattice->grdsiz=2.;
lattice->radall=15.;
lattice->ready(pdb);

for(chn=pdb->chn;chn;chn=chn->next)
for(res_temp=chn->res;res_temp;res_temp=res_temp->next)
{
if(res_temp->flag<-9999) lattice->puton(res_temp,0,1000);
else lattice->puton(res_temp,0,3);
}
//
lattice->putonhetatm();
//
float xyzv[200];
for(chn=pdb->chn;chn;chn=chn->next)
for(res_temp=chn->res;res_temp;res_temp=res_temp->next)
{
   if(res_temp->flag<-9999) continue;
   res_temp->transfer(xyzv,0);
   if(flag%10==0) crtmore(res_temp,nummore,cutbb);
   else crtmore(res_temp,10000,1000.);
   res_temp->transfer(xyzv,1);
   if(includeself) doincludeitself(res_temp);
   cerr<<"finishing assembling rotamer to residue: "<<res_temp->name<<res_temp->oid<<endl;
}

}

void Scap::doincludeitself(Res *r) {

  Res *rr;

  rr=new Res(r);

  rr->configure();
  
  if(includeself==1) { //as the first one
  	Res *rr0=r->more;
  	r->more=rr;
  	rr->more=rr0;
  }
  else if(includeself==2&&r->more) { //as the regular rotamer
	Res *rr0=r->more->more;
        r->more->more=rr;
        rr->more=rr0;
  } 
  else if(includeself==2&&r->more==0) {
        r->more=rr;
  } 

  float e;
  Res *r2=rr;
  if(proring==0) {
  if(strchr(force,'1')&&r2->name!='P'&&ring&&torsionring) e=tormax*r2->clashnoring(2,1)*1.0+clash(lattice,r2,4,100);
  else if(strchr(force,'1')&&r2->name!='P'&&ring) e=tormax*r2->clashnoring(2,0)*1.0+clash(lattice,r2,4,100);
  else if(strchr(force,'1')&&r2->name!='P') e=tormax*r2->clash(2)*1.0+clash(lattice,r2,4,100);
  else e=clash(lattice,r2,4,100);
  }
  else {
  if(strchr(force,'1')&&ring&&torsionring) e=tormax*r2->clashnoring(2,1)*1.0+clash(lattice,r2,4,100);
  else if(strchr(force,'1')&&ring) e=tormax*r2->clashnoring(2,0)*1.0+clash(lattice,r2,4,100);
  else if(strchr(force,'1')) e=tormax*r2->clash(2)*1.0+clash(lattice,r2,4,100);
  else e=clash(lattice,r2,4,100);
  }
  rr->flag=e;
  if(TRES.logg>3)cerr<<"include self..."<<rr->name<<rr->id<<" "<<rr->flag<<endl;
}

int Scap::crtmore(Res *s,int num,float cutoff)
{
  Res *res_temp0,*res_temp;
  Atm *a1,*a2,*a3,*a4;
  float x,y;
  int m,n,i;
  Rotate side;
  float xyz[100];
//find the last more residue for add on
  s->transfer(xyz,0);
  if(rott==0)
  {
    a1=(*s)[4];a3=(*s)[" CD "];a4=(*s)[" N  "];
    res_temp0=0; int nn=0;
    //for(res_temp=TRES.rotamer->next->chn->isres(s->name,0)->more;res_temp;res_temp=res_temp->more)
    for(res_temp=TRES.findrotamer("sidechain")->chn->isres(s->name,0)->more;res_temp;res_temp=res_temp->more)
    {
      if(cbeta==0)side.hook(res_temp,s,4);
      else side.hook(res_temp,s);
      
      a2=(*res_temp)[4];
      x=a1->dihedral(0);
      y=a2->dihedral(0);
      y=y-x;
      side.rotate(a1,y);
      if(TRES.logg==-1&&a3&&a4&&s->name=='P') {
	//s->dihedral(stderr);
	//res_temp->dihedral(stderr);
	cerr<<res_temp->name<<s->oid<<" "<<TRES.distance(a3,a4)<<endl;
      }
      if(a3&&a4&&s->name=='P'&&TRES.distance(a3,a4)>1.8) {
	if(TRES.logg)  cerr<<res_temp->name<<res_temp->oid<<" "<<TRES.distance(a3,a4)<<endl;
	continue;
      }
      if(TRES.logg>3&&a3&&a4&&s->name=='P') {
	cerr<<res_temp->name<<res_temp->oid<<" "<<TRES.distance(a3,a4)<<endl;
      }
      if(res_temp0==0) {
		s->more=0;
		s->more=new Res(s);
		s->more->configure();
		res_temp0=s->more;
      }
      else {
		Res *a=s->more;s->more=0;
		res_temp0->more=new Res(s);
		res_temp0->more->configure();
		res_temp0=res_temp0->more;
		s->more=a;
      }
      //res_temp0->dihedral(stdout);
      nn++;
    }
    if(layermode==0) remove(s,num,cutoff);
    else removelayermode(s,num,cutoff);
    return 1;
  }

  n=rotnum[s->id0];
  if(n==0) return 0;
  
  m=0;

  for(res_temp0=s;res_temp0->more;res_temp0=res_temp0->more);

  res_temp=s;

  for(i=0;i<n;i++)
  {
    res_temp0->more=crtmore(res_temp,i);
    res_temp0=res_temp0->more;
    m=1;
  }
  remove(res_temp,num,cutoff);
  return m;
}
 
Res *Scap::crtmore(Res *s,int n)
{
  Res *res_temp0;
  Atm *atm_temp;
  float x,y,z;   
  int k;
  Rotate side;

  if(rotamer[s->id0][n]==0) return 0;
  res_temp0=new Res(s);
  res_temp0->configure();
  k=0; 
  for(atm_temp=res_temp0->atm;atm_temp;atm_temp=atm_temp->next)
  {
    if(atm_temp->tatm->id<4) continue;
    if(atm_temp->tatm->rotate==1)
    {
      x=rotamer[s->id0][n][k];
      y=x-atm_temp->chi;
      z=fabs(y);
      if(z>0.1)side.rotate(atm_temp,y);
      atm_temp->chi=x;
      k++;
    }
  }
  return res_temp0;
}

void Scap::removefast(Res *r1,int num,float cutoff) {
Qsort cc;
Res *r2;
int m,i,u,j,temp1[20000];
float e,temp[20000];
Res *all[20000],*all1[20000];
float *all2[20000];
float take;
u=0;
Atm *c;
for(r2=r1->more;r2;r2=r2->more) {
	c=findfaratom(r2);
	if(c) e=clashfast(lattice,c);
	else  e=0;
	all1[u]=r2;
  	temp[u]=e;
  	u++;
}

cc.sort(temp,u,temp1);

  m=0;
  for(i=0;i<u;i++)
  {
    j=temp1[i];
    r2=all1[j];
    all[i]=r2;
    r2->more=0;
    r2->next=0;
    if(temp[i]>=1&&i>20)
    {
      delete r2;
      all[i]=0;
    }
  }

  j=0;
  for(i=0;i<u;i++)if(all[i]) j++;
  if(TRES.logg==-2) cerr<<"total remaining..."<<u<<" "<<j<<endl;
  u=j;

  r2=r1;
  for(i=0;i<u;i++)
  {
    r2->more=all[i];
    r2=r2->more;
    if(rott)if(TRES.logg)cerr<<r2->name<<r2->id<<" "<<i<<" "<<rotnum[r1->id0]<<" "<<r2->flag<<endl;
    else if(TRES.logg)cerr<<r2->name<<r2->id<<" "<<i<<" "<<r2->flag<<endl;
  }

  r2->more=0;
}

Atm *Scap::findfaratom(Res *r) {

Atm *a;
Atm *b;
b=r->isatmid(4);
if(b==0) return 0;
float d,d0;
Atm *c=0;
d0=0;
for(a=r->atm;a;a=a->next) {
        if(a->tatm->name[1]=='H') continue;
        if(a->tatm->id<=4) continue;
        d=TRES.distsqr(a->xyz,b->xyz);
        if(d>d0) {
                d0=d;
                c=a;
        }
}
return c;
}


void Scap::remove(Res *r1,int num,float cutoff)
{
Qsort cc;
Res *r2;
int m,i,u,j,temp1[20000];
float e,temp[20000];
Res *all[20000],*all1[20000];
float *all2[20000];
float take;

if((flag%10000)/1000>=2) take=2;
else take=2;

if(r1->flag<-9999) return;
if(fastmode==1) removefast(r1,num,cutoff);
u=0;

//calculate the energy interact with backbone
for(r2=r1->more;r2;r2=r2->more)
{
  //cerr<<r2->id<<" "<<r2->name<<" "<<u<<endl;
  if(proring==0) {
     if(strchr(force,'1')&&r2->name!='P'&&ring&&torsionring) e=tormax*r2->clashnoring(2,1)*1.0+clash(lattice,r2,4,100); 
     else if(strchr(force,'1')&&r2->name!='P'&&ring) e=tormax*r2->clashnoring(2,0)*1.0+clash(lattice,r2,4,100);
     else if(strchr(force,'1')&&r2->name!='P') e=tormax*r2->clash(2)*1.0+clash(lattice,r2,4,100);
     else e=clash(lattice,r2,4,100);
  }
  else {
     if(strchr(force,'1')&&ring&&torsionring) e=tormax*r2->clashnoring(2,1)*1.0+clash(lattice,r2,4,100);
     else if(strchr(force,'1')&&ring) e=tormax*r2->clashnoring(2,0)*1.0+clash(lattice,r2,4,100);
     else if(strchr(force,'1')) e=tormax*r2->clash(2)*1.0+clash(lattice,r2,4,100);
     else e=clash(lattice,r2,4,100);
  }
  all1[u]=r2;
  r2->flag=e;
  temp[u]=e;
  i=r1->tres->numrot-2;
  if(rott)rotamer[r1->id0][u][i]=e; 
  u++;
}

/*
if(colony>1) {
    r1->giveid();
    float re=r1->colonycoeff();
    for(r2=r1->more;r2;r2=r2->more) {
       float ent=0;
       for(Res *r3=r1->more;r3;r3=r3->more) {
              e=r2->directrmsd(r3,4,1000);
              ent+=exp(-e*e*e*re-r3->flag/8.31/0.3);
       }
       if(ent<1) ent=r2->flag;
       else ent=-8.31*0.3*log(ent);
       temp[r2->id0]=ent;
    }
}
*/


//sort the energy

  int num1;
  if(colony>1) num1=num;
  else         num1=num;

  cc.sort(temp,u,temp1);

  m=0;
  for(i=0;i<u;i++)
  {
    j=temp1[i];
    r2=all1[j];
    all[i]=r2;
    r2->more=0;
    r2->next=0;
    if((i>num1||r2->flag>cutoff)&&i>20)
    {
      if(r2->flag<take*cutoff) if(rott)all2[m++]=rotamer[r1->id0][j];
      else
      {
        if(rott){delete [] rotamer[r1->id0][j];rotamer[r1->id0][j]=0;}
      }
      
      delete r2;
      all[i]=0;
    }
    else
    {
      if(rott)
      {
        delete [] rotamer[r1->id0][j];
        rotamer[r1->id0][j]=0;
      }
    }
  } 

//recalculate the rotamer library for this residue

  for(i=0;i<m;i++)
  {
    if(rott)rotamer[r1->id0][i]=all2[i];
  }
  if(rott)rotnum[r1->id0]=m;

  j=0;
  for(i=0;i<u;i++)if(all[i]) j++;
  u=j;

  r2=r1;
  for(i=0;i<u;i++)
  {
    r2->more=all[i];
    r2=r2->more;
    if(rott)if(TRES.logg)cerr<<r2->name<<r2->id<<" "<<i<<" "<<rotnum[r1->id0]<<" "<<r2->flag<<endl;
    else if(TRES.logg)cerr<<r2->name<<r2->id<<" "<<i<<" "<<r2->flag<<endl;
  }

  r2->more=0;

  if(colony>1) {

    //calculate colony energy
    r1->giveid();
    float re=colonycoeff(r1);
    //float re=r1->colonycoeff();
    float ce=1;
    if(coleng) ce=colengcoeff(r1);
     
    u=0;
    for(r2=r1->more;r2;r2=r2->more) {
       float ent=0;
       Res *r3;
       
       for(r3=r1->more;r3;r3=r3->more) {
              e=r2->directrmsd(r3,4,1000);
	      e=pow(e,ncolony);
	      ent+=exp(-e*re-ce*r3->flag/8.31/0.3);
              //ent+=exp(-e*e*e*re-r3->flag/8.31/0.3);
       }
       if(ent<1) ent=r2->flag;
       else ent=-8.31*0.3*log(ent);
       temp[u]=ent;       
       all1[u]=r2;
       u++;
    }

    cc.sort(temp,u,temp1);
    
    m=0;
    for(i=0;i<u;i++)
    {
      j=temp1[i];
      r2=all1[j];
      all[i]=r2;
      r2->more=0;
      r2->next=0;
    }
    r2=r1;
    for(i=0;i<u;i++)
    {
      if(i>num) {
	 delete all[i]; continue;	
      }
      r2->more=all[i];
      r2=r2->more;
      //if(rott)cerr<<r2->name<<r2->id<<" "<<u<<" "<<rotnum[r1->id0]<<" "<<r2->flag<<endl;
      //else cerr<<r2->name<<r2->id<<" "<<u<<" "<<r2->flag<<endl;
    }
    r2->more=0;
  }
  //cerr<<r1->id<<r1->name<<" "<<u<<endl;
}

void Scap::removelayermode(Res *r1,int num,float cutoff)
{
Qsort cc;
Res *r2;
int m,i,u,j,temp1[20000];
float e,temp[20000];
Res *all[20000],*all1[20000];
float *all2[20000];
float take;

if((flag%10000)/1000>=2) take=2;
else take=2;

if(r1->flag<-9999) return;
if(fastmode==1) removefast(r1,num,cutoff);
u=0;

//calculate the energy interact with backbone
for(r2=r1->more;r2;r2=r2->more)
{
  //cerr<<r2->id<<" "<<r2->name<<" "<<u<<endl;
  if(proring==0) {
     if(strchr(force,'1')&&r2->name!='P'&&ring&&torsionring) e=tormax*r2->clashnoring(2,1)*1.0;//+clash(lattice,r2,4,100); 
     else if(strchr(force,'1')&&r2->name!='P'&&ring) e=tormax*r2->clashnoring(2,0)*1.0;//+clash(lattice,r2,4,100);
     else if(strchr(force,'1')&&r2->name!='P') e=tormax*r2->clash(2)*1.0;//+clash(lattice,r2,4,100);
     else e=0;//clash(lattice,r2,4,100);
  }
  else {
     if(strchr(force,'1')&&ring&&torsionring) e=tormax*r2->clashnoring(2,1)*1.0;//+clash(lattice,r2,4,100);
     else if(strchr(force,'1')&&ring) e=tormax*r2->clashnoring(2,0)*1.0;//+clash(lattice,r2,4,100);
     else if(strchr(force,'1')) e=tormax*r2->clash(2)*1.0;//+clash(lattice,r2,4,100);
     else e=0;//clash(lattice,r2,4,100);
  }
  all1[u]=r2;
  r2->flag=e;
  temp[u]=e;
  i=r1->tres->numrot-2;
  if(rott)rotamer[r1->id0][u][i]=e; 
  u++;
}

/*
if(colony>1) {
    r1->giveid();
    float re=r1->colonycoeff();
    for(r2=r1->more;r2;r2=r2->more) {
       float ent=0;
       for(Res *r3=r1->more;r3;r3=r3->more) {
              e=r2->directrmsd(r3,4,1000);
              ent+=exp(-e*e*e*re-r3->flag/8.31/0.3);
       }
       if(ent<1) ent=r2->flag;
       else ent=-8.31*0.3*log(ent);
       temp[r2->id0]=ent;
    }
}
*/


//sort the energy

  int num1;
  if(colony>1) num1=num;
  else         num1=num;

  cc.sort(temp,u,temp1);

  m=0;
  for(i=0;i<u;i++)
  {
    j=temp1[i];
    r2=all1[j];
    all[i]=r2;
    r2->more=0;
    r2->next=0;
    if((i>num1||r2->flag>cutoff)&&i>20)
    {
      if(r2->flag<take*cutoff) if(rott)all2[m++]=rotamer[r1->id0][j];
      else
      {
        if(rott){delete [] rotamer[r1->id0][j];rotamer[r1->id0][j]=0;}
      }
      
      delete r2;
      all[i]=0;
    }
    else
    {
      if(rott)
      {
        delete [] rotamer[r1->id0][j];
        rotamer[r1->id0][j]=0;
      }
    }
  } 

//recalculate the rotamer library for this residue

  for(i=0;i<m;i++)
  {
    if(rott)rotamer[r1->id0][i]=all2[i];
  }
  if(rott)rotnum[r1->id0]=m;

  j=0;
  for(i=0;i<u;i++)if(all[i]) j++;
  u=j;

  r2=r1;
  for(i=0;i<u;i++)
  {
    r2->more=all[i];
    r2=r2->more;
    if(rott)if(TRES.logg)cerr<<r2->name<<r2->id<<" "<<i<<" "<<rotnum[r1->id0]<<" "<<r2->flag<<endl;
    else if(TRES.logg)cerr<<r2->name<<r2->id<<" "<<i<<" "<<r2->flag<<endl;
  }

  r2->more=0;

  if(colony>1) {

    //calculate colony energy
    r1->giveid();
    float re=colonycoeff(r1);
    //float re=r1->colonycoeff();
    float ce=1;
    if(coleng) ce=colengcoeff(r1);
     
    u=0;
    for(r2=r1->more;r2;r2=r2->more) {
       float ent=0;
       Res *r3;
       
       for(r3=r1->more;r3;r3=r3->more) {
              e=r2->directrmsd(r3,4,1000);
	      e=pow(e,ncolony);
	      ent+=exp(-e*re-ce*r3->flag/8.31/0.3);
              //ent+=exp(-e*e*e*re-r3->flag/8.31/0.3);
       }
       if(ent<1) ent=r2->flag;
       else ent=-8.31*0.3*log(ent);
       temp[u]=ent;       
       all1[u]=r2;
       u++;
    }

    cc.sort(temp,u,temp1);
    
    m=0;
    for(i=0;i<u;i++)
    {
      j=temp1[i];
      r2=all1[j];
      all[i]=r2;
      r2->more=0;
      r2->next=0;
    }
    r2=r1;
    for(i=0;i<u;i++)
    {
      if(i>num) {
	 delete all[i]; continue;	
      }
      r2->more=all[i];
      r2=r2->more;
      //if(rott)cerr<<r2->name<<r2->id<<" "<<u<<" "<<rotnum[r1->id0]<<" "<<r2->flag<<endl;
      //else cerr<<r2->name<<r2->id<<" "<<u<<" "<<r2->flag<<endl;
    }
    r2->more=0;
  }
  //cerr<<r1->id<<r1->name<<" "<<u<<endl;
}


float Scap::colonycoeff(Res *tr) {

	if(tr==0) return 0; 

        int n=0;float d=0;
	int m=0;int m1=0;
	Res *r2,*r3;
        for(r2=tr->more;r2;r2=r2->more) {
		m++;
		if(m>300) break;
		m1=0;
                for(r3=r2->more;r3;r3=r3->more) {
                        if(r3==r2) continue;
			m1++;
			if(m1>300) break;
                        d+=r2->directrmsd(r3,4,100);
                        n++;
                }
        }

        if(n==0) return 0;
        d=d/n;
        d=log(1/coleffect)/pow(d,nncolony);
        return d;
}

float Scap::colengcoeff(Res *tr) {

        if(tr==0) return 0;

        int n=0;float d=0;
        int m=0;int m1=0;
        Res *r2,*r3;

	float re=engrmsd(tr);

        for(r2=tr->more;r2;r2=r2->more) {
                m++;
                if(m>300) break;
                m1=0;
                for(r3=r2->more;r3;r3=r3->more) {
                        if(r3==r2) continue;
                        m1++;
                        if(m1>300) break;
			float dd=fabs(r3->flag-r2->flag);
			if(dd>re) continue;
                        d+=dd;
                        n++;
                }
        }

        if(n==0) return 1;
        d=d/n;
        d=log(2)*8.31*0.298/d;
        return d;
}

float Scap::colengcoeff(Res *tr,float *entemp) {

        if(tr==0) return 0;

        int n=0;float d=0;
        int m=0;int m1=0;
        Res *r2,*r3;

        float re=engrmsd(tr,entemp);

        for(r2=tr->more;r2;r2=r2->more) {
                m++;
                if(m>300) break;
                m1=0;
                for(r3=r2->more;r3;r3=r3->more) {
                        if(r3==r2) continue;
                        m1++;
                        if(m1>300) break;
                        float dd=fabs(entemp[r3->id0]-entemp[r2->id0]);
                        if(dd>re) continue;
                        d+=dd;
                        n++;
                }
        }

        if(n==0) return 1;
        d=d/n;
        d=log(2)*8.31*0.298/d;
        return d;
}


float Scap::averageeng(Res *tr) {

	Res *r2;

	//average energy
	float e=0;
	int n=0;
	for(r2=tr->more;r2;r2=r2->more) {
		e+=r2->flag;
		n++;
	}

	if(n==0) return 1;

	e= e/n;


	//average rmsd of energy

	float re=0;

	n=0;
	for(r2=tr->more;r2;r2=r2->more) {
                re+=(r2->flag-e)*(r2->flag-e);
                n++;
        }

	re=sqrt(re/n);


	//average coeff

	float d=0;
	float av=0;
	n=0;
	for(r2=tr->more;r2;r2=r2->more) {
		if(r2->flag-e>re) continue;

		d=fabs(r2->flag-e);
		
		av+=d;
		n++;
	}

	if(n==0) return 1;

	av=av/n;
	if(av==0) return 1;

	return log(2)*8.31*0.298/av;
}

float Scap::engrmsd(Res *tr) {

        Res *r2;

        //average energy
        float e=0;
        int n=0;
        for(r2=tr->more;r2;r2=r2->more) {
                e+=r2->flag;
                n++;
        }

        if(n==0) return 1;

        e= e/n;


        //average rmsd of energy

        float re=0;

        n=0;
        for(r2=tr->more;r2;r2=r2->more) {
                re+=(r2->flag-e)*(r2->flag-e);
                n++;
        }

        re=sqrt(re/n);


	return re;
}

float Scap::engrmsd(Res *tr,float *entemp) {

        Res *r2;

        //average energy

        float e=0;
        int n=0;
        for(r2=tr->more;r2;r2=r2->more) {
                e+=entemp[r2->id0];
                n++;
        }

        if(n==0) return 1;

        e= e/n;


        //average rmsd of energy

        float re=0;

        n=0;
        for(r2=tr->more;r2;r2=r2->more) {
                re+=(entemp[r2->id0]-e)*(entemp[r2->id0]-e);
                n++;
        }

        re=sqrt(re/n);


        return re;
}


float Scap::clash(Lattice *lattice,Res *r,int begin,int end)
{
float e;
Atm *a,*b,*side[20];
int t,i;

t=0;

for(a=r->atm;a;a=a->next)
{
  if(a->tatm->name[1]=='H') 
  {
    if(strchr(force,'4')==0) break;
    if(fastmodeH==1) continue; //no H-H atom interaction.
    b=a->bond[0];
    if(b->tatm->id>=begin&&b->tatm->id<=end)side[t++]=a;
  }
  else if(a->tatm->id>=begin&&a->tatm->id<=end)side[t++]=a;
}

e=0;

for(i=t-1;i>=0;i--)
{
  e+=clash(lattice,side[i]);
  if(e>99999) return e;
}
return e;
}

float Scap::clashfast(Lattice *lattice,Atm *s) {

if(s->res->name=='C'||s->res->name=='P') return 0;
if(s->tatm->name[1]=='H') return 0;

lattice->getcell(s,1.0);
cerr<<s->name<<" "<<s->res->name<<" "<<lattice->nget<<endl;
float e=0;
int i,j;
float d,r;
r=0;
for(i=lattice->nget-1;i>=0;i--) {
	Atm *a=lattice->obtain[i]->atm;
	if(a==0) continue;
	if(a->res==s->res) continue; //if(a->res->id==s->res->id) continue;//same residue
	if(a->tatm->name[1]=='H') r=s->tatm->eng->radius;
	else r=a->tatm->eng->radius+s->tatm->eng->radius;
	j=s->tatm->hbond*a->tatm->hbond;
	if(j) r=3;
	d=TRES.distsqr(a->xyz,s->xyz);
	if(d<0.25*r*r) e+=1;
}
return e;
}
float Scap::clash(Lattice *lattice,Atm *s)
{
Atm *a;
int i,j;
float w,e,d,r,c;
float dis2;
e=0;
w=4*cut*cut;
lattice->getcell(s,cut);
lattice->cutoff(s,cut);
if(strchr(force,'a')) lattice->getrescutcharge(s);
else if(strchr(force,'A')) lattice->getrescut(s);
//lattice->getrescutcharge(s);
for(i=lattice->nget-1;i>=0;i--)
{
 a=lattice->obtain[i]->atm; 
 if(a->res==s->res) continue;//if(a->res->id==s->res->id) continue;//same residue
 
 if(s->res->name=='P'&&a->res->id-s->res->id==-1&&a->res->chn==s->res->chn) {//same chain of P
    if(a->tatm->id==1||a->tatm->id==2||a->tatm->id==3) continue;
 }
 else if(s->res->name=='P'&&a->res->chn->ishet==1) {//in case of unknown chained residues
    if(strcmp(s->name," CD ")==0) continue;
 }

 d=TRES.distsqr(a->xyz,s->xyz);
 dis2=d;

 if(emilcharge) {
	if((a->mask&TRES.constant->surfacecalculated)==0) {
		if(a->res->chn->ishet==0) latsurfv->calcarea(a);
		else   a->area=12.6*(probe+a->tatm->eng->radius)*(probe+a->tatm->eng->radius);                 
	}
	if((s->mask&TRES.constant->surfacecalculated)==0) {
                latsurfv->calcarea(s);
        }
 }

 if(emilcharge&&emilnewcharge&&lattice==this->lattice0&&strchr(force,'2')) {
        if(emilnewcharge==1) {
		float probe=latsurfv->probe_radius;
        	float g2=a->area/12.6/(probe+a->tatm->eng->radius)/(a->tatm->eng->radius+probe);
        	if(g2>1) g2=1;
        	float g1=s->area/12.6/(probe+s->tatm->eng->radius)/(s->tatm->eng->radius+probe);
        	if(g1>1) g1=1;
        	float cr2=a->tatm->eng->charge;
        	float cr1=s->tatm->eng->charge;
        	float rd=sqrt(dis2);
        	if(rd<3.5) rd=3.5;
        	//float aa=cr2*cr1*40*0.59/rd*(1.25-g1)*(1.25-g2);
        	//aa+=-cr1*cr1*40*0.59/(s->tatm->eng->radius+probe)*g1*g1;//(1.25-g1)*(1.25-g1);
		float aa=332.3*cr2*cr1/80/rd*sqrt((1+19*(1-g1))*(1+19*(1-g2)));
		aa+=-cr1*cr1*0.5*332.3/4/(s->tatm->eng->radius+probe)*g1;
        	e+=aa;
		continue;
	}
 }
 //if(TRES.logg>3) { 
 // 	if(d<2.0) cerr<<a->res->name<<a->res->id0<<" "<<a->name<<":"<<s->res->name<<s->id0<<" "<<s->name<<" "<<d<<endl;
 //}
 if(d<0.1) return 9999999;
 else if(d>w&&strchr(force,'2')==0) continue;  
 r=a->tatm->eng->radius+s->tatm->eng->radius;
 j=s->tatm->hbond*a->tatm->hbond;


 if(j==2||j==3||j==6)
 {
  r=3.;
  if(strchr(force,'3')&&dis2<25) {
	float ee=s->ishbond(a,1.5,40,hbondeng);  
	//cerr<<ee<<" "<<sqrt(d)<<endl;
	e+=ee;
  }
 }
 else 
 {
   if(a->tatm->name[1]=='H'&&s->tatm->name[1]=='H'&&r>1.){
	if(adaptvdw==110) {
		r=r/2;
	}
	else if(adaptvdw>=100) {
           r=r; 
           //cerr<<a->name<<" "<<s->name<<r<<endl; 
        }
	else {
	      r=1.;
        }        
   }
 }

 if(a->res->name=='C'&&s->res->name=='C')
 {
    if(a->tatm->id==5&&s->tatm->id==5) {r=2.0; e+=-10;}
    else if(a->tatm->id==5&&s->tatm->id==4) r=3.;
    else if(a->tatm->id==4&&s->tatm->id==5) r=3.;
 }
 
 if(adaptvdw==101) {
	r=r*1.10;
 }
 else if(adaptvdw==102) {
	if(a->tatm->name[1]=='H'&&s->tatm->name[1]=='H'&&r>1.) r=r*1.1;
	else r=r*0.9;
 }
 
 r=r*r;
 c=sqrt(s->tatm->eng->epslon*a->tatm->eng->epslon);
 r=r/d;

 //detect outenforce: expanding surface vdw
 float edone=1;
 if(outenforce==10) edone=sqrt(r);
 //end

 r=r*r*r; 
 float ed=c*(r*r-2*r);

 if(outenforce==10) {
    if(outforcemag==0) {
    	if(solventbury[s->res->id0]<0.5) ed=c*(r*r*edone*edone-2*r*edone);
    }
    else if(outforcemag==1) {
	if(solventbury[s->res->id0]>=0.9) ed=c*(r*r/edone/edone-2*r/edone);
	else if(solventbury[s->res->id0]<0.5) ed=c*(r*r*edone*edone-2*r*edone);
    }
    else if(outforcemag==2) {
	if(solventbury[s->res->id0]<0.5) ed=c*(r*r/edone/edone-2*r/edone);
    }
    else if(outforcemag==3) {
        if(solventbury[s->res->id0]<0.5) ed=c*(r*r/edone/edone/edone/edone-2*r/edone/edone);
    }
    else if(outforcemag==4) {
        if(solventbury[s->res->id0]<0.5) ed=c*(sqrt(r*r)-2*sqrt(r));
    }
 } 

 if(iweight) {
	if(solventbury&&iterateweight<1) {
		ed=ed*iterateweight*solventbury[a->res->id0];
	}
 }

 if(adaptvdw==1) {
    int mr=0;
    j=a->tatm->hbond*s->tatm->hbond;  
    if(j==2||j==3||j==6||j==9) mr=1;
    int isp=a->tatm->ispolar*s->tatm->ispolar; 
    if(ed<0) {
	if(mr==1)       ed= ed*1.50;
	else if(isp==1) ed= ed*1.50;
	else if(isp==2) ed= ed*1.25;
	else if(isp==3) ed= ed*0.75;
	else            ed= ed;
    }
    else {
	if(mr==1)       ed= ed*.50;
        else if(isp==1) ed= ed*.50;
        else if(isp==2) ed= ed*0.75;
        else if(isp==3) ed= ed*1.25;
	else            ed= ed;
    } 
 }
 else if(adaptvdw==2) {
    int mr=0;
    j=a->tatm->hbond*s->tatm->hbond;
    if(j==2||j==3||j==6||j==9) mr=1;
    int isp=a->tatm->ispolar*s->tatm->ispolar;
    if(ed<0) {
        if(mr==1)       ed= ed*1.250;
        else if(isp==1) ed= ed*1.250;
        else if(isp==2) ed= ed*1.125;
        else if(isp==3) ed= ed*0.875;
	else            ed= ed;
    }
    else {
        if(mr==1)       ed= ed*.250;
        else if(isp==1) ed= ed*.250;
        else if(isp==2) ed= ed*0.875;
        else if(isp==3) ed= ed*1.125;
	else            ed= ed;
    }
 }

 e+=ed;

 //e+=c*(r*r-2*r);

 if(emilnewcharge&&emilcharge&&strchr(force,'2')) {
	
	if(emilnewcharge==1) {
                float probe=latsurfv->probe_radius;
                float g2=a->area/12.6/(probe+a->tatm->eng->radius)/(a->tatm->eng->radius+probe);
                if(g2>1) g2=1;
                float g1=s->area/12.6/(probe+s->tatm->eng->radius)/(s->tatm->eng->radius+probe);
                if(g1>1) g1=1;
                float cr2=a->tatm->eng->charge;
                float cr1=s->tatm->eng->charge;
                float rd=sqrt(dis2);
                if(rd<3.5) rd=3.5;
                //float aa=cr2*cr1*40*0.59/rd*(1.25-g1)*(1.25-g2);
                //aa+=-cr1*cr1*40*0.59/(s->tatm->eng->radius+probe)*g1*g1;//(1.25-g1)*(1.25-g1);
                float aa=332.3*cr2*cr1/80/rd*sqrt((1+19*(1-g1))*(1+19*(1-g2)));
                aa+=-cr1*cr1*0.5*332.3/4/(s->tatm->eng->radius+probe)*g1;
                e+=aa;
                continue;
        }
 } 

 if(strchr(force,'2')&&s->tatm->eng->charge!=0&&a->tatm->eng->charge!=0&&emilnewcharge==0&&emilcharge==0) {
	if(dsqrt==1) {
		d=sqrt(d);
		if(d<3.5) d=3.5;
	}
	else if(dsqrt==2) {
		d=sqrt(d);
                if(d<3.5) d=d;
		else      d=d*d;
	}
	else {
  		if(d<3.5*3.5) d=3.5*3.5;
	}
	//ed=s->tatm->eng->charge*a->tatm->eng->charge/d*166.15*2;
	ed=s->tatm->eng->charge*a->tatm->eng->charge/d*166.15*2/dielectric;
	if(TRES.logg>3) cerr<<ed<<endl;
	if(iweight) {
        	if(solventbury&&iterateweight<1) {
			ed=ed*iterateweight*solventbury[a->res->id0];
		}
 	}
        if(outcharge==10) {
		ed=s->tatm->eng->charge*a->tatm->eng->charge/d*166.15*2/dielectric*solventbury[s->res->id0]*solventbury[a->res->id0];
        }
	//if(ed>100) cerr<<ed<<endl;
	e+=ed;
  	//e+=s->tatm->eng->charge*a->tatm->eng->charge/d*166.15*2;
	//cerr<<s->tatm->name<<" "<<a->tatm->name<<" "<<e<<endl;
 }
 
 
}
//cerr<<e<<endl;
return e;
}


void Scap::setsolventorder() {
	
	 
	int i,j;
	Chn *ch;
	Res *r;
	
	if(solventorder) delete [] solventorder;
	solventorder=0;

	int n=0;
	for(ch=pdb->chn;ch;ch=ch->next)
	for(r=ch->res;r;r=r->next) n++;

	if(n==0) return ;

	solventorder=new Res*[n+100];
	for(i=0;i<n+100;i++) solventorder[i]=0;
	solventbury=new float[n+1000];
	for(i=0;i<n+1000;i++) solventbury[i]=0;
	float *tmp=new float[n+100];
	Res **tt=new Res*[n+100]; 
	int *order=new int[n+100];

	float dd=0;
	n=0;
	for(ch=pdb->chn;ch;ch=ch->next)
	for(r=ch->res;r;r=r->next) {
		dd=r->bury(4,100);
		dd=fabs(dd);
		solventbury[r->id0]=dd;
		if(solventorderrandom==0) {
			tmp[n]=-dd;
 		}
		else if(solventorderrandom==1) {
			tmp[n]=dd;
		}
		else if(solventorderrandom==2) {
			tmp[n]=random()%1000;
		}
		//make it random
		//tmp[n]=random()&100000;
		//
		order[n]=n;
		tt[n]=r;
		n++;
	}

	Qsort cc;
	cc.sort(tmp,n,order);
	for(i=0;i<n;i++) {
		j=order[i];
		solventorder[i]=tt[j];
		if(tt[j]) {
			if(TRES.logg)cerr<<"residue: "<<solventorder[i]->name<<solventorder[i]->id<<" percentage buried: "<<solventbury[tt[j]->id0]<<endl;
		}
	}
	delete [] tmp;
	delete [] tt;
	delete [] order;
}

void Scap::setsolventorderanyway() {
	
	 
	int i;
	Chn *ch;
	Res *r;
	
	if(solventorder) delete [] solventorder;
	solventorder=0;

	int n=0;
	for(ch=pdb->chn;ch;ch=ch->next)
	for(r=ch->res;r;r=r->next) n++;

	if(n==0) return;

	solventorder=new Res*[n+100];
	solventbury=new float[n+1000];
	for(i=0;i<n+100;i++) {
		solventorder[i]=0;
		solventbury[i]=0;
	}
	 

	n=0;
	for(ch=pdb->chn;ch;ch=ch->next)
	for(r=ch->res;r;r=r->next) {
		 
		solventorder[n]=r;
		solventbury[r->id0]=fabs(r->bury(4,100));
		//cerr<<n<<" "<<r->id0<<" "<<r->name<<" "<<solventbury[r->id0]<<endl;
		n++;
	}
}

float Scap::scpred(int rannum){

	return scpred("output",0,rannum);
}

float Scap::scpred(char *output,char *filnam,int rannum)
{
Strhandler cc;
FILE *fp;
Res *r,*r1,*r2;
float d,e,f,*xyz;
int n,i,j,k,flg,m;
char line[100],*poit;
float ee[1000],*ener0,*flagt;
int   eeran[2000];
Chn *chn;
char **outname=0;
outname=new char*[rannum*2+1000];
int ir=0;
for(ir=0;ir<rannum*2+1000;ir++) outname[ir]=0;

if(pdb==0) return 0;

/*
if(pdborg) delete pdborg;pdborg=0;
pdborg=new Pdb(pdb);
pdborg->configure();
*/

if(nncolony<0.0001&&ncolony>0.001) nncolony=ncolony;

cerr<<"calculate solvent accessible surface...."<<endl;

cerr<<"building multiple rotamers for each side-chain..."<<endl;

pdb->dihedral(0);

if(delphi) {
char nn[1000];
sprintf(nn,"%s.siz",delfile);
TRES.writedelphisize(nn);
}

if(emilcharge) {
	if(latsurfv==0) {
		latsurfv=new LatSurfv();
		latsurfv->ready(pdb,2,1.4);
	}	
}

int en=emilcharge;
emilcharge=0;
sidelib(filnam);
emilcharge=en;
if(disulfide) adddisulfide();
n=pdb->getresnum();
ener=new float[n+100];
ener0=new float[n+100];
flagt=new float[n+100];
i=100;
for(chn=pdb->chn;chn;chn=chn->next)
for(r=chn->res;r;r=r->next)
{
 i+=r->tres->number*3;
 ener[r->id0]=10000;
}
xyz=new float[i];

pdb->transfer(xyz,0);
e=100000;
 
//calculating the initial energy of the system
lattice0->putoff();
lattice0->flag=0;
lattice0->grdsiz=2.;
lattice0->radall=15.;
lattice0->ready(pdb);
lattice0->puton(pdb,0,1000); 
e=0;
for(chn=pdb->chn;chn;chn=chn->next)
for(r=chn->res;r;r=r->next)
{
 if(r->flag<-9999)continue;
 if(proring==0) {
 if(strchr(force,'1')&&r->name!='P'&&ring&&torsionring)e+=tormax*r->clashnoring(2,1);
 else if(strchr(force,'1')&&r->name!='P'&&ring)e+=tormax*r->clashnoring(2,0);
 else if(strchr(force,'1')&&r->name!='P')e+=tormax*r->clash(2);
 }
 else {
 if(strchr(force,'1')&&ring&&torsionring)e+=tormax*r->clashnoring(2,1);
 else if(strchr(force,'1')&&ring)e+=tormax*r->clashnoring(2,0);
 else if(strchr(force,'1'))e+=tormax*r->clash(2);
 }
 e+=clash(lattice0,r,4,100);
}
if(TRES.logg>2) cerr<<"the total initial energy of "<<output<<" is:"<<e<<endl;

if(rannum<=0)
{
 for(chn=pdb->chn;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next)
 {
    if(r->flag<-9999)continue;
    d=10000;
    r2=chn->pdb->next->chn->isres(r->id);
    for(r1=r->more;r1;r1=r1->more)
    {
      n=(int)r2->rmsd(r1,5,1000)/1000000;
      e=r2->rmsd(r1,5,1000)-n*1000000;
      if(e<d) 
      {
        d=e;
        r1->transfer(xyz,0);
        f=r1->flag;      
      }
    }
    
    for(i=0;i<rotnum[r->id0];i++)
    {
      j=r->tres->numrot-2;
      r->dihedral(0);
      r1=crtmore(r,i);r1->next=0;r1->more=0;
      r1->flag=rotamer[r->id0][i][j];
      if(r1==0) continue;
      n=(int)r2->rmsd(r1,5,1000)/1000000;
      e=r2->rmsd(r1,5,1000)-n*1000000;
      if(e<d)
      {
        d=e;
        r1->transfer(xyz,0);
        f=r1->flag;
      }
      delete r1;
    }
    r->transfer(xyz,1);
    r->flag=f;
  }
}

flg=0;
if(rannum==0)
{
 flg=1;
 rannum=1;
}
else if(rannum==-1) {write(); chn->transfer(xyz,0);}

lattice0->putoff();
lattice0->flag=0;
lattice0->grdsiz=2.;
lattice0->radall=15.;
lattice0->ready(pdb);
//lattice0->puton(chn,0,3);

for(chn=pdb->chn;chn;chn=chn->next)
for(r=chn->res;r;r=r->next)
{
if(r->flag<-9999) lattice0->puton(r,0,1000);
else lattice0->puton(r,0,3);
}

cerr<<"building random initial conformation..."<<endl;

poit=strstr(output,".pdb");
if(poit) *poit='\0';
sprintf(line,"%s_random",output);
//fpp=fopen(line,"w");
/*
for(chn=pdb->chn;chn;chn=chn->next)
for(r=chn->res;r;r=r->next)
{
 if(r->flag<-9999){continue;}//r->write(fpp);continue;}
 for(r1=r->more;r1;r1=r1->more) k++;
 if(k==0) k=1;
 r1=r->ismore(k-1);
 if(r1==0) r1=r;
 r1->chn=chn;
 //r1->write(fpp);
}
*/
//fflush(fpp);

e=1000000;

for(chn=pdb->chn;chn;chn=chn->next)
for(r=chn->res;r;r=r->next)
ener0[r->id0]=100000;

averagermsd();
pdb->giveresmoreid();
int totout=0;ir=0;
for(i=0;i<rannum;i++)
{
 j=random();
 if(i==0) j=0;
 if(rannum>1000) j=rannum;
 for(chn=pdb->chn;chn;chn=chn->next)
 for(r=chn->res;r&&flg==0;r=r->next)
 {
   if(r->flag<-9999)continue;
   if(i>arbt)
   {
     if(ener0[r->id0]<0) m=(j+r->name+r->id0)%4;
     else if(ener0[r->id0]<-3) m=(j+r->name+r->id0)%3;
     else m=(j+r->name+r->id0)%2; 
   }
   else m=0;
   
   if(m) continue;

   k=0;
   for(r1=r->more;r1;r1=r1->more) 
   { 
    if(k<5) k++;
    else if(r1->flag<0) k++;
   }
   if(k==0) k=1;
   n=j%k;  
   if(includeself==1) {
	if(i==1) {
        	if(r->more&&r->more->more) r1=r->more->more;
        	else r1=findminimalmore(r);
   	}
   	else if(colony>1&&i==2) {
        	r1=findminimaltorsionmore(r);
   	}
   	else if(singletorsion==1&&i==3) {
        	r1=findminimaltorsionmore(r);
   	}
   	else if(bordertorsion==1&&i==4) {
        	r1=findminimalbordermore(r);
   	}
   	else r1=r->ismore(n);
   }
   else {
	if(colony>1&&i==1) {
        	//r1=findminimaltorsionmore(r);
		if(r->more&&r->more->more) r1=r->more->more;
                else r1=findminimalmore(r);
   	}
   	else if(singletorsion==1&&i==2) {
        	r1=findminimaltorsionmore(r);
   	}
   	else if(bordertorsion==1&&i==3) {
        	r1=findminimalbordermore(r);
   	}
   	else r1=r->ismore(n);
   }
   if(r1==0) {r1=r->more;if(r1==0) continue;}
   r->transfer(r1);
   r->flag=r1->flag;
 }
 seqn=i;
 if(i==0) {
  	if(solventeff) {	
		TRES.surface(10,1.4,-1);
		for(chn=pdb->chn;chn;chn=chn->next)chn->surface(10,1.4);
		setsolventorder();
  	}
  	else if(solventcolony) {  
		TRES.surface(10,1.4,-1);
		for(chn=pdb->chn;chn;chn=chn->next)chn->surface(10,1.4);
		setsolventorderanyway();
  	}
	else if(outenforce) {
		TRES.surface(10,1.4,-1);
                for(chn=pdb->chn;chn;chn=chn->next)chn->surface(10,1.4);
                setsolventorderanyway();       
		outenforce=10;
	}
        else if(outcharge) {
                TRES.surface(10,1.4,-1);
                for(chn=pdb->chn;chn;chn=chn->next)chn->surface(10,1.4);
                setsolventorderanyway();
		outcharge=10;
        }
  	else {
		setsolventorderanyway();
  	}
 }

 if(i<outinitial) {
    char namee[100];
    sprintf(namee,"%s_initial_%i.pdb",output,i);
    pdb->write(namee);
 }
 cerr<<endl;
 cerr<<"sidechain prediciton on "<<i<<"th candidates.."<<endl;
 cerr<<"this may take about 30min depending on the rotamer library and protein size.."<<endl;
 cerr<<"please be patient..."<<endl;
 cerr<<endl; 
 if(i==rannum-1) d=iterate(1);
 else            d=iterate(0);

/*
 for(r=chn->res;r;r=r->next) 
 {
   if(r->flag>10||r->flag<-9999) continue;
   r->dihedral(0);
   insert(r);
 }
*/ 

 ee[i]=d;
 if(TRES.logg>3) cerr<<i<<"...."<<d<<endl;
 if(nout&&resultout) {
    int ii=0;
    int jj=0;
    for(ii=0;ii<i;ii++){
        if(ee[ii]==ee[i]) {
                jj=1;
                break;
        }
    }
    if(jj==0) {
        char nf[1000];

        sprintf(nf,"%s_scap.%i.pdb",output,totout);
	outname[ir]=strdup(nf);
	//cerr<<"write down scap prediction for diffeent initial conformations:"<<nf<<endl;
        pdb->write(nf);
	eeran[ir]=ee[i];
	ir++;
        totout++;
    }
 }
 if(d<e) 
 {
   e=d;mine=d;
   pdb->transfer(xyz,0);
   for(chn=pdb->chn;chn;chn=chn->next)
   for(r=chn->res;r;r=r->next)
   {flagt[r->id0]=r->flag;ener0[r->id0]=ener[r->id0];}
 }
 if(TRES.logg) cerr<<"calculation..."<<i<<"..."<<j<<"minimized energy:"<<d<<endl;
 //write();
 pdb->transfer(xyz,1);
 for(chn=pdb->chn;chn;chn=chn->next)
 for(r=chn->res;r;r=r->next) r->flag=flagt[r->id0];
 if(rannum>1000) break;
}

pdb->transfer(xyz,1);
if(rannum>1&&rannum<=1000)
{
cout<<"minimized energy:"<<e<<endl; 
//write();
}

if(resultout==0) goto reout;
sprintf(line,"%s_scap.pdb",output);
//cerr<<"write down the final scap prediction: "<<line<<endl;
fp=fopen(line,"w");
pdb->setatmoid();
pdb->write(fp);
char *lineout;lineout=strdup(line);
//pdb->writeold(fp);
for(j=0;j<i;j++)
{
//fprintf(fp,"the total energy:%d,%f\n",eeran[j],ee[j]);
}
float eout;eout=e;
//fprintf(fp,"the total sidechain energy:%f\n",e);
fflush(fp);
fclose(fp);

int ire;ire=0;
cerr<<endl;
while(outname[ire]) {
cerr<<"output intermediate "<<ire<<"th structure..."<<outname[ire]<<" energy:"<<eeran[ire]<<endl;
outname[ire]=cc.strdel(outname[ire]);
ire++;
}
if(outname) delete [] outname;
cerr<<endl;
if(nout==2) {
sprintf(line,"%s_scap_rtm.pdb",output);
cerr<<"write out scap out with all rotamers..."<<line<<endl;
pdb->setallatmid0();
pdb->writeoldmore(line);
}
cerr<<endl;
cerr<<"write down the final scap prediction: "<<lineout<<endl;
cerr<<endl;
lineout=cc.strdel(lineout);
printf("the total sidechain energy:%f\n",eout);
if(eout>100) 
{
cerr<<"there is unresolved clashes"<<endl;
cerr<<"you may remove clashes by using larger rotamer library and larger number of initial conformations.."<<endl; 
}
reout:
for(i=0;i<pdb->getresnum();i++)
{
 for(j=0;rotnum&&j<rotnum[i];j++)
 {
   if(rotamer[i][j]) {delete [] rotamer[i][j];rotamer[i][j]=0;}
 }
 if(rotamer[i]) delete [] rotamer[i];rotamer[i]=0;
}
if(rotamer) delete [] rotamer; rotamer=0;
if(rotnum) delete [] rotnum;rotnum=0;
if(xyz) delete [] xyz;xyz=0;
if(ener) delete [] ener;ener=0;
if(ener0) delete [] ener0;ener0=0;
if(flagt) delete [] flagt;flagt=0;
if(armsd) delete [] armsd;armsd=0;
return e;
}

float Scap::mutate()
{
singletorsion=1;
colonyline=1;
colony=2;  
ncolony=1;
nncolony=1;
nummore=100;
bmax=2;
tormax=2;
ring=1;
arbt=0;
includeself=0;
resultout=0;

float e=0;
sidelib(0);
for(Chn *c=pdb->chn;c;c=c->next)
for(Res *r=c->res;r;r=r->next) {
	if(r->flag<=-9999) continue;

	if(r->more==0) continue;
	float d=99999;
	/*
	for(Res *rr=r->more;rr;rr=rr->more)  {
		if(rr==r->more) d=r->more->flag;
		if(rr->flag<d) d=rr->flag;
	}
	*/	
	d=r->more->flag;
	r->transfer(r->more,1);
	e+=d;
}
return e;
}



float Scap::calcdelphi(Res *r) {

	float totc=0;
	int nc=0;
	for(Tatm *ta=r->tres->tatm;ta;ta=ta->next) {
		totc+=ta->eng->charge;
		if(fabs(ta->eng->charge)>0.01) nc++;
	}

	if(nc==0) return 0;

	if(fabs(totc)<0.01&&allcharges==0) return 0;

	
	if(TRES.logg)cerr<<r->id0<<r->name<<"  test..."<<endl;
	pdb->writedelphi(delfile,r,r->id0+1,4,100);
	writedelphiparm();
	char line[200];
	sprintf(line,"./delphi %s.prm|grep energy:>%s.tmp",delfile,delfile);
	system(line);
	
	float tot=0; 
	
        FILE *fp=0;
        if(delphiself) {
	sprintf(line,"%s.tmp",delfile);
	fp=fopen(line,"r");
	int io=0;
	while(fgets(line,200,fp)!=NULL) {
		char *s=strstr(line,"corrected reaction field energy:");
		if(s) {
			tot+=atof(s+strlen("corrected reaction field energy:"));		
			io++;
			break;
		}
	}
	if(TRES.logg)cerr<<"the reaction field energy is: "<<tot*0.59<<endl;
	if(io==0)cerr<<"the reaction field energy is: "<<tot*0.59<<endl;
	fclose(fp);
        }
	
	int nn=pdb->manyatm();

	Atm **myatm=new Atm*[nn+100];

	int i=0;
	for(i=0;i<nn+100;i++) myatm[i]=0; 
	nn=0;
	Chn *c;
	
	for(c=pdb->chn;c;c=c->next)
	for(Res *r=c->res;r;r=r->next)
	for(Atm *a=r->atm;a;a=a->next) {
		myatm[nn++]=a;
	}	
 	if(TRES.logg==-10) cerr<<nn<<endl;
	
	for(c=pdb->hetatm;c;c=c->next) 
	for(Res *r=c->res;r;r=r->next) 
	for(Res *r0=r;r0;r0=r0->more)
	for(Atm *atm_temp=r0->atm;atm_temp;atm_temp=atm_temp->next){
		Atm *s=r->isatm(atm_temp->name);
		if(s==0) continue;
		myatm[nn++]=atm_temp;
  	}

 	if(TRES.logg==-10) cerr<<nn<<endl;

	sprintf(line,"%s.out",delfile);	
	fp=fopen(line,"r");
	i=0;nn=0;

	float x,y,z,ci,p,t;

        while(fgets(line,200,fp)!=NULL) {
                char *s=strstr(line,"ATOM COORDINATES");
                if(s&&i==0) {
                	i=1;
			continue;        
                }
		else if(s==0&&i==0) {
			continue;
		}
		else {
			sscanf(line,"%f %f %f %f %f",&x,&y,&z,&ci,&p);
			t=myatm[nn]->xyz[0]-x+myatm[nn]->xyz[1]-y+myatm[nn]->xyz[2]-z;
			if(fabs(t)>5) {
				if(TRES.logg)cerr<<"the atom does not match "<<myatm[nn]->name<<"  "<<myatm[nn]->id0<<myatm[nn]->name<<"  "<<t<<endl;
				nn++;
				continue;
			}	
			if(myatm[nn]->res!=r) {
				tot+=myatm[nn]->tatm->eng->charge*p;
			}
			else if(myatm[nn]->res==r) {
				if(delphionself) {
					if(strcmp(myatm[nn]->tatm->name," HN ")==0) {
						tot+=myatm[nn]->tatm->eng->charge*p;
					}
					else if(myatm[nn]->tatm->id==0) {
						tot+=(myatm[nn]->tatm->eng->charge+myatm[nn]->tatm->next->eng->charge)*p;
					}
					else if(myatm[nn]->tatm->id==3||myatm[nn]->tatm->id==2) {
						tot+=myatm[nn]->tatm->eng->charge*p;
					}
				}
			}
			nn++;			
			if(myatm[nn]==0) break;
		}
        }
 	if(TRES.logg==-10) cerr<<nn<<endl;
        fclose(fp);
	delete [] myatm;

	return tot*0.59;
}

void Scap::writedelphiparm() {
	
	char line[100];
	sprintf(line,"%s.prm",delfile);

	FILE *fp=fopen(line,"w");

	fprintf(fp,"scale=%f\n",delphiscale);
	//fprintf(fp,"gsize=%d\n",delphigsize);
	fprintf(fp,"linit=%d\n",delphilinit);
	fprintf(fp,"perfil=%f\n",perfil);
	fprintf(fp,"indi=%f\n",indi);
	fprintf(fp,"exdi=%f\n",exdi);
	fprintf(fp,"prbrad=%f\n",prbrad);
	fprintf(fp,"nonit=0.000000\n");
	fprintf(fp,"in(pdb,file=\"%s.pdb\")\n",delfile);
	fprintf(fp,"in(crg,file=\"%s.crg\")\n",delfile);
	fprintf(fp,"in(siz,file=\"%s.siz\")\n",delfile);
	fprintf(fp,"in(frc,file=\"self\")\n");
	fprintf(fp,"energy(s)\n");
	fprintf(fp,"out(frc,file=\"%s.out\")\n",delfile);
	fclose(fp);
}


float Scap::iterate(int fff)
{
Res *r1,*r2,*r3,*r4;
float xyz[200],d,e,tot,pre,cflg,tot0;
int i,n,j,rotm,rotm1;
int *comb,mm[5],rotm2;
int step;
step=1;
comb=combination(1);
for(i=0;i<5;i++)mm[i]=step;
for(i=0;i<5;i++) for(j=0;j<5;j++) if(j<=i) mm[4-j]*=2*step+1;

lattice->putoff();
lattice->flag=0;
lattice->ishet=0;//modified to take off hetatms
lattice->grdsiz=2.;
lattice->radall=15.;
lattice->ready(pdb);

if(emilcharge) {
lattice1->putoff();
lattice1->flag=0;
lattice1->grdsiz=2.;
lattice1->radall=15.;
lattice1->ready(pdb);
}

Chn *chn;
for(chn=pdb->chn;chn;chn=chn->next)
for(r1=chn->res;r1;r1=r1->next)
{
if(r1->flag>-9999) lattice->puton(r1,4,1000);
if(emilcharge) {
	lattice1->puton(r1,0,3);
}
}


for(chn=pdb->chn;chn;chn=chn->next)
for(r1=chn->res;r1;r1=r1->next) ener[r1->id0]=1000;

float enetemp[2000];
float ent,ent0,pre0;
//int ite;

n=0;rotm1=1;rotm2=1;tot0=-1000000;
int nite=0;
while(solventorder[nite])nite++;
int itt=0;

float areaorg=0;
if(hydrophobic==1) {
  pdb->setflg(0);
}
else if(hydrophobic==2) {
  pdb->setflg(1);
  pdb->chn->surface(10,1.4);
  areaorg=pdb->chn->area; 
}
else if(hydrophobic==3) {
  /*
  pdb->setflg(1);
  pdb->chn->surface(10,1.4);
  areaorg=pdb->chn->area;
  */
  areaorg=calcularea(pdb->chn);
  if(TRES.logg)cerr<<"the total inital area is..."<<areaorg<<endl;
}
else if(hydrophobic==4) {
  //pdb->setflg(1);
  //pdb->chn->surface(10,0);
  //areaorg=pdb->chn->area;
  areaorg=calcularea(pdb->chn,0);
  if(TRES.logg)cerr<<"the total inital area is..."<<areaorg<<endl;
}
else if(hydrophobic==5) {
  //pdb->setflg(1);
  //pdb->chn->surface(10,0);
  //areaorg=pdb->chn->area;
  //areaorg=calcularea(pdb->chn,0);
  //cerr<<"the total inital area is..."<<areaorg<<endl;
  //TRES.surface(10,1.4);
  TRES.surface(10,1.4,-1);
  pdb->setflg(1);
  pdb->chn->surface(10,1.4);
  areaorg=pdb->chn->area;
}
else if(hydrophobic==6||hydrophobic==7) {
  //pdb->setflg(1);
  //pdb->chn->surface(10,0);
  //areaorg=pdb->chn->area;
  //areaorg=calcularea(pdb->chn,0);
  //cerr<<"the total inital area is..."<<areaorg<<endl;
  //TRES.surface(10,1.4);
  TRES.surface(10,1.4,-1);
  pdb->setflg(1);
  pdb->chn->surface(10,1.4);
  areaorg=pdb->chn->area;
}
else if(hydrophobic==8) {
  TRES.surface(10,1.4,-1);
  pdb->setflg(1);
  pdb->chn->surface(10,1.4);
  areaorg=pdb->chn->calcareadiff(hydroflg);
}

else if(emilcharge) {
  latsurfv->gb->putoff();
  latsurfv->update(pdb);
  latsurfv->gb->puton(pdb);
  areaorg=latsurfv->calcarea(pdb);
}

for(; ;)
{
 if(mutationonly) break;
 tot=0;n++;rotm=0;
 
 if(emilcharge) {
	latsurfv->clearsurfvmask(pdb);
	latsurfv->gb->putoff();
	latsurfv->update(pdb);
	latsurfv->gb->puton(pdb);
	areaorg=latsurfv->calcarea(pdb);
 }

 if(hydrophobic==6||hydrophobic==7) {
      pdb->chn->surface(8,1.4);
      areaorg=pdb->chn->area;
 }
 else if(hydrophobic==8) {
      pdb->chn->surface(8,1.4);
      areaorg=pdb->chn->calcareadiff(hydroflg);
 }

 //for(chn=pdb->chn;chn;chn=chn->next)
 //for(r1=chn->res;r1;r1=r1->next)
 for(itt=0;itt<nite;itt++)
 {
  r1=solventorder[itt];
  chn=r1->chn;
  if(r1->flag<-9999) continue;
  if(iweight) {
  	if(n==0) iterateweight=0;
  	else if(n==1) iterateweight=0.5;
     	else if(n==2) iterateweight=1.0;
	else iweight=0;
  }
//putoff this residue
 
  if(emilcharge) {
	latsurfv->clearsurfvmask(pdb);
	latsurfv->gb->putoff();
        latsurfv->update(pdb);
	latsurfv->gb->puton(pdb);
	areaorg=latsurfv->calcarea(pdb);
  }

  lattice->putoff(r1,4,100);
  if((flag%10000)/1000>1&&rotm2==0&&fff&&r1->name!='P') d=minimize(r1,comb,mm,0);
  else d=bmax*r1->flag+clash(lattice,r1,4,100);
  
  if(emilcharge&&emilnewcharge) {
        float ee=clash(lattice1,r1,4,100);
        d+=ee;
  }

  pre=d;cflg=r1->flag;
  r1->transfer(xyz,0);
  //int nit=0;
  int isbury=0;

  float buy=0;

  if(hydrophobic==6||hydrophobic==7||hydrophobic==8) buy=r1->bury(4,100);

  for(r2=r1->more;r2;r2=r2->more)
  {
    if(emilcharge) {
	latsurfv->clearsurfvmask(pdb);	
	latsurfv->gb->putoff(r1);
    }

    r1->transfer(r2); r1->flag=r2->flag;
   
    if(emilcharge) {
        latsurfv->update(r1);
	latsurfv->gb->puton(r1);
    }	

    if((flag%10000)/1000>2&&rotm2==0&&fff&&r1->name!='P') e=minimize(r1,comb,mm,1);
    else e=bmax*r2->flag+clash(lattice,r1,4,100);
 
    if(emilcharge) {
	float ee;
	if(emilnewcharge) ee=clash(lattice1,r1,4,100);
	if(emilhydro==1) {
		float nare=pdb->getarea();
		ee+=(nare-areaorg)*0.025;
        	if(TRES.logg)cerr<<"the total area..."<<r1->name<<r1->id<<" "<<nare<<"  "<<areaorg<<endl;
	}
	//latsurfv->gb->putoff();
	//cerr<<r1->name<<r1->id<<" "<<ee<<"  "<<e<<endl;
	//cerr<<"the total area..."<<r1->name<<r1->id<<" "<<latsurfv->calcarea(pdb)<<endl;
	e+=ee;
    }  

    if(delphi==1) {
	if(allcharges==2) {
		float ee=calcdelphi(r1);
                if(TRES.logg)cerr<<r1->name<<r1->id0<<" solvation energy: "<<ee<<endl;
                e+=ee;
	}
	else if(allcharges==1&&strchr("RKHEQDNTS",r1->name)) {
		float ee=calcdelphi(r1);
                if(TRES.logg)cerr<<r1->name<<r1->id0<<" solvation energy: "<<ee<<endl;
                e+=ee;
	}
	else if(strchr("RKHEQDNTS",r1->name)) {
		float ee=calcdelphi(r1);
		if(TRES.logg)cerr<<r1->name<<r1->id0<<" solvation energy: "<<ee<<endl;
		e+=ee;
	}
    }

    /*
    if(fabs(ermsd>0.01) {
	e=e;
    }
    */
    if(hydrophobic==1) {
	r1->setflg(1);
	chn->surface(10,1.4);
	e+=r1->area*0.025;
	r1->setflg(0);
	//cerr<<"area:"<<r1->name<<r1->id0<<" "<<r1->area<<endl;
    }
    else if(hydrophobic==2) { 
	if(strchr("RKEHQYWFML",r1->name)) {
		chn->surface(10,1.4);
		e+=(chn->area-areaorg)*0.05;
		if(TRES.logg)cerr<<"the total area is :" <<r1->name<<" "<<chn->area<<"  "<<areaorg<<" "<<chn->area-areaorg<<endl;
	}
    }
    else if(hydrophobic==3) {
	float are=calcularea(chn);
	if(are==0) are=areaorg;
	if(TRES.logg)cerr<<"the total area is :" <<are<<"  "<<areaorg<<" "<<are-areaorg<<endl;
	e+=(are-areaorg)*0.025;
    }
    else if(hydrophobic==4) {
        float are=calcularea(chn,0);
        if(are==0) are=areaorg;
	if(TRES.logg)cerr<<"the total area is :" <<are<<"  "<<areaorg<<" "<<are-areaorg<<endl;
        e+=(are-areaorg)*0.05;
    }
    else if(hydrophobic==5) {
		
	if(strchr("RKEHQYWFML",r1->name)&&isbury==0) {
                chn->surface(10,1.4);
		if(r1->bury(4,100)>0.5&&r2==r1->more) {
			isbury=1;
			if(TRES.logg)cerr<<r1->name<<r1->id<<" is buryed "<<r1->bury(4,100)<<endl;
		}
                if(isbury==0) e+=(chn->area-areaorg)*0.05;
                if(TRES.logg)cerr<<"the total area is :" <<r1->name<<r1->id<<" "<<chn->area<<"  "<<areaorg<<" "<<chn->area-areaorg<<endl;
        }
    } 
    else if(hydrophobic==6) {

	//if(buy<0.5&&strchr("RKEHQYWFMLNDI",r1->name)) {
	if(buy<0.5) {
		chn->surface(8,1.4);
                e+=(chn->area-areaorg)*0.05;
                if(TRES.logg)cerr<<"the total area is :" <<r1->name<<" "<<chn->area<<"  "<<areaorg<<" "<<chn->area-areaorg<<endl;
	}
	else {
		e+=(chn->area-areaorg)*0.05;
		if(TRES.logg)cerr<<r1->name<<r1->id<<" is buryed "<<buy<<endl;
	}
    }
    else if(hydrophobic==7) {
        //if(buy<0.5&&strchr("RKEHQYWFMLNDI",r1->name)) {
        if(buy<0.9) {
                chn->surface(8,1.4);
                e+=(chn->area-areaorg)*0.05;
                if(TRES.logg)cerr<<"the total area is :" <<r1->name<<" "<<chn->area<<"  "<<areaorg<<" "<<chn->area-areaorg<<endl;
        }
        else {
                e+=(chn->area-areaorg)*0.05;
                if(TRES.logg)cerr<<r1->name<<r1->id<<" is buryed "<<buy<<endl;
        }
    }
    else if(hydrophobic==8) {
	if(buy<0.5) {
                chn->surface(8,1.4);
                e+=(chn->calcareadiff(hydroflg)-areaorg)*0.025;
                if(TRES.logg)cerr<<"the total area is :" <<r1->name<<" "<<chn->area<<"  "<<areaorg<<" "<<chn->area-areaorg<<endl;
        }
        else {
                e+=(pdb->chn->calcareadiff(hydroflg)-areaorg)*0.025;
                if(TRES.logg)cerr<<r1->name<<r1->id<<" is buryed "<<buy<<endl;
        }
    }


    enetemp[r2->id0]=e;
    if(e<d) 
    {
      r1->transfer(xyz,0);d=e;cflg=r1->flag;rotm=1;
      ener[r1->id0]=d;
    }
  }

  if(colony&&armsd[r1->id0]>0) {
    r4=0;pre0=ent0;ent0=0;
    float ce=1;
    if(coleng) ce=colengcoeff(r1,enetemp);			
    int coly=1;
    if(solventcolony) {
	if(solventbury[r1->id0]>0.8) coly=0;
    }
    for(r2=r1->more;r2;r2=r2->more) {
       ent=0;
       for(r3=r1->more;r3&&coly;r3=r3->more) {
	      e=r2->directrmsd(r3,4,1000);
	      e=pow(e,ncolony);
	      if(TRES.logg==-1001) cerr<<r3->id0<<endl;	
	      ent+=exp(-e*armsd[r1->id0]-ce*enetemp[r3->id0]/8.31/0.3);
       }
       if(printouteng) {
		Res *rr0=pdborg->isres(r1->id0);
		cout<<r2->id<<" "<<r1->name<<r1->id0<<" "<<rr0->directrmsd(r2,4,100)<<" "<<enetemp[r2->id0]<<" "<<-8.31*0.3*log(ent)<<endl;
       }
       if(ent<1) ent=enetemp[r2->id0];
       else ent=-8.31*0.3*log(ent);
       if(r2->id0==0) {
	      ent0=ent;
	      r4=r2;		
       }
       else if(ent<ent0){
	      ent0=ent;
              r4=r2;        
       }	
    } 
    if(r4) {
	r1->transfer(r4);r1->flag=r4->flag;
    	r1->transfer(xyz,0);cflg=r1->flag;
    	pre=pre0;
    	d=ent0;
    }
  }
  if(rotm1==0||n==2)
  {
   if(ener[r1->id0]>-1)
   {
     for(i=0;i<rotnum[r1->id0];i++)
     {     
       j=r1->tres->numrot-2;        
       if(rotamer[r1->id0][i][j]>5) break;
       r1->dihedral(0);
       r2=crtmore(r1,i);       
       if(r2==0) break;
       r2->next=0;r2->more=0;
       r1->transfer(r2);r1->flag=rotamer[r1->id0][i][j];
       if((flag%10000)/1000>3&&rotm2==0&&fff&&r1->name!='P') e=minimize(r1,comb,mm,1);
       else e=bmax*rotamer[r1->id0][i][j]+clash(lattice,r1,4,100);
       if(e<d) 
       {
         r1->transfer(xyz,0);d=e;cflg=r1->flag;rotm=1;
         ener[r1->id0]=d;
       }
       delete r2;
     }
   }
  }

  r1->transfer(xyz,1);r1->flag=cflg;
  lattice->puton(r1,4,1000);
  tot+=d;
  e=fabs(pre-d);
  if(e>0.01) rotm=1;
  if(fff==0&&e>5.0) rotm=1;
  if(TRES.logg)cerr<<seqn<<"..."<<n<<" "<<r1->id<<r1->name<<" "<<d;
  if(TRES.logg)cerr<<" "<<d-pre<<" "<<rotm<<endl;
}
//cerr<<"the total energy of this round: "<<tot<<endl;
if(colony) if(TRES.logg)cerr<<seqn<<"...the total colony energy after "<<n<<" iterations is: "<<tot<<endl;
else if(TRES.logg)cerr<<seqn<<"...the total energy after "<<n<<" iterations is: "<<tot<<endl;
if(n>10||fabs(tot-tot0)<0.1)rotm=0;
if(rotm==0&&rotm1==0&&rotm2==0) break;
if(rotm==0&&rotm1==0) {rotm2=0;if((flag%10000)/1000<=1)break;}
if(rotm==0) {rotm1=0;}
if(n>2&&tot-fabs(tot-tot0)-mine>5.) break;
tot0=tot;
cerr<<"rounds of finding the best rotamers:"<<n<<" with energy: "<<tot0<<endl;
//if(n>5) break;
if(n>niterate) break;
}
if((flag%10000)/1000>=1&&fff)tot=minimize();
else
{
  if(colonyline==0) tot=0;
  for(chn=pdb->chn;chn;chn=chn->next)
  for(r1=chn->res;r1;r1=r1->next)
  {
    if(r1->flag<-9999) continue;
    ener[r1->id0]=r1->flag+clash(lattice,r1,4,100);
    if(colonyline==0)tot+=ener[r1->id0];
  }
}
return tot;
}

float Scap::minimize()
{
Res *r1;
int m,*comb,i,j,m1,mm[5];
float xyz[200],xyz1[200],co[5];
Atm *a,*rote[20];
float tot,f,e,d,pre,cflg;
int rotm,step,nnn;
Rotate side;
comb=combination(1);
cerr<<"refine the sidechain conformation using 2-5 degree rotation..."<<endl;
for(i=0;i<5;i++)mm[i]=1;
for(i=0;i<5;i++) for(j=0;j<5;j++) if(j<=i) mm[4-j]*=3;

lattice->putoff();
lattice->flag=0;
lattice->grdsiz=2.;
lattice->radall=15.;
lattice->ready(pdb);
Chn *chn;
for(chn=pdb->chn;chn;chn=chn->next)
for(r1=chn->res;r1;r1=r1->next)
{
if(r1->flag>-9999) lattice->puton(r1,4,1000);
}


nnn=0;step=5;
for(; ;)
{
  rotm=0;nnn++;tot=0;
  for(chn=pdb->chn;chn;chn=chn->next)
  for(r1=chn->res;r1;r1=r1->next)
  {
    if(r1->flag<-9999) continue;
    if(r1->name=='P') continue;
    m=r1->rotable(rote,4,100);
    if(m>=2)m=2;
    if(m>0)m1=mm[m-1];
    else m1=0;
    lattice->putoff(r1,4,100);
    r1->transfer(xyz,0);
    r1->transfer(xyz1,0);
    d=bmax*r1->flag+clash(lattice,r1,4,100);
    pre=d;
    if(TRES.logg)cerr<<"original..."<<r1->flag<<" "<<d<<endl;
    for(i=0;i<m1;i++)
    {
      r1->transfer(xyz1,1);
      for(j=0;j<m;j++) co[j]=comb[i*5+j+5-m];
      for(j=0;j<m;j++)
      {
        a=rote[j];
        f=step*co[j];
        if(f==0) continue;
        side.rotate(a,f);
      }
      if(proring==0) {
      if(strchr(force,'1')&&r1->name!='P'&&ring&&torsionring) cflg=tormax*r1->clashnoring(2,1)+clash(lattice0,r1,4,100);
      else if(strchr(force,'1')&&r1->name!='P'&&ring) cflg=tormax*r1->clashnoring(2,0)+clash(lattice0,r1,4,100); 
      else if(strchr(force,'1')&&r1->name!='P') cflg=tormax*r1->clash(2)+clash(lattice0,r1,4,100);
      else cflg=clash(lattice0,r1,4,100);
      }
      else {
      if(strchr(force,'1')&&ring&&torsionring) cflg=tormax*r1->clashnoring(2,1)+clash(lattice0,r1,4,100);
      else if(strchr(force,'1')&&ring) cflg=tormax*r1->clashnoring(2,0)+clash(lattice0,r1,4,100);
      else if(strchr(force,'1')) cflg=tormax*r1->clash(2)+clash(lattice0,r1,4,100);
      else cflg=clash(lattice0,r1,4,100);
      }
      e=bmax*cflg+clash(lattice,r1,4,100);
      if(TRES.logg)cerr<<i<<" "<<cflg<<endl;
      if(e<d)
      {
        r1->transfer(xyz,0);
        r1->flag=cflg;
        d=e;
      }
    }
    r1->transfer(xyz,1);
    lattice->puton(r1,4,100);
    tot+=d; 
    ener[r1->id0]=d;
    if(fabs(pre-d)>0.01) 
    {
      rotm=1; 
      if(TRES.logg)cerr<<seqn<<"..."<<nnn<<"  minimized..."<<r1->id<<r1->name<<" "<<d<<" ";
      if(TRES.logg)cerr<<pre-d<<" "<<step<<endl;
    }
  }

  if(TRES.logg)cerr<<seqn<<"...the total energy after minimization is: "<<tot<<endl;
  if((nnn>5||rotm==0)&&step==2) break; 
  if(rotm==0||nnn>5){step=step/2;nnn=0;}
}
delete [] comb;

tot=0;
for(chn=pdb->chn;chn;chn=chn->next)
for(r1=chn->res;r1;r1=r1->next)
{
 if(r1->flag<-9999) continue;
 if(proring==0) {
 if(strchr(force,'1')&&r1->name!='P'&&ring&&torsionring) ener[r1->id0]=tormax*r1->clashnoring(2,1)+clash(lattice0,r1,4,100);
 else if(strchr(force,'1')&&r1->name!='P'&&ring) ener[r1->id0]=tormax*r1->clashnoring(2,0)+clash(lattice0,r1,4,100);
 else if(strchr(force,'1')&&r1->name!='P') ener[r1->id0]=tormax*r1->clash(2)+clash(lattice0,r1,4,100);
 else ener[r1->id0]=clash(lattice0,r1,4,100);
 }
 else {
 if(strchr(force,'1')&&ring&&torsionring) ener[r1->id0]=tormax*r1->clashnoring(2,1)+clash(lattice0,r1,4,100);
 else if(strchr(force,'1')&&ring) ener[r1->id0]=tormax*r1->clashnoring(2,0)+clash(lattice0,r1,4,100);
 else if(strchr(force,'1')) ener[r1->id0]=tormax*r1->clash(2)+clash(lattice0,r1,4,100);
 else ener[r1->id0]=clash(lattice0,r1,4,100);
 }
 tot+=ener[r1->id0]*bmax+clash(lattice,r1,4,100);
}
return tot;
}

float Scap::minimize(Res *r1,int *comb,int *mm,int flg)
{
int m,i,j,m1,k;
float xyz[200],xyz1[200],co[5];
Atm *a,*rote[20];
float pre,f,e,d;
int step,nnn;
Rotate side;
step=5;

m=r1->rotable(rote,4,100);
if(m>=5)m=5;
if(m>0)m1=mm[m-1];
else m1=0;

lattice->putoff(r1,4,100);
r1->transfer(xyz,0);
r1->transfer(xyz1,0);
d=r1->flag+clash(lattice,r1,4,100);
pre=d;
nnn=0;
for(; ;)
{
nnn=0;
for(i=0;i<m1;i++)
{
  r1->transfer(xyz1,1);
  for(j=0;j<m;j++) co[j]=comb[i*5+j+5-m];
  k=0;
  for(j=0;j<m;j++) if(co[j]!=0) k=1;
  if(k==0) continue; 
  for(j=0;j<m;j++)
  {
    a=rote[j];
    f=step*co[j];
    if(f==0) continue;
    side.rotate(a,f);
  }
  if(proring==0) {
  if(strchr(force,'1')&&r1->name!='P'&&ring&&torsionring) f=tormax*r1->clashnoring(2,1)+clash(lattice0,r1,4,100);
  else if(strchr(force,'1')&&r1->name!='P'&&ring) f=tormax*r1->clashnoring(2,0)+clash(lattice0,r1,4,100);
  else if(strchr(force,'1')&&r1->name!='P') f=tormax*r1->clash(2)+clash(lattice0,r1,4,100);
  else f=clash(lattice0,r1,4,100);
  }
  else {
  if(strchr(force,'1')&&ring&&torsionring) f=tormax*r1->clashnoring(2,1)+clash(lattice0,r1,4,100);
  else if(strchr(force,'1')&&ring) f=tormax*r1->clashnoring(2,0)+clash(lattice0,r1,4,100);
  else if(strchr(force,'1')) f=tormax*r1->clash(2)+clash(lattice0,r1,4,100);
  else f=clash(lattice0,r1,4,100);
  }
  e=bmax*f+clash(lattice,r1,4,100);
   
  if(e<d)
  {
    r1->transfer(xyz,0);r1->flag=f;
    d=e;
    nnn=1;
  }
}
  TRES.copy(xyz,xyz1,r1->tres->number*3+10);
  if(nnn==0&&step<=2) break;
  step=step/2;  
  if(flg)break;
}
r1->transfer(xyz,1);
lattice->puton(r1,4,100);
if(fabs(ener[r1->id0]-d)>0.01) 
{
 if(TRES.logg)cerr<<seqn<<"..."<<nnn<<"  minimized..."<<r1->id<<r1->name<<" "<<d<<" ";
 if(TRES.logg)cerr<<pre-d<<endl;
}
return d;
}

int *Scap::combination(int m)
{

int a[5];
int n,i;
int *s;

n=2*m+1;
i=n*n*n*n*n*5;

s=new int[i];

for(n=0;n<i;n++) s[n]=-999;

n=0;

for(a[0]=-m;a[0]<=m;a[0]++)
for(a[1]=-m;a[1]<=m;a[1]++)
for(a[2]=-m;a[2]<=m;a[2]++)
for(a[3]=-m;a[3]<=m;a[3]++)
for(a[4]=-m;a[4]<=m;a[4]++)
{
 for(i=0;i<5;i++) s[n*5+i]=a[i];
 n++; 
}
return s;
}

void Scap::bbdep(char *name)
{
FILE *fp;
char line[256],line0[256];
char c;
int i;
int bbnum;
float bbdep[1000][8];
int k,j;
float ii[8],x;
Res *r;
Atm *a;

bbnum=0;
line0[0]='\0';
c='0';
fp=fopen(name,"r");
pdb->dihedral(0);
while(fgets(line,256,fp)!=NULL)
{
  if(line0[0]!='\0'&&strncmp(line,line0,15)!=0) //return bbnum;
  {
    if(TRES.logg)cerr<<c<<" "<<bbdep[0][0]<<" "<<bbdep[0][1]<<endl;
    x=0;
    j=0;
    for(i=0;i<bbnum;i++)
    {
      x+=bbdep[i][2];
      j++;
      if(x>1) break;
    }
    bbnum=j;
    for(Chn *chn=pdb->chn;chn;chn=chn->next)
    for(r=chn->res;r;r=r->next)
    {
      if(r->name!=c)continue;
      if(r->flag<-9999) continue;
      a=r->atm->next;
      x=a->chi;
      if(x==999)x=-60;
      x=fabs(x-bbdep[0][0]);
      if(x>180)x=fabs(360-x);
      if(x>5) continue;
      x=a->next->chi;
      if(x==999)x=60;
      x=fabs(x-bbdep[0][1]);
      if(x>180)x=fabs(360-x);
      if(x>5) continue;
      if(rotnum[r->id0]!=0) continue;
      for(i=0;i<bbnum;i++)
      {
        j=rotnum[r->id0];
        for(k=0;k<5;k++)
        rotamer[r->id0][j][k]=bbdep[i][3+k];
        rotnum[r->id0]++;  
      }
    }
    bbnum=0;
    line0[0]='\0';
  }
  if(bbnum>=999) continue;
  c=TRES.swap(line);
  if(c=='?') continue;
  sscanf(line+3,"%f%f",ii,ii+1);
  sscanf(line+32,"%f%f%f%f%f",ii+2,ii+3,ii+4,ii+5,ii+6);
  for(i=0;i<2;i++) bbdep[bbnum][i]=ii[i]+5;
  for(i=2;i<7;i++) bbdep[bbnum][i]=ii[i];
  bbdep[bbnum][7]=180;
  bbnum++;
  strncpy(line0,line,16);
}
fclose(fp);
}

void Scap::bbind(char *name) 
{
FILE *fp; 
int i,m,k,n;
Tres *tres;
float ii[7];
Res *r,*rr[30][500];
char line[256];
int  mr,*comb,mm[5],nn,nnn,h,co[5],nrr[30];

for(i=0;i<30;i++) nrr[i]=0;

for(Chn *chn=pdb->chn;chn;chn=chn->next)
for(r=chn->res;r;r=r->next)
{
  mr=r->name-'A';
  i=nrr[mr];
  rr[mr][i]=r;
  nrr[mr]++;
}

nn=flag/10000;

comb=combination(1);
for(i=0;i<5;i++)mm[i]=1;
for(i=0;i<5;i++) for(k=0;k<5;k++) if(k<=i) mm[4-k]*=3;

fp=fopen(name,"r");

if(fp==0){cerr<<"could open rotamer library\n";exit(0);}

while(fgets(line,256,fp)!=NULL)
{
tres=TRES[line];
if(tres==0) continue;
nnn=tres->numrot-2;
if(tres==0) continue;
for(i=0;i<5;i++) ii[i]=180;
if(nnn==1)sscanf(line+3,"%f",ii+0);
else if(nnn==2)sscanf(line+3,"%f%f",ii+0,ii+1);
else if(nnn==3)sscanf(line+3,"%f%f%f",ii+0,ii+1,ii+2);
else if(nnn==4)sscanf(line+3,"%f%f%f%f",ii+0,ii+1,ii+2,ii+3);
else if(nnn==5)sscanf(line+3,"%f%f%f%f%f",ii+0,ii+1,ii+2,ii+3,ii+4);
if(nnn>nn&&nn>0)nnn=nn;
for(mr=0;mr<nrr[tres->name-'A'];mr++)
{
r=rr[tres->name-'A'][mr];
//for(r=chn->res;r;r=r->next)
//{
 if(r->flag<-9999)continue;
 if(r->name!=tres->name) continue;
 n=r->id0;
 if(nn>0&&nnn>0)
 {
   for(h=0;h<mm[nnn-1];h++) 
   {
     m=rotnum[n];
     for(k=0;k<tres->numrot-2;k++)co[k]=0;
     for(k=0;k<nnn;k++)co[k]=comb[h*5+k+5-nnn];
     for(k=0;k<tres->numrot-2;k++)
     {
       if(k==4) rotamer[n][m][k]=ii[k]+10*co[k]; 
       else rotamer[n][m][k]=ii[k]+20*co[k];
     }
     rotnum[n]++;
   }
 }
 else 
 {
   m=rotnum[n];
   for(k=0;k<tres->numrot-2;k++)
   {
     rotamer[n][m][k]=ii[k];
   }
   rotnum[n]++;
 }
}
}
fclose(fp);
return;
}

void Scap::rmsd(Chn *s2,float angle,FILE *fp)
{
float d;
Res *r1,*r2;
Atm *a1,*a2;
int i,j,k,n,m;
float x1,x2,x3[256];
int nn1[256],xx1[256],xx2[256],nn2[256],xx3[256],nn3[256];
if(fp) fprintf(fp,"%s %d\n\n",s2->pdb->name,s2->number);
Rotate side;

Chn *chn=pdb->chn;
i=chn->number;
for(j=0;j<256;j++) {xx1[j]=0;nn1[j]=0;xx2[j]=0;nn2[j]=0;xx3[j]=0;nn3[j]=0;}
for(j=0;j<i;j++)
{
  k=0;m=0;
  r1=(*chn)[j];
  r2=(*s2)[j];
  //if(strchr("PCGA",r1->name)) continue;
  if(r1->flag<-9999) continue;
  a1=(*r1)[4];
  a2=(*r2)[4];
  if(a1==0||a2==0) continue;
 
  x1=a1->dihedral(0); 
  x2=a2->dihedral(0);
  d=fabs(x1-x2);  
  if(d>180) d=360-d;
  if(a1->tatm->balance==1&&d>90) 
  { 
     side.rotate(a1,180);
     d=180-d;
  }
  if(d<=angle)  {xx1[r1->name]++;m=1;}
  nn1[r1->name]++;

  a1=a1->next;  
  a2=a2->next;
  if(a1==0||a2==0) continue;
  x1=a1->dihedral(0);
  x2=a2->dihedral(0);
  d=fabs(x1-x2);
  if(d>180) d=360-d;
  if(a1->tatm->balance==1&&d>90) 
  {
    side.rotate(a1,180);
    d=180-d;
  }
  if(d<=angle)  {xx2[r1->name]++;k=1;}
  nn2[r1->name]++;

  if(k==1&&m==1)
  { 
    xx3[r1->name]++;
  }
  nn3[r1->name]++;
}

for(j=0;j<256;j++)
{
if(nn1[j]==0) continue;
x1= xx1[j];
x1=x1/nn1[j];
fprintf(fp,"%c %f %d ",char(j),x1,nn1[j]);
if(nn2[j]==0) {/*fprintf(fp,"\n");*/continue;}
x1=xx2[j];
x1=x1/nn2[j];
fprintf(fp,"%f %d ",x1,nn2[j]);
x1=xx3[j];
if(nn3[j]==0) 
{//fprintf(fp,"\n");
 continue;
}
x1=x1/nn3[j];
fprintf(fp,"%f %d\n",x1,nn3[j]);
}

i=pdb->chn->number;
for(j=0;j<256;j++) {x3[j]=0;xx1[j]=0;nn1[j]=0;xx2[j]=0;nn2[j]=0;xx3[j]=0;nn3[j]=0;}

d=0;m=0;
for(j=0;j<i;j++)
{
 r1=(*chn)[j];
 r2=(*s2)[j];
 //if(strchr("PCGA",r1->name)) continue; 
 if(r1->flag<-9999)continue;
 n=r1->number;
 for(k=5;k<n;k++)
 {
  a1=(*r1)[k];
  a2=(*r2)[k];
  if(a1==0||a2==0) continue;
  if(a1->good==0||a2->good==0) continue;
  if(a1->name[1]=='H'||a2->name[1]=='H') continue;
  x3[r1->name]+=TRES.distsqr(a1->xyz,a2->xyz);
  nn1[r1->name]++;
  d+=TRES.distsqr(a1->xyz,a2->xyz);
  m++;
 }
}

for(j=0;j<256;j++)
{
if(nn1[j]==0) continue;
x1= x3[j];
x1=x1/nn1[j];
x1=sqrt(x1);
fprintf(fp,"%c %f %d\n ",char(j),x1,nn1[j]);
}

d=sqrt(d/m);
fprintf(fp,"the rmsd for %s is %f\n",s2->pdb->name,d);
fprintf(fp,"the above torsional angle is based on:%f\n",angle); 
}

void Scap::printhelp()
{
printf("scap -f -h -m -o -i -c pdbfile rotamerfile residuelist\n");
printf("-f     topology flag. choices are -f=100 meaning all atom or -f=102 meaning heavy atom.optional.default heavy atom, \n");
printf("-h     print out help file.optional\n");
printf("-m     minimization flag.optional.default is off\n");
printf("-o     all results from the initial conformations written to disk.optional.default is off");
printf("-i     number of initial conformations: format: -i=100, -i=10, -i=1, etc....default is 5\n");
printf("-c     colony energy turned off.default is off. optional\n");
printf("file   pdb file\n");
printf("rotamerfile  rotamer filename, optional. \n");
printf("residuelist  residue list whose sidechain is to be predicted or mutated.optional\n");
}
void Scap::printnewhelp()
{
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"scap is a protein side-chain program with the following capabilities:\n");
fprintf(stderr,"A. predict side-chain conformations of a whole protein\n");
fprintf(stderr,"B. predict specified side-chains in a protein\n");
fprintf(stderr,"C. mutate specified residues in a protein\n");
fprintf(stderr,"D. sidechain modeling with a finite difference PB calculation using DELPHI\n");
fprintf(stderr,"\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage: scap -prm num -rtm num -seed num -min num -out num -ini num -self num  file.pdb res_scap.list\n");
fprintf(stderr,"-prm   	force parameter flag. default is 1\n");
fprintf(stderr,"	-prm 1: CHARMM22 with all atom model\n");
fprintf(stderr,"	-prm 2: AMBER with all atom model\n");
fprintf(stderr,"	-prm 3: CHARMM with heavy atom model\n");
fprintf(stderr,"	-prm 4: AMBER with heavy atom model\n");
fprintf(stderr,"-min   	side-chain refinement with 2 degree steps. default is 0.\n");
fprintf(stderr,"	num ranges from 0-4 with 4 most thorough refinement.\n");
fprintf(stderr,"-seed   random seed number\n");
fprintf(stderr,"-out   	output model. num from 0-2. default is 0.\n");
fprintf(stderr,"	-out 0: only the final structure outputted.\n");
fprintf(stderr,"	-out 1: output the final and intermediate structure.\n");
fprintf(stderr,"	-out 2: output all residue conformations of different rotamers plus -out 1\n");
fprintf(stderr,"-ini   	the number of initial structures tried.  default is 3.\n");
fprintf(stderr,"-self	option to retain original sidechain as rotamer.default is 0\n");
fprintf(stderr,"	-self 0: original sidechain conformations not used as rotamer.\n"); 
fprintf(stderr,"	-self 1: original sidechain conformations used as a regular rotamer\n"); 
fprintf(stderr,"	-self 2: original sidechain conformations used as rotamer and used as the first initial conformations.\n"); 
fprintf(stderr,"-force  force field flag. default is vdw and torsional energy.\n");
fprintf(stderr,"        -s: force field is vdw, torsional and electrostatics from DELPHI\n");
fprintf(stderr,"        -h: force field is vdw, torsional and hydrogen bond from STICKLE method\n");
fprintf(stderr,"        -c: force field is vdw, torsional and electrostatics based on distance-dependent dielectric\n");
fprintf(stderr,"-scale  grid scale for delphi program. default is 1.0\n");
fprintf(stderr,"-prbrad probe radius for delphi program. default is 1.6\n");
fprintf(stderr,"-perfil percent of box filled for delphi program. default is 70\n");
fprintf(stderr,"-indi   protein interior dielectric constant for delphi program. default is 4.0\n");
fprintf(stderr,"-exdi   protein exterior dielectric constant for delphi program. default is 80\n");
/*
fprintf(stderr,"-heta   atoms of non-standard amino acids.  default is 1.\n");
fprintf(stderr,"        -heta 1: consider atoms of  non-standard amino acids.\n");
fprintf(stderr,"        -heta 0: not consider atoms of  non-standard amino acids.\n");
*/
fprintf(stderr,"-rtm    side-chain rotamer library.  default is 1.\n");
fprintf(stderr,"	-rtm 1: large side-chain rotamer library\n");
fprintf(stderr,"        -rtm 2: mix side-chain rotamer library\n");
fprintf(stderr,"        -rtm 3: medium side-chain rotamer library\n");
fprintf(stderr,"        -rtm 4: small side-chain rotamer library\n");
fprintf(stderr,"file.pdb   	pdb file\n");
fprintf(stderr,"res_scap.list  	residue list whose sidechain is to be predicted or mutated.optional. default is none\n");
fprintf(stderr,"		the residuelist must be ended with _scap.list\n");

}

void Scap::printeznewhelp()
{
printf("ezscap -prm num -res str -out num -cid char file.pdb\n");
printf("-prm   	force parameter flag. default is 1\n");
printf("	-prm 1: CHARMM22 with all atom model\n");
printf("	-prm 2: AMBER with all atom model\n");
printf("	-prm 3: CHARMM with heavy atom model\n");
printf("	-prm 4: AMBER with heavy atom model\n");
printf("-res   	residue to be predicted or mutated.with format as A3N:\n");
printf("	mutate A3 to N.\n");
printf("-out   	output model. num from 0-2. default is 0.\n");
printf("	-out 0: only the final structure outputted.\n");
printf("	-out 1: output the final and intermediate structure.\n");
printf("	-out 2: output all residue conformations of different rotamers plus -out 1\n");
printf("-cid   	chain id.  default is the first chain\n");
printf("-rtm    side-chain rotamer library.  default is 1.\n");
printf("	-rtm 1: large side-chain rotamer library\n");
printf("        -rtm 2: mix side-chain rotamer library\n");
printf("        -rtm 3: medium side-chain rotamer library\n");
printf("        -rtm 4: small side-chain rotamer library\n");
printf("file.pdb   	pdb file\n");

}
void Scap::write()
{
Res *r1;
float d,low[20],high[20];
char sec[20];
int i;
low[0]=1.0;high[0]=1.1; sec[0]='0';
low[1]=0.9;high[1]=1.0; sec[1]='0';
low[2]=0.8;high[2]=0.9;sec[2]='0';
low[3]=0.7;high[3]=0.8;sec[3]='0';
low[4]=0.6;high[4]=0.7;sec[4]='0';
low[5]=0.5;high[5]=0.6;sec[5]='0';
low[6]=0.4;high[6]=0.5;sec[6]='0';
low[7]=0.3;high[7]=0.4;sec[7]='0';
low[8]=0.2;high[8]=0.3;sec[8]='0';
low[9]=0.1;high[9]=0.2;sec[9]='0';
low[10]=0.;high[10]=0.1;sec[10]='0';
low[11]=0;high[11]=1.0;sec[11]='h';
low[12]=0;high[12]=1.0;sec[12]='e';
low[13]=0;high[13]=1.0;sec[13]='-';
low[14]=0;high[14]=1;sec[14]='0';

for(i=0;i<15;i++)
{
  cout<<"****"<<i<<"** "<<low[i]<<" "<<high[i]<<" "<<sec[i]<<"****"<<endl;
  for(r1=pdb->chn->res;r1;r1=r1->next)
  {
    r1->flag=0;
    if(strchr("PGAC",r1->name)) r1->flag=-999999;
    d=r1->bury(4,100);
    if(d<low[i]||d>high[i]) r1->flag=-99999;
    if(sec[i]!='0') if(r1->sec!=sec[i]) r1->flag=-99999;
  }
  rmsd(pdb->next->chn,20,stdout);
  fflush(stdout);
}
}

void Scap::insert(Res *r)
{
  Res *r1,*r2;
  Atm *a[10],*b[10];
  int i,j,k,n;
  float x,y;
  i=r->rotable(a,4,100);
  n=0;
  for(r1=r->more;r1;r1=r1->more)
  { 
    j=r1->rotable(b,4,100);
    if(i!=j) return;
    x=0;
    for(k=0;k<i;k++)
    { 
      y=fabs(a[k]->chi-b[k]->chi);
      if(y>180) y=360-y;
      if(y>10) {x=1000;break;}
      x+=y*y;
    }
    if(x<25) return;
    n++;
    r2=r1;
  }
  r2->more=new Res(r);
  r2->more->configure();
  if(TRES.logg)cerr<<seqn<<"....insert: "<< r->name<<r->id0<<"..."<<n+1<<endl;
  r2->more->more=0;
  r2->more->next=0;
}

void Scap::setlist(char chnid,int resid,char resold,char resnew) {
  
  Chn *chn;
  Res *rr;

  pdb->setflgr(-99999);
 
  chn=(*pdb)[chnid];   
  if(chn==0) {
	cerr<<"there is no chain with id:"<<chnid<<endl;
	cerr<<"check commands: -cid"<<endl;
	exit(0);
  }

  rr=chn->isresoid(resid);
  if(rr==0) {  
	cerr<<"there is no residue with id:"<<resid<<endl;
	cerr<<"check commands: -res"<<endl;
	exit(0);
  }

  if(rr->name!=resold) {
	cerr<<"the residue with id:"<<resid<<" has different name:"<<rr->name<<" "<<resold<<endl;
	cerr<<"check commands: -res"<<endl;
	exit(0);
  }
  
  if(resold==resnew) cerr<<"****residue :"<<rr->oid<<rr->name<<"  will be predicted***"<<endl;
  else	   	     cerr<<"****residue :"<<rr->oid<<rr->name<<"  will be mutated to:"<<resnew<<endl;
  if(resold!=resnew) rr->mutateResidue(resnew);
  
  if(rr->name=='G') rr->flag=-99999;
  else rr->flag=10000;
}

void Scap::readlist(char *f) {
  
  FILE *fp; 

  if(f==0||strlen(f)<1) return;
  fp=fopen(f,"r");

  if(fp==0) {
	if(TRES.logg)cerr<<"could not open file:"<<f<<endl;	
  	return;
  }

  char line[1000];
  Chn *chn;
  Res *rr;

  pdb->setflgr(-99999);
 
  Strhandler cc;
  char **tt;
  char c,r;
  int n;

  //read list such as: A,123,B; 123,B; A,123;
  while(fgets(line,256,fp)!=NULL) {

	if(line[0]=='!') continue;
	if(line[0]=='\n') continue;
	if(strlen(line)==0) continue;
	char *nu=strdup(line);
	nu=cc.clearendchar(nu,"\n\t\r ");
	if(nu==0) continue;
	strcpy(line,nu);nu=cc.strdel(nu);
	
	tt=cc.pairbytoken(line,",");
        tt=cc.clearendchar(tt,"\n\t, ");
	int nn=cc.gettotnum(tt);
 	if(nn==3) {
		if(strlen(tt[0])!=0) c=tt[0][0];
		else 		     c=' ';
		if(strlen(tt[1])!=0) n=atoi(tt[1]);
		else {
			cerr<<"error in scap list file"<<endl;
			cerr<<line<<endl;
			cerr<<"the residue id does not exist"<<endl;
			exit(0);
		}
		if(strlen(tt[2])!=0) r=tt[2][0];
		else {
			r='0';
		}		
	}
	else if(nn==2) {
		if(strlen(tt[0])==0&&strlen(tt[1])==0) {
			cerr<<"do not understand ..."<<line<<endl;
                        exit(0);
		}
		
		if(strlen(tt[1])>0&&tt[1][0]>='0'&&tt[1][0]<='9') {
			if(strlen(tt[0])==0) c=' ';
			else                c=tt[0][0];
                        n=atoi(tt[1]);
                        r='0';
		}
		else {
			c=' ';
			if(strlen(tt[0])==0) {
				cerr<<"error in scap list file"<<endl;
				cerr<<line<<endl;
				cerr<<"the residue id not exist"<<endl;
				exit(0);
			}
                        else n=atoi(tt[0]);
			if(strlen(tt[1])>0) r=tt[1][0];
			else r='0';
		}
	}
	else if(nn==1) {
		c=' ';
		r='0';
		if(strlen(tt[0])==0) {
                         cerr<<"error in scap list file"<<endl;
                         cerr<<line<<endl;
                         cerr<<"the residue id not exist"<<endl;
                         exit(0);
                }
		else n=atoi(tt[0]);
	}
	else {
		cerr<<"do not understand ..."<<line<<endl;
		continue;
	}

	tt=cc.strdel(tt);
        chn=(*pdb)[c];
	if(chn==0) {
		cerr<<"Warning! chain: "<<line<<" does not exist!"<<endl;
		continue;
   	}
	
	//rr=chn->isres(n-chn->start);
	rr=chn->isresoid(n);

	if(rr==0) {
		cerr<<"Warning! residue: "<<line<<" does not exist!"<<endl;
		continue;
	}
	
	if(r=='0'||rr->name==r) cerr<<"****residue :"<<rr->oid<<rr->name<<" in chain: "<<rr->chn->id<<"  will be predicted***"<<endl;
	else	   cerr<<"****residue :"<<rr->oid<<rr->name<<" in chain: "<<rr->chn->id<<"  will be mutated to:"<<r<<endl;
	Tres *tr=TRES[r];
        if(tr) {
		if(r!=rr->name) rr->mutateResidue(r);
	}
	else if(r!=' '&&r!='0'){
		cerr<<"there is no residue with id:"<<r<<" in chain: "<<c<<endl;
		cerr<<"mutation is not performed"<<endl;
	}
	if(rr->name=='G') rr->flag=-99999;
	else rr->flag=10000;
	n++;
  } 	

  if(n==0) {
 	for(chn=pdb->chn;chn;chn=chn->next)
	for(rr=chn->res;rr;rr=rr->next) {
		if(rr->name=='A'||rr->name=='G') rr->flag=-99999;
		else rr->flag=10000;
	}
  }
}

void Scap::hookside() {

  for(Chn *chn=pdb->chn;chn;chn=chn->next)
  for(Res *r=chn->res;r;r=r->next)
  {
      if(r->flag<-9999) continue;
      if(r->name=='G') continue;
      r->addsidechain();
  }
  pdb->configure();
}

void Scap::averagermsd() {
 
  Res *r1,*r2,*r3; 
  
  int n=pdb->getresnum(); 
  if(armsd) delete [] armsd;
  armsd=new float[n];
  for(int i=0;i<n;i++) armsd[i]=0;

  for(Chn *chn=pdb->chn;chn;chn=chn->next)
  for(r1=chn->res;r1;r1=r1->next) {
	if(r1->flag<-9999) continue;
	armsd[r1->id0]=0;
	n=0;int mm=0,mm1=0;
	for(r2=r1->more;r2;r2=r2->more) {
		mm++;mm1=0;
		if(mm>300) break;
		for(r3=r2->more;r3;r3=r3->more) {
			if(r3==r2) continue;
			mm1++;
			if(mm1>300) break;
			armsd[r1->id0]+=r2->directrmsd(r3,4,100);
			n++;
		}
	}
	if(n==0) continue;
	armsd[r1->id0]=armsd[r1->id0]/n;
	if(TRES.logg)cerr<<"the average rmsd of sidechain rotamers of residue:"<<r1->id<<r1->name<<" is:"<<armsd[r1->id0]<<endl;
	armsd[r1->id0]=log(1/coleffect)/pow(armsd[r1->id0],ncolony);
	//armsd[r1->id0]=log(2)/pow(armsd[r1->id0],ncolony);
	//armsd[r1->id0]=log(2)/pow(armsd[r1->id0],3);
  }
}

Res *Scap::findminimaltorsionmore(Res *s) {

  Res *r1,*r2;

  float d=100000;

  r2=0;
  int n=0;
  for(r1=s->more;r1;r1=r1->more) {
    
    float e;
    if(proring==0) {
    if(r1->name!='P'&&ring&&torsionring) e=r1->clashnoring(2,1);
    else if(r1->name!='P'&&ring) e=r1->clashnoring(2,0);
    else if(r1->name!='P') e=r1->clash(2);
    else e=r1->flag;
    }
    else {
    if(ring&&torsionring) e=r1->clashnoring(2,1);
    else if(ring) e=r1->clashnoring(2,0);
    else e=r1->clash(2);
    }
    if(n==0) {
      d=e;
      r2=r1;
    }
    else if(e<d) {
      d=e;
      r2=r1;
    }
    n++;
  }
  return r2;
}

Res *Scap::findminimalbordermore(Res *s)
{

  Res *r1,*r2;

  float d=100000;

  r2=0;
  int n=0;
  for(r1=s->more;r1;r1=r1->more) {
    
    float e;
    if(proring==0) {
    if(r1->name!='P'&&ring&&torsionring)e=r1->clashnoring(3,1);
    else if(r1->name!='P'&&ring)e=r1->clashnoring(3,0);
    else if(r1->name!='P') e=r1->clash(3);
    else e=r1->flag;
    }
    else {
    if(ring&&torsionring)e=r1->clashnoring(3,1);
    else if(ring)e=r1->clashnoring(3,0);
    else e=r1->clash(3);
    }
    if(n==0) {
      d=e;
      r2=r1;
    }
    else if(e<d) {
      d=e;
      r2=r1;
    }
    n++;
  }
  return r2;
}


Res *Scap::findminimalmore(Res *s) {

  Res *r1,*r2;

  float d=100000;

  r2=0;
  int n=0;
  for(r1=s->more;r1;r1=r1->more) {
    if(n==0) {
      d=r1->flag;
      r2=r1;
    }
    else if(r1->flag<d) {
      d=r1->flag;
      r2=r1;
    }
    n++;
  }
  return r2;
}

void Scap::addbackboneh() {

  Chn *chn;
  Res *r;
  Atm *a;

  for(chn=pdb->chn;chn;chn=chn->next)chn->header();

  
  for(chn=pdb->chn;chn;chn=chn->next)
  {
   Res *r1=0;
  for(r=chn->res;r;r=r->next)
  {
	if(r->isatm(" HN ")==0&&r->tres->isatm(" HN ")) {
		if(r1&&r->id-r1->id==1) {
		   a=createnewatm(r," HN ");
		   if(a) r->attachAtm(a);
		   r->configure();
		}
	}
	if(r->name=='G') {
		if(r->isatm("1HA ")==0&&r->tres->isatm("1HA ")) {
			a=createnewatm(r,"1HA ");
			if(a) r->attachAtm(a);
			r->configure();
		}
		if(r->isatm("2HA ")==0&&r->tres->isatm("2HA ")) {
			a=createnewatm(r,"2HA ");
			if(a) r->attachAtm(a);
			r->configure();
		}
	}
	else {
		if(r->isatm(" HA ")==0&&r->tres->isatm(" HA ")) {
			a=createnewatm(r," HA ");
			if(a) r->attachAtm(a);
			r->configure();
		}
	}
	r1=r;
  }
  }
  pdb->configure();
}

Atm *Scap::createnewatm(Res *r,char *s) {

  AtmGeom g;
  Atm *a;
  Tatm *t;//,*t2;
  float d;
  int yes=0;

  t=r->tres->isatm(s);  //hn atom
  if(t==0) return 0;
  if(strcmp(s," HN ")==0) {
	float *n,*ca,*c;
        float *th,*tn;//,*to;

        //h
        //h=aa->xyz;
        th=t->xyz;

        //n
        a=(*r)[0];
        if(a==0) return 0;
        n=a->xyz;
        tn=a->tatm->xyz;

        //ca

        a=(*r)[1];
        if(a==0) return 0;
        ca=a->xyz;
        //tca=a->tatm->xyz;


        //o
	/*
        o=r->temp+6;
        to=r->tres->head+6;
	*/
        //c
        c=r->temp+3;
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
        //float de=d;
	
	//add o

        //ange=TRES.angle(th,tn,tc)+TRES.angle(c,n,o);
        //ange=TRES.angle(th,tn,to);//+angee;
        //d=dd*dd+TRES.distsqr(n,o)-2*dd*TRES.distance(n,o)*cos(ange*3.14/180);
        //d=sqrt(d);
        //g.addbounds(o,d,0.001);
        //g.findcoo();
        g.findmincoo(2,0.0005);
        //cerr<<"nh.."<<TRES.angle(th,tn,to)<<"  "<<TRES.angle(g.coo,n,o)<<endl;
        //cerr<<"nhd..."<<de<<" "<<TRES.distance(g.coo,c)<<" "<<g.evaluate()<<"  "<<g.evaluaterange()<<endl;

 
  }
  else if(strcmp(s," HA ")==0) {


	float *n,*ca,*cb,*c;
	float *th,*tn,*tca,*tcb,*tc;

	//h
	//h=aa->xyz;
	th=t->xyz;

	//n
	a=(*r)[0];
	if(a==0) return 0;
	n=a->xyz;
	tn=a->tatm->xyz;

	//ca
	a=(*r)[1];
        if(a==0) return 0;
        ca=a->xyz;
        tca=a->tatm->xyz;
	
	//cb
	a=(*r)[4];
	if(a==0) return 0;
	cb=a->xyz;
	tcb=a->tatm->xyz;
 
	//c
	a=(*r)[2];
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
  else if(strcmp(s,"1HA ")==0) {

		
	float *n,*ca,*c;
	float *th,*tn,*tca,*tc;

	//h
	//h=aa->xyz;
	th=t->xyz;

	//n
	a=(*r)[0];
	if(a==0) return 0;
	n=a->xyz;
	tn=a->tatm->xyz;

	//ca
	a=(*r)[1];
        if(a==0) return 0;
        ca=a->xyz;
        tca=a->tatm->xyz;
	
	/*
	//cb
	a=(*r)[4];
	if(a==0) return 0;
	cb=a->xyz;
	tcb=a->tatm->xyz;
 	*/

	//c
	a=(*r)[2];
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

	/*//add cb
	ange=TRES.angle(th,tca,tcb);
        d=dd*dd+TRES.distsqr(ca,cb)-2*dd*TRES.distance(ca,cb)*cos(ange*3.14/180);
        d=sqrt(d);
        g.addbounds(cb,d,0.3);
 	*/
  	g.findmincoo(0.06);

 
  }
  else if(strcmp(s,"2HA ")==0) {

	float *n,*ca,*cb,*c;
	float *th,*tn,*tca,*tcb,*tc;

	//h
	//h=aa->xyz;
	th=t->xyz;

	//n
	a=(*r)[0];
	if(a==0) return 0;
	n=a->xyz;
	tn=a->tatm->xyz;

	//ca
	a=(*r)[1];
        if(a==0) return 0;
        ca=a->xyz;
        tca=a->tatm->xyz;
	
	//cb
	a=(*r)["1HA "];
	if(a==0) return 0;
	cb=a->xyz;
	tcb=a->tatm->xyz;
 
	//c
	a=(*r)[2];
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

  
  if(yes) g.num--;
  d=g.evaluate();
  if(d>0.2) {
	if(TRES.logg)cerr<<"h atoms "<<s<<" in residue "<<(r->id+r->chn->start)<<r->name<<" is not accurated predicted "<<d<<endl;
  }
  if(r->name!='G'&&d>0.2) return 0;
  else if(r->name=='G'&&d>0.5) return 0;  
  a=new Atm(t);
  a->res=r;
  for(int i=0;i<3;i++) a->xyz[i]=g.coo[i];
  return a;
}


Atm *Scap::createnewatmold(Res *r,char *s) {

  AtmGeom g;
  Atm *a;
  Tatm *t,*t2;
  float d;
  int yes=0;

  t=r->tres->isatm(s);  //hn atom
  if(t==0) return 0;
  if(strcmp(s," HN ")==0) {

	//N atom
	t2=r->tres->isatm(0);
	d=TRES.distance(t->xyz,t2->xyz); //hn--n
	a=(*r)[t2->id];
 	if(a==0) return 0;
	g.addbounds(a->xyz,d,0.01);

	//CA 
	t2=r->tres->isatm(1);
        d=TRES.distance(t->xyz,t2->xyz); //hn--CA
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.02);

	//O
	d=TRES.distance(t->xyz,r->tres->head+6); //hn--CA
        g.addbounds(r->temp+6,d,0.02);

	//C
        d=TRES.distance(t->xyz,r->tres->head+3); //hn--CA
        g.addbounds(r->temp+3,d,0.02);
	
  }
  else if(strcmp(s," HA ")==0) {
	//N atom
        t2=r->tres->isatm(0);
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.02);

        //CA 
        t2=r->tres->isatm(1);
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--CA
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.01);

        //C
	t2=r->tres->isatm(2);
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.02);

	//CB
        t2=r->tres->isatm(" CB ");
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.02);

  }
  else if(strcmp(s,"1HA ")==0) {
	//N atom
        t2=r->tres->isatm(0);
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.02);

        //CA 
        t2=r->tres->isatm(1);
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--CA
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.01);

        //C
        t2=r->tres->isatm(2);
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.02);

	/*
	//2HA
        t2=r->tres->isatm("2HA ");
        if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a) g.addbounds(a->xyz,d,0.05);
	*/
  }
  else if(strcmp(s,"2HA ")==0) {
        //N atom
        t2=r->tres->isatm(0);
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.02);

        //CA
        t2=r->tres->isatm(1);
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--CA
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.01);

	//C
        t2=r->tres->isatm(2);
        if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a==0) return 0;
        g.addbounds(a->xyz,d,0.02);

        //1HA
        t2=r->tres->isatm("1HA ");
	if(t2==0) return 0;
        d=TRES.distance(t->xyz,t2->xyz); //hn--n
        a=(*r)[t2->id];
        if(a) g.addbounds(a->xyz,1.6,0.8);
	if(a) yes=1;
  }

  g.findcoo();
  if(yes) g.num--;
  d=g.evaluate();
  if(d>0.2) {
	if(TRES.logg)cerr<<"h atoms "<<s<<" in residue "<<(r->id+r->chn->start)<<r->name<<" is not accurated predicted "<<d<<endl;
  }
  if(r->name!='G'&&d>0.2) return 0;
  else if(r->name=='G'&&d>0.5) return 0;  
  a=new Atm(t);
  a->res=r;
  for(int i=0;i<3;i++) a->xyz[i]=g.coo[i];
  return a;
}

void Scap::initialside(int flg) {

  if(initside==0) return;

  Res *r,*rr;
  Atm *a,*b;
  int i;
 
  Chn *chn;

  for(chn=pdb->chn;chn;chn=chn->next)
  for(rr=chn->res;rr;rr=rr->next)
  for(r=rr;r;r=r->more)
  for(a=r->atm;a;a=a->next)
  a->flag=1;

  for(chn=pdb->chn;chn;chn=chn->next)
  for(rr=chn->res;rr;rr=rr->next)
  for(r=rr;r;r=r->more)
  {
    if(r->flag<-9999) continue;
    for(a=r->atm;a;a=a->next)
    {
      a->flag=1;
      if(flg==1)  //not hydrogen
      {
         if(a->tatm->name[1]=='H') a->flag=0;
      }
      else if(flg==2) //no H on sidechain
      {
         if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>3) a->flag=0;
      }
      else if(flg==3)
      {
         if(!strchr("PGA",a->tatm->tres->name))
         {
           if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>3) a->flag=0;
           else if(a->tatm->name[1]!='H'&&a->tatm->id>4) a->flag=0; //no backbone
         }
      }
      else if(flg==4)
      {
        if(a->tatm->id>4) a->flag=0;  //not beyond cb
      }
      else if(flg==5) {  //no sidechain and hcb but including cb
	  if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>3) a->flag=0;
	  else if(a->tatm->name[1]!='H'&&a->tatm->id>4) a->flag=0;
      }
      else if(flg==6) { //only backbone atoms
          if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>3) a->flag=0;
          else if(a->tatm->name[1]!='H'&&a->tatm->id>3) a->flag=0;
      }
    }
  }

  //rebonding
  for(chn=pdb->chn;chn;chn=chn->next)
  for(rr=chn->res;rr;rr=rr->next)
  for(r=rr;r;r=r->more)
  {
    for(a=r->atm;a;a=a->next)
    {
      for(i=0;i<a->tatm->nbond;i++)
      {
        b=a->bond[i];
        if(b&&b->flag==0) a->bond[i]=0;
      }
    }
  }
  
  for(chn=pdb->chn;chn;chn=chn->next)
  for(rr=chn->res;rr;rr=rr->next)
  for(r=rr;r;r=r->more)
  {
    b=r->atm;a=r->atm;
    while(a)
    {
      if(a->flag==0)
      {
        b->next=a->next;
        a->next=0;
        delete a;    
        a=b->next;
        continue;    
      }
      b=a;   
      a=a->next;
    }
  }
}

void Scap::adddisulfide(){

Chn *chn;
Res *r,*r1,*r2,*r3;
Atm *a1,*a2,*a3,*a4;
float d,d1,d2,d3;

for(chn=pdb->chn;chn;chn=chn->next) {
  for(r=chn->res;r;r=r->next) {
      if(r->name!='C') continue;
      for(r1=r->next;r1;r1=r1->next) {
	  if(r1->name!='C') continue;
	  d=TRES.distance(r->atm->next,r1->atm->next);
	  if(d>7||d<5) continue;
	  float d4=10000;
          Res *rr1,*rr2;
	  for(r2=r;r2;r2=r2->more) {	 
             if(r->flag<=-9999&&r!=r2) continue;
             if(fabs(r->flag-10000)<0.1&&r2==r) continue;
             a1=(*r2)[" SG "];
	     a2=(*r2)[" CB "];
             if(a1==0||a2==0) continue;
	     for(r3=r1;r3;r3=r3->more) {
	          if(r1->flag<=-9999&&r1!=r3) continue;
		  if(fabs(r1->flag-10000)<0.1&&r3==r1) continue;
		  a3=(*r3)[" SG "];
 	          a4=(*r3)[" CB "];
		  if(a3==0||a4==0) continue;
		  d1=fabs(TRES.distance(a1,a3)-2.0);
	          d2=fabs(TRES.distance(a1,a4)-3.0);
		  d3=fabs(TRES.distance(a2,a3)-3.0);
		  d=d1*d1+d2*d2+d3*d3;
	          d=sqrt(d/3);
		  if(d<d4) {
			/*
			if(r2!=r) {
				if(r2->flag>0) r2->flag+=-3;
				else 	       r2->flag+=-3;
			}
			if(r3!=r1) {
				if(r3->flag>0) r3->flag+=-3;
                                else           r3->flag+=-3;
			}
			*/
			d4=d;
			rr1=r2;
			rr2=r3;
	          }		
	     }
	  }
	  if(d4<0.3) {
		r2=rr1;
		r3=rr2;
		if(r2!=r) {
                    if(r2->flag>0) r2->flag+=-10;
                    else           r2->flag+=-10;
                }
                if(r3!=r1) {
                    if(r3->flag>0) r3->flag+=-10;
                    else           r3->flag+=-10;
                }
	  }
      }
  }
}

}

void Scap::writeburylist(char *s,float a) {

	Chn *cn;
	Res *r;
	for(cn=pdb->chn;cn;cn=cn->next) cn->surface(10,1.4);

	FILE *fp;

	fp=fopen(s,"w");

	for(cn=pdb->chn;cn;cn=cn->next) {
		for(r=cn->res;r;r=r->next){
			float b=r->bury(0,100);
			if(b>a) continue;
			char c=cn->id;
			//if(b<a) {
			if(c!=' '){
				fprintf(fp,"%c %i\n",c,r->id+cn->start);
			}
			else {
				fprintf(fp,"%i\n",r->id+cn->start); 
			} 
			//}
			if(c!=' '){
                                fprintf(stderr,"%c %i %c %f\n",c,r->id+cn->start,r->name,b);
                        }
                        else {
                                fprintf(stderr,"%i %c %f\n",r->id+cn->start,r->name,b);
                        } 

		}
	}
	fclose(fp);

}

float Scap::calcularea(Chn *c) {

	
	c->write("temp.pdb");

	system("surfv1 1 1.4 rsh.siz temp.pdb>area.out");
	
	FILE *fp=fopen("area.out","r");

	char line[1000];

	while(fgets(line,256,fp)!=NULL) {

		char *s=strstr(line,"accessible area =");

		if(s) {
		
			s+=strlen("accessible area =");	

			float e=atof(s);
			
			//cerr<<"the total area is :" <<e<<endl;
			fclose(fp);
			return e;
		}
	}

	fclose(fp);

	return 0;

}

float Scap::calcularea(Chn *c,float ff ) {


        c->write("temp.pdb");

        char line[1000];
	sprintf(line,"surfv1 1 %f rsh.siz temp.pdb>area.out",ff);
        //system("surfv1 1 1.4 rsh.siz temp.pdb>area.out");
	//cerr<<line<<endl;
	system(line);

        FILE *fp=fopen("area.out","r");


        while(fgets(line,256,fp)!=NULL) {

                char *s=strstr(line,"accessible area =");

                if(s) {

                        s+=strlen("accessible area =");

                        float e=atof(s);

                        //cerr<<"the total area is :" <<e<<endl;
                        fclose(fp);
                        return e;
                }
        }

        fclose(fp);

        return 0;

}


void Scap::setcbetaflg() {
        Chn *c;
        Res *r;
	int n=pdb->maxresid0();
	if(cbetaflg) delete [] cbetaflg;cbetaflg=0;
	cbetaflg=new int[n*2+100];
	int i=0;
	for(i=0;i<n*2+100;i++) cbetaflg[i]=0;

        for(c=pdb->chn;c;c=c->next) {
                for(r=c->res;r;r=r->next) {
                        if(r->name=='G') continue;
                        Atm *a=r->isatmid(4);
                        if(a&&r->name!='P') continue;
                        else if(a==0) {
				cbetaflg[r->id0]=1;
                        }
                        else if(r->name=='P') {
                                Atm *a3=(*r)[" CD "];
                                Atm *a4=(*r)[" N  "];
                                if(a3&&a4&&TRES.distance(a3,a4)<1.8) continue;
				cbetaflg[r->id0]=1;
                        }
                        else {
				cbetaflg[r->id0]=1;
                        }
                }
        }
}

void Scap::calcbeta() {
	Chn *c;
	Res *r;
	cerr<<"calculating the average position of CB atoms..."<<endl;
	for(c=pdb->chn;c;c=c->next) {
		for(r=c->res;r;r=r->next) {
			if(r->name=='G') continue;
			if(cbetaflg[r->id0]==0) continue;
			if(r->flag<=-99990) continue;
			cerr<<"calculate the average position of CB of residue:"<<r->name<<r->oid<<endl;
			calcbeta(r);
		}
	}
}

void Scap::calcbeta(Res *r) {

  Res *res_temp;
  Atm *a1,*a2,*a3,*a4;
  float x,y;
  int i;
  Rotate side;
  float xyz[200];
  //float xyz0[200];
  float xyzcb[3];
  int   numcb=0;
  r->transfer(xyz,0);
  for(i=0;i<3;i++) xyzcb[i]=0;
  if(rott==0)
  {
    a1=(*r)[4];a3=(*r)[" CD "];a4=(*r)[" N  "];
    int nt=0; 
    for(res_temp=TRES.findrotamer("sidechain")->chn->isres(r->name,0)->more;res_temp;res_temp=res_temp->more)
    {
      //if(cbeta==0)side.hook(res_temp,r,4);
      //else side.hook(res_temp,r);
      if(nt>30) break;
      //res_temp->transfer(xyz0,0);
      side.hook(res_temp,r);
      if(r->name=='P') {	
      	a2=(*res_temp)[4];
      	x=a1->dihedral(0);
      	y=a2->dihedral(0);
      	y=y-x;
      	side.rotate(a1,y);
      	if(a3&&a4&&r->name=='P'&&TRES.distance(a3,a4)>1.8) {
		//res_temp->transfer(xyz0,1);
		continue;
	}
      }
      for(i=0;i<3;i++) xyzcb[i]+=a1->xyz[i];
      numcb++;
      r->transfer(xyz,1);
      //res_temp->transfer(xyz0,1);
      nt++;
    }
  }  
  r->transfer(xyz,1);
  if(numcb==0) return;
  for(i=0;i<3;i++) xyzcb[i]=xyzcb[i]/numcb;
  if(TRES.logg>3)  cerr<<r->name<<r->id0<<" "<<TRES.distance(a1->xyz,xyzcb)<<endl;
  for(i=0;i<3;i++) a1->xyz[i]=xyzcb[i]; 
  if(TRES.logg>3)  pdb->write("cb.pdb");
  return;
}
