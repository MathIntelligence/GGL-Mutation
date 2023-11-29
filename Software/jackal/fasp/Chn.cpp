#include"source.h"

void Chn::initial() {
  ishet=0;
  more=0;
  res=0;
  id=' ';
  number=0;
  next=0;
  pdb=0;
  area=0;
  start=0;
  nemp=0;
  temp=0;
  energy=0;
  seqcard=0;
}

Chn::Chn()
{
  initial();
}

Chn::Chn(Chn *s)
{
 initial();
 number=s->number;
 start=s->start;
 id=s->id;
 ishet=s->ishet;
 if(s->seqcard) seqcard=strdup(s->seqcard);
 if(s->res) {
	res=s->res->cloneres();
 }
 Res *r,*s1;
 for(r=res;r;r=r->next) for(s1=r;s1;s1=s1->more) s1->chn=this;
}

Chn* Chn::clonechain() 
{
 Chn *t=new Chn(this);
 if(next) t->next=next->clonechain();
 else t->next=0;
 return t;
} 

Chn::~Chn() 
{
if(res)  {delete res;res=0;}
if(next) {delete next;next=0;}
if(seqcard) {delete [] seqcard;seqcard=0;}
if(temp) delete [] temp;temp=0;
if(more) {delete more;more=0;}
}
 
void Chn::write(FILE *fp, Res *from, int to)
{
  Res *r;
  for(r=from;r;r=r->next)
  {
    if(r->id>to) break;
    r->write(fp);
  }
  fflush(fp);
} 

void Chn::write(FILE *fp)
{
  Res *res_temp;
  for(res_temp=res;res_temp;res_temp=res_temp->next)
  res_temp->write(fp);
}

void Chn::writeold(char *f)
{
	if(f==0) return;
	FILE *fp=fopen(f,"w");
	if(fp==0) {
		cerr<<"could not write to the disk"<<endl;
		return;
	}
	writeold(fp);
	fclose(fp);
}

void Chn::writeold(FILE *fp)
{
  Res *res_temp;
  for(res_temp=res;res_temp;res_temp=res_temp->next)
  res_temp->writeold(fp);
}
void Chn::writerescard(FILE *fp)
{
  Res *res_temp;
  for(res_temp=res;res_temp;res_temp=res_temp->next)
  res_temp->writerescard(fp);
}
void Chn::writeoldmore(FILE *fp)
{
  Res *res_temp;
  for(res_temp=res;res_temp;res_temp=res_temp->next)
  res_temp->writeoldmore(fp);
}

Atm *Chn::getatmbyoid(int n) {

  Res *res_temp;
  for(res_temp=res;res_temp;res_temp=res_temp->next)
  for(Atm *a=res_temp->atm;a;a=a->next) {
	if(a->oid==n) return a;
  }
  return 0;
}

void Chn::setoid(int n) {
  Res *res_temp;
  for(res_temp=res;res_temp;res_temp=res_temp->next) res_temp->setoid(n);
}

void Chn::write(char *s, Res *from, int to)
{
  FILE *fp;
  fp=fopen(s,"w");
  write(fp,from,to);
  fflush(fp);
  fclose(fp);

}

void Chn::write(char *s,Res *from,int to,int flg)
{
  FILE *fp;
  Res *r;
  Atm *a;

  fp=fopen(s,"w");
  for(r=from;r;r=r->next)
  {
    if(r->id0>to) break;
    for(a=r->atm;a;a=a->next)
    {
      if(flg==0) {a->write(fp);continue;}
      if(flg==1)
      {
         if(a->tatm->isbackbone()||a->tatm->id==4)
         {
           if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id==1) continue;
           a->write(fp);
           continue;
         }
      }
    }
  }

  fflush(fp);
  fclose(fp);
} 

void Chn::write(char *s)
{
  FILE *fp;
  fp=fopen(s,"w");
  write(fp);
  fflush(fp);
  fclose(fp);
}

void Chn::write(FILE *fp,char c)
{
  Res *s;
  for(s=res;s;s=s->next)
  if(s->name==c) s->write(fp);
  fflush(fp);
}

Res *Chn::operator[](int resn)
{
 Res *res_temp;
 for(res_temp=res;res_temp;res_temp=res_temp->next)
 if(res_temp->id0==resn) return res_temp; 
 return 0;
}

void Chn::dihedral() {

        FILE *fp=0;
        dihedral(fp);

}

void Chn::dihedral(char *s) {

	FILE *fp=fopen(s,"w");
	dihedral(fp);

}
void Chn::dihedral(FILE *fp)
{
Res *res_temp;
for(res_temp=res;res_temp;res_temp=res_temp->next)
res_temp->dihedral(fp);
}

void Chn::dihedraloption(FILE *fp,int m)
{
Res *res_temp;
for(res_temp=res;res_temp;res_temp=res_temp->next)
res_temp->dihedraloption(fp,m);
}

void Chn::dihedral(Res *r,int n) {

Res *res_temp;
for(res_temp=r;res_temp&&res_temp->id0<=n;res_temp=res_temp->next)
res_temp->dihedral((FILE *)0);	

}
void Chn::setbackbonetorsion(Res *s,int n) {

	Res *r;
	for(r=s;r&&r->id0<=n;r=r->next) r->setbackbonetorsion();

}
void Chn::dihedral(FILE *fp,int m)
{
Res *res_temp;
for(res_temp=res;res_temp;res_temp=res_temp->next)
res_temp->dihedral(fp,m);
}


void Chn::dihedral(FILE *fp,char c)
{
Res *res_temp;
for(res_temp=res;res_temp;res_temp=res_temp->next)
if(res_temp->name==c)res_temp->dihedral(fp);
}

void Chn::dihedral(FILE *fp,char c,int m)
{
Res *res_temp;
for(res_temp=res;res_temp;res_temp=res_temp->next)
if(res_temp->name==c)res_temp->dihedral(fp,m);
}

Atm *Chn::isatm(int i)
{
Res *r1;
Atm *a1;
for(r1=res;r1;r1=r1->next)
for(a1=r1->atm;a1;a1=a1->next)
if(a1->id0==i) return a1;
return 0;
}

void Chn::transfer(Res *s1,Res *s2,int n)
{
  Res *a1,*a2;
  int m;

  if(s1->name!=s2->name) 
  {
    cerr<<"could not copy between these two segments"<<endl;
  }

  a1=s1;a2=s2;m=0;
  while(a1&&a2)
  {
    if(a1->name!=a2->name) 
    {
     cerr<<"could not copy between these two segments"<<endl;
    }
    else
    {
       a1->transfer(a2); 
    }
    a1=a1->next;a2=a2->next;
    m++;
    if(m==n) break;
  }
}

void Chn::transfer(Chn *s)
{
Res *r1,*r2;
for(r1=res;r1;r1=r1->next)
{
  r2=(*s)[r1->id0];
  if(r1==0||r2==0) continue;
  if(r1->name!=r2->name) continue;
  r1->transfer(r2); 
}
}

int Chn::good(int b,int e,int m)
{
int i;
Res *s;
i=0;
for(s=res;s;s=s->next)
{
i+=s->good(b,e,m);
}
return i;
}


int Chn::ssbond()
{
Res *r,*r1;
Atm *a,*a1,*a2;
int i;
float d;
i=0;
for(r=res;r;r=r->next) 
{
if(r->name!='C') continue;
a=(*r)[" SG "];
if(a==0) continue;
for(r1=r->next;r1;r1=r1->next)
{
  if(r1->name!='C') continue;
  if(r1->id-r->id<2) continue;
  a1=(*r1)[" SG "]; 
  if(a1==0) continue;
  d=TRES.distance(a->xyz,a1->xyz);
  if(d>3||d<1) continue; 
  a2=a->bond[0];
  if(a2==0) continue;
  d=TRES.distance(a2->xyz,a1->xyz);
  if(d>4||d<2) continue;
  a2=a1->bond[0];
  if(a2==0) continue; 
  d=TRES.distance(a2->xyz,a->xyz);
  if(d>4||d<2) continue;
  a->bond[2]=a1;
  i++;
}
}
return i;
}

Res *Chn::isres(int i)
{
Res *r;
for(r=res;r;r=r->next)if(r->id==i) return r;
return 0;
}

Res *Chn::isres0(int i)
{
Res *r;
for(r=res;r;r=r->next)if(r->id0==i) return r;
return 0;
}

Res *Chn::isresoid(int i)
{
Res *r;
for(r=res;r;r=r->next)if(r->oid==i) return r;
return 0;
}

Res *Chn::findsmallres(int i)
{
int j;
Res *r;
for(j=i;j>=0;j--) {
	r=isresoid(j);
	if(r) return r;
}
return 0;
}



Res *Chn::isres(char c,int i)
{
  Res *r;
  int j;
  j=0;
  for(r=res;r;r=r->next)
  {
    if(r->name==c) 
    {
       if(i==j) return r;
       j++;
    }
  }
  return 0;
}

void Chn::surface(int m,float probe)
{
Tres *tr;
Tatm *ta;
Res *r;
Atm *a;
float *arer,*polar,xyz[3];
float t1,t2,f1,f2,c3,*temp,d;
int n,i,j,k,kk,ii,*covered,*order;
Map map;
Lattice lat;
Qsort cc;
float coo[3];

if(m==0) m=10;
polar=new float[3*m*2*m]; 
arer=new float[m*2*m];
covered=new int[m*2*m];

i=0;
for(r=res;r;r=r->next)
for(a=r->atm;a;a=a->next)
i++;

temp=new float[i];
order=new int[i];

area=0;

//add the probe to atoms radius;

for(tr=&TRES;tr;tr=tr->next)
for(ta=tr->tatm;ta;ta=ta->next)
{
  ta->eng->radius=ta->eng->radius+probe;
}

t1=3.1415926/m;
t2=2*3.1415926/(2*m);
c3=t1*t2;

k=0;
for(i=0;i<m;i++)
for(j=0;j<2*m;j++)
{
   f1=(i+0.5)*t1;f2=j*t2;
   polar[k*3+0]=sin(f1)*cos(f2);
   polar[k*3+1]=sin(f1)*sin(f2);
   polar[k*3+2]=cos(f1);
   arer[k]=c3*sin(f1);
   k++;
}
ii=k;

//set up the lattice
lat.putoff();
lat.grdsiz=2.0;
lat.ready(pdb);

lat.grdsiz=2;
lat.flag=1; 
lat.puton();//lat.puton(this,0,100);

map.grdsiz=2*probe;
map.clear();
map.ready(pdb);
map.puton();
for(r=res;r;r=r->next)
{
 r->area=0;
 for(a=r->atm;a;a=a->next)
 {
   if(a->temp){delete [] a->temp;a->temp=0;a->nemp=0;}
   a->area=0;
   if(a->flag==0) {
	//cerr<<"flag is zero..."<<a->name<<" "<<a->res->name<<a->res->id<<endl;
	continue;
   }
   if(map.cover(a)) continue;
   lat.getcell(a,a->tatm->eng->radius+lat.radmax); 

   for(i=0;i<lat.nget;i++)
   {
     temp[i]=TRES.distance(a->xyz,lat.obtain[i]->atm->xyz)-
             lat.obtain[i]->atm->tatm->eng->radius; 
   } 
   cc.sort(temp,lat.nget,order);
   for(i=0;i<lat.nget;i++) 
   if(temp[i]>a->tatm->eng->radius)  
   {
     lat.nget=i;break;
   }

   for(i=0;i<ii;i++)covered[i]=0;
   for(i=0;i<3;i++) coo[i]=0;

   kk=-1;
   for(i=0;i<m;i++)
   for(j=0;j<2*m;j++)
   {  
       kk++;
       xyz[0]=polar[kk*3+0]*a->tatm->eng->radius+a->xyz[0];
       xyz[1]=polar[kk*3+1]*a->tatm->eng->radius+a->xyz[1];
       xyz[2]=polar[kk*3+2]*a->tatm->eng->radius+a->xyz[2];
       
       for(k=0;k<lat.nget;k++)
       {
         n=order[k];
         d=lat.obtain[n]->atm->tatm->eng->radius;
         d=d*d;
         if(TRES.distsqr(xyz,lat.obtain[n]->atm->xyz)<d)
         {
           covered[kk]=1;break;     
         }
       }
       if(covered[kk]==1) continue;
       d=arer[kk]*a->tatm->eng->radius*a->tatm->eng->radius;
       coo[0]+=d*xyz[0];
       coo[1]+=d*xyz[1];
       coo[2]+=d*xyz[2];
       a->area+=d;
       //cerr<<a->name<<a->area<<endl;
   }
   if(a->area>0)
   {
       for(i=0;i<3;i++) coo[i]=coo[i]/a->area;
       d=TRES.distance(coo,a->xyz)/a->tatm->eng->radius;
       for(i=0;i<3;i++) coo[i]=a->xyz[i]+(coo[i]-a->xyz[i])/d;
       a->nemp=7;a->temp=new float[a->nemp]; 
       for(i=0;i<a->nemp;i++) a->temp[i]=0;
       for(i=0;i<3;i++) a->temp[i]=coo[i];
       for(i=0;i<3;i++) a->temp[i+3]=(coo[i]-a->xyz[i])/a->tatm->eng->radius;
   }
   r->area+=a->area;
 }
 area+=r->area;
 
}

//set back the radius
 
lat.flag=0;

for(tr=&TRES;tr;tr=tr->next)
for(ta=tr->tatm;ta;ta=ta->next)
ta->eng->radius-=probe;
delete [] arer;
delete [] polar;
delete [] temp;
delete [] order;
delete [] covered;
}
void Chn::surfaceold(int m,float probe)
{
Tres *tr;
Tatm *ta;
Res *r;
Atm *a;
float *arer,*polar,xyz[3];
float t1,t2,f1,f2,c3,*temp,d;
int n,i,j,k,kk,ii,*covered,*order;
Map map;
Lattice lat;
Qsort cc;
float coo[3];

if(m==0) m=10;
polar=new float[3*m*2*m]; 
arer=new float[m*2*m];
covered=new int[m*2*m];

i=0;
for(r=res;r;r=r->next)
for(a=r->atm;a;a=a->next)
i++;

temp=new float[i];
order=new int[i];

area=0;

//add the probe to atoms radius;

for(tr=&TRES;tr;tr=tr->next)
for(ta=tr->tatm;ta;ta=ta->next)
{
  ta->eng->radius=ta->eng->radius+probe;
}

t1=3.1415926/m;
t2=2*3.1415926/(2*m);
c3=t1*t2;

k=0;
for(i=0;i<m;i++)
for(j=0;j<2*m;j++)
{
   f1=(i+0.5)*t1;f2=j*t2;
   polar[k*3+0]=sin(f1)*cos(f2);
   polar[k*3+1]=sin(f1)*sin(f2);
   polar[k*3+2]=cos(f1);
   arer[k]=c3*sin(f1);
   k++;
}
ii=k;

//set up the lattice
lat.putoff();
lat.grdsiz=2.0;
lat.ready(pdb);

lat.grdsiz=2;
lat.flag=1; 
lat.puton(this,0,100);

map.grdsiz=2*probe;
map.clear();
map.ready(pdb);
map.puton();
for(r=res;r;r=r->next)
{
 r->area=0;
 for(a=r->atm;a;a=a->next)
 {
   if(a->temp){delete [] a->temp;a->temp=0;a->nemp=0;}
   a->area=0;
   if(a->flag==0) {
	//cerr<<"flag is zero..."<<a->name<<" "<<a->res->name<<a->res->id<<endl;
	continue;
   }
   if(map.cover(a)) continue;
   lat.getcell(a,a->tatm->eng->radius+lat.radmax); 

   for(i=0;i<lat.nget;i++)
   {
     temp[i]=TRES.distance(a->xyz,lat.obtain[i]->atm->xyz)-
             lat.obtain[i]->atm->tatm->eng->radius; 
   } 
   cc.sort(temp,lat.nget,order);
   for(i=0;i<lat.nget;i++) 
   if(temp[i]>a->tatm->eng->radius)  
   {
     lat.nget=i;break;
   }

   for(i=0;i<ii;i++)covered[i]=0;
   for(i=0;i<3;i++) coo[i]=0;

   kk=-1;
   for(i=0;i<m;i++)
   for(j=0;j<2*m;j++)
   {  
       kk++;
       xyz[0]=polar[kk*3+0]*a->tatm->eng->radius+a->xyz[0];
       xyz[1]=polar[kk*3+1]*a->tatm->eng->radius+a->xyz[1];
       xyz[2]=polar[kk*3+2]*a->tatm->eng->radius+a->xyz[2];
       
       for(k=0;k<lat.nget;k++)
       {
         n=order[k];
         d=lat.obtain[n]->atm->tatm->eng->radius;
         d=d*d;
         if(TRES.distsqr(xyz,lat.obtain[n]->atm->xyz)<d)
         {
           covered[kk]=1;break;     
         }
       }
       if(covered[kk]==1) continue;
       d=arer[kk]*a->tatm->eng->radius*a->tatm->eng->radius;
       coo[0]+=d*xyz[0];
       coo[1]+=d*xyz[1];
       coo[2]+=d*xyz[2];
       a->area+=d;
       //cerr<<a->name<<a->area<<endl;
   }
   if(a->area>0)
   {
       for(i=0;i<3;i++) coo[i]=coo[i]/a->area;
       d=TRES.distance(coo,a->xyz)/a->tatm->eng->radius;
       for(i=0;i<3;i++) coo[i]=a->xyz[i]+(coo[i]-a->xyz[i])/d;
       a->nemp=7;a->temp=new float[a->nemp]; 
       for(i=0;i<a->nemp;i++) a->temp[i]=0;
       for(i=0;i<3;i++) a->temp[i]=coo[i];
       for(i=0;i<3;i++) a->temp[i+3]=(coo[i]-a->xyz[i])/a->tatm->eng->radius;
   }
   r->area+=a->area;
 }
 area+=r->area;
 
}

//set back the radius
 
lat.flag=0;

for(tr=&TRES;tr;tr=tr->next)
for(ta=tr->tatm;ta;ta=ta->next)
ta->eng->radius-=probe;
delete [] arer;
delete [] polar;
delete [] temp;
delete [] order;
delete [] covered;
}

void Chn::surface(char *ss,int m,float probe)
{
Res *r;
Atm *a;
FILE *fp;
fp=fopen(ss,"w");
surface(m,probe);
for(r=res;r;r=r->next) 
for(a=r->atm;a;a=a->next)
fprintf(fp,"%s %5i %s %5i %8.3f %8.3f\n",r->tres->name3,r->id,
           a->name,a->id0,a->area,r->area);
fprintf(fp,"total area:  %8.3f\n",area);
fflush(fp);
}

void Chn::clearhbond() {

	Res *r;
	for(r=res;r;r=r->next){
		r->sec='-';
		if(r->hbond) delete r->hbond;
		r->hbond=0;
	}

}

void Chn::buildhbond(){

	//set up lattice model


	Lattice *lattice;

	lattice=new Lattice;
	lattice->flag=0;
	lattice->grdsiz=2.0;
	lattice->radall=15.;
	lattice->putoff();
	lattice->ready(pdb);
	lattice->puton(this);
		
	//find hbond list
	Res *r;
	Atm *a,*s;
	int i;
	float e;

	for(r=res;r;r=r->next)
	for(a=r->atm;a;a=a->next){
		if(a->name[1]!='N'&&a->name[1]!='O') continue;
		lattice->getcell(a,5.);
		for(i=lattice->nget-1;i>=0;i--) {
			s=lattice->obtain[i]->atm;
			if(s->name[1]!='N'&&s->name[1]!='O') continue;
			if(s->res->id0<=a->res->id0) continue;
			if(a->dsspdefinedhbond(s)==0) continue;
		
			e=0;	
			r->addhbondlist(a,s,e,1);
			s->res->addhbondlist(a,s,e,1);		 
			/*
			if(a->tatm->hbond==1) {
				r->addhbondlist(a,s,e,1);
				s->res->addhbondlist(a,s,e,1);
			}
			else if(s->tatm->hbond==1) {
				r->addhbondlist(s,a,e,1);
				s->res->addhbondlist(s,a,e,1);
			}
			else if(a->tatm->hbond==2&&s->tatm->hbond==3) {
				r->addhbondlist(s,a,e,1);
				s->res->addhbondlist(s,a,e,1);
			}		
			else if(s->tatm->hbond==2&&a->tatm->hbond==3){
				r->addhbondlist(a,s,e,1);
				s->res->addhbondlist(a,s,e,1);
			}
			else if(s->tatm->hbond==3&&a->tatm->hbond==3){
                                r->addhbondlist(a,s,e,1);
                                s->res->addhbondlist(a,s,e,1);
                        }
			*/
		}
	}
	delete lattice;lattice=0;
}

void Chn::buildhbond(int from,int end){

	//set up lattice model


	Lattice *lattice;

	lattice=new Lattice;
	lattice->flag=0;
	lattice->grdsiz=2.0;
	lattice->radall=15.;
	lattice->putoff();
	lattice->ready(pdb);
	lattice->puton(this);
		
	//find hbond list
	Res *r;
	Atm *a,*s;
	int i;
	float e;

	for(r=res;r;r=r->next)
	for(a=r->atm;a;a=a->next){
		if(a->name[1]!='N'&&a->name[1]!='O') continue;
		lattice->getcell(a,5.);
		for(i=lattice->nget-1;i>=0;i--) {
			s=lattice->obtain[i]->atm;
			if(s->name[1]!='N'&&s->name[1]!='O') continue;
			if(s->res->id0<=a->res->id0) continue;
			if(a->dsspdefinedhbond(s)==0) continue;
			if(a->tatm->id<from||a->tatm->id>end) continue;
			if(s->tatm->id<from||s->tatm->id>end) continue;
			e=0;	
			r->addhbondlist(a,s,e,1);
			s->res->addhbondlist(a,s,e,1);		 
			/*
			if(a->tatm->hbond==1) {
				r->addhbondlist(a,s,e,1);
				s->res->addhbondlist(a,s,e,1);
			}
			else if(s->tatm->hbond==1) {
				r->addhbondlist(s,a,e,1);
				s->res->addhbondlist(s,a,e,1);
			}
			else if(a->tatm->hbond==2&&s->tatm->hbond==3) {
				r->addhbondlist(s,a,e,1);
				s->res->addhbondlist(s,a,e,1);
			}		
			else if(s->tatm->hbond==2&&a->tatm->hbond==3){
				r->addhbondlist(a,s,e,1);
				s->res->addhbondlist(a,s,e,1);
			}
			else if(s->tatm->hbond==3&&a->tatm->hbond==3){
                                r->addhbondlist(a,s,e,1);
                                s->res->addhbondlist(a,s,e,1);
                        }
			*/
		}
	}
	delete lattice;lattice=0;
}
void Chn::buildcatracehbond(){

	//set up lattice model


	Lattice *lattice;

	lattice=new Lattice;
	lattice->flag=0;
	lattice->grdsiz=2.0;
	lattice->radall=15.;
	lattice->putoff();
	lattice->ready(pdb);
	lattice->puton(this);
		
	//find hbond list
	Res *r;
	Atm *a,*s;
	int i;
	
	Atm *a1,*b1;
	DistPopular *hdist=TRES.popbin->gethbond();
	if(hdist==0) {
		return;
	}
	for(r=res;r;r=r->next) 
	for(a=r->atm;a;a=a->next){
		if(a->name[1]!='N'&&a->name[1]!='O') continue;
		lattice->getcell(a,6.);
		if(a->name[1]=='O') a1=a->bond[0];
		else		    a1=a;
		for(i=lattice->nget-1;i>=0;i--) {
			s=lattice->obtain[i]->atm;
			if(s->name[1]!='N'&&s->name[1]!='O') continue;
			if(s->res->id0<=a->res->id0) continue;
			//if(a->dsspdefinedhbond(s)==0) continue;
			if(s->name[1]=='O') b1=s->bond[0];
			else		    b1=s;
			
			DistPopular * p=hdist->getDistPopular(a1->res->tres,a1->name,b1->name,b1->res->tres->name);			 
			if(p==0) continue; 
			float e=TRES.distance(a1,b1);
			if(e>=p->dist[0]-0.5&&e<=p->dist[1]+0.5) {
				e=0;
				r->addhbondlist(a,s,e,1);
				s->res->addhbondlist(a,s,e,1);	
			}	 
		}
	}
}

void Chn::writehbondparm(){

        //set up lattice model

        Lattice *lattice;

        lattice=new Lattice;
        lattice->flag=0;
        lattice->grdsiz=2.0;
        lattice->radall=15.;
        lattice->putoff();
        lattice->ready(pdb);
        lattice->puton(this);
               
        //find hbond list
        Res *r;
        Atm *a,*s;
        int i;
        //float e;

        for(r=res;r;r=r->next)
        for(a=r->atm;a;a=a->next){
                if(a->name[1]!='N'&&a->name[1]!='O') continue;
                lattice->getcell(a,5.);
                for(i=lattice->nget-1;i>=0;i--) {
                        s=lattice->obtain[i]->atm;
                        if(s->name[1]!='N'&&s->name[1]!='O') continue;
                        if(s->res->id0<=a->res->id0) continue;
			if(a->dsspdefinedhbond(s)==0) continue;
                        a->writehbondparm(s);
                }
        }
}


void Chn::buildssbond(){

        //find hbond list
        Res *r,*r0;
        Atm *a,*s;
        //int j;
        //float e;

        for(r=res;r;r=r->next) {
		if(r->name!='C') continue;	
		a=r->isatm(" SG ");
		if(a==0) continue;
		for(r0=r->next;r0;r0=r0->next) {
			if(r0->name!='C') continue;
			
			s=r0->isatm(" SG ");
			if(s==0) continue;		
			//e=TRES.distance(a,s);
			//if(e<3.0) continue;
			if(isssbond(r,r0)==0) continue;
			r->addhbondlist(a,s,0,2);
			r0->addhbondlist(a,s,0,2);
		}
	}
}

void Chn::buildssbond(int from,int end){

        //find hbond list
        Res *r,*r0;
        Atm *a,*s;
        //int j;
        //float e;

        for(r=res;r;r=r->next) {
		if(r->name!='C') continue;	
		a=r->isatm(" SG ");
		if(a==0) continue;
		for(r0=r->next;r0;r0=r0->next) {
			if(r0->name!='C') continue;
			
			s=r0->isatm(" SG ");
			if(s==0) continue;		
			//e=TRES.distance(a,s);
			//if(e<3.0) continue;
			if(isssbond(r,r0)==0) continue;
			if(a->tatm->id<from||a->tatm->id>end) continue;
			if(s->tatm->id<from||s->tatm->id>end) continue;
			r->addhbondlist(a,s,0,2);
			r0->addhbondlist(a,s,0,2);
		}
	}
}


void Chn::secstr(FILE *fp)
{
Lattice *lattice;
Res *r,*r0,*r1;
Atm *a,*s;
int i,j;
lattice=new Lattice;
lattice->putoff();
lattice->flag=0;
lattice->grdsiz=2.;
lattice->radall=15.;
lattice->ready(pdb);

//put chain on lattice

for(r=res;r;r=r->next)
{
a=(*r)[3];
lattice->puton(a);
}

for(r=res;r;r=r->next)
{
 r->sec='-';
 for(a=r->atm;a;a=a->next)
 {
  a->hbond[0]=0;
  a->hbond[1]=0;
 }
}

for(r=res;r;r=r->next)
{ 
   a=r->atm;
   lattice->getcell(a,5.); 
   for(i=lattice->nget-1;i>=0;i--)
   {
     s=lattice->obtain[i]->atm;
     if(abs(s->res->id-a->res->id)<2) continue;
     if(a->ishbond(s,1.2,60,-1)>-0.001) j=0;
     else j=1;
     if(j==1)
     {
      if(a->hbond[0]==0) a->hbond[0]=s;
      else if(a->hbond[1]==0) a->hbond[1]=s;                            
     }
   } 
}

//decide the secondary structure

//detect helix
for(r=res;r;r=r->next)
{
  a=r->atm;
  if(a==0) continue;

  i=0;j=0;
  if(a->hbond[0])i=abs(a->hbond[0]->res->id-r->id);
  if(a->hbond[1])j=abs(a->hbond[1]->res->id-r->id);
  if(i==4||j==4) r->sec='h';
}

//connect break helix

r0=0;
for(r=res;r;r=r->next)
{
 if(r->next&&r0) 
 {
  if(r->sec=='-'&&r0->sec=='h'&&r->next->sec=='h') r->sec='h';
 }
 r0=r;
}

//kick out fake helix
r0=0;
i=0;
for(r=res;r;r=r->next)
{
 if(r->sec=='h') i++;
 else 
 {
   if(i<3)
   {   
     for(r1=r0;r1;r1=r1->next)
     {
       if(r1->id>r->id) break;
       r1->sec='-';
     }
   }
   i=0;
   r0=0;
 }
 if(i==1) r0=r;
}

//detecting beta sheet

for(r=res;r;r=r->next)
{
  if(r->sec=='h') continue;
  a=r->atm;
  if(a==0||a->tatm->id!=0) continue;//not the N atom
  for(i=0;i<2;i++)
  {
  if(a->hbond[i]) 
  {
    r1=a->hbond[i]->res;//bonded with N atom residues
    s=r1->atm;
    if(s&&r1->sec!='h'&&s->tatm->id==0)
    {
       if(s->hbond[0]&&s->hbond[0]->res->id-r->id==2)
       {
         r->sec='e';
         r1->sec='e';
       }
       if(s->hbond[1]&&s->hbond[1]->res->id-r->id==2)
       {
         r->sec='e';
         r1->sec='e';
       }
       if(s->hbond[0]&&s->hbond[0]->res->id==r->id)
       { 
         r->sec='e';
         r1->sec='e';
       }
       if(s->hbond[1]&&s->hbond[1]->res->id==r->id)
       { 
         r->sec='e';
         r1->sec='e';
       }
    }
    j=a->hbond[i]->res->id-2;
    r1=(*this)[j]; 
    if(r1&&r1->sec!='h'&&r1->atm&&r1->atm->tatm->id==0)
    {
      s=r1->atm->hbond[0];
      if(s&&s->res->id==r->id+2)
      {
        r->sec='e';r1->sec='e';
      }
      s=r1->atm->hbond[1];
      if(s&&s->res->id==r->id+2) 
      {
        r->sec='e';r1->sec='e';
      }
    }
  }
  }
} 

 
r0=0;
for(r=res;r;r=r->next)
{
 if(r->next&&r0)
 {
  if(r->sec=='-'&&r0->sec=='e'&&r->next->sec=='e') r->sec='e';
 }
 r0=r;
}

//kick out fake helix
r0=0;
i=0;
res->sec='-';
for(r=res;r;r=r->next)
{
 if(r->sec=='e') i++;
 else
 {
   if(i<3)
   {
     for(r1=r0;r1;r1=r1->next)
     {
       if(r1->id>r->id) break;
       r1->sec='-';
     }
   }
   i=0;
   r0=0;
 }
 if(i==1) r0=r;
 if(r->next==0) r->sec='-';
}


if(fp)
{
  for(r=res;r;r=r->next)
  {
    fprintf(fp,"%c %5i %c\n",r->name,r->id,r->sec);
  }
}
fflush(fp);
delete lattice;
}

char *Chn::getseq()
{
int j;
Res *r;
int num;
char *chain_seq;
int xx;

//determine the number of residues including nonexist.

setflg(res,10000000,1);
surface(10,1.4);

r=res;
num=0;
for(r=res;r;r=r->next) num=r->id;
num++;
chain_seq=new char[num*3+3];

for(j=0;j<num*3+3;j++) chain_seq[j]='?';

for(r=res;r;r=r->next)
{
chain_seq[r->id]=r->name;
chain_seq[num+r->id+1]=r->sec;
if(chain_seq[num+r->id+1]==' ')chain_seq[num+r->id+1]='-';
xx=(int)(r->area/r->tres->area*100);
xx=xx/10;
chain_seq[num*2+2+r->id]='0'+xx;
}
chain_seq[num]='\0';
chain_seq[2*num+1]='\0';
chain_seq[3*num+2]='\0';
return chain_seq;
}

int Chn::manyres(){
	int n=0;
	for(Res *r=res;r;r=r->next)n++;
	return n;
}
char *Chn::get2d()
{
	
	int n=manyres();
	
	char *out=new char[n+1];
	n=0;
	for(Res *r=res;r;r=r->next) {
		out[n++]=r->sec;
	}
	out[n]='\0';
	return out;
}

int Chn::getseq(char *s)
{
Res *r;
int i; 
i=0;
for(r=res;r;r=r->next)
{
s[i++]=r->name;
}
s[i]='\0';
return i;
}

char *Chn::getseqn(Res *tt,int n){
Res *r;
int i;
char *s;
s=new char[n-tt->id0+100];
i=0;
for(r=tt;r;r=r->next)
{
if(r->id0>n) break;
s[i++]=r->name;
}
s[i]='\0';
return s;
}

char *Chn::getseqn() {
Res *r;
int i;
char *s;
r=lastres();
if(r==0) return 0;
s=new char[r->id0+100];
i=0;
for(r=res;r;r=r->next)
{
s[i++]=r->name;
}
s[i]='\0';
return s;
}
void Chn::create(char *seqn) //create phsuedo pdb
{
int i,k,j;
Res *res_temp,*r1;
Tres *tres_temp;
char ch;
Rotate rot;

if(seqn==0) return;

k=strlen(seqn);
if(k==0) return;
if(res) {delete res;res=0;}
number=0;
for(i=0;i<k;i++)
{
  ch=seqn[i];
  tres_temp=TRES[ch]; 
  if(tres_temp==0) cerr<<"Warning! no suitable residue topology found for:"<< ch<<endl;
  
  if(!res){res=new Res(tres_temp); res_temp=res;}
  else  {res_temp->next=new Res(tres_temp); res_temp=res_temp->next;}

  res_temp->temp=new float[9];
  for(j=0;j<9;j++) res_temp->temp[j]=tres_temp->head[j];
  res_temp->nemp=9;
  res_temp->chn=this;
  res_temp->id=i;
}
configure();
// connext the chain
  for(res_temp=res;res_temp;res_temp=res_temp->next)
  {
     r1=res_temp->next;
     if(r1)
     {
       rot.link(res_temp,r1,r1->id0,1);
     }
  }

}

void Chn::create(Tres *tres)
{
 Tres *a;
 Res *r;

 r=0;
 for(a=tres;a;a=a->next)
 {
   if(r==0) { res=new Res(a); r=res; }
   else { r->next=new Res(a); r=r->next;} 

   r->chn=this;
   r->id=a->id;
   r->temp=new float[9];
   TRES.copy(a->head,r->temp,9);
   r->nemp=9;
   r->chn=this;
 }
}

void Chn::create(Res *r,int n2)
{ 
Res *r1,*r2;
if(res){delete res;res=0;}
r2=0;
for(r1=r;r1;r1=r1->next)
{
  //if(r1->id0>n2) continue;
  if(r1->id0>n2) break;
  if(!r2)
  {
    r2=new Res(r1);     
    r2->chn=this;
    res=r2;
  }
  else
  {
    r2->next=new Res(r1);
    r2=r2->next;    
    r2->chn=this;
  }
}

}

void Chn::create(Chn *s)
{
create(s->res,10000);
}

void Chn::transfer(float s)
{
Res *r;
for(r=res;r;r=r->next)r->transfer(s);
}

void Chn::transfer(float *s,int f)
{
Res *r;
int i;
i=0;
for(r=res;r;r=r->next)
{
 r->transfer(s+i,f);
 i+=r->tres->number*3;
}

}

void Chn::transfer(Res *rr,float *s,int f)
{

Res *r;
int i;
i=0;
for(r=res;r;r=r->next)
{
 if(rr->id0==r->id0) {r->transfer(s+i,f);return;}
 i+=r->tres->number*3;
}
} 

void Chn::transfer(float *s,Res *from,int to, int f)
{
Res *r;
int i;
i=0;
for(r=from;r;r=r->next)
{
 if(r->id0>to) break;
 r->transfer(s+i,f);
 i+=r->tres->number*3;
}
}

void Chn::transfertemp(float *s,Res *from,int to, int f)
{
Res *r;
int i;
i=0;
for(r=from;r;r=r->next)
{
 if(r->id0>to) break;
 r->transfertemp(s+i,f);
 i+=9;
}
}

float *Chn::gettransfertemp(Res *from,int to)
{
float *out=new float[(to-from->id0+100)*9];
Res *r;
int i;
i=0;
for(r=from;r;r=r->next)
{
 if(r->id0>to) break;
 r->transfertemp(out+i,0);
 i+=9;
}
return out;
}


void Chn::copytemp(float *s,int f)
{
Res *r;
int i,j;
i=0;j=0;
for(r=res;r;r=r->next)
{
  for(i=0;i<r->nemp;i++)
  {
    if(f==1)r->temp[i]=s[j+i];
    else s[j+i]=r->temp[i];
  }
  j+=r->nemp; 
}
}

void Chn::bondleng()
{
  Res *r;
  for(r=res;r;r=r->next)
  {
    r->bondleng();
  }
}

void Chn::transform(int flg)
{
  transform(res,10000000,flg);
}

void Chn::transform(Res *ss,int nn,int flg)
{
  Res *r,*rr;
  Atm *a,*b;
  int i;
 
  setflg(res,10000000,1);

  for(rr=ss;rr&&rr->id0<=nn;rr=rr->next)
  for(r=rr;r;r=r->more)
  {
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
      else if(flg==7) { //only backbone atoms
	  if(strcmp(a->tatm->name," CA ")!=0) a->flag=0;
      }
      else if(flg==8)
      {
          if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>0) a->flag=0;
          else if(a->tatm->name[1]!='H'&&a->tatm->id>4) a->flag=0; //no backbone
      }

    }
  }

  //rebonding

  for(rr=ss;rr&&rr->id0<=nn;rr=rr->next)
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
  
 
  for(rr=ss;rr&&rr->id0<=nn;rr=rr->next)
  for(r=rr;r;r=r->more)
  {
    b=r->atm;a=r->atm;
    while(a)
    {
      if(a->flag==0&&a!=r->atm)
      {
        b->next=a->next;
        a->next=0;
        delete a;    
        a=b->next;
        continue;    
      }
      else if(a->flag==0&&a==r->atm)
      {
        r->atm=a->next;
	b=r->atm;
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


void Chn::transform(int flg,int g)
{
  Res *r,*rr;
  Atm *a,*b;
  int i;
 
  int nn=1000000;
  setflg(res,10000000,1);
  
  for(rr=res;rr&&rr->id0<=nn;rr=rr->next)
  for(r=rr;r;r=r->more)
  {
    if(r->flag==g) continue;
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
      else if(flg==7) {
	 if(a->tatm->id!=1) a->flag=0;
      }
    }
  }

  //rebonding

  for(rr=res;rr&&rr->id0<=nn;rr=rr->next)
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
  
 
  for(rr=res;rr&&rr->id0<=nn;rr=rr->next)
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




void Chn::header()
{
  header(res,100000);
}


void Chn::header(int n1,int n2){
  Res *r1;

  Res *r;

  r1=0;
  for(r=res;r;r=r->next) {
	if(r->id0>=n1) {
		r1=r;
		break;
	}

  }
  if(r1) header(r1,n2);
}

void Chn::header(Res *s,int n)
{
Res *r1,*r2;
Atm *atm,*aa;
float d;
int i;
Rotate rot;
//copy the N and HN file

if(s==0) return;

for(r2=s;r2;r2=r2->next)
{
 if(r2->id0>n) break;
 if(r2->temp) {delete [] r2->temp; r2->temp=0;}
 r2->temp=new float[9];
 for(i=0;i<9;i++) r2->temp[i]=999999;
 r2->nemp=9;
}


r1=s->chn->isres0(s->id0-1);
for(r2=s;r2;r2=r2->next)
{ 
  if(r2->id0>n) break;
  if(r1==0) {
	rot.hookheader(r2);
	r1=r2;
	continue;
  }
  aa=(*r2)[" N  "];
  if(aa==0) 
  {
    //cerr<<"could not header at atom N"<<r2->name<<r2->id<<" ";
    //cerr<<r2->chn->pdb->name<<endl;
    r1=r2;
    rot.hookheader(r2);
    continue;
  }
  if(r1->id-r2->id!=-1) 
  {
    //cerr<<"could not header discontinuous"<<r2->name<<r2->id<<" ";
    //cerr<<r2->chn->pdb->name<<endl;
    rot.hookheader(r2);
    r1=r2;
    continue;
  }
  atm=(*r1)[" CA "];
  if(atm==0) 
  {
   //cerr<<"could not header at atom CA"<<r1->name<<r1->id<<" ";
   //cerr<<r1->chn->pdb->name<<" "<<endl;
   rot.hookheader(r2);
   r1=r2;
   continue;
  }
  TRES.copy(atm->xyz,r2->temp,3);
  atm=(*r1)[" C  "];
  if(atm==0) 
  {
    //cerr<<"could not header at atom N"<<r2->name<<r2->id<<" ";
    //cerr<<r2->chn->pdb->name<<endl;
    rot.hookheader(r2);
    r1=r2;
    continue;
  }
  TRES.copy(atm->xyz,r2->temp+3,3);
  atm=(*r1)[" O  "];
  if(atm==0) 
  {
    //cerr<<"could not header at HN"<<r2->name<<r2->id<<" ";
    //cerr<<r2->chn->pdb->name<<endl;
    rot.hookheader(r2);
    r1=r2;
    continue;
  }
  TRES.copy(atm->xyz,r2->temp+6,3);
  d=TRES.distance(aa->xyz,r2->temp+3);
  if(d>1.7)
  {
    //cerr<<"could not header at distance too large"<<r2->name<<r2->id<<" ";
    //cerr<<r2->chn->pdb->name<<endl;
    rot.hookheader(r2);
    r1=r2;
    continue;
  } 
  r1=r2;
}
}

 
void Chn::headerpdbfix(Res *s,int n)
{
  header(s,n);
  Res *r2;
  for(r2=s;r2;r2=r2->next) {
	int i,n;
	n=0;
	for(i=0;i<4;i++) {
		Atm *a=r2->isatmid(i);
		if(a) n++;
	}
	if(n!=4) {//backbone not complete
		if(r2->temp) delete [] r2->temp;
		r2->temp=0;r2->nemp=0;
	}
  }
}
 


void Chn::setflg(int flg) {
 setflg(res,1000000,flg);
}
void Chn::setflg(Res *s,int n,int flg)
{
  Res *r;
  Atm *a; 
  for(r=s;r;r=r->next)
  {
    if(r->id0>n) return;
    for(a=r->atm;a;a=a->next)
    {
      a->flag=flg;
    }
  }
}

void Chn::setflg(int *give,Res *s,int n,int flg)
{
  Res *r;
  Atm *a;
  int i;
  i=0;
  for(r=s;r;r=r->next)
  {
    if(r->id0>n) return;
    for(a=r->atm;a;a=a->next)
    { 
      if(flg==0) give[i]=a->flag;
      else       a->flag=give[i];
      i++;
    }
  }
}

void Chn::setflgr(int flg) {

	setflgr(res,1000000,flg);

}
void Chn::setflgr(Res *s,int n,int flg)
{
  Res *r;
  for(r=s;r;r=r->next)
  {
    if(r->id0>n) return;
    r->flag=flg;
  }
}


void Chn::setflgr(char c,int flg)
{
  Res *r;
  for(r=res;r;r=r->next)
  {
    if(r->name==c) r->flag=flg;
  }
}

void Chn::setflgr(int *give,Res *s,int n,int flg)
{
  Res *r;
  int i;
  i=0;
  for(r=s;r;r=r->next)
  {
    if(r->id0>n) return;
    if(flg==0) give[i]=(int)(r->flag);
    else       r->flag=give[i];
    i++;
  }
}

int Chn::setflg(int *give,Res *s,int n,int val,int flg)
{
  Res *r;
  Atm *a;
  int i;
  for(r=s;r;r=r->next)
  {
    if(r->id0>n) return i;
    if(flg==0) 
    {
      give[r->id0]=val;i++;
    }
    else
    {
      for(a=r->atm;a;a=a->next) {give[a->id0]=val;i++;}
    }
  }
  return i;
}

void Chn::onlybc()
{
  Res *r;
  Atm *a,*a1,*aa[50];
  Tatm *ta;
  int i,n;

  for(r=res;r;r=r->next)
  {
    if(r->flag<-9999) continue;
    for(i=0;i<50;i++) aa[i]=0;
    i=0;
    for(ta=r->tres->tatm;ta;ta=ta->next)
    {
       a=(*r)[ta->id];
       if(a==0)
       {
         if(ta->isbackbone()) continue;
         a1=new Atm;
         strcpy(a1->name,ta->name);
         a1->res=r;
         a1->tatm=ta;
         aa[i++]=a1; 
       } 
       else 
       {
         aa[i++]=a;
       }
    }

    for(n=0;n<i;n++)
    {
      if(n==0) r->atm=aa[n];
      else  aa[n-1]->next=aa[n];
    }
    aa[n-1]->next=0;
  }

  //configure();
}

void Chn::setoff(float *xyz,float cut) 
{
  Res *r;
  Atm *a;
  float d;
  for(r=res;r;r=r->next)
  for(a=r->atm;a;a=a->next)
  {
    d=TRES.distance(a->xyz,xyz);
    if(d<cut) a->flag=1;
  }
}
void Chn::setoff(Res *r0,int n,float cut)
{
 Res *r;
 Atm *a;
 Lattice lat;
 int i;
 float d;

 lat.putoff();
 lat.grdsiz=2.0;
 lat.ready(pdb);
 lat.grdsiz=2;
 lat.radall=15;
 lat.flag=0;
 lat.puton(this);
 lat.putoff(r0,n);
 setflg(res,1000000,0);
 setflg(r0,n,1);

 for(r=r0;r;r=r->next)
 {
   if(r->id0>n) break;
   for(a=r->atm;a;a=a->next)
   {
     lat.getcell(a,cut);
     for(i=0;i<lat.nget;i++)
     {
       d=TRES.distance(a->xyz,lat.obtain[i]->atm->xyz);
       if(d>cut) continue;
       lat.putoff(lat.obtain[i]->atm);
       lat.obtain[i]->atm->flag=1;
     }
   }
 }
}

void Chn::configure()
{
int i,k,n,j,m,nn;
Chn *chn_temp;
Res *res_temp,*res_temp1;
Tatm *tatm_temp;
Atm *atm_temp,*atm_temp1;
Atm *atm[30];
 
//put the atom on order; re-order it according to tres
chn_temp=this;
for(res_temp=chn_temp->res;res_temp;res_temp=res_temp->next)
{
  for(i=0;i<30;i++)atm[i]=0;     
  for(atm_temp=res_temp->atm;atm_temp;atm_temp=atm_temp->next)
  { 
    atm[atm_temp->tatm->id]=atm_temp;
  }
  j=0;
  for(i=0;i<res_temp->tres->number;i++)
  {
   if(atm[i]==0) continue;
   atm[j++]=atm[i];
  }
  res_temp->atm=atm[0];
  for(i=0;i<j-1;i++)
  {
    atm[i]->next=atm[i+1];
  }
  atm[j-1]->next=0;
}

// detect the number of residue and atoms
// number and residue id0 and atom id0 and atom id and residue number

n=0;m=0;nn=0;
chn_temp->number=0;
start=res->id;
for(res_temp=chn_temp->res;res_temp;res_temp=res_temp->next)
{
  res_temp->chn=this;
  res_temp->id0=m++;
  chn_temp->number++;

  //atom id
  res_temp->id=res_temp->id-start;
  res_temp->number=0;res_temp->nhydr=0;
  for(atm_temp=res_temp->atm;atm_temp;atm_temp=atm_temp->next)
  {
    res_temp->number++;
    atm_temp->id0=n++;
    if(atm_temp->name[1]=='H') res_temp->nhydr++;
  }

  //atom id0 with tres complete supposed 
  for(tatm_temp=res_temp->tres->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  {
    atm_temp=(*res_temp)[tatm_temp->id];
    if(atm_temp)atm_temp->id=nn;
    nn++;
  }
}

// detect the standard residue

for(res_temp=chn_temp->res;res_temp;res_temp=res_temp->next)
{
  n=res_temp->number-res_temp->nhydr;
  k=res_temp->tres->number-res_temp->tres->nhydr;
  if(n!=k) 
  {
    res_temp->regular=0;
    continue;
  }
}

// assign topology to the pdb read and connectivity

res_temp1=0;
for(res_temp=chn_temp->res;res_temp;res_temp=res_temp->next)
{
   for(atm_temp=res_temp->atm;atm_temp;atm_temp=atm_temp->next)
   {
     for(i=0;i<4;i++)
     {
       tatm_temp=atm_temp->tatm->bond[i];
       if(!tatm_temp) continue;
       atm_temp->bond[i]=(*res_temp)[tatm_temp->id];
     }
     if(atm_temp->tatm->id==0)
     {
       atm_temp->bond[0]=0;
       if(res_temp1) atm_temp1=(*res_temp1)[2]; 
       if(res_temp1==0||atm_temp1==0) {
	 atm_temp->bond[0]=0;
       }
       /*
       else if(TRES.distance(atm_temp->xyz,atm_temp1->xyz)>2) {
	 atm_temp->bond[0]=0;
       }
       */
       else if(res_temp1&&res_temp->id-res_temp1->id==1)
       {
         atm_temp->bond[0]=atm_temp1;
       }
       else
       {
         atm_temp->bond[0]=0;
       }
     }
     if(atm_temp->tatm->id==2)
     {
       atm_temp->bond[1]=0;
       if(res_temp->next)atm_temp1=(*res_temp->next)[0];
       if(res_temp->next==0||atm_temp1==0) {
	 atm_temp->bond[1]=0;
       }
       /*
       else if(TRES.distance(atm_temp->xyz,atm_temp1->xyz)>2) {
	 atm_temp->bond[1]=0;
       }	
       */  		
       else if(res_temp->next&&res_temp->next->id-res_temp->id==1) 
       {
         atm_temp->bond[1]=(*res_temp->next)[0];
       }
       else 
       {
        atm_temp->bond[1]=0;
       } 
     }
   } 
   res_temp1=res_temp;   
}
}

int Chn::manyatm()
{
int i;
Res *r;
Atm *a;

i=0;
for(r=res;r;r=r->next)
for(a=r->atm;a;a=a->next)i++;

return i;
}

void Chn::seqout(FILE *fp)
{
int n,m,k,j,i;
char *s;
int num;

i=-1;
m=-1;
n=-1;

s=getseq();
num=strlen(s);

for(k=0;k<num/50+1;k++)
{
   for(j=0;j<50;j++){
    i++;
    if(i>=num) break;
    if(i%10==0) fprintf(fp," ");
    fprintf(fp,"%c",s[i]);
   }
   fprintf(fp,"\n");
}
fprintf(fp,"\n\n");

i=-1;
m=-1;
n=-1;

for(k=0;k<num/50+1;k++) 
{
   for(j=0;j<50;j++){
    i++;
    if(i>=num) break;
    if(i%10==0) fprintf(fp," ");
    fprintf(fp,"%c",s[i]);
   }
   fprintf(fp,"\n");
   for(j=0;j<50;j++){
     m++;
     if(m>=num) break;
     if(m%10==0) fprintf(fp," ");
     fprintf(fp,"%c",s[1+num+m]);
   }
   fprintf(fp,"\n");
   for(j=0;j<50;j++){
     n++;
     if(n>=num) break;
     if(n%10==0) fprintf(fp," ");
     fprintf(fp,"%c",s[2+2*num+n]);
   }
   fprintf(fp,"\n\n");if(i>=num-1)  {delete [] s; return;}
}
}

float *Chn::center(Res *s,int end,int fa,int fb)
{
  Res *r;
  Atm *am;
  float *coo;
  int i,j;

  coo=new float[3];
  coo[0]=coo[1]=coo[2]=0;
  i=0;
  for(r=s;r;r=r->next)
  {
    if(r->id0>end) break;
    for(am=r->atm;am;am=am->next)
    {
       if(fa==0&&am->tatm->id!=1) continue;
       if(fa==1&&am->tatm->id>=4) continue;
       if(fa==2&&am->tatm->name[1]=='H') continue;
       for(j=0;j<3;j++) coo[j]+=am->xyz[j];
       i++;
    }
  }

  for(j=0;j<3;j++)  coo[j]= coo[j]/i;
  if(fb==0) return coo;
  for(r=s->chn->res;r;r=r->next)
  for(am=r->atm;am;am=am->next)
  {
     for(j=0;j<3;j++) am->xyz[j]-=coo[j];
  }
  return coo;
}

int Chn::lastid(int flg)
{
  Res *r;
  Atm *a;
  int i;
  for(r=res;r;r=r->next)
  {
    if(flg==0)
    {
      i=r->id0;
    }
    else
    {
      for(a=r->atm;a;a=a->next)
      {
        i=a->id0;
      }
    }
  }
  return i;
}

Res *Chn::lastres() {

	Res *r;
	if(res==0) return 0;
	for(r=res;r->next;r=r->next);
	return r;

}

void Chn::addresiduesonly(int from,int end,char *seqn) {

	Res *r;
	Res *r0;
	Res *r1; 
	Chn temp;
	Rotate rot;
	//addresid(start); 
	start=0;
	temp.create(seqn);
	
	header();
	//temp.configure();

	r0=isresoid(from-1);
        r=isresoid(end+1);
	r1=temp.lastres();
	if(r==0&&r0==0) {
		cerr<<"the residues at both end not exist:"<<from<<" "<<end<<endl;	
		exit(0);
	}
	if(r0) {
		Res *ta=isresoid(from);	   
		Res *tb=isresoid(end);
		
		int ii=r0->id+1;
        	int nn=r0->id0+1;
		//int ee=r0->oid+1;
		Res *tr;
        	for(tr=temp.res;tr!=NULL;tr=tr->next) {
                	tr->id=ii;ii++;
                	tr->id0=nn;nn++;
			//tr->oid=ee;ee++;
        	}
		Res *rr;
		int n=r1->id+1;
		int n0=r1->id0+1;
		for(rr=r;rr;rr=rr->next) {
			rr->id=n;n++;
			rr->id0=n0;n0++;
		}
		temp.res->oid=ta->oid;
		r1->oid=tb->oid;
	}	
	else if(r) {
		Res *ta=isresoid(from);	   
		Res *tb=isresoid(end);

		int ii=r->id-1;
        	int nn=r->id0-1;
		//int ee=r->oid-1;
		Res *tr;
        	for(tr=temp.res;tr!=NULL;tr=tr->next) {
                	tr->id=ii;ii--;
                	tr->id0=nn;nn--;
			//tr->oid=ee;ee--;
        	}
		temp.res->oid=ta->oid;
		r1->oid=tb->oid;
	}
  
	if(r0==NULL) {
		res=temp.res;			
		rot.link(r1,r,r1->id0,0);
		r1->next=r;	
	}
	else if(r==0) {
		r0->next=temp.res;
		rot.link(r0,temp.res,temp.res->id0,1);
	}
	else {
		rot.link(r0,temp.res,temp.res->id0,1);
		if(strlen(seqn)>=2) rot.link(r1,r,r1->id0,0);
	        else if(temp.res->regular==1) {
			Res *re=isresoid(from+1);
			float dn=TRES.distance(re->temp+3,re->temp+6)-
			         TRES.distance(temp.res->atm->next->next,temp.res->atm->next->next->next);
			dn=fabs(dn);
 			float ter[1000];
			float ddr=1000000;
			for(int ii=-180;ii<180;ii+=5) {
				rot.rotate(temp.res->atm->next,5,1);
				for(int ii1=-180;ii1;ii1+=5) {
					rot.rotate(temp.res->atm->next->next,5,1);
					float d=TRES.distance(temp.res->atm->next->next->xyz,re->temp+3)+
						TRES.distance(temp.res->atm->next->next->next->xyz,re->temp+6);
					if(d<ddr) {
						ddr=d;
						temp.transfer(ter,temp.res,temp.res->id0,0);
					}	
					//cout<<d<<" "<<ii<<" "<<ii1<<endl;
					//if(d<dn+0.5) goto rebre;
				}
			}
			temp.transfer(ter,temp.res,temp.res->id0,1);
			float d=TRES.distance(temp.res->atm->next->next->xyz,re->temp+3)+
                                                TRES.distance(temp.res->atm->next->next->next->xyz,re->temp+6);
			if(TRES.logg>3) cerr<<d<<endl;
		}
		rebre:
		r0->next=temp.res;
		r1->next=r;
		
	}

	temp.res=0;
	temp.pdb=0;
	temp.next=0;
	pdb->configure();
	pdb->setlastresorder();
	if(TRES.logg>3) pdb->write("ss");
	pdb->header();
}	


void Chn::addresid(int d) {
 for(Res *r=res;r;r=r->next) r->id+=d;
}


void Chn::addresid0(int d) {
 for(Res *r=res;r;r=r->next) r->id0+=d;
}


void Chn::addresoid(int d) {
 for(Res *r=res;r;r=r->next) r->oid+=d;
}




void Chn::addatmid0(int d) {

Res *r;
Atm *a; 
for(r=res;r;r=r->next)
for(a=r->atm;a;a=a->next) a->id0+=d;

}

void Chn::idequal(Res *r1,Res *r2,int n) {

int i=0;

Res *r3,*r4;

r3=r1;
r4=r2;

while(r3&&r4) {

if(i>=n) return;  

r1->id=r2->id;
r1->id0=r2->id0;

r3=r1->next;
r4=r2->next;

i++;
}	

}

void Chn::addsidechain() {

  for(Res *r=res;r;r=r->next) r->addsidechain();
}

void Chn::addhatoms() {
  for(Res *r=res;r;r=r->next) r->addhatoms(0);
}

void Chn::addmissingatms(Res *rr,int n){
  Res *r;
  for(r=rr;r;r=r->next) {
	if(r->id0>n) break;
	r->addmissingatms();
  }
}
void Chn::addhatoms(int f) {
  for(Res *r=res;r;r=r->next) r->addhatoms(f);
}

void Chn::buryangle(char *s,int from,int to) {
  
  FILE *fp=fopen(s,"w");
  for(Res *r=res;r;r=r->next) r->buryangle(fp,from,to);
}

void Chn::buryangle(FILE *fp,int from,int to) {
  
  for(Res *r=res;r;r=r->next) r->buryangle(fp,from,to);
}

void Chn::setseqcard() {

  if(seqcard) delete [] seqcard;

  seqcard=0;

  Res *r=lastres();

  int n=r->id;

  seqcard=new char[n+2];

  for(int i=0;i<n+2;i++) {
	seqcard[i]='X';
  }
  seqcard[n+1]='\0'; 

  for(r=res;r;r=r->next) {
       seqcard[r->id]=r->name;
  }
}

void Chn::setseqcardnogap() {

  if(seqcard) delete [] seqcard;

  seqcard=0;

  Res *r=lastres();

  int n=r->id;

  seqcard=new char[n+2];

  for(int i=0;i<n+2;i++) {
        seqcard[i]='X';
  }
  seqcard[n+1]='\0';

  int m=0;
  for(r=res;r;r=r->next) {
       seqcard[m++]=r->name;
  }
  seqcard[m]='\0';
}

void Chn::setchnonly(){

	next=0;	
	pdb=0;
	nemp=0;
	if(temp) {delete [] temp;temp=0;}
	area=0;
}

void Chn::deletehbondlist(){

	Res *r;

	for(r=res;r;r=r->next) {

		if(r->hbond) {delete r->hbond;r->hbond=0;}
	}

	return;
}

void Chn::setsecstr() {

//decide the secondary structure

//detect helix
Res *r;

for(r=res;r;r=r->next)
{
  if(r->ishbondedwithres(r->id+4)) isres(r->id+4)->sec='h';
  if(r->ishbondedwithres(r->id+3)) isres(r->id+3)->sec='h';
  if(r->ishbondedwithres(r->id+5)) isres(r->id+5)->sec='h';
}

//connect break helix

Res *r0=0;
for(r=res;r;r=r->next)
{
 if(r->next&&r0)
 {
  if(r->sec=='-'&&r0->sec=='h'&&r->next->sec=='h') r->sec='h';
 }
 r0=r;
}
//return;
//kick out fake helix
r0=0;
int i=0;
Res *r1;
for(r=res;r;r=r->next)
{
 if(r->sec=='h') i++;
 else
 {
   if(i<3)
   {  
     for(r1=r0;r1;r1=r1->next)
     {
       if(r1->id>r->id) break;
       r1->sec='-';
     }
   }
   i=0;
   r0=0;
 }
 if(i==1) r0=r;
}

//detecting beta sheet

for(r=res;r;r=r->next)
{
    Atm *n=r->isatmid(0);
    Atm *o=r->isatmid(3);
    for(r0=res;r0;r0=r0->next) {
	Atm *n0=r0->isatmid(0);
	Atm *o0=r0->isatmid(3);
	if(r->ishbondedwithatm(n,o0)&&r->ishbondedwithatm(o,n0)) {
		r->sec='e';	
		r0->sec='e';
	}
    }
}

r0=0;
for(r=res;r;r=r->next)
{
 if(r->next&&r0)
 {
  if(r->sec=='-'&&r0->sec=='e'&&r->next->sec=='e') r->sec='e';
 }
 r0=r;
}

//kick out fake sheet
r0=0;
i=0;
res->sec='-';
for(r=res;r;r=r->next)
{
 if(r->sec=='e') i++;
 else
 {
   if(i<3)
   {
     for(r1=r0;r1;r1=r1->next)
     {
       if(r1->id>r->id) break;
       r1->sec='-';
     }
   }
   i=0;
   r0=0;
 }
 if(i==1) r0=r;
 if(r->next==0) r->sec='-';
}

}


	void Chn::findalphahelix( int pitch) {

    		//HBondList *hbond;
    		Chn *chain;
    		Res *group;
    		Res *first;
    		Res *ptr;
    		Atm *srcCA;
    		Atm *dstCA;
    		int res,dist,prev;

   	 	/* Protein chains only! */
    		//hbond = 0;//molecule.hBond;
    		//for( chain=molecule.strand; chain!=0; chain=chain.next )
		chain=this;
    		if( (first=chain->res)!=0) {   
			prev = 0; dist = 0;
        		for( group=chain->res; group!=0; group=group->next ) {   
				//if( MolStruct.IsAmino(group.refno) && (srcCA=group.getAtomInResidue(1))!=0)
				srcCA=group->isatmid(0);
				if(srcCA)
            			{  
				 	if( dist==pitch ) {   
						res = 0;
                    				dstCA=first->isatmid(3);//first.getAtomInResidue(1);
						if(group->ishbondedwithatm(dstCA,srcCA)) res=1;
						
                    				/*while( hbond!=0 && hbond.srcCA == srcCA ) {  
				 			if( hbond.dstCA==dstCA ) res=1;
                        				hbond = hbond.next;
                    				}*/
						

                    				if( res !=0) {   
							Res *las=isres(first->id-1);
							las=0;							
							if( prev!=0||(las&&las->sec=='h')) {   
								//if((first.struc&MolStruct.HelixFlag)==0) InfoHelixCount++;
	
                            					ptr = first;
                           			 		do {
                                					ptr->sec = 'h';
                                					ptr = ptr->next;
                            					} while( ptr != group );
								//group->sec='h';
                        				} else prev = 1;
                    				} else prev = 0;
                			}  
					/*else {
						while((hbond!=0) &&(hbond.srcCA==srcCA)) hbond = hbond.next;
					}*/
            			} else prev = 0;

            			if( group->sec=='h') {   
					first = group; prev = 0; dist = 1;
            			} 
				else if( dist==pitch ) {   
					first = first->next;
            			} else dist++;
        		}
    		}/* 
		else if( first!=0 && MolStruct.IsNucleo(first.refno)){
        		while( hbond!=0&&!MolStruct.IsAminoBackbone(hbond.src.refno))
           		 hbond = hbond.next;
		}*/

		
		
	}

void Chn::writesecondary(FILE *fp) {
        for(Res *r=res;r;r=r->next) fprintf(fp,"%c%i %c\n",r->name,r->id,r->sec);
}


void Chn::writehbondlist(FILE *fp) {
	for(Res *r=res;r;r=r->next) r->writehbondlist(fp);
}

void Chn::writessbondlist(FILE *fp) {
        for(Res *r=res;r;r=r->next) r->writessbondlist(fp);
}


	void Chn::findbetasheets() {
   		Chn *chain;
    		int ladder;
    		int count;
		Res *nexti,*curri,*previ; 
		//Atm *cprevi,*ccurri,*cnexti;
		//HBondList *hcurri,*hnexti;
		chain=this;
    		//hnexti = molecule.hBond;
    		//for( chain=molecule.strand; chain!=0; chain=chain.next )
        		if( (nexti = chain->res)!=0)
            			if(1) {//if( MolStruct.IsProtein(nexti.refno) ) {   
					count = 1;
                			ladder = 0;
                			do {
                    				//cnexti = nexti->isatmid(0); //nexti.getAtomInResidue(1);

                    				if( count == 3 ) {   
							testladder(previ,curri,nexti);
                        				if(curri->sec=='e') {   
								if(  ladder==0) {   
									//InfoLadderCount++;
                                					ladder = 1;
                            					}
                        				} else ladder = 0;
                    				} else count++;

                    				//cprevi = ccurri; ccurri = cnexti;
						previ=curri;
                    				curri = nexti;   //hcurri = hnexti;
                    				//while( (hnexti!=0) && (hnexti.srcCA==cnexti) )
                        				//hnexti = hnexti.next;
                			} while((nexti = nexti->next)!=0);

            			} 
				/*else if( MolStruct.IsNucleo(nexti.refno) ){
                			while( hnexti!=0 && !MolStruct.IsAminoBackbone(hnexti.src.refno) )
                    				hnexti = hnexti.next;
				}*/				
	}



 
	void Chn::testladder(Res *previ,Res *curri,Res *nexti)
	{
		//Chn* chain;
    		//Atm *cprevj=0, *ccurrj=0, *cnextj=0;
		Res *prevj=0;
    		//HBondList *hcurrj=0,*hnextj=0;
    		Res *currj=0, *nextj=0;
    		int count, result, found;

		//chain=this;
    		/* Already part of atleast one ladder */
    		//found = curri.flag & MolStruct.SheetFlag;
		if(curri->flag=='e') found=1;
		else found=0;

		if(nexti==0) return;
    		nextj = nexti->next;
		if( nextj==0 )return;
    		//hnextj = hnexti;
		
    		//while( hnextj!=0 && hnextj.srcCA==cnexti )
        	//	hnextj = hnextj.next;

    		//while( 1 )
    		//{   if( nextj!=0 )
            		//if( MolStruct.IsProtein(chain.residue.refno) )
			if(1)
            		{   count = 1;
                		do {
                    			//cnextj = nextj->isatmid(0);//nextj.getAtomInResidue(1);
                    			if( count == 3 )
                    			{   
						if( nexti->ishbondedsource(currj)&&
                            			    currj->ishbondedsource(previ))
                        			{   
							result = 1;
                        			} 
						else if(nextj->ishbondedsource(curri)&&
                                                        curri->ishbondedsource(prevj))
                        			{   
							result = 1;
                        			} 
						else if(nexti->ishbondedsource(prevj) &&
                                   			nextj->ishbondedsource(previ))
                        			{   
							result = -1;
                        			} 
						else if(currj->ishbondedsource(curri) &&
                                   			curri->ishbondedsource(currj))
                        			{   
							result = -1;
                        			} 
						else result =0;

                        			if( result!=0 )
                        			{   
							curri->sec='e';//struc |= MolStruct.SheetFlag;
                            				currj->sec='e';//.struc |= MolStruct.SheetFlag;
                            				if( found!=0 ) return;
                            				found = 1;
                        			}
                    			} else count++;

                    			//cprevj = ccurrj; ccurrj = cnextj;
					prevj = currj;
                    			currj = nextj;   //hcurrj = hnextj;
					
                    			//while( hnextj!=0 && hnextj.srcCA==cnextj )
                        			//hnextj = hnextj.next;
                		} while( nextj&&(nextj = nextj->next)!=0 );

            		} /*else if( MolStruct.IsNucleo(chain.residue.refno) )
                		while( hnextj!=0 && !MolStruct.IsAminoBackbone(hnextj.src.refno) )
                    			hnextj = hnextj.next;*/

        	    /*if( (chain = chain.next)!=0 )
      		    {   nextj = chain.residue;
      		    } else return;*/
  		//}
	}







	void Chn::findturnstructure() {
    		Atm *aptr[5];
    		Res* group;
   		Res *prev=0;
    		Atm *ptr;
    		float ux,uy,uz,mu;
    		float vx,vy,vz,mv;
    		int i,len;
   		double CosKappa;

    		//for( chain=molecule.strand; chain!=null; chain=chain.next )
		Chn *chain=this;
        		if( chain->res)
       			 {   
				len = 0;  //found = 0;
            			for( group=chain->res; group!=0; group=group->next )
            			{    
					ptr = group->isatmid(1); //group.getAtomInResidue(1);
                 			if(ptr!=0&&(chain->isres(group->id-1)==0||chain->isres(group->id+1)==0))
                 			{   
						//found = 0;
                     				len = 0;
                 			} 
					else if( len==5 )
                 			{   
						for( i=0; i<4; i++ ) aptr[i] = aptr[i+1];
                     				len = 4;
                 			} 
					else if( len==2 ) prev = group;

                			aptr[len++] = ptr;
                 			if( len==5 )
                 			{   
						if( (prev->sec!='e'&&prev->sec!='h') &&
                         				aptr[0]!=0 && aptr[2]!=0 && aptr[4]!=0 )
                     				{   
							ux = ((aptr[2]->xyz[0] - aptr[0]->xyz[0]));
                         				uy = ((aptr[2]->xyz[1] - aptr[0]->xyz[1]));
                         				uz = ((aptr[2]->xyz[2] - aptr[0]->xyz[2]));

                         				vx = ((aptr[4]->xyz[0] - aptr[2]->xyz[0]));
                         				vy = ((aptr[4]->xyz[1] - aptr[2]->xyz[1]));
                         				vz = ((aptr[4]->xyz[2] - aptr[2]->xyz[2]));

                         				mu = (ux*ux + uy*uy + uz*uz);
                         				mv = (vx*vx + vz*vz + vy*vy);
                         				if( mu!=0 &&mv!=0)
                         				{   
								CosKappa =  (ux*vx + uy*vy + uz*vz);
                             					CosKappa /= sqrt(mu*mv );
                             					if( CosKappa<0.34202014332567 ) {   
									//if(  found==0 )  InfoTurnCount++;
                                 					prev->sec = 't'; //MolStruct.TurnFlag;
                             					}
                         				}
                     				}
                     				//if(prev->sec=='t') found=1;
						//else               found=0; //found = prev.struc&MolStruct.TurnFlag;
                     				prev = prev->next;
                	 		} /* len==5 */
            			}
        		}
	}

void Chn::setthreestatesec()
{

for(Res *r=res;r;r=r->next) {
	if(r->sec!='e'&&r->sec!='h') r->sec='-';
}
}
void Chn::setdsspstr() {
findalphahelix(4);
findalphahelix(3);
findalphahelix(5);
findbetasheets();
findturnstructure();
}

void Chn::setdsspstr(int n) {
findalphahelix(n);
findbetasheets();
findturnstructure();
}


float Chn::getarea() {
float e=0;
for(Res *r=res;r;r=r->next)e+=r->getarea();
return e;
}

int Chn::isnearnext(Res *r0,Res *r1) {

	Res *a,*b;

	//out of possible
	if(r0==r1) return 1;

	if(r0==0||r1==0) return 0;
	else if(r0->chn!=r1->chn) return 0;
	else if(abs(r0->id-r1->id)!=1) return 0;

	if(r0->id<r1->id) {a=r0;b=r1;}
	else if(r0->id>r1->id) {a=r1;b=r0;}
	else return 0;

	//check 
	
	for(int i=1;i<=3;i++) {

		Atm *t=a->isatmid(i);
		//Tatm *tt=a->tres->isatm(i);
		if(t==0) return 0;
		for(int j=0;j<=1;j++) {
			Atm *t0=b->isatmid(j);
			Tatm *tt0=b->tres->isatm(j);
			if(t0==0||tt0==0) return 0;
			float d=TRES.distance(t->xyz,t0->xyz);
			
			float dd=TRES.distance(b->tres->head+3*i-3,tt0->xyz);		
			if(fabs(d-dd)>1) return 0;
		}
	}
	
	return 1;
}

int Chn::isssbond(Res *r0,Res *r1) {

	if(r0==r1) return 1;

	if(r0==0||r1==0) return 0;
	else if(r0->chn!=r1->chn) return 0;
        else if(abs(r0->id-r1->id)<=1) return 0;

	Atm *cb0,*s0;
	Atm *cb1,*s1;

	cb0=r0->isatm(" CB ");
	s0=r0->isatm(" SG ");

	cb1=r1->isatm(" CB ");
        s1=r1->isatm(" SG ");

	if(cb0==0||s0==0||cb1==0||s1==0) return 0;

	float d1=TRES.distance(cb0->xyz,s1->xyz);
	float d2=TRES.distance(cb1->xyz,s0->xyz);
	float d3=TRES.distance(s0->xyz,s1->xyz);

	if(d3<3&&d1<4&&d2<4) return 1;

	return 0;
}

void Chn::setreversethread() {

	Res *r0=0;

	for(Res *r=res;r;r=r->next){
		r->last=r0;
		r0=r;
	}
}

void Chn::setlastresorder() {

	setreversethread();
}

void Chn::setreversezero() {


        for(Res *r=res;r;r=r->next){
                r->last=0;
        }
}

void Chn::writeoutrotamer(FILE *fpp,char c) {


        Res *r=isres(c,0);

        if(r==0) return;

        if(fpp==0) return;

        int n=0;
        for(Res *rr=r->more;rr;rr=rr->more) {

                rr->write(fpp,4);
                fflush(fpp);
                Res *rtt=rr;
                fprintf(fpp,"HEAD 1 %8.3f %8.3f %8.3f\n",
                rtt->temp[0],rtt->temp[1],rtt->temp[2]);
                fprintf(fpp,"HEAD 2 %8.3f %8.3f %8.3f\n",
                rtt->temp[3],rtt->temp[4],rtt->temp[5]);
                fprintf(fpp,"HEAD 3 %8.3f %8.3f %8.3f\n",
                rtt->temp[6],rtt->temp[7],rtt->temp[8]);
                fflush(fpp);
		n++;
        }

        cerr<<"total for residue :"<<r->name<<" "<<"is "<<n<<endl;
}

void Chn::insertresidues(char *seqn,int posn,int direct) {
 
	Res *r=isres(posn); 

	Res *begin,*end;

	if(r==0) {
		cerr<<"the residue at position:"<<posn<<" does not exist"<<endl;
	}

	//prepare header
	//header();

	//create the sequence according to seqn

	Chn temp;
	temp.create(seqn);
 	
	begin=temp.res;
	end=temp.lastres();

	//insert residue
	if(direct==1) {
		//begin=r;
		Res *r1=r->next;
		r->next=temp.res;
		int sd=r->id;
		Res *r2,*r3; 
		for(r2=temp.res;r2;r2=r2->next) {
			sd++;
			r2->id=sd;
			r2->chn=this;
			r3=r2;
		}
		//temp.res=0;
		//end=r1;
		r3->next=r1;
		
		int ddo=1;
		if(r1==0) ddo=0;
		else if(r1->id>r3->id) ddo=0; 
		if(ddo) {
 			int sdd=r3->id-r1->id+1;		
			for(r2=r1;r2&&ddo;r2=r2->next) {			 
				r2->id+=sdd;	
			}
		}	
	}
	else if(direct==-1) {
		
		//end=r;
		int n=temp.manyres();
		Res *r1=isres0(r->id0-1);
		Res *r2=temp.lastres();
		r2->next=r; 
		int sd=r->id-n-1;
		
		//Res *r3; 
		for(r2=temp.res;r2;r2=r2->next) {
			if(r2==r) break;
			sd++;
			r2->id=sd;
			r2->chn=this;
			//r3=r2;
		}
		//begin=r1;	
		if(r1) r1->next=temp.res;
		//temp.res=0;
		 
		
		int ddo=1;
		if(r1==0) {
			ddo=0;
			res=temp.res;
		}
		else if(r1->id<temp.res->id) ddo=0; 
		if(ddo) {
 			int sdd=temp.res->id-r1->id-1;		
			for(r2=res;r2&&ddo;r2=r2->next) {							 
				r2->id+=sdd;	
				if(r2==r1) break;
			}		
		}
			 
	}
	 
	int s=start;
	pdb->configure();
	start=s;
 	temp.write("temp");
	temp.res=0;
	temp.pdb=0;
	temp.next=0;
	pdb->write("ss");
	 
	header();
	Rotate rot;

	if(direct==-1) {
		Res *rr=isres0(end->id0+1);
		rot.link(end,rr,begin->id0,0);
		//rot.link(end,rr,1000,1);
	}	
	else if(direct==1) {
		Res *rr=isres0(begin->id0-1);
                rot.link(rr,begin,end->id0,1);
	}
	/*
	Bound *bound=new Bound;
	//
	//temp.create(seqn);
	//temp.configure();
	//pdb->chn=&temp;
	//
	bound->model=pdb;
	bound->flag|=TRES.constant->pdbundeletable;
	bound->ready();
 
	//
	//bound->setpredt(1);
	// 
	
	bound->setpredt(0);	
	bound->setpredt(begin,end->id0+1,1);		
	 
	bound->setcloseresidue();
	//bound->delsidechainbound();
	bound->maxnit=500;
	bound->adjust();
	bound->reorder();
	bound->clearredundancy();
	bound->setimproper();
	bound->settag();	
	bound->boundsmooth(5);
	bound->writeoutdistbound("out0");
	bound->setbuffertag();
	bound->writeoutdistbound("out");
	bound->setrange();
	bound->setbond();
	bound->setcurve(1);
	bound->writeatmbound("out_1");
	bound->boundsquare(2);
	bound->printoutimproper();
	//bound->setbuddy();
	bound->searchstructure(10);
	*/
}

void Chn::allnearbond(int n) {
	
	for(Res *r=res;r;r=r->next)
	for(Atm *a=r->atm;a;a=a->next)
	a->allnear(n,0);

}

void Chn::replace(char *seqn,int a,int b) {
	
	Res *ra=isres0(a);
	Res *rb=isres0(b);

	if(ra==0||rb==0) {
		cerr<<"the residue at: "<<a<<" and/or "<<b<<" does not exist"<<endl;
		return;
	}

	Res *ra0=isres0(a-1);
	Res *rb0=isres0(b+1);

	int n1=strlen(seqn);
	int n2=b-a+1;

	if(ra0&&rb0) { 
		ra0->next=rb0;
		rb->next=0;
		delete ra;
		pdb->configure();
		insertresidues(seqn,ra0->id0,1);
	}	
	else if(ra0) {
		ra0->next=0;
		rb->next=0;
                delete ra;
                pdb->configure();
                insertresidues(seqn,ra0->id0,1);
	}
	else if(rb0) {
		res=rb0;
		rb->next=0;
		delete ra;
		pdb->configure();       
                insertresidues(seqn,rb0->id0,-1);   
	}
	else {
		if(n1==0) {
			if(res) delete res;
			res=0;
			pdb->configure();
			return;
		}
		else if(n2<=0) {
			Chn temp;
        		temp.create(seqn);
			res=temp.res;
			temp.res=0;
			pdb->configure();
                        return;
		}
	}
}

int Chn::islinked(Res *r0,Res *r){

	if(r0==0||r==0) return 0;
	Atm *a=r0->isatmid(1);
	Atm *c=r0->isatmid(2);
	Atm *o=r0->isatmid(3);

	if(a==0||c==0||o==0) return 0;

	float d1=TRES.distance(a->xyz,r->temp);
	float d2=TRES.distsqr(c->xyz,r->temp+3);
	float d3=TRES.distsqr(o->xyz,r->temp+6);

	if(d1>0.5||d2>0.5||d3>0.5) return 0;

	return 1;
}

int Chn::isdatabaselinked(Res *r0,Res *r){

        if(r0==0||r==0) return 0;

	DistPopular *dst=TRES.popbin->getDistPopular("internal");
	if(dst==0) {
		return islinked(r0,r);
        }

	Atm *aa0[3],*aa[3];
        aa0[0]=r0->isatmid(1);
        aa0[1]=r0->isatmid(2);
        aa0[2]=r0->isatmid(3);

        aa[0]=r->isatmid(0);
	aa[1]=r->isatmid(1);


	int i,j;

	for(i=0;i<2;i++) {
		if(aa[i]==0) return 0;
		for(j=0;j<3;j++) {
			if(aa0[j]==0) return 0;
			DistPopular *d=dst->getDistPopular(r->tres,aa[i]->name,aa0[j]->name,r0->id-r->id);
			if(d==0) return 0;
			float e=TRES.distance(aa[i],aa0[j]);
			if(e<d->dist[0]-0.1) return 0;
			if(e>d->dist[1]+0.1) return 0;			
		}
	}

        return 1;
}


void Chn::deletecontact() {

	for(Res *r=res;r;r=r->next) r->deletecontact();

}

void Chn::buildcontact(Lattice *lat,Res *r0,int n,float dcut) {

	deletecontact();

	if(lat==0||r0==0) return;

	for(Res *r=r0;r;r=r->next) {
		if(r->id0>n) break;	

		for(Atm *a=r->atm;a;a=a->next) {
			lat->getcell(a,dcut);
			lat->cutoff(a,dcut);		
			int m=lat->nget;
			a->contact=new Contact();			
			a->contact->setatmarray(m);	
			m=0;			
			for(int i=0;i<lat->nget;i++) {
				Atm *a1=lat->obtain[i]->atm;
				if(a->isnear(a1)) continue;
				a->contact->atm[m++]=a1;
			}
			a->contact->num=m;
		}
	}
}

void Chn::setfreenextonmore(){

	for(Res *r=res;r;r=r->next) {
		r->setfreenextonmore();
		if(r->more) delete r->more;r->more=0;
	}
}

void Chn::setallnear(Res *r,int n) {

	Res *t;
	Atm *a;

	for(t=r;t;t=t->next) {
		if(t->id0>n) break; 

		for(a=t->atm;a;a=a->next) {
			if(a->near) delete [] a->near;a->near=0;
			if(a->nnear) delete [] a->nnear;a->nnear=0;
			a->allnear(3,0);
		}
	}
}

void Chn::printenergy(FILE *fp,Res *r,int n) {
	
	Res *t;

        for(t=r;t;t=t->next) {
                if(t->id0>n) break;

		fprintf(fp,"energy: %c%5i %f\n",t->name,t->id0,r->energy);
	}

}

void Chn::setresenergy(Res *r,int n) {
	Res *t;
	for(t=r;t;t=t->next) {
		if(t->id0>n) break;
		t->energy=0;
	}
}
float Chn::getresenergy(Res *r,int n,int f) {
	
	Res *t;
        Atm *a;

	float ee=0;
        for(t=r;t;t=t->next) {

                if(t->id0>n) break;

		t->energy=0;
		for(a=t->atm;a;a=a->next) {
			if(f==3) { //only backbone
				if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>3) continue;
				else if(a->tatm->id>3) continue;
			}
			else if(f==4) {
				if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id<=3) continue;
                                else if(a->tatm->id<=3) continue;
			}
			if(a->contact==0) continue;
			int i;
			a->energy=0;
			for(i=0;i<a->contact->num;i++) {
				int j=a->contact->atm[i]->res->id0;
				if(f==0) if(j>=r->id0&&j<=n) continue; //not the same segment
				if(f==1) if(j==t->id0) continue;  //not the same residue
				if(f==2) if(abs(j-t->id0)<=1) continue;	//not the neighboring residue
				a->energy+=a->contact->energy[i];
			}
			t->energy+=a->energy;
		}
		ee+=t->energy;
                printf("energy: %c%5i    %f\n",t->name,t->id0,t->energy);
        }

	return ee;
}
void Chn::setresenergy(Res *r,int n,int f) {
	
	Res *t;
        Atm *a;

        for(t=r;t;t=t->next) {

                if(t->id0>n) break;

		t->energy=0;
		for(a=t->atm;a;a=a->next) {
			if(f==3) { //only backbone
				if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>3) continue;
				else if(a->tatm->id>3) continue;
			}
			else if(f==4) {
				if(a->tatm->name[1]=='H'&&a->tatm->bond[0]->id<=3) continue;
                                else if(a->tatm->id<=3) continue;
			}
			if(a->contact==0) continue;
			int i;
			a->energy=0;
			for(i=0;i<a->contact->num;i++) {
				int j=a->contact->atm[i]->res->id0;
				if(f==0) if(j>=r->id0&&j<=n) continue; //not the same segment
				if(f==1) if(j==t->id0) continue;  //not the same residue
				if(f==2) if(abs(j-t->id0)<=1) continue;	//not the neighboring residue
				a->energy+=a->contact->energy[i];
			}
			t->energy+=a->energy;
		}
                if(TRES.logg) printf("energy: %c%5i    %f\n",t->name,t->id0,t->energy);
        }


}

Res *Chn::findmaxclash(Res *r1,int n) {

	Res *r;
	Res *s0=0;
	float d=-999999;
	
	for(r=r1;r;r=r->next) {
		if(r->id0>n) break;
		if(r==r1||r->energy>d) {
			s0=r;
			d=r->energy;
		}
	}
	return s0;
}
 
int Chn::isalllinked(Res *s,int n){

	Res *r;
	
	for(r=s;r;r=r->next) {
		if(r->id0>n||r->next==0) break;
		if(islinked(r,r->next)==0) return 0;
	}

	Res *t=isres0(n);
	if(t==0) t=lastres();
	for(r=t;r;r=r->last) {
		if(r->id0<s->id0||r->last==0) break;
                if(islinked(r->last,r)==0) return 0;
	}
	return 1;
}

int Chn::isalllinked(int n0,int n){

        Res *r,*s;

	s=isres0(n0);
	if(s==0) s=res;

        for(r=s;r;r=r->next) {
                if(r->id0>n||r->next==0) break;
                if(islinked(r,r->next)==0&&r->next->id-r->id>1) return 0;
        }

        Res *t=isres0(n);
        if(t==0) t=lastres();
        for(r=t;r;r=r->last) {
                if(r->id0<s->id0||r->last==0) break;
                if(islinked(r->last,r)==0&&r->id-r->last->id>1) return 0;
        }
        return 1;
}

float *Chn::gettransfer(Res *r0,int n){

	Res *r;
	int nseq=0;
  	for(r=r0;r&&r->id0<=n;r=r->next) nseq+=r->tres->number*3;

	float *xyz=new float[nseq];

	transfer(xyz,r0,n,0);

	return xyz;
}

int Chn::calcnumber(Res *r0,int end) {

	int nseq=0;
	Res *r;
  	for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
	return nseq;
}

void Chn::increaseid(Res *r0,int end,int nn) {

	Res *r;
        for(r=r0;r&&r->id0<=end;r=r->next)r->id+=nn;
}

int Chn::checkpdb(int flg) {

	Res *r;

	header();
	int nerr=0;
	//check id sequence
	for(r=pdb->chn->res;r;r=r->next) {
		if(r->next==0) continue;
		if(r->next->id-r->id>1&&pdb->chn->islinked(r,r->next)==1) {	
			cerr<<"warning!"<<endl;
			cerr<<"the residue "<<r->next->name<<r->next->oid<<" "<<r->name<<r->oid<<" is linked!"<<endl;
			cerr<<"they should have sequential id numbers."<<endl;
			cerr<<"the program has fixed it accordingly."<<endl;
			int n=r->next->id-r->id-1;
			Res *t;
			for(t=r->next;t;t=t->next) t->id-=n;	
			nerr++;
		}		 
	}

	//check breaker
	for(r=pdb->chn->res;r;r=r->next) {
		if(r->next==0) continue;
		if(r->next->id-r->id>1)  {
			cerr<<endl<<endl<<endl<<"Warning......"<<endl;
			cerr<<"the pdb file:"<<pdb->name<<" has breaker at:";
			cerr<<r->name<<r->oid<<" "<<r->next->name<<r->next->oid<<endl;
			cerr<<endl<<endl;
			nerr++;
		}
	}
	
	//check completeness
	for(r=pdb->chn->res;r;r=r->next) {
		Tatm *t;
		for(t=r->tres->tatm;t;t=t->next) {
			Atm *a=r->isatm(t->name);
			if(a==0) {
				cerr<<"atom missing in residue "<<r->name<<r->oid<<" : "<<t->name<<endl;
				nerr++;
			}
		}
	} 	

	return nerr;
}

float Chn::calcareadiff(int flg) {

	Res *r;
	Atm *a;

	float e=0;

	for(r=res;r;r=r->next) {  
	 	if(strchr("RKEDNQHTS",r->name)) e+=-r->area;
		else				 e+=r->area;
	}

	return e;
} 
