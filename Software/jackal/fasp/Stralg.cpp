#include"source.h"

Stralg::Stralg()
{
  xyza=0;xyzb=0;
  alga=0;algb=0;
  onlyaligned=0;
}

Stralg::~Stralg()
{
  if(xyza) {delete [] xyza;xyza=0;}
  if(xyzb) {delete [] xyzb;xyzb=0;}
  if(alga) {delete [] alga;alga=0;}
  if(algb) {delete [] algb;algb=0;}
}

float Stralg::superimposestill(int flg)
{
  Res *r;
  Atm *a,*a0;
  int nseq;
  int n,k,i,j; 
  float *coo;
  
  if(xyza) {delete [] xyza;xyza=0;}
  if(xyzb) {delete [] xyzb;xyzb=0;}

  Chn *chn;
  nseq=0; 
  if(flag==1) {
     for(chn=alga[0]->chn->pdb->chn;chn;chn=chn->next)
     for(r=chn->res;r;r=r->next) nseq+=r->tres->number*3; 
  }
  else if(flag==2){
     for(r=alga[0]->chn->res;r;r=r->next) nseq+=r->tres->number*3;
  }
  if(nseq==0) {cerr<<" no structure segment"<<endl; return -1;}

  xyza=new float[nseq+10];

  if(flag==1) {
    alga[0]->chn->pdb->transfer(xyza,0); 
  }
  else if(flag==2){
    alga[0]->chn->transfer(xyza,0);
  }

  coo=center(alga,flg); delete [] coo;

  nseq=0; 
  if(flag==1) {

     for(chn=algb[0]->chn->pdb->chn;chn;chn=chn->next)
     for(r=chn->res;r;r=r->next) nseq+=r->tres->number*3;
  }
  else if(flag==2) {
     for(r=alga[0]->chn->res;r;r=r->next) nseq+=r->tres->number*3;
  }

  if(nseq==0) {cerr<<" no structure segment"<<endl;return -1;}

  xyzb=new float[nseq+10];
   
  if(flag==1) algb[0]->chn->pdb->transfer(xyzb,0);
  else if(flag==2) algb[0]->chn->transfer(xyzb,0);

  coo=center(algb,flg); delete [] coo;
  
  i=0;while(alga[i])i++;
  j=0;while(algb[j])j++;
  if(i!=j)
  {
    cerr<<"the number of equivalence is not match:"<<i<<" "<<j<<endl;
    int nnn=0; 
    for(nnn=0;nnn<max(i,j);nnn++) {
        Res *ri,*ri0;
        ri=alga[nnn];ri0=algb[nnn];
	if(ri==0) cerr<<"none  ";
        else      cerr<<ri->name<<ri->id; 
        if(ri0==0) cerr<<" none"<<endl;
        else       cerr<<" "<<ri0->name<<ri0->id<<endl;
    }
 
    return -1;
  }
  k=i;

  float f0,fx0x,fx0y,fx0z,fy0x,fy0y,fy0z,fz0x,fz0y,fz0z;


  f0=fx0x=fx0y=fx0z=fy0x=fy0y=fy0z=fz0x=fz0y=fz0z=0;
  n=0;
  for(i=0;i<k;i++)
  {
    for(a=alga[i]->atm;a;a=a->next)
    {
      if(flg==0&&a->tatm->id!=1) continue;
      if(flg==1&&a->tatm->name[1]!='H'&&a->tatm->id>=4) continue;
      if(flg==1&&a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>=3) continue;
      if(flg==2&&a->tatm->name[1]=='H') continue;
      if(flg==3&&a->tatm->id>2) continue;
      a0=(*algb[i])[a->tatm->name];
      if(a0==0) continue;
      f0+=a->xyz[0]*a->xyz[0]+a->xyz[1]*a->xyz[1]+a->xyz[2]*a->xyz[2];  
      f0+=a0->xyz[0]*a0->xyz[0]+a0->xyz[1]*a0->xyz[1]+a0->xyz[2]*a0->xyz[2];
      fx0x+=a0->xyz[0]*a->xyz[0];fx0y+=a0->xyz[0]*a->xyz[1];fx0z+=a0->xyz[0]*a->xyz[2];
      fy0x+=a0->xyz[1]*a->xyz[0];fy0y+=a0->xyz[1]*a->xyz[1];fy0z+=a0->xyz[1]*a->xyz[2];
      fz0x+=a0->xyz[2]*a->xyz[0];fz0y+=a0->xyz[2]*a->xyz[1];fz0z+=a0->xyz[2]*a->xyz[2];
      n++;
    }
  }

  f0=f0/2;   

  float l1,m1,n1,l2,m2,n2,l3,m3,n3;
  float dd0[3],ss[3],cc[3];
  int rec,ff;
  float step;
  float dmin,d,e;
  int i1,i2,i3;

  for(i=0;i<3;i++) dd0[i]=0;
  step=0.17;
  dmin=f0-fx0x-fy0y-fz0z;
  if(dmin<0) dmin=0;
  cerr<<" the initial rmsd: "<<sqrt(dmin*2/n)<<endl;

  e=3.14/3;
  for(i1=0;i1<3;i1++)
  {
    //break;
    d=i1*e;
    ss[0]=sin(d);cc[0]=cos(d);
    for(i2=0;i2<6;i2++)
    {
      d=i2*e;
      ss[1]=sin(d);cc[1]=cos(d);
      for(i3=0;i3<6;i3++)
      {
        d=i3*e;
        ss[2]=sin(d);cc[2]=cos(d);
        l1=cc[1]*cc[2]-cc[0]*ss[1]*ss[2];l2=-cc[1]*ss[2]-cc[0]*ss[1]*cc[2];l3=ss[0]*ss[1];
        m1=ss[1]*cc[2]+cc[0]*cc[1]*ss[2];m2=-ss[1]*ss[2]+cc[0]*cc[1]*cc[2];m3=-ss[0]*cc[1];
        n1=ss[0]*ss[2];n2=ss[0]*cc[2];n3=cc[0];
        d=f0;
        d+=-l1*fx0x-l2*fy0x-l3*fz0x;
        d+=-m1*fx0y-m2*fy0y-m3*fz0y;
        d+=-n1*fx0z-n2*fy0z-n3*fz0z;
        if(d<dmin)
        {
           dd0[0]=i1*e;dd0[1]=i2*e;dd0[2]=i3*e;
           //cerr<<" "<<d<<" "<<step<<" "<<dd0[0]<<" "<<dd0[1]<<"  "<<dd0[2]<<endl;
           dmin=d;
        }
      }
    }
  }

  for(i=0;i<3;i++) dd0[i]=0;dmin=f0-fx0x-fy0y-fz0z;step=0.17;
  for(rec=0;rec<200;rec++)
  {
    ff=0; 

    for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++) {ss[j]=sin(dd0[j]);cc[j]=cos(dd0[j]);}
      for(j=-1;j<=1;j+=1)
      {
        if(j==0) continue;
        e=dd0[i]+j*step;
        ss[i]=sin(e);cc[i]=cos(e);
        l1=cc[1]*cc[2]-cc[0]*ss[1]*ss[2];l2=-cc[1]*ss[2]-cc[0]*ss[1]*cc[2];l3=ss[0]*ss[1];
        m1=ss[1]*cc[2]+cc[0]*cc[1]*ss[2];m2=-ss[1]*ss[2]+cc[0]*cc[1]*cc[2];m3=-ss[0]*cc[1];
        n1=ss[0]*ss[2];n2=ss[0]*cc[2];n3=cc[0];
        d=f0;
        d+=-l1*fx0x-l2*fy0x-l3*fz0x;
        d+=-m1*fx0y-m2*fy0y-m3*fz0y;
        d+=-n1*fx0z-n2*fy0z-n3*fz0z; 
        if(d<dmin)
        {
           //cerr<<rec<<" "<<i<<" "<<j<<" "<<step<<" "<<d<<endl;
           dd0[i]=e;
           dmin=d;
           ff=1;
           
        }
      }
    }
    if(ff==0) step=step/2;
    if(step<0.17/100) break;
    if(dmin<=0) {dmin=0;break;}
  }
  
  dmin*=2;
  //rotate 
  ss[0]=sin(dd0[0]);cc[0]=cos(dd0[0]);
  ss[1]=sin(dd0[1]);cc[1]=cos(dd0[1]);
  ss[2]=sin(dd0[2]);cc[2]=cos(dd0[2]);

  l1=cc[1]*cc[2]-cc[0]*ss[1]*ss[2];l2=-cc[1]*ss[2]-cc[0]*ss[1]*cc[2];l3=ss[0]*ss[1];
  m1=ss[1]*cc[2]+cc[0]*cc[1]*ss[2];m2=-ss[1]*ss[2]+cc[0]*cc[1]*cc[2];m3=-ss[0]*cc[1];
  n1=ss[0]*ss[2];n2=ss[0]*cc[2];n3=cc[0];

  if(flag==1) {
   for(chn=alga[0]->chn->pdb->chn;chn;chn=chn->next)
   for(r=chn->res;r;r=r->next)
   for(a=r->atm;a;a=a->next)
   {
      ss[0]=a->xyz[0]*l1+a->xyz[1]*m1+a->xyz[2]*n1;
      ss[1]=a->xyz[0]*l2+a->xyz[1]*m2+a->xyz[2]*n2;
      ss[2]=a->xyz[0]*l3+a->xyz[1]*m3+a->xyz[2]*n3;
      for(i=0;i<3;i++) a->xyz[i]=ss[i];
    }
  }
  else if(flag==2) {
    for(r=alga[0]->chn->res;r;r=r->next)
    for(a=r->atm;a;a=a->next)
    {
      ss[0]=a->xyz[0]*l1+a->xyz[1]*m1+a->xyz[2]*n1;
      ss[1]=a->xyz[0]*l2+a->xyz[1]*m2+a->xyz[2]*n2;
      ss[2]=a->xyz[0]*l3+a->xyz[1]*m3+a->xyz[2]*n3;
      for(i=0;i<3;i++) a->xyz[i]=ss[i];
    }
  }
  if(TRES.logg>1) cerr<<"final rmsd "<<sqrt(dmin/n)<<endl;
  return sqrt(dmin/n);
}

float Stralg::superimpose(int flg)
{
  Res *r;
  Atm *a,*a0;
  int nseq;
  int n,k,i,j; 
  float *cooa,*coob;
  
  if(xyza) {delete [] xyza;xyza=0;}
  if(xyzb) {delete [] xyzb;xyzb=0;}

  Chn *chn;
  nseq=0; 
  if(flag==1) {
     for(chn=alga[0]->chn->pdb->chn;chn;chn=chn->next)
     for(r=chn->res;r;r=r->next) nseq+=r->tres->number*3; 
  }
  else if(flag==2){
     for(r=alga[0]->chn->res;r;r=r->next) nseq+=r->tres->number*3;
  }
  if(nseq==0) {cerr<<" no structure segment"<<endl; return -1;}

  xyza=new float[nseq+10];

  if(flag==1) {
    alga[0]->chn->pdb->transfer(xyza,0); 
  }
  else if(flag==2){
    alga[0]->chn->transfer(xyza,0);
  }

  cooa=center(alga,flg); //delete [] coo;

  nseq=0; 
  if(flag==1) {

     for(chn=algb[0]->chn->pdb->chn;chn;chn=chn->next)
     for(r=chn->res;r;r=r->next) nseq+=r->tres->number*3;
  }
  else if(flag==2) {
     for(r=alga[0]->chn->res;r;r=r->next) nseq+=r->tres->number*3;
  }

  if(nseq==0) {cerr<<" no structure segment"<<endl;return -1;}

  xyzb=new float[nseq+10];
   
  if(flag==1) algb[0]->chn->pdb->transfer(xyzb,0);
  else if(flag==2) algb[0]->chn->transfer(xyzb,0);

  coob=center(algb,flg); //delete [] coo;
  
  
  i=0;while(alga[i])i++;
  j=0;while(algb[j])j++;
  if(i!=j)
  {
    cerr<<"the number of equivalence is not match:"<<i<<" "<<j<<endl;
    int nnn=0; 
    for(nnn=0;nnn<max(i,j);nnn++) {
        Res *ri,*ri0;
        ri=alga[nnn];ri0=algb[nnn];
	if(ri==0) cerr<<"none  ";
        else      cerr<<ri->name<<ri->id; 
        if(ri0==0) cerr<<" none"<<endl;
        else       cerr<<" "<<ri0->name<<ri0->id<<endl;
    }
 
    return -1;
  }
  k=i;

  float f0,fx0x,fx0y,fx0z,fy0x,fy0y,fy0z,fz0x,fz0y,fz0z;


  f0=fx0x=fx0y=fx0z=fy0x=fy0y=fy0z=fz0x=fz0y=fz0z=0;
  n=0;
  for(i=0;i<k;i++)
  {
    for(a=alga[i]->atm;a;a=a->next)
    {
      if(flg==0&&a->tatm->id!=1) continue;
      if(flg==1&&a->tatm->name[1]!='H'&&a->tatm->id>=4) continue;
      if(flg==1&&a->tatm->name[1]=='H'&&a->tatm->bond[0]->id>3) continue;
      if(flg==2&&a->tatm->name[1]=='H') continue;
      if(flg==3&&a->tatm->id>2) continue;
      a0=(*algb[i])[a->tatm->name];
      if(a0==0) continue;
      f0+=a->xyz[0]*a->xyz[0]+a->xyz[1]*a->xyz[1]+a->xyz[2]*a->xyz[2];  
      f0+=a0->xyz[0]*a0->xyz[0]+a0->xyz[1]*a0->xyz[1]+a0->xyz[2]*a0->xyz[2];
      fx0x+=a0->xyz[0]*a->xyz[0];fx0y+=a0->xyz[0]*a->xyz[1];fx0z+=a0->xyz[0]*a->xyz[2];
      fy0x+=a0->xyz[1]*a->xyz[0];fy0y+=a0->xyz[1]*a->xyz[1];fy0z+=a0->xyz[1]*a->xyz[2];
      fz0x+=a0->xyz[2]*a->xyz[0];fz0y+=a0->xyz[2]*a->xyz[1];fz0z+=a0->xyz[2]*a->xyz[2];
      n++;
    }
  }

  f0=f0/2;   

  float l1,m1,n1,l2,m2,n2,l3,m3,n3;
  float dd0[3],ss[3],cc[3];
  int rec,ff;
  float step;
  float dmin,d,e;
  int i1,i2,i3;

  for(i=0;i<3;i++) dd0[i]=0;
  step=0.17;
  dmin=f0-fx0x-fy0y-fz0z;
  if(dmin<0) dmin=0;
  if(TRES.logg>1) cerr<<" the initial rmsd: "<<sqrt(dmin*2/n)<<endl;

  e=3.14/3;
  for(i1=0;i1<3;i1++)
  {
    //break;
    d=i1*e;
    ss[0]=sin(d);cc[0]=cos(d);
    for(i2=0;i2<6;i2++)
    {
      d=i2*e;
      ss[1]=sin(d);cc[1]=cos(d);
      for(i3=0;i3<6;i3++)
      {
        d=i3*e;
        ss[2]=sin(d);cc[2]=cos(d);
        l1=cc[1]*cc[2]-cc[0]*ss[1]*ss[2];l2=-cc[1]*ss[2]-cc[0]*ss[1]*cc[2];l3=ss[0]*ss[1];
        m1=ss[1]*cc[2]+cc[0]*cc[1]*ss[2];m2=-ss[1]*ss[2]+cc[0]*cc[1]*cc[2];m3=-ss[0]*cc[1];
        n1=ss[0]*ss[2];n2=ss[0]*cc[2];n3=cc[0];
        d=f0;
        d+=-l1*fx0x-l2*fy0x-l3*fz0x;
        d+=-m1*fx0y-m2*fy0y-m3*fz0y;
        d+=-n1*fx0z-n2*fy0z-n3*fz0z;
        if(d<dmin)
        {
           dd0[0]=i1*e;dd0[1]=i2*e;dd0[2]=i3*e;
           //cerr<<" "<<d<<" "<<step<<" "<<dd0[0]<<" "<<dd0[1]<<"  "<<dd0[2]<<endl;
           dmin=d;
        }
      }
    }
  }

  for(i=0;i<3;i++) dd0[i]=0;dmin=f0-fx0x-fy0y-fz0z;step=0.17;
  for(rec=0;rec<200;rec++)
  {
    ff=0; 

    for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++) {ss[j]=sin(dd0[j]);cc[j]=cos(dd0[j]);}
      for(j=-1;j<=1;j+=1)
      {
        if(j==0) continue;
        e=dd0[i]+j*step;
        ss[i]=sin(e);cc[i]=cos(e);
        l1=cc[1]*cc[2]-cc[0]*ss[1]*ss[2];l2=-cc[1]*ss[2]-cc[0]*ss[1]*cc[2];l3=ss[0]*ss[1];
        m1=ss[1]*cc[2]+cc[0]*cc[1]*ss[2];m2=-ss[1]*ss[2]+cc[0]*cc[1]*cc[2];m3=-ss[0]*cc[1];
        n1=ss[0]*ss[2];n2=ss[0]*cc[2];n3=cc[0];
        d=f0;
        d+=-l1*fx0x-l2*fy0x-l3*fz0x;
        d+=-m1*fx0y-m2*fy0y-m3*fz0y;
        d+=-n1*fx0z-n2*fy0z-n3*fz0z; 
        if(d<dmin)
        {
           //cerr<<rec<<" "<<i<<" "<<j<<" "<<step<<" "<<d<<endl;
           dd0[i]=e;
           dmin=d;
           ff=1;
           
        }
      }
    }
    if(ff==0) step=step/2;
    if(step<0.17/100) break;
    if(dmin<=0) {dmin=0;break;}
  }
  
  dmin*=2;
  //rotate 
  ss[0]=sin(dd0[0]);cc[0]=cos(dd0[0]);
  ss[1]=sin(dd0[1]);cc[1]=cos(dd0[1]);
  ss[2]=sin(dd0[2]);cc[2]=cos(dd0[2]);

  l1=cc[1]*cc[2]-cc[0]*ss[1]*ss[2];l2=-cc[1]*ss[2]-cc[0]*ss[1]*cc[2];l3=ss[0]*ss[1];
  m1=ss[1]*cc[2]+cc[0]*cc[1]*ss[2];m2=-ss[1]*ss[2]+cc[0]*cc[1]*cc[2];m3=-ss[0]*cc[1];
  n1=ss[0]*ss[2];n2=ss[0]*cc[2];n3=cc[0];
  int ii=0;while(alga[ii])ii++;
  if(flag==1&&onlyaligned==0) {
   for(chn=alga[0]->chn->pdb->chn;chn;chn=chn->next)
   for(r=chn->res;r;r=r->next) {
   	for(a=r->atm;a;a=a->next)
   	{
      		ss[0]=a->xyz[0]*l1+a->xyz[1]*m1+a->xyz[2]*n1;
      		ss[1]=a->xyz[0]*l2+a->xyz[1]*m2+a->xyz[2]*n2;
      		ss[2]=a->xyz[0]*l3+a->xyz[1]*m3+a->xyz[2]*n3;
      		for(i=0;i<3;i++) a->xyz[i]=ss[i];
    	}
        
	if(r->temp) {
                int i,j;
                for(i=0;i<r->nemp/3;i++) {
		   ss[0]=r->temp[i*3]*l1+r->temp[1+i*3]*m1+r->temp[2+i*3]*n1;
      		   ss[1]=r->temp[i*3]*l2+r->temp[1+i*3]*m2+r->temp[2+i*3]*n2;
      		   ss[2]=r->temp[i*3]*l3+r->temp[1+i*3]*m3+r->temp[2+i*3]*n3;
                   for(j=0;j<3;j++) r->temp[j+i*3]=ss[j];
		}
        }
        
   }	
  }
  else if(flag==2&&onlyaligned==0) {
    for(r=alga[0]->chn->res;r;r=r->next) {
    	for(a=r->atm;a;a=a->next)
    	{
      		ss[0]=a->xyz[0]*l1+a->xyz[1]*m1+a->xyz[2]*n1;
      		ss[1]=a->xyz[0]*l2+a->xyz[1]*m2+a->xyz[2]*n2;
      		ss[2]=a->xyz[0]*l3+a->xyz[1]*m3+a->xyz[2]*n3;
      		for(i=0;i<3;i++) a->xyz[i]=ss[i];
    	}
	if(r->temp) {
                int i,j;
                for(i=0;i<r->nemp/3;i++) {
		   ss[0]=r->temp[i*3]*l1+r->temp[1+i*3]*m1+r->temp[2+i*3]*n1;
      		   ss[1]=r->temp[i*3]*l2+r->temp[1+i*3]*m2+r->temp[2+i*3]*n2;
      		   ss[2]=r->temp[i*3]*l3+r->temp[1+i*3]*m3+r->temp[2+i*3]*n3;
                   for(j=0;j<3;j++) r->temp[j+i*3]=ss[j];
		}
        }
    }
  }
  if(flag==1&&onlyaligned==1) {
   int k;
   for(k=0;k<ii;k++) {
	r=alga[k];
   	for(a=r->atm;a;a=a->next)
   	{
      		ss[0]=a->xyz[0]*l1+a->xyz[1]*m1+a->xyz[2]*n1;
      		ss[1]=a->xyz[0]*l2+a->xyz[1]*m2+a->xyz[2]*n2;
      		ss[2]=a->xyz[0]*l3+a->xyz[1]*m3+a->xyz[2]*n3;
      		for(i=0;i<3;i++) a->xyz[i]=ss[i];
    	}
        
	if(r->temp) {
                int i,j;
                for(i=0;i<r->nemp/3;i++) {
		   ss[0]=r->temp[i*3]*l1+r->temp[1+i*3]*m1+r->temp[2+i*3]*n1;
      		   ss[1]=r->temp[i*3]*l2+r->temp[1+i*3]*m2+r->temp[2+i*3]*n2;
      		   ss[2]=r->temp[i*3]*l3+r->temp[1+i*3]*m3+r->temp[2+i*3]*n3;
                   for(j=0;j<3;j++) r->temp[j+i*3]=ss[j];
		}
        }
        
   }	
  }
  else if(flag==2&&onlyaligned==1) {
    for(k=0;k<ii;k++) {
   	r=alga[k];
    	for(a=r->atm;a;a=a->next)
    	{
      		ss[0]=a->xyz[0]*l1+a->xyz[1]*m1+a->xyz[2]*n1;
      		ss[1]=a->xyz[0]*l2+a->xyz[1]*m2+a->xyz[2]*n2;
      		ss[2]=a->xyz[0]*l3+a->xyz[1]*m3+a->xyz[2]*n3;
      		for(i=0;i<3;i++) a->xyz[i]=ss[i];
    	}
	if(r->temp) {
                int i,j;
                for(i=0;i<r->nemp/3;i++) {
		   ss[0]=r->temp[i*3]*l1+r->temp[1+i*3]*m1+r->temp[2+i*3]*n1;
      		   ss[1]=r->temp[i*3]*l2+r->temp[1+i*3]*m2+r->temp[2+i*3]*n2;
      		   ss[2]=r->temp[i*3]*l3+r->temp[1+i*3]*m3+r->temp[2+i*3]*n3;
                   for(j=0;j<3;j++) r->temp[j+i*3]=ss[j];
		}
        }
    }
  }
  if(TRES.logg>1) cerr<<"final rmsd "<<sqrt(dmin/n)<<endl;
  setcenter(coob,algb,flg);
  setcenter(coob,alga,flg);
  delete [] coob;delete [] cooa;
  return sqrt(dmin/n);
}


  
float Stralg::superimpose(Chn *chna,Chn *chnb,int flg)
{
  Res *r;
  int n,na,nb;
  na=0; for(r=chna->res;r;r=r->next) na++;
  nb=0; for(r=chna->res;r;r=r->next) nb++;
  int nre=na*nb;
  n=max(na,nb);
  alga=new Res*[2*n];
  algb=new Res*[2*n];
  for(na=0;na<2*n;na++) {alga[na]=0;algb[na]=0;}

  Algn algn;

  chna->setseqcard();
  chnb->setseqcard();
  algn.setsequence(chna->seqcard,chnb->seqcard);
  algn.defalgn();
  
  int n1,n2;
  n1=0;
  n2=0;
  int slen=algn.getroutelength();
  char **result=algn.output(stderr);
  na=0;
  int i;
  for(i=0;i<slen;i++) {
	  Tres *a=TRES[result[0][i]];
	  Tres *b=TRES[result[1][i]];
          if(a) n1++;
          if(b) n2++;
          if(a&&b) {
	      alga[na]=(*chna)[n1-1];
              algb[na]=(*chnb)[n2-1];
	      //cerr<<na<<" "<<alga[na]->name<<alga[na]->id0<<" "<<algb[na]->name<<algb[na]->id0<<endl;
	      na++;
          }
  }
  for(i=0;i<3;i++) delete [] result[i];
  delete [] result;
  /*
  for(na=0;na<n;na++)
  {
    alga[na]=(*chna)[na];algb[na]=(*chnb)[na];
  }
  */
  cerr<<"the aligned residue: "<<na<<" "<<sqrt(na*na*1.0/nre)*100<<"%"<<endl;
  return superimpose(flg);
}


float Stralg::superimposesimple(Res *rr,int nn,Res *tt,int ee,int flg)
{
  Chn *chna,*chnb;
  Res *r,*t;
  int n,na,nb;
  chna=rr->chn;
  chnb=tt->chn;
  na=0; for(r=chna->res;r;r=r->next) na++;
  nb=0; for(r=chna->res;r;r=r->next) nb++;
  int nre=na*nb;
  n=max(na,nb);
  alga=new Res*[2*n];
  algb=new Res*[2*n];
  for(na=0;na<2*n;na++) {alga[na]=0;algb[na]=0;}

   
  na=0;
  r=rr;t=tt;
  while(r&&t&&r->id0<=nn&&t->id0<=ee) {
	 alga[na]=r;
	 algb[na]=t;
	 na++;	
	 r=r->next;
	 t=t->next;
  }
 
  cerr<<"the aligned residue: "<<na<<endl;
  return superimpose(flg);
}



float Stralg::simplermsd(Chn *chna,Chn *chnb,int flg) {

  float d=0;
  Res *r,*r1;

  int n=0;
  for(r=chna->res;r;r=r->next) {

	r1=(*chnb)[r->id0];
	if(r1==0) continue;
	if(r->atm==0||r->atm->next==0||r1->atm==0||r1->atm->next==0) continue;
	d+=TRES.distsqr(r->atm->next->xyz,r1->atm->next->xyz);
	n++;
  }
  return sqrt(d/n);
}

float Stralg::superimpose(Pdb *chna,Pdb *chnb,int flg)
{
  Res *r;
  Chn *chn;
  int n,na,nb;

  na=0; 
  for(chn=chna->chn;chn;chn=chn->next)
  for(r=chn->res;r;r=r->next) na++;

  nb=0; 
  for(chn=chnb->chn;chn;chn=chn->next)
  for(r=chn->res;r;r=r->next) nb++;

  n=min(na,nb);
  alga=new Res*[2*n];
  algb=new Res*[2*n];
  for(na=0;na<2*n;na++) {alga[na]=0;algb[na]=0;}
  for(na=0;na<n;na++)
  {
    alga[na]=chna->isres(na);algb[na]=chnb->isres(na);
  }
  return superimpose(flg);
}


void Stralg::randmz(Chn *chn,float x,float y,float z)
{
  Res *r;
  Atm *a;
  float l1,m1,n1,l2,m2,n2,l3,m3,n3;
  float s1,c1,s2,c2,s3,c3;

 //rotate
  cout<<x<<" "<<y<<" "<<z<<endl;
  s1=sin(x);c1=cos(x);
  s2=sin(y);c2=cos(y);
  s3=sin(z);c3=cos(z);
  l1=c2*c3-c1*s2*s3;l2=-c2*s3-c1*s2*c3;l3=s1*s2;
  m1=s2*c3+c1*c2*s3;m2=-s2*s3+c1*c2*c3;m3=-s1*c2;
  n1=s1*s3;n2=s1*c3;n3=c1;

  for(r=chn->res;r;r=r->next)
  for(a=r->atm;a;a=a->next)
  {
    s1=a->xyz[0]*l1+a->xyz[1]*m1+a->xyz[2]*n1;
    s2=a->xyz[0]*l2+a->xyz[1]*m2+a->xyz[2]*n2;
    s3=a->xyz[0]*l3+a->xyz[1]*m3+a->xyz[2]*n3;
    a->xyz[0]=s1;a->xyz[1]=s2;a->xyz[2]=s3;
  }
}


float *Stralg::center(Res **s,int flg)
{
  Res *r;
  Atm *am;
  float *coo;
  int k,i,j;

  coo=new float[3];
  coo[0]=coo[1]=coo[2]=0;
  i=0;
  k=0;
  while(s[k])
  {
     r=s[k];
     for(am=r->atm;am;am=am->next) 
     {
       if(flg==0&&am->tatm->id!=1) continue;  //only ca
       if(flg==1&&am->tatm->name[1]!='H'&&am->tatm->id>=4) continue;  //only backbone
       if(flg==1&&am->tatm->name[1]=='H'&&am->tatm->bond[0]->id>=3) continue;
       if(flg==2&&am->tatm->name[1]=='H') continue; //only heavy atoms
       for(j=0;j<3;j++) coo[j]+=am->xyz[j];
       i++;
     }
     k++;
  }

  if(i==0||k==0) {cerr<<" no structure segment"<<endl;return 0;}
  for(j=0;j<3;j++)  coo[j]= coo[j]/i;
  
  Chn *chn;
  Pdb *pdb;

  if(flag==1) {
    pdb=s[0]->chn->pdb;

    for(chn=pdb->chn;chn;chn=chn->next)
    for(r=chn->res;r;r=r->next) {
	if(r->temp) {
		int i,j;
		for(i=0;i<r->nemp/3;i++)
		 for(j=0;j<3;j++) r->temp[j+i*3]-=coo[j]; 
	}
    	for(am=r->atm;am;am=am->next)
    	{
       		for(j=0;j<3;j++) am->xyz[j]-=coo[j];
    	}
    }
  }
  else {

    chn=s[0]->chn;
    for(r=chn->res;r;r=r->next) {
	if(r->temp) {
		int i,j;
		for(i=0;i<r->nemp/3;i++)
		 for(j=0;j<3;j++) r->temp[j+i*3]-=coo[j]; 
	}
    	for(am=r->atm;am;am=am->next)
    	{
       		for(j=0;j<3;j++) am->xyz[j]-=coo[j];
    	}
    }
  }

  	
  return coo;
}

float *Stralg::getcenter(Res **s,int flg)
{
  Res *r;
  Atm *am;
  float *coo;
  int k,i,j;

  coo=new float[3];
  coo[0]=coo[1]=coo[2]=0;
  i=0;
  k=0;
  while(s[k])
  {
     r=s[k];
     for(am=r->atm;am;am=am->next) 
     {
       if(flg==0&&am->tatm->id!=1) continue;  //only ca
       if(flg==1&&am->tatm->name[1]!='H'&&am->tatm->id>=4) continue;  //only backbone
       if(flg==1&&am->tatm->name[1]=='H'&&am->tatm->bond[0]->id>=3) continue;
       if(flg==2&&am->tatm->name[1]=='H') continue; //only heavy atoms
       for(j=0;j<3;j++) coo[j]+=am->xyz[j];
       i++;
     }
     k++;
  }

  if(i==0||k==0) {cerr<<" no structure segment"<<endl;return 0;}
  for(j=0;j<3;j++)  coo[j]= coo[j]/i;
  
 
  return coo;
}

void Stralg::setcenter(float *coo, Res **s,int flg)
{
  Res *r;
  Atm *am;
 
  int k,i,j;

  Chn *chn;
  Pdb *pdb;
  int ii=0;while(alga[ii])ii++;
  if(flag==1&&onlyaligned==0) {
    pdb=s[0]->chn->pdb;

    for(chn=pdb->chn;chn;chn=chn->next)
    for(r=chn->res;r;r=r->next) {
	if(r->temp) {
		int i,j;
		for(i=0;i<r->nemp/3;i++)
		 for(j=0;j<3;j++) r->temp[j+i*3]+=coo[j]; 
	}
    	for(am=r->atm;am;am=am->next)
    	{
       		for(j=0;j<3;j++) am->xyz[j]+=coo[j];
    	}
    }
  }
  else if(onlyaligned==0) {

    chn=s[0]->chn;
    for(r=chn->res;r;r=r->next) {
	if(r->temp) {
		int i,j;
		for(i=0;i<r->nemp/3;i++)
		 for(j=0;j<3;j++) r->temp[j+i*3]+=coo[j]; 
	}
    	for(am=r->atm;am;am=am->next)
    	{
       		for(j=0;j<3;j++) am->xyz[j]+=coo[j];
    	}
    }
  }
  else if(flag==1&&onlyaligned==1) {
 
    int k;
    for(k=0;k<ii;k++) {
	r=alga[k];
	if(r->temp) {
		int i,j;
		for(i=0;i<r->nemp/3;i++)
		 for(j=0;j<3;j++) r->temp[j+i*3]+=coo[j]; 
	}
    	for(am=r->atm;am;am=am->next)
    	{
       		for(j=0;j<3;j++) am->xyz[j]+=coo[j];
    	}
    }
  }
  else if(onlyaligned==1) {
    int k;
    //chn=s[0]->chn;
    //for(r=chn->res;r;r=r->next) {
    for(k=0;k<ii;k++) {
	r=alga[k];
	if(r->temp) {
		int i,j;
		for(i=0;i<r->nemp/3;i++)
		 for(j=0;j<3;j++) r->temp[j+i*3]+=coo[j]; 
	}
    	for(am=r->atm;am;am=am->next)
    	{
       		for(j=0;j<3;j++) am->xyz[j]+=coo[j];
    	}
    }
  }
}
