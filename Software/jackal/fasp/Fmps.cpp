#include"source.h"

void Fmps::initial() {
 pdb=0;
 end=start=-1;
 segment=0;
 lat=new Lattice;
 step=0.17;
 cutoff=6.;
 direct=1;
 pdbout=1;
 outlet=0.5;
 flat=1;
 cycle=1000;
 cid='1';
 force=strdup("v");
 side=1;
}

Fmps::Fmps()
{
 initial();
}

 
Fmps::~Fmps()
{
 if(segment) {delete segment;segment=0;}
 if(lat)     {delete lat;lat=0;}
 if(force)   {delete [] force;force=0;}
}

void Fmps::ready()
{
  Chn *chn;
  Res *rr;
  Atm *aa;

  if(cid=='1')  cid=pdb->chn->id;

  //create another one segment class
	
  chn=(*pdb)[cid];
  chn->header();
  
  //find near residues
  for(chn=pdb->chn;chn;chn=chn->next)
  for(rr=chn->res;rr;rr=rr->next)
  for(aa=rr->atm;aa;aa=aa->next)
  aa->allnear(3,0);
}

int Fmps::create()
{
  Chn *chn; 
  Res *r,*r1,*r3,*r4;
  int n;

  chn=(*pdb)[cid];
  if(chn==0) {cerr<<"could not create segment"<<endl;exit(0);}

  if(start==-1) start=chn->res->id0;
  if(end==-1)   end=chn->lastres()->id0;

  r=(*chn)[start];
  r1=(*chn)[end];
  r3=(*chn)[start-1];
  r4=(*chn)[end+1];
 
  if(r3==0&&r4==0) 
  { 
     start=chn->res->id0;
     for(r1=chn->res;r1->next;r1=r1->next);
     end=r1->id0;direct=1;
     n=0;
  } 
  else if(r3==0) 
  {
    r=chn->res; 
    start=r->id0;
    direct=0;
    n=1;
  }
  else if(r4==0)  
  { 
    for(r1=chn->res;r1->next;r1=r1->next); 
    end=r1->id0;direct=1;
    n=0;
  } 
  else if(r3&&r4) 
  {
    n=-1; direct=1;
  }

  create(r,end);
  
  return n;
}

void Fmps::create(Res *s,int n)
{
  if(segment) {delete segment;segment=0;} 
  segment=new Pdb;
  segment->chn=new Chn;
  segment->chn->pdb=segment;
  segment->chn->create(s,n);
  segment->chn->id=s->chn->id;
  segment->chn->start=s->chn->start;
  segment->name=new char[100];
  segment->configure();
  Res *r1;
  for(Res *r=segment->chn->res;r;r=r->next) {
	r->id0+=s->id0;
	r->id+=s->id;
	//cerr<<"aaa  "<<r->name<<r->id0<<" "<<s->id0<<endl;
	//cerr<<"aaa  "<<r->name<<r->id<<" "<<s->id<<endl;
        if(r->temp) delete [] r->temp;
	r->nemp=9;
	r->temp=new float[9];
	r1=(*s->chn)[r->id0];
	TRES.copy(r1->temp,r->temp,9);
  }
  strcpy("unknown",segment->name);
}
 
float Fmps::rmsd(int flg)
{
Res *r,*rr;
Atm *a,*aa;
int i,j;
float d;
Chn *chn=(*pdb)[cid];

d=0;
j=0;
for(i=start;i<=end;i++)
{
  r=(*chn)[i];
  rr=(*segment->chn)[i];
  if(r==0||rr==0) continue;
  if(flg==2)
  {
   d+=TRES.distsqr(r->atm->xyz,rr->atm->xyz);
   d+=TRES.distsqr(r->atm->next->xyz,rr->atm->next->xyz);
   d+=TRES.distsqr(r->atm->next->next->xyz,rr->atm->next->next->xyz);
   j+=3;
  }
  else if(flg==1)
  {
   d+=TRES.distsqr(r->atm->next->xyz,rr->atm->next->xyz);
   j+=1; 
  }
  else if(flg==0)
  {
    for(a=r->atm;a;a=a->next)
    { 
      if(a->tatm->name[1]=='H') break;
      aa=(*rr)[a->tatm->id];
      if(aa==0) continue;
      d+=TRES.distsqr(a->xyz,aa->xyz);
      j++;
    }
  }
}
d=sqrt(d/j);
return d;
}

void Fmps::printhelp()
{

}

float Fmps::clashall() {

  lat->putoff();
  lat->flag=1;
  lat->grdsiz=2.;
  lat->radall=15.;
  lat->ready(pdb);
  lat->puton(pdb);

  pdb->setflg(1);
  
  float e=0;
  Chn *cchn;
  for(cchn=pdb->chn;cchn;cchn=cchn->next)
  {
    e+=clash(cchn->res,100000);
  }
  return e;
}
 
float Fmps::clashall(Res *from,int to) {

  lat->putoff();
  lat->flag=1;
  lat->grdsiz=2.;
  lat->radall=15.;
  lat->ready(pdb);
  lat->puton(pdb);
  pdb->setflg(-1);
  from->chn->setflg(from,to,1);

  float e=0;

  e=clash(from,to);

  return e;
}


float Fmps::clash(Res *s,int to)
{
  Res *r;
  float e;
  e=0;
  for(r=s;r;r=r->next)
  {
    if(r->id0>to) break;  
    e+=clash(r); 
  }
  return e;
}

float Fmps::clash(Res *s) 
{ 
  Atm *a; 
  float e=0;
  for(a=s->atm;a;a=a->next)
  {
     if(a->flag<0) continue;
     e+=clash(a);
  } 
  return e;
}

float Fmps::clash(Atm *aa0)
{
  Atm *aa1,*all[50];
  float dr,dn,e,d,elon,xr[3];
  int i,j,m,n,isp;
  float ta,tb,tc,tn;
  int mr; 
  i=4*TRES.smt;
  ta=TRES.smooth[i+0];
  tb=TRES.smooth[i+1];
  tc=TRES.smooth[i+2];
  tn=TRES.smooth[i+3];

  j=0;n=0;
  while(aa0->near[j])
  {
    if(lat->exist(aa0->near[j]))
    {  
       all[n++]=aa0->near[j];
       lat->putoff(aa0->near[j]);
    }
    j++;  
  }

  lat->getcell(aa0,cutoff);
  lat->cutoff(aa0,cutoff);
  e=0;
  
  for(m=0;m<lat->nget;m++)
  {
     if(strchr(force,'v')==0&&strchr(force,'V')==0) break;
     if(strchr(force,'V')) if(aa0->tatm->id!=1) continue;
     aa1=lat->obtain[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
     mr=0;
     d=aa0->tatm->eng->radius+aa1->tatm->eng->radius;
     j=aa0->tatm->hbond*aa1->tatm->hbond;
     if(j==2||j==3||j==6||j==9){d=3.0;mr=1;}
     if(j==2&&aa0->tatm->id<4&&aa1->tatm->id<4) mr=3;
     if(aa0->tatm->hbond==5) d-=aa0->tatm->eng->radius;
     if(aa1->tatm->hbond==5) d-=aa1->tatm->eng->radius;
     if(d<=0) continue;
     if(aa0->res->name=='C'&&aa1->res->name=='C')
     {
        if(aa0->tatm->id==5&&aa1->tatm->id==5) {d=2.0;mr=2;}
        else if(aa0->tatm->id==5&&aa1->tatm->id==4) d=3.;
        else if(aa0->tatm->id==4&&aa1->tatm->id==5) d=3.;
     }
     if(aa0->tatm->name[1]=='H'&&aa1->tatm->name[1]=='H') d=1.0;
     for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j])/d;
     dr=TRES.distsqr(xr);
     if(dr<0.1) dr=0.1; 
     else if(dr>4) continue;
     if(aa0->res->name=='C'&&aa1->res->name=='C')
     {
        if(aa0->tatm->id==5&&aa1->tatm->id==5)e+=-10.*exp(-(dr-1)*(dr-1)); 
     }
     elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     elon=ta*elon*exp(-dr*tb);
     dr=1/dr;
     dn=pow(dr,tn);  
     dn=elon*dn*(dn-tc);     
     isp=aa0->tatm->ispolar*aa1->tatm->ispolar;
     if(dn<0)
     {
      if(mr==1)       e+= dn*1.50;
      else if(mr==2)  e+= dn*5;
      else if(mr==3)  e+= dn*3.0;
      else if(isp==1) e+= dn*1.50;  
      else if(isp==2) e+= dn*1.25; 
      else if(isp==3) e+= dn*0.75; 
      else            e+= dn;
     }
     else
     {
      if(mr==1)       e+= dn*0.50;
      else if(mr==2)  e+= dn*0.20;
      else if(isp==1) e+= dn*0.50; 
      else if(isp==2) e+= dn*0.75;
      else if(isp==3) e+= dn*1.25;
      else            e+= dn;
     }
  }
  
  if(strchr(force,'E')||strchr(force,'e')||strchr(force,'H')||\
     strchr(force,'h')||strchr(force,'i')) 
  {
    lat->resonly();
    for(m=0;m<lat->nget;m++)
    {
      if(aa0->tatm->eng->charge==0) break;
      if(strchr(force,'E')==0&&strchr(force,'e')==0&&strchr(force,'i')==0) 
      {
        if(aa0->tatm->name[1]!='O'&&aa0->tatm->name[1]!='N') break;
      }
      for(aa1=lat->obtain[m]->atm->res->atm;aa1;aa1=aa1->next)
      {
         if(aa1->tatm->eng->charge==0) continue;
         if(strchr(force,'E')==0&&strchr(force,'e')==0||strchr(force,'i')==0) 
         {
           if(aa1->tatm->name[1]!='O'&&aa1->tatm->name[1]!='N') continue;
         }

         if(aa1->res==aa0->res)continue;
         if(aa1->res->id-aa0->res->id==1)
         {
           if((aa0->tatm->name[1]=='H'&&aa0->bond[0]->tatm->id<=4)||
              (aa0->tatm->name[1]!='H'&&aa0->tatm->id<=4))
           {
             continue;
             //if(aa1->tatm->id==0||aa1->tatm->id==1) continue;
             //if(aa1->tatm->name[1]=='H'&&aa1->bond[0]->tatm->id==0) continue;
             //if(aa1->tatm->name[1]=='H'&&aa1->bond[0]->tatm->id==1) continue;
           }
         }
         if(aa1->res->id-aa0->res->id==-1)
         {
           if((aa0->tatm->name[1]=='H'&&aa0->bond[0]->tatm->id<=4)||
              (aa0->tatm->name[1]!='H'&&aa0->tatm->id<=4))
           {
             continue;
             //if(aa1->tatm->id==2||aa1->tatm->id==3) continue;
           }
         }
         d=TRES.distance(aa1->xyz,aa0->xyz);
         dr=aa0->tatm->eng->radius+aa1->tatm->eng->radius;
         if(d<0.7*dr) d=0.7*dr;
         if(strchr(force,'E')||strchr(force,'e')||strchr(force,'i'))
         {
           e+=aa0->tatm->eng->charge*aa1->tatm->eng->charge/d/d*332.3/10;
         }
         if(strchr(force,'h')||strchr(force,'H'))
         {
           if(aa0->tatm->id<=4&&aa1->tatm->id<=4)      e+=aa0->ishbond(aa1,1.5,50,-0.5);
           else if(aa0->tatm->id<=4&&aa1->tatm->id>=4) e+=aa0->ishbond(aa1,1.5,50,-0.35);
           else if(aa0->tatm->id>=4&&aa1->tatm->id<=4) e+=aa0->ishbond(aa1,1.5,50,-0.35);
           else                                        e+=aa0->ishbond(aa1,1.5,50,-0.25);
         }
      }
    }
  }
 
  if(strchr(force,'e')||strchr(force,'h')||strchr(force,'i'))
  {
    i=0;
    for(m=0;m<lat->nget;m++) if(lat->obtain[m]->atm->res==aa0->res) i=1;
    if(i==0) {lat->obtain[lat->nget]=lat->atom[aa0->id0-lat->offset];lat->nget++;}
    if(aa0->area>0&&(aa0->tatm->name[1]=='O'||aa0->tatm->name[1]=='N')) 
    {
      dr=(aa0->tatm->eng->radius+1.4)*(aa0->tatm->eng->radius+1.4)*12.56;
      dr=dr*0.6;
      dr=aa0->area/dr;
      if(aa0->tatm->id<4)
      {
         if(dr>0.9) e+=-0.15;   
         else e+=-0.15*exp(2*dr-2);
      }
      else
      {
         if(dr>0.9) e+=-0.25;
         else e+=-0.25*exp(2*dr-2);
      }
    }
    for(m=0;m<lat->nget;m++)
    {
      if(strchr(force,'e')==0) break; 
      if(aa0->tatm->eng->charge==0) break;
      for(aa1=lat->obtain[m]->atm->res->atm;aa1;aa1=aa1->next)
      {
        if(aa1->area==0) continue;
        d=TRES.distance(aa1->temp,aa0->xyz);
        dr=aa0->tatm->eng->radius+aa1->tatm->eng->radius;
        if(d<0.7*dr) d=0.7*dr;
        e+=aa0->tatm->eng->charge*aa1->temp[6]/d/d*332.3/10;
      }
    }
  }

  for(i=0;i<n;i++)lat->puton(all[i]);

  return e;  
}

float Fmps::segmin(int fast)
{
  Res *r,*r0,*r1,*last,*last0;
  Atm *aa2,*aa0,*aa1,*all[50];
  float dn,dca,d,dmin;
  float xr[3],xr1[3],xr2[3];
  float *coeff[6],*unknown,offset;
  float **effco,*exp1st,*xyz1,*xyz,*xyz0,lamba[6],*exp2st;
  int i,j,k,m,mseq,nseq,*mcsc,mdel; 
  int rec,rec1,rec2,*order,nuse,endcon;
  int **affect,n;
  Qsort qsort;
  Rotate rot;
  float ta,tb,tc,tn;
  float dcut,di; 
  int isp;
  float *dist;
  int qut;
  int mr;
  int only,onlysum;
  Chn *chn=(*pdb)[cid];
  Chn *cchn;

  if(end-start+1<=3) return 0;

  //the first residue
  r0=(*chn)[start];
  //the number of backbone + sidechain torsional angle
  mseq=0; 
  for(r1=r0;r1;r1=r1->next) 
  {
    if(r1->id0>end) break;
    for(aa0=r1->atm;aa0;aa0=aa0->next)
    {
      if(aa0->tatm->rotate==0) continue;
      
      j=0;
      for(i=0;i<aa0->tatm->nbond;i++) if(aa0->bond[i]) j++;
      if(j<=1) continue;
      mseq++;
    }
  }

  mcsc=new int[mseq]; //if the rotate belong to mainchain or sidechain
  dist=new float[mseq];
  for(i=0;i<mseq;i++) dist[i]=1; 
  mseq=0;
  for(r1=r0;r1;r1=r1->next)
  {
    if(r1->id0>end) break;
    for(aa0=r1->atm;aa0;aa0=aa0->next)
    {
      if(aa0->tatm->rotate==0) continue;
      j=0;
      for(i=0;i<aa0->tatm->nbond;i++) if(aa0->bond[i]) j++;
      if(j<=1) continue;
      if(aa0->tatm->id<3)  mcsc[mseq++]=1;
      else                 mcsc[mseq++]=0;
    }
  }

  //inherent end constrain distance and end constrain existed or not

  endcon=1;

  
  if(direct) 
  {
     last=(*segment->chn)[end]; last0=(*chn)[end];
     if(last0->next==0) endcon=0;
     else
     {
        dn=TRES.distance((*last)[1]->xyz,(*last)[2]->xyz);
        dca=TRES.distance((*last0)[1]->xyz,(*last0)[2]->xyz);
        offset=fabs(dn-dca);
     }
  }
  else 
  {
     last=(*segment->chn)[start];last0=(*chn)[start];
     if(last0->id0==chn->res->id0) endcon=0;
     else
     {
        dn=TRES.distance((*last)[1]->xyz,(*last)[0]->xyz);
        dca=TRES.distance((*last0)[1]->xyz,(*last0)[0]->xyz);
        offset=fabs(dn-dca);
     }
  }
  
  //assing space 

  exp1st=new float[mseq+1];
  exp2st=new float[mseq];
  unknown=new float[mseq];
  order=new int[mseq];
  for(i=0;i<6;i++) coeff[i]=new float[mseq+1];

//preparing lattice

  nseq=0;
  for(cchn=pdb->chn;cchn;cchn=cchn->next)
  for(r=cchn->res;r;r=r->next)
  {
    for(aa0=r->atm;aa0;aa0=aa0->next)
    {
       if(r->id0<start||r->id0>end||cchn!=chn) {aa0->flag=-1;continue;}
       if(r->id0==start&&start>chn->res->id0)
       {
         if(aa0->tatm->id==0) {aa0->flag=-1;continue;}
         if(aa0->tatm->id==1) {aa0->flag=-1;continue;}
         if(strncmp(aa0->tatm->name," HN ",4)==0) {aa0->flag=-1;continue;}
       }
       if(r->id0==end&&r->next)
       {
         if(aa0->tatm->id>=1&&aa0->tatm->id<=3){aa0->flag=-1;continue;}
       }
       aa0->flag=nseq++; 
    }
  }

  nseq=nseq*3;


  affect=new int*[nseq/3]; 
  for(i=0;i<nseq/3;i++) affect[i]=new int[mseq];
  for(i=0;i<nseq/3;i++)
  for(j=0;j<mseq;j++) affect[i][j]=-1;

  
  effco=new float*[nseq];
  for(i=0;i<nseq;i++) effco[i]=new float[mseq]; 


  //store the chain conformation
  m=0;
  for(r=r0;r;r=r->next) if(r->id0<=end) m+=r->tres->number*3;
  xyz0=new float[m];
  xyz=new float[m];
  xyz1=new float[m];
  chn->transfer(xyz0,r0,end,0); 
  segment->chn->transfer(xyz,0);

  //preparing lattice
  lat->putoff();
  lat->flag=1;
  lat->grdsiz=2.;
  lat->radall=15.;
  lat->ready(pdb);
  for(cchn=pdb->chn;cchn;cchn=cchn->next)
  for(r=cchn->res;r;r=r->next)
  {
    if(r->id0>=start&&r->id0<=end&&cchn==chn) continue;
    lat->puton(r);
  }
 
  ta=TRES.smooth[TRES.smt*4+0];
  tb=TRES.smooth[TRES.smt*4+1];
  tc=TRES.smooth[TRES.smt*4+2];
  tn=TRES.smooth[TRES.smt*4+3];

  dmin=1000000; rec=0;rec1=0;rec2=0;qut=0;
  only=0;onlysum=0;
  for(rec=0;rec<cycle;rec++)
  {

    for(i=0;i<6;i++)
    {
      for(j=0;j<mseq+1;j++) coeff[i][j]=0; 
    }

    for(j=0;j<mseq;j++) {unknown[j]=0;exp2st[j]=0;dist[j]=0.1;}
   
    for(j=0;j<mseq+1;j++) exp1st[j]=0;  

    for(i=0;i<nseq;i++)
    for(j=0;j<mseq;j++)
    {
       effco[i][j]=0;
    }

    for(i=0;i<nseq/3;i++)
    for(j=0;j<mseq;j++)
    affect[i][j]=-1;

    chn->transfer(xyz0,r0,end,1);

    if(endcon)
    {
      if(direct) {aa0=last->atm->next;aa1=last0->atm->next; }
      else { aa0=last->atm;aa1=last0->atm; } 
      dn=TRES.distance(aa0->xyz,aa1->xyz);
      dca=TRES.distance(aa0->next->xyz,aa1->next->xyz);
      dcut=dn+dca-offset;
      
      if(dcut>40) dcut=40;
      if(dcut>0.4*(end-start)&&onlysum>10)  only=1;
      else                                  only=0;
      if(only) onlysum=0;
      else     onlysum++;
      for(i=0;i<3;i++)
      { 
        coeff[i][mseq]=aa0->xyz[i]-aa1->xyz[i];
        coeff[i+3][mseq]=aa0->next->xyz[i]-aa1->next->xyz[i];
      }
    }
     
    segment->chn->transfer(xyz1,0);
    chn->transfer(xyz1,r0,end,1);

    for(r=r0;r;r=r->next)
    {
      if(r->id0>end) break;
      lat->putoff(r);
      lat->puton(r);
    }
     
    m=0; 
    for(r=r0;r;r=r->next)
    {
      if(r->id0>end) break;

      for(aa2=r->atm;aa2;aa2=aa2->next)
      {
	if(aa2->tatm->rotate==0) continue;
        j=0;
        for(i=0;i<aa2->tatm->nbond;i++) if(aa2->bond[i]) j++;       
        if(j<=1) continue; //kick out no children's case torsion angle

        aa0=aa2->bond[0];
        for(j=0;j<3;j++)xr[j]=aa2->xyz[j]-aa0->xyz[j];
        d=sqrt(xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]);
        if(d==0) {m++;continue;} 
        for(j=0;j<3;j++)xr[j]=xr[j]/d; //unit vector of axis

        // create equation satifying end constraint
        rec1=0;
        if(aa2->tatm->id>3) goto re300;
        for(k=0;k<3;k++)
        {
          if(endcon==0) break;
          if(direct&&k==0)continue;
          else if(!direct&&k==2) continue; 

          for(j=0;j<3;j++)xr1[j]=(*last)[k]->xyz[j]-aa0->xyz[j]; //aa0-N
          d=xr1[0]*xr[0]+xr1[1]*xr[1]+xr1[2]*xr[2];
          for(j=0;j<3;j++)xr1[j]=xr1[j]-d*xr[j];
          dist[m]+=TRES.distance(xr1);
          xr2[0]=xr[1]*xr1[2]-xr[2]*xr1[1];
          xr2[1]=xr[2]*xr1[0]-xr[0]*xr1[2];
          xr2[2]=xr[0]*xr1[1]-xr[1]*xr1[0];
          for(j=0;j<3;j++) coeff[rec1*3+j][m]=xr2[j];
          rec1++;
        }
   
        re300: 
        //create coefficients due to angle change

        i=aa2->tatm->id;
        for(r1=r0;r1;r1=r1->next)
        { 
          if(only) break;
          if(r1->id0>end) break;
          if(mcsc[m]==0&&r1!=r) continue;
          if(direct&&r1->id0<r->id0) continue;
          else if(!direct&&r1->id0>r->id0)continue;

          for(aa1=r1->atm;aa1;aa1=aa1->next)
          {         
            if(aa1->flag<0) continue;
            if(mcsc[m]==0)
            {
              if(aa1->tatm->name[1]=='H')
              {
                if(aa1->tatm->bond[0]->id<i) 
                continue;
              }
              else
              {
                if(aa1->tatm->id<i) 
                continue;
              }
            }

            if(r1==r&&(i==1||i==2))
            {
              if(aa1->tatm->id==1) continue;

              if(i==1){if(aa1->tatm->id==0) continue;}
              else {if(aa1->tatm->id==2) continue;}

              if(direct)
              {
               if(i==1){if(strncmp(aa1->tatm->name," HN ",4)==0) continue;}
               else {if(aa1->tatm->id!=3) continue;}
              }
              else
              {
               if(i==1){if(strncmp(aa1->tatm->name," HN ",4)!=0) continue;}
               else {if(aa1->tatm->id==3) continue;}
              }
            }

            for(j=0;j<3;j++)xr1[j]=aa1->xyz[j]-aa0->xyz[j]; //aa0-N
            d=xr1[0]*xr[0]+xr1[1]*xr[1]+xr1[2]*xr[2];
            for(j=0;j<3;j++)xr1[j]=xr1[j]-d*xr[j];
            xr2[0]=xr[1]*xr1[2]-xr[2]*xr1[1];
            xr2[1]=xr[2]*xr1[0]-xr[0]*xr1[2];
            xr2[2]=xr[0]*xr1[1]-xr[1]*xr1[0];
            k=aa1->flag; 
            for(j=0;j<3;j++) effco[k*3+j][m]=xr2[j];
          }
        }
	m++;
      }        
    }
    
    //calculate the constrain equation
    if(endcon) 
    {
      d=0;
      for(i=0;i<6;i++)
      {
         lamba[i]=0;
         for(j=0;j<mseq;j++)
         {
           lamba[i]+=coeff[i][j]*coeff[i][j]; 
         }
         if(lamba[i]>d) d=lamba[i];
      }

      for(i=0;i<6;i++)
      {
         dn=sqrt(d/lamba[i]);
         for(j=0;j<mseq+1;j++)
         {
           coeff[i][j]=coeff[i][j]*dn;
         }
      }
      nuse=qsort.slnpd(coeff,order,6,mseq+1);
    }
    else
    {
      nuse=0;
      for(j=0;j<mseq;j++) order[j]=j;
    }

  
    //check effects;
    for(i=0;i<nseq/3;i++)
    {
       m=0;
       for(j=0;j<mseq;j++)
       {
         k=order[j];
         if(effco[i*3][k]==0&&effco[i*3+1][k]==0&&effco[i*3+2][k]==0) continue;
         affect[i][m]=j;
         m++;
       }
    }

    //expanding energy terms 
    for(r=r0;r;r=r->next)
    {
      if(only) break;
      if(r->id0>end) break; 
      for(aa0=r->atm;aa0;aa0=aa0->next)
      {
        if(aa0->flag<0) continue;
        j=0; mdel=0;
        while(aa0->near[j]) //put off neighboring atoms
        {
          if(lat->exist(aa0->near[j]))
          {
            all[mdel++]=aa0->near[j];
            lat->putoff(aa0->near[j]);  
          }
          j++;
        }


        lat->getcell(aa0,cutoff);
        lat->cutoff(aa0,cutoff);
        for(m=0;m<lat->nget;m++)
        {
          aa1=lat->obtain[m]->atm;
          if(aa0->flag<0&&aa1->flag<0) continue;
          if(aa0->flag>=0&&aa1->flag>=0){if(aa1->flag<aa0->flag) continue;}
          mr=0; 
          d=aa0->tatm->eng->radius+aa1->tatm->eng->radius;
          j=aa0->tatm->hbond*aa1->tatm->hbond;
          if(j==2||j==3||j==6||j==9) {d=3.0;mr=1;}
          if(aa0->tatm->hbond==5) d-=aa0->tatm->eng->radius;
          if(aa1->tatm->hbond==5) d-=aa1->tatm->eng->radius;
          if(d<=0) continue;  
          if(aa0->res->name=='C'&&aa1->res->name=='C')
          {
             if(aa0->tatm->id==5&&aa1->tatm->id==5) {d=2.0;mr=2;}
             else if(aa0->tatm->id==5&&aa1->tatm->id==4) d=3.;
             else if(aa0->tatm->id==4&&aa1->tatm->id==5) d=3.;
          }
          for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j])/d;
          dn=TRES.distsqr(xr);

          if(dn<0.1) dn=0.1; 
          else if(dn>4) continue;
          xr2[0]=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
          xr2[0]=ta*xr2[0]*exp(-dn*tb);
          dn=1/dn;
          xr2[1]=pow(dn,tn);  
          xr2[2]=xr2[1]*xr2[1];
          dca=tn*xr2[1]*dn*(tc-2*xr2[1])
             +tb*xr2[1]*(tc-xr2[1]);
          dca=dca*xr2[0];
          di=xr2[0]*xr2[1]*(xr2[1]-tc);
          //if(di>1000) di=1000;
          //if(dca>1000) dca=1000;
          isp=aa0->tatm->ispolar*aa1->tatm->ispolar;
          if(di<0)
          {
             if(mr==1)       {di=1.50*di;dca=1.50*dca;}
             else if(mr==2)  {di=5*di;dca=5*dca;}
             else if(isp==1) {di=1.50*di;dca=1.50*dca;}
             else if(isp==2) {di=1.25*di;dca=1.25*dca;}
             else if(isp==3) {di=0.75*di;dca=0.75*dca;}
          }
          else
          {
             if(mr==1)       {di=0.50*di;dca=0.50*dca;}
             else if(mr==2)  {di=0.20*di;dca=0.20*dca;}
             else if(isp==1) {di=0.50*di;dca=0.50*dca;}
             else if(isp==2) {di=0.75*di;dca=0.75*dca;}
             else if(isp==3) {di=1.25*di;dca=1.25*dca;}
          }
          exp1st[mseq]+=di; 
	  //cerr<<m<<" "<<aa0->res->name<<" "<<aa0->name<<" "<<aa0->xyz[1]<<" "<<aa1->res->name<<" "<<aa1->name<<" "<<aa1->xyz[1]<<" "<<di<<endl;	
          for(j=0;j<3;j++) xr1[j]=2*xr[j]*dca/d; 
          
          //calculating electrostatics
          /*
          if(fabs(aa0->tatm->eng->charge)>0.0001&&fabs(aa1->tatm->eng->charge)>0.001)
          {
            di=sqrt(dn/d/d);
            dca=aa0->tatm->eng->charge*aa1->tatm->eng->charge*di*332.3/20;
            di=-0.5*dca*di*di;
            exp1st[mseq]+=dca;
            for(j=0;j<3;j++) xr1[j]+=2*xr[j]*d*di;
          }            
          */         
          for(n=0;n<mseq;n++)
          {
             i=affect[aa0->flag][n];
             if(i==-1) break;
             k=order[i];
             for(j=0;j<3;j++)    xr2[j]=effco[aa0->flag*3+j][k];
             exp1st[i]+=xr1[0]*xr2[0]+xr1[1]*xr2[1]+xr1[2]*xr2[2];
          }
          for(n=0;n<mseq;n++)
          {
             if(aa1->flag<0) break;
             i=affect[aa1->flag][n];
             if(i==-1) break;
             k=order[i];
             for(j=0;j<3;j++) xr2[j]=-effco[aa1->flag*3+j][k];
             exp1st[i]+=xr1[0]*xr2[0]+xr1[1]*xr2[1]+xr1[2]*xr2[2];
          }
        }
        for(j=0;j<mdel;j++) lat->puton(all[j]);
      }
    }
    //exit(0);

    if(fast==0&&rec2>8&&qut) break;
    if(fast==1&&rec2>1&&qut) break; 
    if(qut&&fast==1&&fabs(exp1st[mseq]-dmin)<0.1&&dcut<outlet) break;

    di=1;
    for(i=0;i<mseq;i++) if(fabs(exp1st[i])>di) di=fabs(exp1st[i]);

    for(i=0;i<mseq;i++) 
    {
      d=fabs(exp1st[i]);
      if(d<1) d=1;
      if(flat%10) exp2st[i]=di;
      else          exp2st[i]=d;
      if((flat%100)/10) exp2st[i]=exp2st[i]*dist[order[i]]/((flat%100)/10);
      if((flat%1000)/100) exp2st[i]=exp2st[i]*exp(dcut/((flat%1000)/100)); 
      if(flat/1000)      exp2st[i]= exp2st[i]*(flat/1000);
      if(mcsc[order[i]]==0) exp2st[i]=fabs(exp1st[i]);
    }

    for(i=0;i<nuse;i++) lamba[i]=coeff[i][mseq];

    for(rec1=0;rec1<30;rec1++)
    {
       m=0;
       for(j=nuse;j<mseq;j++)
       {
         if(mcsc[order[j]]==0&&only) continue;
         d=exp1st[j];dn=2*exp2st[j];
         for(i=0;i<nuse;i++)
         {
           lamba[i]-=coeff[i][j]*unknown[j];
           dn+=coeff[i][j]*coeff[i][j]*exp2st[i];
           d+=exp1st[i]*coeff[i][j];
           d+=exp2st[i]*lamba[i]*coeff[i][j];
         }
         if(dn<0.001) d=unknown[j];
         else d=(-d/dn)*0.7+0.3*unknown[j];
         if(d>step)d=step;
         else if(d<-step) d=-step;
         for(i=0;i<nuse;i++) lamba[i]=lamba[i]+coeff[i][j]*d;
         if(fabs(d-unknown[j])>0.01*fabs(d))m=1;
         unknown[j]=d;
       }
       if(m==0) break;
    }
    
    for(i=0;i<nuse;i++) 
    {
       d=lamba[i];
       if(d>step)d=step;
       else if(d<-step) d=-step;
       unknown[i]=d;
    }
    
    d=0; //calculating decreased energy 
    for(i=0;i<mseq;i++)d+=exp1st[i]*unknown[i];

    if(dmin>exp1st[mseq]&&dcut<outlet&&only==0) //store lower energy conformation
    {
      segment->chn->transfer(xyz,0);dmin=exp1st[mseq];rec2=0;qut=1;
    }
    else if(dmin<exp1st[mseq]&&dcut<outlet&&only==0)
    {
      rec2++;
    }

    if(dcut<outlet&&rec==1&&only==0) {dmin=exp1st[mseq];qut=1;}

    for(i=0;i<mseq;i++)
    {
       k=order[i];
       exp2st[k]=unknown[i];
    }
    for(i=0;i<mseq;i++) unknown[i]=exp2st[i];

    cerr<<rec<<"energy:  "<<only<<"  "<<exp1st[mseq]<<"  "<<dcut<<"  "<<d<<endl;


    m=0; 
    for(r=r0;r;r=r->next)
    {
      if(r->id0>end) continue;
      r1=(*segment->chn)[r->id0];

      for(aa2=r->atm;aa2;aa2=aa2->next)
      {
        if(aa2->tatm->rotate==0) continue;
        j=0;
        for(i=0;i<aa2->tatm->nbond;i++) if(aa2->bond[i]) j++;
        if(j<=1) continue; //kick out no children's case torsion angle
        aa1=(*r1)[aa2->tatm->id]; 
        d=unknown[m]*57.29578;
        if(d==0) {m++;continue;}
        if(aa1->tatm->id<3)rot.rotate(aa1,d,last->id0,direct);
        else rot.rotate(aa1,d);
        m++;
      }
    }
    if(rec>cycle/2&qut==0) break;
  }
   
  chn->transfer(xyz0,r0,end,1); 
  if(qut)segment->chn->transfer(xyz,1);
  if(qut)setend(0);
  for(i=0;i<6;i++) delete [] coeff[i];
  delete [] unknown;
  delete [] order;
  delete [] exp2st;
  delete [] exp1st;
  delete [] xyz;
  delete [] xyz0;
  for(i=0;i<nseq;i++) delete [] effco[i];
  delete [] effco;
  for(i=0;i<nseq/3;i++) delete [] affect[i];
  delete [] affect;
  delete [] xyz1;delete [] mcsc;
  delete [] dist;
  if(qut==0) dmin=pow(10,25);
  return dmin;
}

int Fmps::setend(int flg)
{
//take care of O or HN atoms

  Res *last, *last0;
  Atm *aa0,*aa1;
  float d;
  Chn *chn=(*pdb)[cid];

  if(direct==1) {last0=(*chn)[end];last=(*segment->chn)[end];}
  else          {last0=(*chn)[start];last=(*segment->chn)[start];}

  if(direct==1&&last0->next==0) return 1;
  if(direct==0&&start>=chn->res->id0) return 1;

  if(direct==1)
  {
    d=TRES.distance((*last)[1]->xyz,(*last0)[1]->xyz)+
      TRES.distance((*last)[2]->xyz,(*last0)[2]->xyz);
    if(d>outlet) return 0;
   
    if(flg)
    {
       aa1=(*last)[1]; aa0=(*last0)[1];
       if(aa1&&aa0) aa1->transfer(aa0);
       aa1=(*last)[2]; aa0=(*last0)[2];
       if(aa1&&aa0) aa1->transfer(aa0);
       aa1=(*last)[3]; aa0=(*last0)[3];
       if(aa1&&aa0) aa1->transfer(aa0);
    }
    else
    {
       aa1=(*last)[3]; aa0=(*last0)[3];
       if(aa1&&aa0) aa1->transfer(aa0);
    }
  }
  else
  {
    d=TRES.distance((*last)[1]->xyz,(*last0)[1]->xyz)+
      TRES.distance((*last)[0]->xyz,(*last0)[0]->xyz);
    if(d>outlet) return 0;
   
    if(flg)
    {
       aa1=(*last)[0]; aa0=(*last0)[0];
       if(aa1&&aa0) aa1->transfer(aa0);
       aa1=(*last)[1]; aa0=(*last0)[1];
       if(aa1&&aa0) aa1->transfer(aa0);
       aa1=(*last)[" HN "]; aa0=(*last0)[" HN "];
       if(aa1&&aa0) aa1->transfer(aa0);
    }
    else
    {
       aa1=(*last)[" HN "]; aa0=(*last0)[" HN "];
       if(aa1&&aa0) aa1->transfer(aa0);
    }
  }
  return 1;
}

 
void Fmps::minseg(int flg)
{
  float xyz[10000];
  create();
  segmin(flg);
  Res *r0;
  Chn *chn=(*pdb)[cid];
  r0=(*chn)[start];
  segment->chn->transfer(xyz,0);
  cerr<<"rmsd   "<<rmsd(0)<<endl;
  chn->transfer(xyz,r0,end,1);
  chn->write("test");
}

void Fmps::scap() {

  if(side==0) return;
  Scap scap;
  scap.pdb=pdb;
  scap.colony=2;
  scap.ncolony=2;
  scap.nummore=100;
  scap.tormax=2;
  scap.arbt=0;
  scap.resultout=0;
  scap.pdb->setflgr(10000);
  scap.pdb->setflgr('A',-99999);
  scap.pdb->setflgr('P',-99999);
  scap.pdb->setflgr('G',-99999);
  scap.scpred(pdb->name,"",1);
  Chn *chn;
  Res *r,*r1;
  for(chn=pdb->chn;chn;chn=chn->next)
  for(r=chn->res;r;r=r->next)  { 
     for(r1=r->more;r1;r1=r1->more) {
	r1->next=0;
        r1->chn=0;
     }
     if(r->more) {
	delete r->more;
	r->more=0;
     }
  }
}

void Fmps::minall(int nn)
{
  int j;
  float xyz[30000];
  rusage begin, ended;

  start=0;end=10000;
  Chn *chn=(*pdb)[cid];
  chn->transfer(xyz,0);
  cerr<<"original total energy..."<<clashall()<<endl;

  cycle=10; 
  cutoff=3.0;
  getrusage(0, &begin);
  for(j=0;j<nn;j++)
  {
     tep=5;flag=1;
     scap();
     minimize(1);
     start=0;end=10000;
     create();
     segment->chn->transfer(xyz,1);
     cout<<"the final total energy..."<<clashall()<<endl;
     cout<<" rmsd: "<<rmsd(0)<<endl;
  }
  scap();

  getrusage(0, &ended);
  float ss;
  ss=(ended.ru_stime.tv_usec-begin.ru_stime.tv_usec)/1000000.+(ended.ru_stime.tv_sec-begin.ru_stime.tv_sec);
  ss+=(ended.ru_utime.tv_usec-begin.ru_utime.tv_usec)/1000000.+(ended.ru_utime.tv_sec-begin.ru_utime.tv_sec);
  cerr<<"time spent:"<<ss <<"s"<<endl;

}

void Fmps::minimize(int flg)
{
  Res *r0;
  float d;
  float xyz[10000];
  int i;
 
  Chn *chn=(*pdb)[cid];
  for(i=0;i<chn->number;i++)
  {
    start=i;end=i+tep;
    r0=(*chn)[start];
    if(segment) {delete segment;segment=0;}
    create();
    d=segmin(flg);
    cerr<<"minimize...."<<i<<" "<<end<<"  "<<d<<endl;
    segment->chn->transfer(xyz,0);
    chn->transfer(xyz,r0,end,1);
    if(segment) {delete segment;segment=0;}
    if(flag==1) i+=tep/2;
    if(flag==2) i+=tep-1;
  }
}

