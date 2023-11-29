#include"source.h"

void Loopy::initial() {
 pdb=0;
 end=start=-1;
 segment=0;
 lat=new Lattice;
 step=0.17;
 cutoff=6.;
 last=0;
 next=0;
 direct=1;
 arbt=200;
 id=0;
 near=0;
 init=0;
 part=20;
 revs=1;
 pdbout=1;
 outlet=0.5;
 flat=1;
 cycle=1000;
 secd='-';
 cid='-';
 part0=40;
 discut=3.;
 //
 chiangle=0;
 fapr=0;
 nomin=0;
 databaseonly=0;
 numclose=50;
}

Loopy::Loopy()
{
 initial();
}

Loopy::Loopy(Loopy *ss)
{
 initial();
 pdb=0;
 end=ss->end;
 start=ss->start;
 segment=0;
 step=ss->step; 
 cutoff=ss->cutoff;
 last=ss->last;
 next=0;
 direct=ss->direct;
 arbt=ss->arbt;
 id=ss->id;
 near=0;
 init=0;
 part=ss->part;
 revs=ss->revs;
 pdbout=ss->pdbout;
 outlet=ss->outlet;
 flat=ss->flat;
 cycle=ss->cycle;
 secd=ss->secd;
 cid=ss->cid;
 part0=ss->part0;
 //
 chiangle=ss->chiangle;
 fapr=ss->fapr;
 discut=ss->discut;
 databaseonly=ss->databaseonly;
 numclose=ss->numclose;
}

Loopy::~Loopy()
{
 if(segment) {delete segment;segment=0;}
 if(next)    {delete next;next=0;}
 if(lat)     {delete lat;lat=0;}
}

void Loopy::ready()
{
  Chn *cht,*chn;
  Res *rr;
  Atm *aa;

  if(cid=='1')  cid=pdb->chn->id;

  //create another one segment class
	
  chn=(*pdb)[cid];
  chn->header();
  pdb->next=new Pdb(pdb);
  pdb->next->configure();
  chn=(*pdb->next)[cid]; 
  chn->header();

  chn=(*pdb)[cid];
  next=new Loopy(this);
  next->id=this->id+1;
  next->pdb=pdb->next;
  next->last=this;
  next->next=next;

  //create segments
  create();
  next->create();
  next->direct=direct;

  //kick out not needed atoms
  cht=(*next->pdb)[next->cid];
  cht->transform((*cht)[next->start],next->end,3);
  cht=next->segment->chn;
  cht->transform((*cht)[next->start],next->end,3);

  //find near residues
  for(chn=pdb->chn;chn;chn=chn->next)
  for(rr=chn->res;rr;rr=rr->next)
  for(aa=rr->atm;aa;aa=aa->next)
  aa->allnear(3,0);

  for(chn=pdb->next->chn;chn;chn=chn->next)
  for(rr=chn->res;rr;rr=rr->next)
  for(aa=rr->atm;aa;aa=aa->next)
  aa->allnear(3,0);
}

int Loopy::create()
{
  Chn *chn; 
  Res *r,*r1,*r3,*r4;
  int n;

  chn=(*pdb)[cid];
  if(chn==0||end==-1) {cerr<<"could not create segment"<<endl;exit(0);}

  r=(*chn)[start];
  r1=(*chn)[end];
  r3=(*chn)[start-1];
  r4=(*chn)[end+1];
 
  if(r3==0&&r4==0) 
  { 
     cerr<<"warning!!! this is not a loop prediction problem.\n"; 
     cerr<<"it is a protein structure prediction"<<endl ;
     exit(0); 
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

void Loopy::create(Res *s,int n)
{
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
        if(r->temp) delete [] r->temp;
	r->nemp=9;
	r->temp=new float[9];
	r1=(*s->chn)[r->id0];
	TRES.copy(r1->temp,r->temp,9);
  }
  strcpy(segment->name,"unknown");//,segment->name);
}
 
void Loopy::create(char *seqn)
{
  Res *r,*r1;
  Rotate rot;
  segment=new Pdb;
  segment->chn=new Chn;
  segment->chn->pdb=segment;
  segment->chn->create(seqn);
  for(r=segment->chn->res;r;r=r->next)
  {
     r1=r->next;
     if(r1)
     {
       rot.link(r,r1,r1->id0,1);
     }
  }
  segment->chn->configure();
} 

void Loopy::randmz(int seed,int rand)
{
  Res *r;
  Atm *a;
  int j,i,n;
  Rotate rot;

  srandom(seed); //random initiializaton

  for(r=segment->chn->res;r;r=r->next)
  {
    for(a=r->atm;a;a=a->next)
    {
      if(a->tatm->name[1]=='H') break;
      if(a->tatm->rotate==0) continue;
      if(a->tatm->id>=4) continue;
      n=0;
      for(i=0;i<a->tatm->nbond;i++) if(a->bond[i]) n++; 
      if(n==1) continue;
      j=random()%rand;
      if(j%2==0) j=-j;
      if(direct==1) rot.rotate(a,j,1);
      else          rot.rotate(a,j,0);
    }
  }
}

void Loopy::randmzHelix(int seed,int rand)
{
  Res *r;
  Atm *a,*a0;
  int j,i,n;
  Rotate rot;
  float rt;
  srandom(seed); //random initiializaton

  Chn *chn=(*pdb)[cid];
  Res *last,*last0;
  if(direct==1) {
    last=(*segment->chn)[start];
    last0=(*chn)[start];
  }
  else {
    last=(*segment->chn)[end];
    last0=(*chn)[end];
  }

  
  for(r=segment->chn->res;r;r=r->next)
  {
    for(a=r->atm;a;a=a->next)
    {
      if(a->tatm->name[1]=='H') break;
      if(a->tatm->rotate==0) continue;
      if(a->tatm->id>=4) continue;
      n=0;
      for(i=0;i<a->tatm->nbond;i++) if(a->bond[i]) n++;
      if(n==1) continue;
      j=random()%rand;
      if(j%2==0) j=-j;
     
      if(direct==1&&r==segment->chn->res) {
          a0=(*last0)[a->tatm->id];
          rt=a0->dihedral(0);
      }
      else if(direct==0&&r==last) {
          a0=(*last0)[a->tatm->id];
          rt=a0->dihedral(0);
      }
      else rt=a->dihedral(0);

      float ee=0;
      if(a->tatm->id==1&&direct==1) ee=j-rt-60;
      else if(a->tatm->id==2&&direct==1) ee=j-rt-43;
      else if(a->tatm->id==1&&direct==0) ee=j+rt+60;
      else if(a->tatm->id==2&&direct==0) ee=j+rt+43;
      if(direct==1) rot.rotate(a,ee,1);
      else          rot.rotate(a,ee,0);
    }
  }
  //segment->chn->dihedral(stdout);
}

int Loopy::randmzChiangle(int nu)
{
  Res *r;
  Atm *a,*a0;
  int i,n;
  Rotate rot;
  float rt;
  Chn *chn=(*pdb)[cid];

  Chiangle *chia=chiangle->get(nu);
  if(chia==0) return 0; 
  int nn=0;
 
  Res *last,*last0;

  if(direct==1) {
    last=(*segment->chn)[start];
    last0=(*chn)[start];
  }
  else {
    last=(*segment->chn)[end];
    last0=(*chn)[end];
  }

  segment->chn->dihedral("s4"); 
  for(r=segment->chn->res;r;r=r->next)
  {
    for(a=r->atm;a;a=a->next)
    {
      if(a->tatm->name[1]=='H') break;
      if(a->tatm->rotate==0) continue;
      if(a->tatm->id>=4) continue;
      n=0;
      for(i=0;i<a->tatm->nbond;i++) if(a->bond[i]) n++;
      if(n==1) continue;
      
      if(direct==1&&r==segment->chn->res&&a->tatm->id==1) {
          //a0=(*last0)[a->tatm->id]; 	 	
          //rt=a0->dihedral(0);
	  float *w1,*w2,*w3,*w4;
	  if(r->isatmid(0)==0||r->isatmid(1)==0||r->isatmid(2)==0) continue;
	  w1=last0->temp+3;
	  w2=r->isatmid(0)->xyz;
	  w3=r->isatmid(1)->xyz;
	  w4=r->isatmid(2)->xyz;
          rt=TRES.dihedral(w1,w2,w3,w4);
      }
      else if(direct==1&&r->next==0&&a->tatm->id==2) {
	  
          if(databaseonly==0) continue;
	  float *w1,*w2,*w3;
          if(r->isatmid(0)==0||r->isatmid(1)==0||r->isatmid(2)==0) continue;
          w1=r->isatmid(0)->xyz;
          w2=r->isatmid(1)->xyz;
          w3=r->isatmid(2)->xyz;
          AtmGeom g;
         
	  float d1=TRES.distance(r->temp+3,w1);
	  g.addbounds(w3,d1,0.01);
	  d1=TRES.distance(r->temp,w1);
          g.addbounds(w2,d1,0.01); 
          d1=TRES.distance(r->temp+6,w1);
          g.addbounds(r->isatmid(3)->xyz,d1,0.01);     
          g.findcoo(); 
          rt=TRES.dihedral(w1,w2,w3,g.coo); 
	  if(TRES.logg>3) cerr<<"distance ..."<<TRES.distance(g.coo,r->isatmid(3)->xyz)<<" "<<TRES.distance(g.coo,w2)<<" "<<TRES.distance(g.coo,w3)<<" "<<rt<<endl;
	  //continue;
      }
      else if(direct==0&&r==last&&a->tatm->id==2) {
	  a0=(*last0)[a->tatm->id];
          rt=a0->dihedral(0); 
      }
      else rt=a->dihedral(0);
      
      //cerr<<r->name<<r->id0<<" "<<rt<<endl;
      float ee=0;
      if(a->tatm->id==1) ee=chia->getangle(nn,0);
      else if(a->tatm->id==2) ee=chia->getangle(nn,1);
      //if(r==segment->chn->res) {
 	/*if(a->tatm->id==1||a->tatm->id==2){
		cout<<a->name<<" "<<r->name<<r->id0<<" "<<ee<<" "<<rt<<endl;
	}*/
      //}
      if(direct==1)ee=-rt+ee;
      else         ee=rt-ee;
      if(direct==1) rot.rotate(a,ee,1);
      else          rot.rotate(a,ee,0);


      /*
      if(direct==1&&r->next==0&&a->tatm->id==2) {

          float *w1,*w2,*w3,*w4;
          if(r->isatmid(0)==0||r->isatmid(1)==0||r->isatmid(2)==0) continue;
          w1=r->isatmid(0)->xyz;
          w2=r->isatmid(1)->xyz;
          w3=r->isatmid(2)->xyz;
          AtmGeom g;

          float d1=TRES.distance(r->temp+3,w1);
          g.addbounds(w3,d1,0.01);
          d1=TRES.distance(r->temp,w1);
          g.addbounds(w2,d1,0.01);
          d1=TRES.distance(r->temp+6,w1);
          g.addbounds(r->isatmid(3)->xyz,d1,0.01);
          g.findcoo();
          rt=TRES.dihedral(w1,w2,w3,g.coo);
          if(TRES.logg>3)cerr<<" again distance ..."<<TRES.distance(g.coo,r->isatmid(3)->xyz)<<" "<<TRES.distance(g.coo,w2)<<" "<<TRES.distance(g.coo,w3)<<" "<<rt<<endl;
          //continue;
      }
      */


    }
    /*
    if(r->name!='G') {

		Atm *ca=r->isatm(" CA ");
		Atm *n=r->isatm(" N  ");
		Atm *c=r->isatm(" C  ");
		Atm *cb=r->isatm(" CB ");
		DistBound disb;
		float e=disb.getta(ca,n,c,cb);

		//if(e<0.0||e>70.0) {
			cout<<"the improper angle at "<<r->name<<r->id<<" is unconventional:"<<e<<endl;
		//}
    } 
    */
    nn++;
  } 
  segment->chn->dihedral("s3");
  return 1;
}


void Loopy::randmzSheet(int seed,int rand)
{
  Res *r;
  Atm *a,*a0;
  int j,i,n;
  Rotate rot;
  float rt;
  srandom(seed); //random initiializaton

  Res *last,*last0;
  Chn *chn=(*pdb)[cid];
  if(direct==1) {
    last=(*segment->chn)[start];
    last0=(*chn)[start];
  }
  else {
    last=(*segment->chn)[end];
    last0=(*chn)[end];
  }

  for(r=segment->chn->res;r;r=r->next)
  {
    for(a=r->atm;a;a=a->next)
    {
      if(a->tatm->name[1]=='H') break;
      if(a->tatm->rotate==0) continue;
      if(a->tatm->id>=4) continue;
      n=0;
      for(i=0;i<a->tatm->nbond;i++) if(a->bond[i]) n++;
      if(n==1) continue;
      j=random()%rand;
      if(j%2==0) j=-j;
   
      if(direct==1&&r==segment->chn->res) {
          a0=(*last0)[a->tatm->id];
          rt=a0->dihedral(0);
      }
      else if(direct==0&&r==last) {
          a0=(*last0)[a->tatm->id];
          rt=a0->dihedral(0);
      }
      else rt=a->dihedral(0);
      //cerr<<a->name<<a->id<<" "<<rt<<"  ";
      float ee=0;
      if(a->tatm->id==1&&direct==1) ee=j-rt-100;
      else if(a->tatm->id==2&&direct==1) ee=j-rt+130;
      else if(a->tatm->id==1&&direct==0) ee=j+rt+100;
      else if(a->tatm->id==2&&direct==0) ee=j+rt-130;
      if(direct==1) rot.rotate(a,ee,1);
      else          rot.rotate(a,ee,0);
      //rt=a->dihedral(0);

      //cerr<<a->name<<a->id<<" "<<rt<<endl;
    }
  }
}



float Loopy::closeg(Loopy *seg1,Loopy *seg2)
{ 
  Loopy *segen[2];
  Res *r,*last,*last0,*rr;
  Atm *aa[2],*aa0;
  float d,d1,dpre,cut;
  float xr[3],xr1[3],xr2[3];
  float *coeff[12],*unknown,lamba[12];
  int i,j,m,k,mitr,mseq,meg;
  int rec,rec1;
  Rotate rot; 

  for(i=0;i<12;i++) coeff[i]=0;
  unknown=0;
  last=0; 
  last0=0;
  mseq=0;

  segen[0]=seg1; segen[1]=seg2;
  if(segen[0]->end!=segen[1]->start) 
  {
    cerr<<" the two segments not contigent"<<endl;
    return -1;
  }

  if(segen[0]->segment) 
  {	  
    for(r=segen[0]->segment->chn->res;r;r=r->next)
    {
       if(r->id0<segen[0]->start||r->id0>segen[0]->end) continue;
       mseq+=2;
    }
    last=(*segen[0]->segment->chn)[segen[0]->end];
  }  
 
  if(segen[1]->segment) 
  {
    for(r=segen[1]->segment->chn->res;r;r=r->next)
    {
      if(r->id0<segen[1]->start||r->id0>segen[1]->end) continue;
      mseq+=2;
    }
    last0=(*segen[1]->segment->chn)[segen[1]->start]; 
  }
  
  if(segen[0]->segment==0&&segen[1]->segment) 
  { 
     return segen[1]->closeg(); 
  } 
  else if(segen[0]->segment&&segen[1]->segment==0) 
  { 
     return segen[0]->closeg(); 
  } 
  else if(segen[0]->segment==0&&segen[1]->segment==0) 
  { 
     cerr<<" no existing of two segments for connecting"<<end; 
     return -1; 
  } 

  unknown=new float[mseq];
  for(i=0;i<12;i++) coeff[i]=new float[mseq];

  rec=0;rec1=0;dpre=1000;
  
  for(rec=0;rec<50;rec++)
  {
    for(i=0;i<12;i++)
    {
      for(j=0;j<mseq;j++) coeff[i][j]=0; 
    }

    for(j=0;j<mseq;j++) unknown[j]=0;

    for(i=0;i<4;i++)
    {
      aa[0]=(*last)[i];
      aa[1]=(*last0)[i];
      for(j=0;j<3;j++) lamba[i*3+j]=aa[0]->xyz[j]-aa[1]->xyz[j];
    }

    d=TRES.distance(lamba)+TRES.distance(lamba+3)+
      TRES.distance(lamba+6)+TRES.distance(lamba+9);
    //cerr<<rec<<" "<<rec1<<" the total distance: "<<d<<endl;
     
    if(fabs(d-dpre)<0.01*d||rec>100||d<outlet) {cut=d;break;}
 
    dpre=d; 

    m=0;
    for(meg=0;meg<2;meg++)
    {
      if(meg==0) rr=last;
      else rr=last0;
      
      for(r=segen[meg]->segment->chn->res;r;r=r->next)
      {
        aa[0]=r->atm->next;  // CA
        aa[1]=aa[0]->next;   // C
        for(i=0;i<2;i++)
        {
	 if(aa[i]->tatm->rotate==0) {m++;continue;}
        
	  aa0=aa[i]->bond[0]; 
		  
	  if(aa0==0)  
	  {cerr<<"warning in closeg(int) aa0=0"<<endl;exit(0);} 

          for(j=0;j<3;j++)xr[j]=aa[i]->xyz[j]-aa0->xyz[j];
          d=sqrt(xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]);
          
	  if(d==0) {m++;continue;} 
         
	  for(j=0;j<3;j++)xr[j]=xr[j]/d; //unit vector of axis
          
          for(k=0;k<4;k++)
          {
            if(r==last)
            {
              if(i==0) {if(k==0||k==1) continue;}
              else     {if(k!=3) continue;}
            }
            else if(r==last0)
            {
              if(i==0) {if(k!=0) continue;}
              else     {if(k<=1&&k>=3) continue;} 
            }

            for(j=0;j<3;j++)xr1[j]=(*rr)[k]->xyz[j]-aa0->xyz[j]; //aa0-N
            d=xr1[0]*xr[0]+xr1[1]*xr[1]+xr1[2]*xr[2];
            for(j=0;j<3;j++)xr1[j]=xr1[j]-d*xr[j];
            xr2[0]=xr[1]*xr1[2]-xr[2]*xr1[1];
            xr2[1]=xr[2]*xr1[0]-xr[0]*xr1[2];
            xr2[2]=xr[0]*xr1[1]-xr[1]*xr1[0];
            for(j=0;j<3;j++) 
            { 
             if(rr==last) coeff[k*3+j][m]= xr2[j]; 
             else         coeff[k*3+j][m]=-xr2[j]; 
            }
          }
          m++;
        }        
      }
    }

    //create the lamba equation
    for(rec1=0;rec1<30;rec1++)
    {
      mitr=0;
      for(i=0;i<mseq;i++)
      {
        d=0;d1=1;
        for(j=0;j<12;j++) 
        {
          d1+=coeff[j][i]*coeff[j][i];
          lamba[j]-=unknown[i]*coeff[j][i];
          d+=lamba[j]*coeff[j][i];
        }
        if(d1==0) continue;
        d=-d/d1;
        d=0.7*d+0.3*unknown[i];
        if(d>step) d=step;
        else if(d<-step) d=-step;

        if(fabs(d-unknown[i])>0.2*fabs(d))mitr=1;
        for(j=0;j<12;j++)
        {
          lamba[j]+=d*coeff[j][i];
        }
        unknown[i]=d;
      }
      if(mitr==0) break;
    }

    //rotate the segment

    m=0;
    for(meg=0;meg<2;meg++)
    {
      for(r=segen[meg]->segment->chn->res;r;r=r->next)
      {
        if(r->id0<segen[meg]->start||r->id0>segen[meg]->end) continue;
        aa[0]=r->atm->next;
        aa[1]=aa[0]->next;
        for(i=0;i<2;i++)
        {
          if(aa[i]->tatm->rotate==0) {m++;continue;}
          d=unknown[m]*57.29578;
          if(d==0) {m++;continue;}
          if(meg==0) rot.rotate(aa[i],d,last->id0,1); 
          else       rot.rotate(aa[i],d,last0->id0,0);
          m++;
        }
      }
    }
    if(rec>100) break;
  }
  for(i=0;i<12;i++) if(coeff[i]){delete [] coeff[i]; coeff[i]=0;}
  if(unknown) {delete [] unknown; unknown=0;}
  if(cut>2) cut=-2;
  return cut;
}
   

void Loopy::writeseg(float **xyz_all,char *s,int nn)
{
  int i,all;
  FILE *fp;
  char nam[1000];

  all=0;
  while(xyz_all[all])all++;
  for(i=0;i<all;i++)
  {
     if(i>=nn) continue;
     segment->chn->transfer(xyz_all[i],1);
     sprintf(nam,"%s_%d.pdb",s,i);
     fp=fopen(nam,"w");
     fprintf(fp,"rmsd: %f\n",rmsd(2));
     segment->chn->write(fp);
     fflush(fp);
     fclose(fp);
  }
}

void Loopy::write(float **xyz_all,char *s,int nn)
{
  Res *r0,*r;
  int i,nseq,all;
  FILE *fp;
  float *xyz;
  char nam[1000];
  Chn *chn;

  chn=(*pdb)[cid];
  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;

  xyz=new float[nseq];
  chn->transfer(xyz,r0,end,0);

  all=0;
  while(xyz_all[all])all++;
  for(i=0;i<all;i++)
  {
     if(i>=nn) continue;
     chn->transfer(xyz,r0,end,1);
     segment->chn->transfer(xyz_all[i],1);
     sprintf(nam,"%s_%d.pdb",s,i);
     fp=fopen(nam,"w");
     fprintf(fp,"rmsd: %f\n",rmsd(2));
     chn->transfer(xyz_all[i],r0,end,1);
     pdb->write(fp);
     //chn->write(fp);
     fflush(fp);
     fclose(fp);
  }

  chn->transfer(xyz,r0,end,1);
  delete [] xyz;
}

float **Loopy::predt0()
{ 
  Res *r,*r0; 
  int nseq,k,i,j,kep;
  float **xyz_all,**xyz_out,*ent;
  Chn *chn; 
 
  chn=(*pdb)[cid];

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
  
  float *xyzorg = new float[nseq];
  chn->transfer(xyzorg,r0,end,0);

  xyz_out=new float*[revs*2+100];
  for(i=0;i<revs*2+100;i++) xyz_out[i]=0;
  xyz_out[0]=new float[nseq];
  chn->transfer(xyz_out[0],r0,end,0);
 
  if(end!=start) {
    //generate random initial conformation either based on ab-inito or database method
    //the total number of initial conformation would be arbt
    cerr<<"generating loop candidate conformation..."<<endl;
    if(secd=='h') {
    	xyz_all=next->zipper(0,100*discut*(end-start+1),1);
    }
    else if(secd=='e') {
    	xyz_all=next->zipper(0,100*discut*(end-start+1),1);
    }
    else {
	xyz_all=next->zipper(0,100*discut*(end-start+1),1);
    }
    /*
    Segen segen;
    segen.start=start;
    segen.end=end;
    segen.cid=cid;
    segen.pdb=pdb;
    segen.arbt=arbt;
    segen.part=part;
    segen.part0=part0;
    segen.ready();
    xyz_all=segen.next->zipper(discut*(end-start+1),1);
    segen.pdb=0;
    */	
  }
  else {
    xyz_all=new float*[2];
    xyz_all[0]=new float[nseq];
    chn->transfer(xyz_all[0],r0,end,0);
    xyz_all[1]=0;
  }
  if(end-start+1<=2) {
     ent=scap(xyz_all,"v");delete [] ent;
     if(end!=start) ent=ensort(xyz_all,"v");delete [] ent;
     if(end!=start&&fapr==0) {ent=trim(xyz_out,"shve",arbt);delete [] ent; }
     if(end!=start&&fapr==1) {ent=trim(xyz_out,"hv",arbt);delete [] ent; }
     kep=0;while(xyz_all[kep])kep++;
     for(i=0;i<revs&&i<kep;i++)
     {
       xyz_out[i+1]=xyz_all[i];
     }
     for(i=revs;i<kep;i++) delete [] xyz_all[i];
     delete [] xyz_all;
     xyz_all=copy(xyz_out);
     goto re1;
  }

  chn->transfer(xyzorg,r0,end,1);

  cerr<<"sort candidates according to energy..."<<endl;
  ent=next->trim(xyz_all,"v",next->part0);delete [] ent;
  //ent=next->trim(xyz_all,"v",arbt);delete [] ent;

  if(fapr==0) {
     next->cycle=200;next->outlet=0.3*(end-start);next->flat=300;
     cerr<<"minimizing loop candidates..."<<endl;	
     cerr<<"this step may take minutes to hours depending on the number of loop candidates..."<<endl;
     int ne=TRES.logg;
     TRES.logg=-1;
     ent=next->expand(xyz_all,0); delete [] ent;
     TRES.logg=ne;
  }


  next->cycle=200;next->outlet=0.3;next->flat=111;
  if(fapr==0) {
	cerr<<"minimizing loop candidates again..."<<endl;
	cerr<<"this step may take minutes to hours depending on the number of loop candidates..."<<endl;
	int ne=TRES.logg;
	TRES.logg=-1;
	ent=next->expand(xyz_all,0); delete [] ent;
	TRES.logg=ne;
  }
  else {
	cerr<<"minimizing loop candidates..."<<endl;
	cerr<<"this step may take minutes to hours depending on the number of loop candidates..."<<endl;
	int ne=TRES.logg;
        TRES.logg=-1;
        ent=next->expand(xyz_all,0); delete [] ent;
        TRES.logg=ne;
  }
  init=xyz_all;
  j=1;
  cycle=200;outlet=0.3;flat=111;
  
  kep=0;while(init&&init[kep])kep++;
  if(revs>kep) revs=kep;
  next->revs=revs;
  chn->transfer(xyzorg,r0,end,1); 
  for(i=0;i<revs;i++)
  {
    cerr<<"working on the "<<i<<"th outputs..."<<endl;
    xyz_all=predt(i);
    kep=0;while(xyz_all[kep])kep++;
    if(kep>0) xyz_out[j++]=xyz_all[0];
    for(k=1;k<kep;k++) delete [] xyz_all[k];
    delete [] xyz_all;
  }
  chn->transfer(xyzorg,r0,end,1);  

  xyz_all=copy(xyz_out);

  re1:

  if(pdbout)
  {
      int totfile;
      char namem[100];
      totfile=0;while(xyz_all[totfile])totfile++;
      r0=(*chn)[start];
      for(i=1;i<totfile;i++)
      {
          chn->transfer(xyz_all[i],r0,end,1);
	  if(totfile>2) {
           	sprintf(namem,"%s_loopy.%i.pdb",pdb->name,i);
	  }
	  else {
		sprintf(namem,"%s_loopy.pdb",pdb->name);
	  }
	  cerr<<"output final result..."<<namem<<endl;
          pdb->write(namem);
      }   
  }

  chn->transfer(xyzorg,r0,end,1);  

  if(fapr==0) {
    ent=trim(xyz_out,"v",arbt,1);delete [] ent;
    ent=trim(xyz_out,"sv",arbt,1);delete [] ent;
    ent=trim(xyz_out,"shv",arbt,1);delete [] ent;
    ent=trim(xyz_out,"shve",arbt,1);delete [] ent;
  }
 
  kep=0;while(init&&init[kep])kep++;
  for(k=0;k<kep;k++) delete [] init[k];
  delete [] init;init=0;
  kep=0;while(xyz_out[kep])kep++;
  for(k=0;k<kep;k++) delete [] xyz_out[k];
  delete [] xyz_out;

  chn->transfer(xyzorg,r0,end,1);

  if(xyzorg)delete [] xyzorg;xyzorg=0;

  //delete the first one
  kep=0;while(xyz_all&&xyz_all[kep])kep++;
  if(kep&&xyz_all[0]) {delete [] xyz_all[0];xyz_all[0]=0;}
  int ii; 
  int mm=0;
  for(ii=0;ii<kep;ii++) {
	if(xyz_all[ii]==0) continue;
	xyz_all[mm++]=xyz_all[ii];
  }
  xyz_all[mm]=0;
  return xyz_all;
}

float **Loopy::predt(int flg)
{
Res *r,*r0;
float **xyz_all,*xyz,*ent;
int i,j,kep,nseq;
float close;
int cop;

Chn *chn=(*pdb)[cid]; 

r0=(*chn)[start];
nseq=0;
for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;

xyz=new float[nseq];
chn->transfer(xyz,r0,end,0);

close=0.1*(end-start);

//choose needed colony

kep=0;while(init[kep])kep++;
xyz_all=new float*[2*kep];
for(i=0;i<2*kep;i++)xyz_all[i]=0;
j=0;
cerr<<flg<<"th outputs are based from loop candidates between: "<<flg*kep/revs<<"-"<<(flg+1)*kep/revs<<endl;
for(i=flg*kep/revs;i<(flg+1)*kep/revs;i++)
{
  xyz_all[j]=new float[nseq];
  TRES.copy(init[i],xyz_all[j],nseq);
  j++;
  if(j>=10&&flg>0) break;
}

//side-chain

cerr<<"add side-chains to candidate conformations..."<<endl;
int ne=TRES.logg;
TRES.logg=-1;
ent=scap(xyz_all,"v");delete [] ent;
TRES.logg=ne;

if(fapr==0) cerr<<"minimizing again..."<<endl;
if(fapr==0) {
	ne=TRES.logg;
	TRES.logg=-2;
	ent=minimize(xyz_all,1);delete [] ent;
	TRES.logg=ne;
}
//filter on vw 
cerr<<"discarding bad candidates..."<<endl;
ent=trim(xyz_all,"vw",0);delete [] ent;

if(fapr==1) {
  chn->transfer(xyz,r0,end,1);
  if(xyz) {delete [] xyz; xyz=0;}
  return xyz_all;
}


//filter on next->vw
ent=next->trim(xyz_all,"vw",0); delete [] ent;

kep=0;while(xyz_all[kep])kep++;
if(kep>1.3*part) {ent=trim(xyz_all,"shvw",part); delete [] ent;}

cop=1; 
 
for(i=0;i<2;i++)
{
  ne=TRES.logg;
  TRES.logg=-1;
  cerr<<"loop fusion to generate more candidates..."<<endl;
  xyz_all=generate(xyz_all,cop,close);
  if(TRES.logg>3)cerr<<"***cycle***: "<<i<<"  "<<j<<endl;
  TRES.logg=ne;
}

for(i=0;i<2;i++)
{
  ne=TRES.logg;
  TRES.logg=-1;
  xyz_all=next->generate(xyz_all,cop,close);
  if(TRES.logg>3)cerr<<"***cycle for backbone only***: "<<i<<" "<<j<<endl;
  TRES.logg=ne;
}

//side-chain
cerr<<"do side-chain prediction again..."<<endl;
ne=TRES.logg;
TRES.logg=-2;
ent=scap(xyz_all,"v");delete [] ent;
cerr<<"minimizing again..."<<endl;
ent=minimize(xyz_all,1);delete [] ent;
TRES.logg=ne;
if(part>10) part=part/2;
for(i=0;i<2;i++)
{
  cerr<<"loop fusion to generate more candidates..."<<endl;
  xyz_all=generate(xyz_all,cop,close);
  if(TRES.logg>3)cerr<<"***cycle***: "<<i<<"  "<<j<<endl;
}
part=next->part;

cerr<<"discarding bad candidates..."<<endl;
ent=trim(xyz_all,"vw",arbt);delete [] ent;

chn->transfer(xyz,r0,end,1);
if(xyz) {delete [] xyz; xyz=0;}
return xyz_all;
}

int Loopy::setend(int flg)
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

float Loopy::closeg()
{
  Res *r0,*r,*r1,*last,*last0;
  Atm *aa[2],*aa0;
  float dn,dca,d;
  float xr[3],xr1[3],xr2[3],cut;
  float *coeff[6],*unknown,off,offset;
  float lamba[6];
  int i,j,k,m,mseq;
  int rec,rec1,*order,nuse;
  Qsort qsort;
  Rotate rot;
  Chn *chn=(*pdb)[cid];

  unknown=0;for(i=0;i<6;i++)coeff[i]=0;
//middle exist!    

  r0=(*segment->chn)[start];
  mseq=0;
  
  if(segment==0) {cerr<<"no closing segment"<<endl;return -1;}

  if(direct==1) 
  {
    r1=(*chn)[end];
    if(r1->next==0) 
    {cerr<<"no closing segment"<<endl;return -1;}
  }
  else if(direct==0) 
  {
    if(start==chn->res->id0)
    {cerr<<"no closing segment"<<endl;return -1;}
  }

  for(r1=r0;r1;r1=r1->next) 
  {
    if(r1->id0>end) break;
    mseq+=2;
  }

  if(direct) 
  {
    last=(*segment->chn)[end]; 
    last0=(*chn)[end];
    dn=TRES.distance((*last)[1]->xyz,(*last)[2]->xyz);
    dca=TRES.distance((*last0)[1]->xyz,(*last0)[2]->xyz);
    offset=fabs(dn-dca);
  }
  else 
  {
    last=(*segment->chn)[start];
    last0=(*chn)[start];
    dn=TRES.distance((*last)[1]->xyz,(*last)[0]->xyz);
    dca=TRES.distance((*last0)[1]->xyz,(*last0)[0]->xyz);
    offset=fabs(dn-dca);
  }

  if(last==0||last0==0) 
  {
    cerr<<"no segment close process!"<<endl;return -1;
  }

  unknown=new float[mseq];
  order=new int[mseq];
  for(i=0;i<6;i++) coeff[i]=new float[mseq+1];

  rec=0;rec1=0;off=1000;cut=1000;

  for(rec=0;rec<numclose;rec++)
  {
    for(i=0;i<6;i++)
    {
      for(j=0;j<mseq+1;j++) coeff[i][j]=0; 
    }

    for(j=0;j<mseq;j++) unknown[j]=0;
   
    if(direct) { aa[0]=last->atm->next; aa[1]=last0->atm->next; }
    else       { aa[0]=last->atm;       aa[1]=last0->atm;       } 
    for(i=0;i<3;i++)
    {
       coeff[i][mseq]=aa[0]->xyz[i]-aa[1]->xyz[i];
       coeff[i+3][mseq]=aa[0]->next->xyz[i]-aa[1]->next->xyz[i];
    }
    dn=TRES.distance(aa[0]->xyz,aa[1]->xyz);
    dca=TRES.distance(aa[0]->next->xyz,aa[1]->next->xyz);
    /* 
    fprintf(stderr,"%3i %3i %c%d: %s:%8.3f  %s:%8.3f\n",rec,rec1,
            aa[1]->res->name,aa[1]->res->id,
            aa[1]->name,dn,aa[1]->next->name,dca);
    */
    d=dn+dca;
    //cerr<<rec<<" "<<d<<endl;
    if(fabs(d-off)<0.01*d||d<outlet) {cut=d;break;} 
     
    if(fabs(off-d)<0.001) {cut=d; break;}

    off=d;

    m=0; 
    for(r=r0;r;r=r->next)
    {
      if(r->id0>end) break;
      aa[0]=r->atm->next;  // CA
      aa[1]=aa[0]->next;   // C
      for(i=0;i<2;i++)
      {
	if(aa[i]->tatm->rotate==0) {m++;continue;}

        aa0=aa[i]->bond[0];
        for(j=0;j<3;j++)xr[j]=aa[i]->xyz[j]-aa0->xyz[j];
        d=sqrt(xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]);
        if(d==0) {m++;continue;} 
        for(j=0;j<3;j++)xr[j]=xr[j]/d; //unit vector of axis

        // create equation satifying end constraint
        rec1=0;
        for(k=0;k<3;k++)
        {
          if(direct&&k==0)continue;
          else if(!direct&&k==2) continue; 

          for(j=0;j<3;j++)xr1[j]=(*last)[k]->xyz[j]-aa0->xyz[j]; //aa0-N
          d=xr1[0]*xr[0]+xr1[1]*xr[1]+xr1[2]*xr[2];
          for(j=0;j<3;j++)xr1[j]=xr1[j]-d*xr[j];
          xr2[0]=xr[1]*xr1[2]-xr[2]*xr1[1];
          xr2[1]=xr[2]*xr1[0]-xr[0]*xr1[2];
          xr2[2]=xr[0]*xr1[1]-xr[1]*xr1[0];
          for(j=0;j<3;j++) coeff[rec1*3+j][m]=xr2[j];
          rec1++;
        }
	m++;
      }        
    }
    
    //calculate the constrain equation

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
    
    for(i=0;i<nuse;i++) lamba[i]=coeff[i][mseq];     
    for(rec1=0;rec1<30;rec1++)
    { 
      m=0;
      for(j=nuse;j<mseq;j++)
      {
        k=order[j];d=0;dn=1;
        for(i=0;i<nuse;i++)
        {
          lamba[i]-=coeff[i][j]*unknown[order[j]];
          dn+=coeff[i][j]*coeff[i][j];
          d+=lamba[i]*coeff[i][j];
        }
        d=(-d/dn)*0.7+0.3*unknown[k];
        if(d>step)d=step;
        else if(d<-step) d=-step;
        for(i=0;i<nuse;i++) lamba[i]=lamba[i]+coeff[i][j]*d;
        if(fabs(d-unknown[k])>0.3*fabs(d))m=1; 
        unknown[k]=d;
      }
      if(m==0) break;
    }

    for(i=0;i<nuse;i++)
    {
       d=lamba[i];
       if(d>step) d=step;
       else if(d<-step) d=-step;
       unknown[order[i]]=d;
    }

    m=0;
    for(r=r0;r;r=r->next)
    {
      if(r->id0>end) break;
      aa[0]=r->atm->next;
      aa[1]=aa[0]->next;
      for(i=0;i<2;i++)
      {
        if(aa[i]->tatm->rotate==0) {m++;continue;}
        d=unknown[m]*57.29578;
        if(d==0) {m++;continue;}
        rot.rotate(aa[i],d,last->id0,direct);
        m++;
      }
    }

    if(rec>100) break;
  }

  for(i=0;i<6;i++) if(coeff[i]) {delete [] coeff[i];coeff[i]=0;}
  if(unknown) {delete [] unknown;unknown=0;}
  /*
  if(databaseonly==1) {
	//cout<<"the distance at ends: "<<cut<<endl;
	return 0.1;
  }	
  */
  if(cut<=0) return cut; 
  else if(cut>offset+outlet+0.1) return -2;
  else return cut;
}



float Loopy::rmsd(int flg)
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

void Loopy::printhelp()
{ 
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"loopy is a protein loop conformation prediction with the following capabilities:\n");
fprintf(stderr,"A. predict loop conformation\n");
fprintf(stderr,"B. mutation,insertion and deletion in loop region\n");
fprintf(stderr,"\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"loopy -prm num -obj num-num -seed num -ini num -cid char -out num -sec char -heta num file.pdb\n");
fprintf(stderr,"-prm   force parameters. default 1\n");
fprintf(stderr,"-prm 1: CHARMM22 with all atom model\n");
fprintf(stderr,"-prm 2: AMBER with all atom model\n");
fprintf(stderr,"-prm 3: CHARMM with heavy atom model\n");
fprintf(stderr,"-prm 4: AMBER with heavy atom model\n");
fprintf(stderr,"-obj num-num: start and end id of loops to be predicted\n");
fprintf(stderr,"-seed num: random seed number\n");
fprintf(stderr,"-ini: number of initial conformations generated.default is 1000\n");
fprintf(stderr,"-cid: chain id. default is the chain with empty space\n");
fprintf(stderr,"-out: number of outputs. default is 1\n");
fprintf(stderr,"-res: sequence of the segment. optional.\n");
fprintf(stderr,"-fast: fast mode. default is 1.\n");
fprintf(stderr,"-fast 0: fast mode off.\n");
fprintf(stderr,"-fast 1: fast mode on. \n");
fprintf(stderr,"-heta: non standard amino acid residues discarded or not. default is 1\n");
fprintf(stderr,"-heta 0: non standard residues discarded\n");
fprintf(stderr,"-heta 1: non standard residues considered\n");
fprintf(stderr,"file.pdb   pdb file\n");
//printf("file_loopy.chi  torsion angle of predefined loop conformations. optional.\n");
}

/*
void Loopy::printoldhelp()
{ 
printf("scap -f -h -o -i -c -t -n -s -c -r pdbfile loop_chi_database\n");
printf("-f     fast prediction, default is slow method\n");
printf("-h     print out help file.optional.optional.\n");
printf("-o     the start and end residue number where loop to be predicted. e.g., -o=20-30, predict loop between 20-30\n");
printf("-i     the number of initial conformations generated. e.g., -i=30, generate 30 random initial conformations.\n");
printf("-c     the chain id where the loop to be predicted. default meaning the first chain.optional. default: the first chain\n");
printf("-t     topology. e.g., -t=a, mean charm22 all atom model. -t=b, mean charmm22 heavy atom model.\n");
printf("       -t=A, mean amber94 all atom model, -t=B, mean amber94 heavy atom model.optional.default:-t=b\n");
printf("-n     number of outputs.e.g., -n=3, mean output 3 loop predictions.optional, default: -n=1\n");
printf("-s     the secondary structure of protein segment. -s=c, loop, -s=h, helix, -s=e, sheet.optional.default:-s=c\n");
printf("-r     the new residue sequnece. e.g., -r=KKK, change the segment to residue K,K,K. optional\n");
printf("pdbfile   pdb file\n");
printf("loop_chi_database  segment conformations. optional. \n");
}
*/

void Loopy::randmz(int seed)
{
  Res *r,*r1,*r2,*r3,*r5;
  Res *last;
  Atm *a,*a1;
  int i,j,m;
  Chn *segp;
  Rotate rot;
  float xr[3],xr1[3],dn,dca;
  Chn *chn=(*pdb)[cid];

  //for the first segment
  if(direct==0) goto re200;
  
  srandom(seed+1); //skip random
 
  segp=segment->chn;

  r3=(*segp)[start];
  for(r=r3;r;r=r->next)
  {
    if(r->id0>end) break;
    
    r1=TRES.rotamer->chn->isres(r->name,0);
    j=random()%r1->nummore;
    r2=r1->ismore(j);    
    if(r==r3)
    {
       if(r3->id0==segp->res->id0) r5=(*chn)[start-1];
       else r5=(*segp)[start-1];

       if(r5==0)
       {
          cerr<<"strange in subroutine randmz..."<<endl;
       }
       else
       {
          r->transfer(r2);
          r->copytemp(r2);
          rot.link(r5,r,r->id0,1); 
       }

       if(r3->id0==segp->res->id0) 
       {
          r5=(*chn)[start];
          r->atm->transfer(r5->atm);
          r->atm->next->transfer(r5->atm->next);
          a=r->isatm(" HN ");
          a1=r5->isatm(" HN ");
          if(a&&a1)a->transfer(a1);
         
          //if(TRES.distance(r->atm,r5->atm)<0.5) r->atm->transfer(r5->atm);
          //if(TRES.distance(r->atm->next,r5->atm->next)<0.5) r->atm->next->transfer(r5->atm->next);
          //a=r->isatm(" HN ");
          //a1=r5->isatm(" HN ");
	   
          //if(a&&a1&&TRES.distance(a,a1)<0.5)a->transfer(a1);
       }
    }
    else
    { 
       r->transfer(r2);
       r->copytemp(r2);
       rot.link(r3,r,r->id0,1);       
    }
    r3=r;
  }
   
  r1=(*segp)[end];r2=(*chn)[end];
  a= r1->atm->next; a1=r2->atm->next;
  dn=TRES.distance(a->xyz,a->next->xyz);
  dca=TRES.distance(a1->xyz,a1->next->xyz);
  for(i=0;i<3;i++) xr[i]=(a->next->xyz[i]-a->xyz[i])/dn*(dca-dn);
  for(i=0;i<3;i++) xr1[i]=a->next->next->xyz[i]-a->next->xyz[i];
  for(i=0;i<3;i++) a->next->xyz[i]+=xr[i];
  for(i=0;i<3;i++) a->next->next->xyz[i]=a->next->xyz[i]+xr1[i];


  return;

  re200:

  //for the next segment;

  segp=segment->chn;
  if(segp->res==0) return;

  last=(*segp)[end];
  r3=last;

  for(m=end;m>=start;m--)
  {
    r=(*segp)[m];
    if(r==0) continue;

    r1=TRES.rotamer->chn->isres(r->name,0);
    j=random()%r1->nummore;
    r2=r1->ismore(j);
    if(r==last)
    {
       if(last->next) r5=last->next;
       else           r5=(*chn)[end+1];
       if(r5==0)
       {
           cerr<<"strange in subroutine randmz..."<<endl;
       }
       else
       {
           r->transfer(r2);
           r->copytemp(r2);
           rot.link(r,r5,r->id0,0);
       }
       r5=(*chn)[r->id0]; 
       r->atm->next->transfer(r5->atm->next); 
       r->atm->next->next->transfer(r5->atm->next->next);
       r->atm->next->next->next->transfer(r5->atm->next->next->next);
    }
    else
    {
       r->transfer(r2);
       r->copytemp(r2);
       rot.link(r,r3,r->id0,0);
    }
    r3=r;
  }
  
  r1=(*segp)[start];r2=(*chn)[start];
  dn=TRES.distance(r1->atm->xyz,r1->atm->next->xyz);
  dca=TRES.distance(r2->atm->xyz,r2->atm->next->xyz);
  for(i=0;i<3;i++) xr[i]=(r1->atm->xyz[i]-r1->atm->next->xyz[i])/dn*(dca-dn);
  a=(*r1)[" HN "];
  if(a) { for(i=0;i<3;i++) xr1[i]=a->xyz[i]-r1->atm->xyz[i]; }
  for(i=0;i<3;i++) r1->atm->xyz[i]+=xr[i];
  if(a) { for(i=0;i<3;i++) a->xyz[i]=r1->atm->xyz[i]+xr1[i]; }
}

float Loopy::clash(Res *s,int to,char *force)
{
  Res *r;
  float e;
  e=0;
  for(r=s;r;r=r->next)
  {
    if(r->id0>to) break;  
    e+=clash(r,force); 
  }
  return e;
}
float Loopy::clash(Res *s,int from,int to,char *force)
{
  Atm *a;
  float e;
  e=0;
  for(a=s->atm;a;a=a->next)
  {
     if(a->tatm->name[1]=='H') break;
     if(a->tatm->id<from||a->tatm->id>to) continue;
     if(a->flag<0) continue;
     e+=clash(a,force);
  }
  return e;
}

float Loopy::clash(Res *s,char *force) 
{ 
  Atm *a; 
  float e;
  e=0;
  for(a=s->atm;a;a=a->next)
  {
     if(a->flag<0) continue;
     e+=clash(a,force);
  } 
  return e;
}

float Loopy::clash(Atm *aa0,char *force)
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
     if(aa1->res->chn->ishet==0&&aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->res->chn->ishet==0&&aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
     mr=0;
     d=aa0->tatm->eng->radius+aa1->tatm->eng->radius;
     j=aa0->tatm->hbond*aa1->tatm->hbond;
     if(j==2||j==3||j==6||j==9){d=3.0;mr=1;}
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
	 if(aa1->res->chn->ishet==1) continue;
         if(strchr(force,'E')==0&&strchr(force,'e')==0||strchr(force,'i')==0) 
         {
           if(aa1->tatm->name[1]!='O'&&aa1->tatm->name[1]!='N') continue;
         }

         if(aa1->res==aa0->res)continue;
         if(aa1->res->chn==aa0->res->chn&&aa1->res->id-aa0->res->id==1)
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
         if(aa1->res->chn==aa0->res->chn&&aa1->res->id-aa0->res->id==-1)
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
         if(aa1->res->chn==aa0->res->chn&&(strchr(force,'h')||strchr(force,'H')))
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
        if(aa1->res->chn->ishet==1) continue;
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

float*  Loopy::scap(float **xyz_all,char *force)
{
 Res *r,*r0;
 int nseq;
 float d,*xyz;
 int all,i;
 float *ent;
 Chn *chn=(*pdb)[cid];

 all=0;
 while(xyz_all[all])all++;
 ent=new float[all];

 r0=(*chn)[start];

 nseq=0;
 for(r=r0;r;r=r->next)
 {
   if(r->id0>end) break;
   nseq+=r->tres->number*3;
 }

 xyz=new float[nseq];
 chn->transfer(xyz,r0,end,0);
 int ne=0;
 if(TRES.logg==-1) cerr<<endl;
 for(i=0;i<all;i++)
 {
   segment->chn->transfer(xyz_all[i],1);
   setend(1);
   d=scap(force);
   ent[i]=d;
   segment->chn->transfer(xyz_all[i],0);
   chn->transfer(xyz,r0,end,1);
   if(i%20==0&&TRES.logg==-1) {
        cerr<<"finished sidechain assembling:"<<i<<endl; 
        ne++;
        //if(ne%60==0) cerr<<endl;
   }
   //cerr<<i<<" "<<"side-chain predicted energy: "<<rmsd(2)<<"  "<<d<<endl;
 }
 if(TRES.logg==-1) cerr<<endl;

 chn->transfer(xyz,r0,end,1);
 
 delete [] xyz;

 return ent;
}


float Loopy::scap(char *force)
{
  Res *r,*r1,*r0;
  Atm *a;
  Rotate rot;
  float d,*xyz,*xyz0;
  int nseq;
  Chn *chn=(*pdb)[cid];

  //return 0;
  xyz=0;
  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r;r=r->next)
  {
    if(r->id0>end) break;
    nseq+=r->tres->number*3;
  }
  xyz=new float[nseq];
  xyz0=new float[nseq];
  chn->transfer(xyz0,r0,end,0);

  //assemble side-chain

  segment->chn->transfer(xyz,0);
  chn->transfer(xyz,r0,end,1);

  for(r=r0;r;r=r->next)
  {
    if(r->id0>end) break;
    if(strchr("PGA",r->tres->name)) continue;
    if(r->more) r1=r->more;
    else        r1=TRES.rotamer->next->chn->isres(r->tres->name,0);
    rot.hook(r1,r,4);
  }

  //preparing lattice
  lat->putoff();
  lat->flag=1;
  lat->ishet=1;
  lat->grdsiz=2.;
  lat->radall=15.;
  lat->ready(pdb);
  lat->putonhetatm(pdb); //modified for hetatm
  for(Chn *cchn=pdb->chn;cchn;cchn=cchn->next)
  for(r=cchn->res;r;r=r->next)
  {
    if(r->id0<start||r->id0>end||cchn!=chn)
    {
       lat->puton(r);continue;
    }
    if(strchr("PGA",r->tres->name)) 
    {
       lat->puton(r);continue;
    }
    lat->puton(r,0,4); 
  }
  chn->setflg(chn->res,10000,-1);

  for(r=chn->res;r;r=r->next) 
  for(a=r->atm;a;a=a->next)
  a->area=0;

  //create more residues of choices

  crtmore(force);

  //puton lattice
  lat->putoff();lat->ishet=0;
  for(r=r0;r;r=r->next)
  {
    if(r->id0>end) break;
    if(strchr("PGA",r->tres->name)) continue;   
    lat->puton(r,5,100);
  }

  //side-chain prediction

  d=scpred(force);

  chn->transfer(xyz,r0,end,0);
  segment->chn->transfer(xyz,1);
  chn->transfer(xyz0,r0,end,1);

  //delete more residue
  for(r=chn->res;r;r=r->next) {
	//if(r->id0>end) break;
	for(r1=r->more;r1;r1=r1->more) {
		r1->next=0;
	}
	if(r->more) delete r->more;
	r->more=0;
  }

  if(xyz){delete [] xyz;xyz=0;}
  if(xyz0){delete [] xyz0;xyz0=0;}
  return d;
}

void Loopy::crtmore(char *force)
{
  Res *r0,*r,*r1,*r2,**rm;
  Rotate rot;
  Atm *a1,*a2;
  float e,x,y;
  float temp[1000];
  Qsort cc;
  int n,order[1000],i,j;
  Chn *chn=(*pdb)[cid];

  i=0;
  for(r=TRES.rotamer->next->chn->res;r;r=r->next)
  {
    j=r->many(); 
    if(j>i) i=j;
  }
  rm=new Res*[end-start+i];

  r0=(*chn)[start];

  //take out more
  /*
  for(r=r0;r;r=r->next) {
    if(r->id0>end) break;
    if(strchr("PGA",r->tres->name)) continue;
    for(r1=r->more;r1;r1=r1->more) {
       r1->next=0;
    }
    delete r->more;r->more=0;
  }
  */
  //add more
  for(r=r0;r;r=r->next)
  {
    if(r->id0>end) break;
    if(strchr("PGA",r->tres->name)) continue;

    chn->setflg(r,r->id0,0);
    a1=(*r)[4];
    n=0;

    if(r->more) r1=r;
    else        r1=TRES.rotamer->next->chn->isres(r->tres->name,0);
    for(r2=r1->more;r2;r2=r2->more)
    {
       if(n==1000) break;
       rot.hook(r2,r,4);     
       a2=(*r2)[4];
       x=a1->dihedral(0);
       y=a2->dihedral(0);
       y=y-x;
       rot.rotate(a1,y);
       e=0;
       for(a2=a1->next;a2;a2=a2->next)
       {
         if(a2->name[1]=='H') break;
         e+=clash(a2,force);
       }
       if(r->more){r2->transfer(r);rm[n]=r2;}
       else       {rm[n]=new Res(r); rm[n]->configure(); }
       temp[n]=e;n++;
    }

    chn->setflg(r,r->id0,-1);
    cc.sort(temp,n,order);
    
    r2=0;
    for(i=0;i<n;i++)
    {     
      j=order[i];
      if(r2==0)
      {
        r->more=rm[j];r2=r->more;r2->flag=temp[i];
        r2->next=0;r2->more=0;rm[j]=0;
      }
      else
      {
        r2->more=rm[j];r2=r2->more;r2->flag=temp[i];   
        r2->next=0;r2->more=0;rm[j]=0;
      }
    }
    r2->more=0;
    r->transfer(r->more);r->flag=r->more->flag;
  }

  delete [] rm;
}


float Loopy::scpred(char *force)
{
 Res *r0,*r,*r1;
 float e,d,dp,dmin,tot;
 float xyz[100];
 int i,n;
 Chn *chn=(*pdb)[cid];

 r0=(*chn)[start];
  
 dp=10000;i=0;
 for(; ;)
 {
  tot=0;i++;
  for(r=r0;r;r=r->next)
  {
    if(r->id0>end) break;
    if(strchr("PGA",r->tres->name)) continue;
    lat->putoff(r,5,100);
    chn->setflg(r,r->id0,0);

    dmin=r->flag+clash(r,5,100,force);
    d=r->flag;
    r->transfer(xyz,0);
    n=0;
    for(r1=r->more;r1;r1=r1->more)
    {   
       n++;
       if(n>10) break; 
       r->transfer(r1);
       e=r1->flag+clash(r,5,100,force);
       if(e<dmin)
       {
         dmin=e;d=r1->flag;
         r->transfer(xyz,0);     
       }
    }
    r->flag=d;
    r->transfer(xyz,1);
    tot+=dmin;
    lat->puton(r,5,100);
    chn->setflg(r,r->id0,-1);
  }
  if(fabs(tot-dp)<0.01*fabs(tot)&&i>1) break;
  //cerr<<i<<"  side-chain minimization: "<<tot<<endl;
  dp=tot;
  if(i>20) break;
 }
 return tot;
}

float Loopy::minimize(float *xyz0,int flg)
{
  float **xyz_all;
  float *ent,d;
  xyz_all=new float*[2];
  xyz_all[0]=xyz0;
  xyz_all[1]=0;
  ent=minimize(xyz_all,1);
  d=ent[0];
  delete [] ent;
  delete [] xyz_all;  
  return d;
}

float* Loopy::minimize(float **xyz_all,int flg)
{
  Res *r0,*r;
  float *xyz;
  int i,nseq,all;
  float *ent,**pot;
  int   *order;
  Qsort cc;
  Chn *chn=(*pdb)[cid];


  xyz=0;
  nseq=0; 

  r0=(*chn)[start];

  all=0;
  while(xyz_all[all])all++;

  ent=new float[all];
  order=new int[all];
  pot=new float*[all];

  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;     
  
  xyz=new float[nseq];  

  chn->transfer(xyz,r0,end,0);

  //rusage tt;
  //getrusage(RUSAGE_SELF,&tt);
  //int ne=0;
  if(TRES.logg==-2) cerr<<endl;
  for(i=0;i<all;i++)
  {      
    segment->chn->transfer(xyz_all[i],1);                          
    ent[i]=segmin(flg); 
    segment->chn->transfer(xyz_all[i],0);          
    if(i%20==0&&TRES.logg==-2) {
	cerr<<"finished minimizing candidates:"<<i<<endl;
        //cerr<<".";
        //ne++;
        //if(ne%60==0) cerr<<endl;
    }
  }
  if(TRES.logg==-2) cerr<<endl; 
  //rusage tt0;
  //getrusage(RUSAGE_SELF,&tt0);

  //cerr<<"\nminimize loos:"<<all<<" with time..."<<tt0.ru_stime.tv_usec-tt.ru_stime.tv_usec<<endl;
  cc.sort(ent,all,order);
  for(i=0;i<all;i++)
  {
    pot[i]=xyz_all[order[i]];
  }

  for(i=0;i<all;i++) xyz_all[i]=pot[i];

  delete [] pot;delete [] order;

  chn->transfer(xyz,r0,end,1);
  delete [] xyz; 

  return ent;
}




float** Loopy::zipper(int seed,float dis,int flg)
{
  //generate loop conformations
  Res  *r1,*r2;
  Res  *r0,*r;  
  Atm  *a,*a1;
  int ne,i,n,nseq;
  float d;
  float *xyz,*xyz0,**xyz_all;
  Chn *chn=(*pdb)[cid];

  xyz=0;xyz_all=0;

  //coordinate buffer

  nseq=0; 
  r0=(*chn)[start];
  for(r=r0;r&&r->id0<=end;r=r->next) 
  nseq+=r->tres->number*3;     
  
  xyz=new float[nseq];  
  xyz0=new float[nseq];
  xyz_all=new float*[arbt+10];
  for(i=0;i<arbt+10;i++) xyz_all[i]=0;

  chn->transfer(xyz,r0,end,0);
  
  r1=(*chn)[end];r2=(*segment->chn)[end]; 
  //find possible hits

  i=0;n=0;
  /*
  Chn datb;
  if(databaseonly==1) {
	char *s=segment->chn->getseqn();
	datb.create(s);
  }
  */
  ne=0;
  for(; ;)
  {
    n++;ne++;
    randmz(ne+seed);
    if(i==0&&ne>40) {
	cerr<<"the segment could not be connected"<<endl;
	exit(0);
    }
    else if(i&&ne/i>=20) {
	cerr<<"the segment could not be connected"<<endl;
	exit(0);
    }
    //segment->chn->transfer(xyz0,0);
    /*
    if(checkomega()==0) {
	n--;continue; 
    }
    */
    //chn->transfer(xyz0,r0,end,1);
    //chn->dihedral(stdout);
    if(chiangle) {
	if(TRES.logg>3)cerr<<"rotate the dihedral to the specified in the database..."<<n-1<<"  "<<databaseonly<<endl;
	//checkomega(stdout);
	if(randmzChiangle(n-1)==0&&databaseonly==1) break;
	//checkomega(stdout);
        //FILE *fp=0;
	if(TRES.logg>3) segment->chn->dihedral(stderr);
    }
    else if(secd=='h'&&end-start+1>3) randmzHelix(seed+n*n,15);
    else if(secd=='e'&&end-start+1>3) randmzSheet(seed+n*n,30);
    if(TRES.logg>3) segment->chn->dihedral(stderr);
    //chn->dihedral(stdout);
    //segment->chn->transfer(xyz0,0);
    //chn->transfer(xyz0,r0,end,1);
    //segment->chn->write("se1");
    //chn->dihedral("s1");
    //chn->write("te1");
    chn->transfer(xyz,r0,end,1);
    //chn->dihedral("s2");
    //chn->write("te2");
    d=TRES.distance(r1->atm->xyz,r2->atm->xyz)+
      TRES.distance(r1->atm->next->xyz,r2->atm->next->xyz);
    /*if(chiangle) {
	cerr<<"database check stem distance...."<<n<<" "<<d<<endl; 
    }*/	
    if(d>dis) continue;
    if(flg==1) {d=closeg();  if(d==-2) continue;  }
    if(flg==2) {d=segmin(1); if(d>1.e20) continue;}
	
    if(direct==1&&r1->next&&end-start+1<=3) {
   	a=(*r1)[2];a1=(*r2)[2];
        d=0;
	if(a&&a1) d+=TRES.distance(a,a1);
	a=(*r1)[3];a1=(*r2)[3];
 	if(a&&a1) d+=TRES.distance(a,a1);
	if(d>1.5) continue;	
    }

    chn->transfer(xyz,r0,end,1);
    if(flg&&outlet<0.5)    {if(setend(1)==0) continue;}
    xyz_all[i]=new float[nseq];
    segment->chn->transfer(xyz_all[i],0);
    if(TRES.logg>3) {
    	pdb->chn->transfer(xyz_all[i],r0,end,1);
    	char sf[100];
    	sprintf(sf,"_%i",i);
    	pdb->write(sf);
    	pdb->chn->transfer(xyz,r0,end,1);
    }
    i++;
    if(i>=arbt&&databaseonly==0) break;
  }
  cerr<<"number of segments generated :"<<n<<endl;
  cerr<<"number of segment successfully connected: "<<i<<endl;
  chn->transfer(xyz,r0,end,1);
  if(xyz) {delete [] xyz;xyz=0;}
  if(xyz0) {delete [] xyz0;xyz=0;}
  return xyz_all;
}

float** Loopy::cross(float **xyz_all,int cot,float cut)
{
  //loop cross
  Chn *cht;
  Res *r,*r0;  
  int all,i,nseq,j,m,m3,m4;
  float e,d,*xyz,**xyz_out;     
  Loopy *seg1,*seg2;    
  Qsort cc;
  int *order,*used;
  float *dist;
  Chn *chn=(*pdb)[cid];

  if(end-start+1<=3) return xyz_all;
 //space
   
  seg1=new Loopy(this); seg2=new Loopy(this);
  seg1->start=start; seg1->end=(start+end)/2;
  seg2->start=(start+end)/2; seg2->end=end;
  seg1->pdb=pdb;seg2->pdb=pdb;
  seg1->create();seg2->create();
  seg1->direct=1;seg2->direct=0;
 
  all=0;
  while(xyz_all[all])all++;
  order=new int[all];
  dist=new float[all];

  r0=(*chn)[start];
  nseq=0;m3=0;
  for(r=r0;r;r=r->next)
  {
    if(r->id0<=end) nseq+=r->tres->number*3;
    if(r->id0<seg2->start) m3+=r->tres->number*3;
  }
  
  m4=m3+seg2->segment->chn->res->tres->number*3;  
  if(cot>all) cot=all;
  used=new int[all*all+all];
  for(i=0;i<all*all+all;i++) used[i]=0;
  xyz_out=new float*[all*(cot+2)];
  for(i=0;i<all*(cot+2);i++) xyz_out[i]=0;

  xyz=new float[nseq];
  chn->transfer(xyz,r0,end,0);

  //create random conformation
  for(i=0;i<all;i++)
  { 
     xyz_out[i]=new float[nseq];
     TRES.copy(xyz_all[i],xyz_out[i],nseq);
  }
  
  int cot0,kit;
  d=0.5*(end-start); m=all;
  for(i=0;i<all;i++)
  {
    for(j=0;j<all;j++)dist[j]=TRES.distance(xyz_all[i]+m3+3,xyz_all[j]+m3+3);
    cc.sort(dist,all,order);
    cot0=0;
    for(kit=0;kit<cot;kit++)
    {  
     j=order[kit];
     if(i==j) continue;
     if(used[i*all+j]==1) continue;
     if(cot0>cot) break;
     e=dist[kit];
     if(e>d||e<cut) continue;
     cht=seg1->segment->chn; r=cht->res;
     cht->transfer(xyz_all[i],r,seg1->end,1);
     cht=seg2->segment->chn; r=cht->res;
     cht->transfer(xyz_all[j]+m3,r,seg2->end,1);                        
     e=closeg(seg1,seg2);
     if(e==-2) continue; 
     xyz_out[m]=new float[nseq];
     cht=seg1->segment->chn; r=cht->res;
     seg1->segment->chn->transfer(xyz_out[m],r,seg1->end,0);
     cht=seg2->segment->chn; r=cht->res->next;
     cht->transfer(xyz_out[m]+m4,r,seg2->end,0);
     m++;cot0++;
     used[i*all+j]=1;used[j*all+i]=1;
     if(m>all+10*all) break;
    }
  }
  chn->transfer(xyz,r0,end,1);
  
  for(i=0;i<all;i++) 
  {
     if(xyz_all[i]) {delete [] xyz_all[i];xyz_all[i]=0;}
  }
  
  delete [] xyz_all;
  seg1->pdb=0;seg2->pdb=0;
  //seg1->rotamer=0;seg2->rotamer=0;
  //seg1->smooth=0; seg2->smooth=0;
  seg1->next=0;seg2->next=0;seg1->lat=0;seg2->lat=0;
  delete seg1;seg1=0;
  delete seg2;seg2=0;
  delete [] order;
  delete [] dist;
  delete [] used;
  delete [] xyz;
  return xyz_out;
}

float Loopy::ensort(float *xyz,char *flg)
{
//v: ver der waals energy 
//V: ver der waals energy only CA
//e: electrostatic energy including solvent
//E: electrostatic energy no solvent
//h: hydrogen bond including solvent hydrogen bond
//H: hydrogen bond not including solvent contribution

  float d;
  float **xyz_all;
  float *ent;
  xyz_all=new float*[2];
  xyz_all[0]=xyz;
  xyz_all[1]=0;
  ent=ensort(xyz_all,flg);
  d=ent[0];
  delete [] ent;
  delete [] xyz_all;
  return d;
}

float* Loopy::ensort(float **xyz_all,char *flg)
{
  Res *r,*r0;
  Atm *a,*a1;
  int i,j,nseq,all,m;
  float e,*xyz,area,area1,area0,d,f;
  float **pot,*ent,*entv,*ente,*ents;
  int  *order;
  Qsort cc;
  int *give,*gives;
  float coo[3];
  Chn *chn=(*pdb)[cid];

  xyz=0;pot=0;order=0;ent=0;
  //space
     
  all=0;
  while(xyz_all[all])all++;
  i=0;
  for(r=chn->res;r;r=r->next)
  for(a=r->atm;a;a=a->next)i++;
  give=new int[i];
  gives=new int[i];

  order =new int[all];
 
  //preparing lattice
  lat->putoff();
  lat->flag=1;
  lat->ishet=1;
  lat->grdsiz=2.;
  lat->radall=15.;
  lat->ready(pdb);
  for(Chn *cchn=pdb->chn;cchn;cchn=cchn->next)
  for(r=cchn->res;r;r=r->next)
  {
    if(r->id0>=start&&r->id0<=end&&cchn==chn) continue;
    lat->puton(r);
  }
  lat->putonhetatm(pdb); //modified to account for hetatm
  //double computing

  r0=(*chn)[start];
  chn->setflg(chn->res,10000,-1);
  chn->setflg(r0,end,1); 
  r0->atm->flag=-1;r0->atm->next->flag=-1;
  a=(*r0)[" HN "]; if(a)a->flag=-1;
  r=(*chn)[end]; a=r->atm->next;
  a->flag=-1;a->next->flag=-1;a->next->next->flag=-1;
  chn->setflg(give,chn->res,100000,0);

  //coordinate buffer

  nseq=0; 
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;         
  xyz=new float[nseq];

  pot=new float*[all];
  for(i=0;i<all;i++) pot[i]=0; 
  
  ent=new float[all];
  entv=new float[all];
  ente=new float[all];
  ents=new float[all];
  for(i=0;i<all;i++) {ent[i]=0;entv[i]=0;ente[i]=0;ents[i]=0;}

  chn->transfer(xyz,r0,end,0);
     
  //  energy
  
  if(strchr(flg,'v')||strchr(flg,'V'))
  {
    char flgv[100];
    int  iiv;
    iiv=0;
    if(strchr(flg,'v'))flgv[iiv++]='v';
    if(strchr(flg,'V'))flgv[iiv++]='V';   
    flgv[iiv++]='\0';

    chn->setflg(give,chn->res,100000,1);
    for(i=0;i<all;i++)
    {     
     chn->transfer(xyz_all[i],r0,end,1);  
     lat->putoff(r0,end);
     lat->puton(r0,end);     
     for(r=r0;r;r=r->next)
     {
       if(r->id0>end) break;
       entv[i]+=clash(r,flgv);
     }      
    }
  }

  if(strchr(flg,'H')||strchr(flg,'E'))
  {
    char flgv[100];
    int  iiv;
    iiv=0;
    if(strchr(flg,'H'))flgv[iiv++]='H';
    if(strchr(flg,'E'))flgv[iiv++]='E';
    flgv[iiv++]='\0';

    chn->setflg(give,chn->res,100000,1);
    for(i=0;i<all;i++)
    {    
     chn->transfer(xyz_all[i],r0,end,1);
     lat->putoff(r0,end);
     lat->puton(r0,end);
     for(r=r0;r;r=r->next)
     {
       if(r->id0>end) break;
       ente[i]+=clash(r,flgv);
     }     
    }
  }



  //calculating farout

  if(strchr(flg,'d'))
  {
    for(i=0;i<all;i++)
    {
      d=0;
      chn->transfer(xyz_all[i],r0,end,1);
      for(j=0;j<all;j++)
      {
        if(j==i) continue;
        segment->chn->transfer(xyz_all[j],1);
        e=rmsd(2);
        if(e<0.01) continue;
        d+=exp(-e*e*e/(end-start)/6);
      }
      ent[i]+=-8.31*0.3*log(d/all);
    }   
  }

  if(strchr(flg,'D'))
  {
    int k1;
    k1=0;
    while(near[k1])k1++;

    for(i=0;i<all;i++)
    {
      d=0;
      chn->transfer(xyz_all[i],r0,end,1);
      for(j=0;j<k1;j++)
      {
        segment->chn->transfer(near[j],1);
        e=rmsd(2);
        d+=exp(-e*e*e/(end-start)/6);
      }
      ent[i]+=-8.31*0.3*log(d/all);
    }
  }

  //calculate surface and electrostatics
  if(strchr(flg,'s')||strchr(flg,'e')||strchr(flg,'h')||strchr(flg,'i'))
  {
    char flgs[100];
    int  iis;
    iis=0;
    if(strchr(flg,'e'))flgs[iis++]='e';
    if(strchr(flg,'h'))flgs[iis++]='h';
    if(strchr(flg,'i'))flgs[iis++]='i';
    flgs[iis++]='\0';

    
    setoff(xyz_all,6.);
    chn->setflg(gives,chn->res,1000000,0);
    for(i=0;i<all;i++)
    {
      chn->setflg(gives,chn->res,1000000,1);
      chn->transfer(xyz_all[i],r0,end,1);
      chn->surface(10,1.4);

      area=0;area0=0;area1=0;
      for(r=chn->res;r;r=r->next)
      for(a=r->atm;a;a=a->next)
      {
        if(a->tatm->ispolar==3)      area0+=a->area;
        else if(a->tatm->ispolar==2) area+=a->area;
        else                         area1+=a->area;
      }
      if(strstr(flg,"s1"))  d=0.025*(area0/4+area/2+area1);
      else                  d=0.05*(area0/4+area/2+area1);
      if(strchr(flg,'s')) ents[i]+=d;

      if(strchr(flg,'e'))
      {
        f=(80-40)/80./12.56;
        for(r=chn->res;r;r=r->next)
        {
          for(a=r->atm;a;a=a->next)
          {
            if(a->area==0||a->nemp==0) continue;
            a->temp[6]=0;
            lat->getcell(a,cutoff);
            lat->resonly();
            for(j=0;j<lat->nget;j++) 
            if(lat->obtain[j]->atm->res==r) {j=-1;break;}
            if(j!=-1)
            {
              lat->obtain[lat->nget]=lat->atom[a->id0-lat->offset];
              lat->nget++;
            }
            for(m=0;m<lat->nget;m++)
            for(a1=lat->obtain[m]->atm->res->atm;a1;a1=a1->next)
            {
              if(a1->tatm->eng->charge==0) continue;
	      if(a1->res->chn->ishet==1) continue; //modifed to accout for hetatm
              for(j=0;j<3;j++) coo[j]=a->temp[j]-a1->xyz[j];
              d=TRES.distance(coo);
              d=d*d*d;
              e=TRES.vectord(coo,a->temp+3);
              e=-e/d;
              a->temp[6]+=f*e*a1->tatm->eng->charge*a->area; 
            } 
          }
        }
      }  
      if(strchr(flg,'e')||strchr(flg,'h')||strchr(flg,'i'))
      {
         chn->setflg(give,chn->res,1000000,1);
         chn->transfer(xyz_all[i],r0,end,1);  
         lat->putoff(r0,end);
         lat->puton(r0,end);     
         for(r=r0;r;r=r->next)
         {
           if(r->id0>end) break;
           ente[i]+=clash(r,flgs);
         }       
      }
      chn->transfer(xyz,r0,end,1);
      segment->chn->transfer(xyz_all[i],1); 
      //cerr<<i<<" surface: "<<rmsd(2)<<"  "<<area+area0<<"  "<<area<<"  "<<area0<<endl;
    }
    
    d=0;
    for(i=0;i<all;i++)d+=ents[i];
    d=d/all;
    for(i=0;i<all;i++)ents[i]=ents[i]-d;
  }

  if(strchr(flg,'n'))
  {
    float dv,de,ds;
    dv=entv[0];
    for(i=0;i<all;i++) if(entv[i]<dv)dv=entv[i];
    ds=ents[0];
    for(i=0;i<all;i++) if(ents[i]<ds)ds=ents[i];
    ds=dv-ds;
    for(i=0;i<all;i++) ents[i]+=ds;
    de=ente[0];
    for(i=0;i<all;i++) if(ente[i]<de)de=ente[i];
    de=dv/3-de;
    for(i=0;i<all;i++) ente[i]+=de; 
  }

  for(i=0;i<all;i++) ent[i]=entv[i]+ents[i]+ente[i];
  //
  float engmax=ent[0];
  for(i=0;i<all;i++) { 
	if(engmax<ent[i])engmax= ent[i];
  }
  //
  if(strchr(flg,'w'))
  {
    float d0;
    for(i=0;i<all;i++) ent[i]=0;
    float et=8.31*0.3;
    for(i=0;i<all;i++)
    {
       
      d=0;d0=0;
      chn->transfer(xyz_all[i],r0,end,1);
      float ex=entv[i]+ente[i]+ents[i];

      for(j=0;j<all;j++)
      {
	float ey=entv[j]+ente[j]+ents[j]; 
	if(ey-ex>20) continue; //modified to speed up calculation
        segment->chn->transfer(xyz_all[j],1);
        e=rmsd(2);
        d=exp(-e*e*e/(end-start)/6);
        d0+=d;
        ent[i]+=-d*exp(-(entv[j]+ente[j]+ents[j])/(et));
      }
      ent[i]=-et*log(-ent[i]);
    }  
  }

  if(strchr(flg,'W'))
  {
    float d0;
    for(i=0;i<all;i++) ent[i]=0;
    float et=8.31*0.3;
    for(i=0;i<all;i++)
    {
      d=0;d0=0;
      chn->transfer(xyz_all[i],r0,end,1);
      for(j=0;j<all;j++)
      {
        segment->chn->transfer(xyz_all[j],1);
        e=rmsd(2);
        d=exp(-e*e*e/(end-start)/6);
        d0+=d;
        ent[i]+=-d*(exp(-entv[j]/(et))+
                    exp(-ents[j]/(et))+
                    exp(-ente[j]/(et)));
      }
      ent[i]=-et*log(-ent[i]);
    } 
  }

  //first part of segment;

  for(i=0;i<all;i++)order[i]=i;

  if(!strstr(flg,"?"))cc.sort(ent,all,order);

  for(i=0;i<all;i++)
  {
     j=order[i];  
     pot[i]=xyz_all[j];
  }   
 
  for(i=0;i<all;i++)
  {  
    xyz_all[i]=pot[i];
    pot[i]=0;
  }

  chn->transfer(xyz,r0,end,1);
  for(r=chn->res;r;r=r->next) 
  for(a=r->atm;a;a=a->next)
  a->area=0;
  
  //delete space!  
  delete [] order; delete [] xyz;   delete [] pot;      
  delete [] entv;  delete [] ente;  delete [] ents;
  delete [] give;  delete [] gives; 
  lat->putoff();lat->ishet=0;
  return ent;
}

float Loopy::noclose(float **xyz_all,float *xyz0) 
{
  Res *r0,*r;
  int nseq,j,n,all;
  float d,*xyz,close;
  Chn *chn=(*pdb)[cid];

  xyz=0;

  all=0;
  while(xyz_all[all])all++;

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
  xyz=new float[nseq];
  chn->transfer(xyz,r0,end,0);
  chn->transfer(xyz0,r0,end,1);

  close=1000;n=0;
  for(j=0;j<all;j++)
  {
    segment->chn->transfer(xyz_all[j],1);
    d=rmsd(2);
    if(d<close)
    {
      close=d;n=j;
    }
  }

  close=n*10000+close;
  chn->transfer(xyz,r0,end,1);
  delete [] xyz;
  return close;
}


int Loopy::noclose(float **xyz_all,float close)
{
  Res *r0,*r;
  int nseq,i,j,m,all;
  float *pot,d,*xyz; 
  int *kepp;
  Chn *chn=(*pdb)[cid];

  pot=0;xyz=0;

  all=0;
  while(xyz_all[all])all++;
  if(all<part+2) return all;
  kepp=new int[all];
  for(i=0;i<all;i++) kepp[i]=1;

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r;r=r->next)
  {
     if(r->id0>end) break;
     nseq+=r->tres->number*3;
  }
  xyz=new float[nseq];
  chn->transfer(xyz,r0,end,0);
  segment->chn->transfer(xyz,1);

  for(i=0;i<all;i++)
  {
    if(kepp[i]==0) continue;
    chn->transfer(xyz_all[i],r0,end,1);

    for(j=i+1;j<all;j++)
    {
      if(xyz_all[j]==0) continue;
      
      segment->chn->transfer(xyz_all[j],1);       
      d=rmsd(2);
      if(d<close) kepp[j]=0;
    }
  }

  m=0;
  for(i=all-1;i>=0;i--)
  {
    if(all-m<part) break;
    if(kepp[i]==0)
    {
      delete [] xyz_all[i];
      xyz_all[i]=0;
      m++;
    }
  }
  
  m=0;
  for(i=0;i<all;i++)
  {
    if(xyz_all[i]==0) continue;
    pot=xyz_all[m];
    xyz_all[m]=xyz_all[i]; 
    xyz_all[i]=pot;
    m++;
  }

  chn->transfer(xyz,r0,end,1);
  segment->chn->transfer(xyz,1);
  if(xyz) {delete [] xyz;xyz=0;}
  if(kepp) {delete [] kepp;kepp=0;}
  return m;
}

float **Loopy::around(float **xyz_all)
{
  float *ent;
  float **xyz_all0,**xyz_all00;
  int kep0,i,j,k,kep;
  Qsort cc;
  float **pot;
  int *order;

  cerr<<"start expanding around..."<<endl;

  kep=0;while(xyz_all[kep])kep++;
  xyz_all0=new float*[kep*12];
  for(i=0;i<kep*12;i++) xyz_all0[i]=0;
  j=0;
  srandom(10000);
  for(i=0;i<kep;i++)
  {
    xyz_all0[j++]=xyz_all[i];
    k=random()%100000;
    xyz_all00=expand(xyz_all[i],k,10,30);
    kep0=0;while(xyz_all00[kep0])kep0++;
    for(k=0;k<kep0;k++) xyz_all0[j++]=xyz_all00[k]; 
    delete [] xyz_all00;
  }
  delete [] xyz_all;
  xyz_all=xyz_all0;
  ent=minimize(xyz_all,1);
  kep=0;while(xyz_all[kep])kep++;
  order=new int[kep];pot=new float*[kep];
  cc.sort(ent,kep,order);
  for(i=0;i<kep;i++) pot[i]=xyz_all[order[i]];
  for(i=0;i<kep;i++) xyz_all[i]=pot[i]; 
  float e,d;
  d=0;
  for(i=0;i<kep;i++)
  {
    segment->chn->transfer(xyz_all[i],1);
    e=rmsd(2);
    d+=e/kep;
    if(TRES.logg>3)cerr<<i<<" after expanding energy of v: "<<rmsd(2)<<" "<<ent[i]<<endl;
  }
  cerr<<" the total number after refill: "<<kep<<" "<<d<<endl;
  delete [] order;delete [] pot; 
  return xyz_all;
}

float *Loopy::expand(float **xyz_all,int flg)
{
  Res *r0,*r;
  int nseq;
  int i,kep;
  float *ent,e,*xyz;
  int kk;
  Qsort cc;
  int *order;
  float **pot;
  char flg0[10];

  Chn *chn=(*pdb)[cid];

  strcpy(flg0,"v"); 

  r0=(*chn)[start];nseq=0;
  for(r=r0;r;r=r->next)if(r->id0<=end) nseq+=r->tres->number*3;
  xyz=new float[nseq];
  chn->transfer(xyz,r0,end,0);

  kep=0;
  while(xyz_all[kep])kep++;
  ent=ensort(xyz_all,flg0);

  kk=part0; 

  float rms,rms0;
  rms=0;
  int no=0;
  if(TRES.logg==-1) cerr<<endl;
  for(i=0;i<kep;i++)
  {
    segment->chn->transfer(xyz_all[i],1);
    rms0=rmsd(2);
    rms+=rms0/kep;
    if(i<kk)
    {
      if(TRES.logg>3)cerr<<i<<" energy of "<<flg0<<": "<<rms0<<" "<<ent[i]<<" ";
      ent[i]=expand(xyz_all[i],flg);
      segment->chn->transfer(xyz_all[i],1);
      if(TRES.logg>3)cerr<<rmsd(2)<<" "<<ent[i]<<endl;
    }
    else    
    {
      if(TRES.logg>3)cerr<<i<<" energy of "<<flg0<<": "<<rms0<<" "<<ent[i]<<endl;
      delete [] xyz_all[i];xyz_all[i]=0;
    }
    if(i%20==0&&TRES.logg==-1) {
	cerr<<"finished minimizing candidates:"<<i<<endl;
	no++;
	//if(no%60==0) cerr<<no<<endl;
    }
  }
  if(TRES.logg==-1) cerr<<endl;
  chn->transfer(xyz,r0,end,1);
  kep=0;
  while(xyz_all[kep])kep++;
  
  if(TRES.logg>3)cerr<<"the total number after "<<flg0<<": "<<kep<<"  "<<rms<<endl;
  
  rms=0;
  pot=new float*[kep];
  order=new int[kep];
  cc.sort(ent,kep,order);
  for(i=0;i<kep;i++) pot[i]=xyz_all[order[i]];
  for(i=0;i<kep;i++) xyz_all[i]=pot[i];
  for(i=0;i<kep;i++)
  {
     if(ent[i]>1.e20) {delete [] xyz_all[i]; xyz_all[i]=0;}
  }
  kep=0;while(xyz_all[kep])kep++;
  for(i=0;i<kep;i++)
  {
    segment->chn->transfer(xyz_all[i],1);
    e=rmsd(2);rms+=e/kep;
    if(TRES.logg>3)cerr<<i<<" energy of "<<flg0<<": "<<e<<" "<<ent[i]<<endl;
  }
  if(TRES.logg>3)cerr<<"the total number after "<<flg0<<": "<<kep<<"  "<<rms<<endl;
  delete [] pot;delete [] order;delete [] xyz;
  return ent;
}

float Loopy::expand(float *xyz0,int flg)
{
  Res *r0,*r;
  float d,*ent,**xyz_all;
  int kep,i,n;
  int nseq;

  Chn *chn=(*pdb)[cid];

  nseq=0;
  r0=(*chn)[start];
  for(r=r0;r;r=r->next)if(r->id0<=end) nseq+=r->tres->number*3;
 
  srandom(100);
  d=minimize(xyz0,1);
  if(d<0||flg==0||fapr==1) return d;

  for(i=0;i<1;i++)
  {
    n=random()%100000;
    xyz_all=expand(xyz0,n,5,10);
    ent=ensort(xyz_all,"v");delete [] ent;
    kep=0;
    while(xyz_all[kep])kep++;
    for(n=1;n<kep;n++) {delete [] xyz_all[n];xyz_all[n]=0;}
    ent=minimize(xyz_all,1);
    if(ent[0]<d) 
    {
      d=ent[0];
      TRES.copy(xyz_all[0],xyz0,nseq);
    } 
    delete [] xyz_all[0];
    delete [] xyz_all;
    delete [] ent;
    if(d<0) break;
  }
  return d;
}

float** Loopy::expand(float *xyz0,int seed,int kep,int ang)
{
//using ranmdz to sample nearby conformations

  Res *r0,*r;
  float d,**xyz_all,*xyz;
  int nseq,i,n;

  Chn *chn=(*pdb)[cid];

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r;r=r->next)
  {
    if(r->id0>end)break;
    nseq+=r->tres->number*3;
  }
  xyz=new float[nseq];
  xyz_all=new float*[kep+10];
  for(i=0;i<kep+10;i++) xyz_all[i]=0;
  //xyz_all[0]=new float[nseq];
  //TRES.copy(xyz0,xyz_all[0],nseq);
  chn->transfer(xyz,r0,end,0); 
  i=-1;n=0;
  for(; ;)
  {
    i++;
    if(n>=kep) break; 
    segment->chn->transfer(xyz0,1);
    randmz(i+seed,ang);
    d=closeg();
    if(d==-2) continue;
    chn->transfer(xyz,r0,end,1);
    if(setend(1)==0) continue;
    xyz_all[n]=new float[nseq];
    segment->chn->transfer(xyz_all[n],0);
    n++;
  }
  chn->transfer(xyz,r0,end,1);
  if(xyz) {delete [] xyz;xyz=0;}
  return xyz_all;
}

float *Loopy::peel(float **xyz_all,char *flg,int got)
{
//recursively peel off
  float *ent;
  int i,kep;
 
  for(i=0;i<10;i++)
  {
    ent=trim(xyz_all,flg,0); 
    kep=0;
    while(xyz_all[kep])kep++;
    if(kep<got+2) break;
  }
  return ent;
}

float **Loopy::refill(float **xyz_all,char *flg,int need)
{
  Res *r,*r0;
  float *ent;
  float **xyz_all0;
  int kep0,i,j,k,kep,nseq;
  float d,e;
  Chn *chn=(*pdb)[cid];
 
  r0=(*chn)[start];

  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
  
  if(TRES.logg>3)cerr<<"start refilling..."<<endl;
  if(!strstr(flg,"?")) {ent=trim(near,flg,3*part);delete [] ent;}
  kep=0;while(xyz_all[kep])kep++;
  kep0=0;while(near[kep0])kep0++;
  xyz_all0=new float*[kep+kep0+10];
  for(i=0;i<kep+kep0+10;i++) xyz_all0[i]=0;
  for(i=0;i<kep;i++) xyz_all0[i]=xyz_all[i];
  k=0;
  for(i=kep;i<need;i++)
  {
    for(j=k;j<kep0;j++)
    {
      d=noclose(xyz_all0,near[j]);
      if(d-int(d)/1000*1000>0.1*(end-start))
      {
        xyz_all0[i]=new float[nseq];
        TRES.copy(near[j],xyz_all0[i],nseq);
        k=j;
        break;
      }
    }
  }
  delete [] xyz_all;
  xyz_all=xyz_all0;
  kep=0;while(xyz_all[kep])kep++;
  d=0;
  for(i=0;i<kep;i++)
  {
    segment->chn->transfer(xyz_all[i],1);
    e=rmsd(2);
    d+=e/kep;
    if(TRES.logg>3)cerr<<i<<" after refill: "<<rmsd(2)<<endl;
  }
  if(TRES.logg>3)cerr<<" the total number after refill: "<<kep<<"  "<<d<<endl;
  return xyz_all;
}

float **Loopy::fill(float **xyz_all,float cut)
{
  Res *r0,*r;
  int i,kep,kep0,nseq;
  float **xyz_all0,*xyz;
  int kk;
  Chn *chn=(*pdb)[cid];

  kep=0;
  while(xyz_all[kep])kep++;

  kep0=0;
  while(near[kep0])kep0++;

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r;r=r->next)if(r->id0<=end)nseq+=r->tres->number*3;
  xyz=new float[nseq];
  chn->transfer(xyz,r0,end,0); 
  
  xyz_all0=new float*[kep+kep0+100];
  for(i=0;i<kep+kep0+100;i++) xyz_all0[i]=0; 

  for(i=0;i<kep;i++) 
  {
    xyz_all0[i]=new float[nseq];
    TRES.copy(xyz_all[i],xyz_all0[i],nseq);
  }

  for(i=0;i<kep0;i++)xyz_all0[i+kep]=near[i];
  delete [] near;
  near=xyz_all0;
  if(id==0) next->near=near;
  else      last->near=near;
  kk=0;
  while(near[kk])kk++;

  kep0=noclose(near,cut);
  kep0=0; while(near[kep0])kep0++;
  for(i=5*part;i<kep0;i++){delete [] near[i];near[i]=0;}
  kep0=0; while(near[kep0])kep0++;
  float e,d;
  d=0;
  for(i=0;i<kep0;i++)
  {
    segment->chn->transfer(near[i],1);
    e=rmsd(2);
    d+=e/kep0;
    //cerr<<"  "<<i<<" fill: "<<rmsd(2)<<endl; 
  }
  if(TRES.logg>3)cerr<<" the total number after fill: "<<kk<<" "<<cut<<"  "<<kep0<<"  "<<d<<endl;
  chn->transfer(xyz,r0,end,1);
  delete [] xyz;
  return xyz_all;
}

float **Loopy::copy(float **xyz_all)
{
  int kep;
  kep=0;
  while(xyz_all[kep])kep++;
  return copy(xyz_all,kep);
}

float **Loopy::copy(float **xyz_all,int all)
{
  Res *r0,*r;
  float **xyz_temp;
  int  i,kep,nseq;
  Chn *chn=(*pdb)[cid];

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next)  nseq+=r->tres->number*3;
  kep=0;
  while(xyz_all[kep])kep++;
  xyz_temp=new float*[kep+10];
  for(i=0;i<kep+10;i++) xyz_temp[i]=0;
  for(i=0;i<all;i++)
  {
    if(i>=kep) break;
    xyz_temp[i]=new float[nseq];
    TRES.copy(xyz_all[i],xyz_temp[i],nseq); 
  }
  return xyz_temp; 
}

float* Loopy::trim(float **xyz_all,float*ent,int haa)
{
  int  i,kep,kk;
  float e,d;
  float **pot,*entt;
  int *kept,*order;
  Qsort cc;
  
  kep=0;
  while(xyz_all[kep])kep++;
  
  kept=new int[kep+10];pot=new float*[kep+10];entt=new float[kep+10];
  order=new int[kep];
  for(i=0;i<kep+10;i++){kept[i]=0;pot[i]=0;}
  for(i=0;i<kep;i++) entt[i]=ent[i];
  cc.sort(entt,kep,order); 

  //sort based on energy
  e=entt[0];
  d=0;
  for(i=0;i<kep;i++) d+=(ent[i]-e)*(ent[i]-e);
  d=sqrt(d/kep)+e;
  e=entt[kep*3/5];
  if(d<e)d=e;
  if(haa<=0)
  {
    for(i=0;i<-haa;i++) kept[i]=1;
  }
  //start
  for(i=0;i<kep;i++)
  {
     kk=order[i];
     if(haa>0)
     {
       if(i>haa) {delete [] xyz_all[kk];xyz_all[kk]=0;}
     }
     else
     {
       if(entt[i]>d&&i>=part&&kept[kk]==0) {delete [] xyz_all[kk];xyz_all[kk]=0;}
     }
  }
  //end
 
  kk=0;
  for(i=0;i<kep;i++)
  {
    if(xyz_all[i]){pot[kk]=xyz_all[i];entt[kk]=ent[i];kk++;}
  }
  for(i=0;i<kep;i++)
  {
    if(i<kk) {xyz_all[i]=pot[i];ent[i]=entt[i];}
    else     xyz_all[i]=0;
  }
  kep=0;
  while(xyz_all[kep])kep++;
  cc.sort(ent,kep,kept);
  for(i=0;i<kep;i++) pot[i]=xyz_all[kept[i]];
  for(i=0;i<kep;i++) xyz_all[i]=pot[i]; 

  //re200:
  xyz_all[kep]=0;
  delete [] kept;delete [] pot;delete [] entt;
  delete [] order;
  return ent;
}

float* Loopy::trim(float **xyz_all,float *ent,float cutf)
{
 int i,kk,kep;
 kep=0;
 while(xyz_all[kep])kep++;
 kk=0;
 for(i=0;i<kep;i++)if(ent[i]<cutf)kk++;
 ent=trim(xyz_all,ent,kk);
 return ent;
}

float* Loopy::trim(float **xyz_all,char *flg,int haa)
{
 int kep,i;
 float e,rms,rms0,*ent,*entt;
 char flg0[100];
 Qsort cc;
 int *order;

 kep=0;
 while(xyz_all[kep])kep++;
 strcpy(flg0,flg); 
 i=strlen(flg);
 flg0[i]='?';flg0[i+1]='\0';
 ent=ensort(xyz_all,flg0);
 order=new int[kep];
 entt=new float[kep];
 for(i=0;i<kep;i++) {entt[i]=ent[i];order[i]=i;}
 if(!strstr(flg,"?"))cc.sort(entt,kep,order); 
 rms=0;
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[order[i]],1);
   e=rmsd(2);
   if(TRES.logg>3)cerr<<"  "<<i<<" energy of "<<flg<<": "<<e<<"  "<<entt[i]<<endl;
   rms+=e/kep;
 }
 ent=trim(xyz_all,ent,haa);
 kep=0;
 while(xyz_all[kep])kep++;
 rms0=0;
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[i],1);
   e=rmsd(2);
   rms0+=e/kep;
 }
 if(TRES.logg>3)cerr<<" the total number after "<<flg<<" "<<kep<<"  "<<rms<<"  "<<rms0<<endl;
 delete [] entt;delete [] order;
 return ent;
}

float* Loopy::trim(float **xyz_all,char *flg,int haa,int fe)
{
 int kep,i;
 float e,rms,rms0,*ent,*entt;
 char flg0[100];
 Qsort cc;
 int *order;

 kep=0;
 while(xyz_all[kep])kep++;
 strcpy(flg0,flg);
 i=strlen(flg);
 flg0[i]='?';flg0[i+1]='\0';
 ent=ensort(xyz_all,flg0);
 order=new int[kep];
 entt=new float[kep];
 for(i=0;i<kep;i++) {entt[i]=ent[i];order[i]=i;}
 if(!strstr(flg,"?"))cc.sort(entt,kep,order);
 rms=0;
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[order[i]],1);
   e=rmsd(2);
   if(fe==0) {
   if(TRES.logg>3)cerr<<"  "<<i<<" energy of "<<flg<<": "<<e<<"  "<<entt[i]<<endl;
   }
   else {
   if(TRES.logg>3)cout<<"  "<<i<<" energy of "<<flg<<": "<<e<<"  "<<entt[i]<<endl; 
   }
   rms+=e/kep;
 }
 ent=trim(xyz_all,ent,haa);
 kep=0;
 while(xyz_all[kep])kep++;
 rms0=0;
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[i],1);
   e=rmsd(2);
   rms0+=e/kep;
 }
 if(TRES.logg>3)cerr<<" the total number after "<<flg<<" "<<kep<<"  "<<rms<<"  "<<rms0<<endl;
 delete [] entt;delete [] order;
 return ent;
}


void Loopy::setoff(float **xyz_all,float ct)
{
    Res *r,*r0;
    Atm *a;
    float d;
    Lattice latt;
    int all,i,j;
    Chn *chn=(*pdb)[cid];

    all=0;
    while(xyz_all[all])all++;

    r0=(*chn)[start];
    latt.putoff();
    latt.grdsiz=2.0;
    latt.ready(pdb);
    latt.grdsiz=2;
    latt.radall=15;
    latt.flag=1;
    latt.puton(pdb);
    latt.putoff(r0,end);
    chn->setflg(chn->res,1000000,0);
    chn->setflg(r0,end,1);

    for(j=0;j<all;j++)
    {
      chn->transfer(xyz_all[j],r0,end,1);
      for(r=r0;r;r=r->next)
      {
         if(r->id0>end) break;
         for(a=r->atm;a;a=a->next)
         {
           latt.getcell(a,ct);
           for(i=0;i<latt.nget;i++)
           {
             if(latt.obtain[i]->atm->flag==1) continue;
             d=TRES.distance(a->xyz,latt.obtain[i]->atm->xyz);
             if(d>ct) continue;
             latt.putoff(latt.obtain[i]->atm);
             latt.obtain[i]->atm->flag=1;
           }
         }
      }
    }
}


float **Loopy::generate(float **xyz_all,int cop,float close)
{
   float *ent;
   int i,kep,kep0;
   kep0=0;while(xyz_all[kep0])kep0++;
   cerr<<"doing loop marriage between different loop candidate pairs..."<<endl;
   xyz_all=cross(xyz_all,cop,close);
   kep=0;while(xyz_all[kep])kep++;
   cerr<<"minimizing the new loops from loop marriage..."<<endl;
   for(i=kep0;i<kep;i++) minimize(xyz_all[i],1);
   cerr<<"discarding bad candidates..."<<endl;
   if(kep>1.3*part) {ent=trim(xyz_all,"vw",0); delete [] ent;}
   kep=noclose(xyz_all,close);
   if(TRES.logg>3)cerr<<" the total number after noclose of "<<close<<" :"<<kep<<endl;
   for(i=(int)(1.3*part);i<kep;i++) {delete [] xyz_all[i];xyz_all[i]=0;}
   return xyz_all;
}
 
int **Loopy::allconf(int nall)
{
  Res *r,*r1;
  int **outid;
  int k,n,m,i,j,jj;
 
  n=1; 
  for(r=segment->chn->res;r;r=r->next)
  {
    r1=TRES.rotamer->chn->isres(r->name,0);
    n=n*r1->nummore;
  }
  if(TRES.logg>3)cerr<<" the total number of conformation is :"<<n<<endl;
  if(nall>n) nall=n;
  n=end-start+1;
  outid=new int*[nall+10];
  for(i=0;i<nall+10;i++) outid[i]=0;
  for(i=0;i<nall;i++) outid[i]=new int[n];
  for(i=0;i<nall;i++) for(j=0;j<n;j++) outid[i][j]=0;

  k=0;n=1;
  for(r=segment->chn->res;r;r=r->next)
  {    
    r1=TRES.rotamer->chn->isres(r->name,0);
    jj=n;
    for(i=0;i<n;i++)
    {
      for(j=1;j<r1->nummore;j++)
      {
        for(m=0;m<end-start+1;m++) outid[jj][m]=outid[i][m];
        outid[jj][k]=j; 
        //for(m=0;m<end-start+1;m++) cerr<<outid[jj][m]<<"  "; 
        //cerr<<endl;
        jj++;
        if(jj>=nall) return outid;
      }
    }
    n=jj;
    k++;
  }
  return outid;
}

float Loopy::segmin(int fast)
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
  lat->ishet=0;
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

    //cerr<<rec<<"energy:  "<<only<<"  "<<exp1st[mseq]<<"  "<<dcut<<"  "<<d<<endl;


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

void Loopy::testclosure()
{
  Res *r,*r0;
  Chn *chn=(*pdb)[cid];

  int nseq,i,j,kep;
  float **xyz_all,**xyz_out,*ent;
  float *xyz;
  float rms,e;
  float d;

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
  xyz=new float[nseq];
  chn->transfer(xyz,r0,end,0);

  xyz_all=zipper(0,30,0);
  xyz_out=copy(xyz_all); 

  kep=0;while(xyz_out[kep])kep++;
  for(i=0;i<kep;i++)
  {
    chn->transfer(xyz,r0,end,1);
    segment->chn->transfer(xyz_all[i],1);
    d=closeg();
    if(d==-2) {xyz_all[i]=0;xyz_out[i]=0;continue;}
    segment->chn->transfer(xyz_all[i],0);
  } 
  j=0;
  for(i=0;i<kep;i++)
  {
    if(xyz_all[i]==0) continue;
    xyz_all[j]=xyz_all[i];
    xyz_out[j]=xyz_out[i];
    j++;
  }
  kep=j;
  if(TRES.logg>3)cerr<<"**********initial**************"<<endl; 
  ent=ensort(xyz_out,"v?");
  kep=0;while(xyz_out[kep])kep++;
  rms=0;
  chn->transfer(xyz,r0,end,1);
  for(i=0;i<kep;i++)
  {
   segment->chn->transfer(xyz_out[i],1);
   e=rmsd(2);
   if(TRES.logg>3)cerr<<"  "<<i<<" energy of "<<e<<"  "<<ent[i]<<endl;
   rms+=e/kep;
  }
  if(TRES.logg>3)cerr<<"average :"<<rms<<endl;
  delete [] ent;

  if(TRES.logg>3)cerr<<"**********random treawk****** **"<<endl;
  ent=ensort(xyz_all,"v?");
  
  rms=0;
  chn->transfer(xyz,r0,end,1);
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[i],1);
   e=rmsd(2);
   if(TRES.logg>3)cerr<<"  "<<i<<" energy of "<<e<<"  "<<ent[i]<<endl;
   rms+=e/kep;
 }
 if(TRES.logg>3)cerr<<"average :"<<rms<<endl;
 delete [] ent;


 for(i=0;i<kep;i++)
 {
    chn->transfer(xyz,r0,end,1);   
    segment->chn->transfer(xyz_all[i],1);
    cycle=200;outlet=0.3;flat=300;
    segmin(1);
    segment->chn->transfer(xyz_all[i],0);
 }
  if(TRES.logg>3)cerr<<"**********after minimization********"<<endl;
  ent=ensort(xyz_all,"v?");
  rms=0;
  chn->transfer(xyz,r0,end,1);
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[i],1);
   e=rmsd(2);
   if(TRES.logg>3)cerr<<"  "<<i<<" energy of "<<e<<"  "<<ent[i]<<endl;
   rms+=e/kep;
 }
 if(TRES.logg>3)cerr<<"average :"<<rms<<endl;
 delete [] ent;


  for(i=0;i<kep;i++)
  {
    chn->transfer(xyz,r0,end,1);
    segment->chn->transfer(xyz_out[i],1);
    cycle=200;outlet=2;flat=300;
    segmin(1);
    cycle=200;outlet=0.3;flat=300;
    segmin(1);
    segment->chn->transfer(xyz_out[i],0);
  }
  if(TRES.logg>3)cerr<<"**********channeled treak*****"<<endl;
  ent=ensort(xyz_out,"v?");
  rms=0;
  chn->transfer(xyz,r0,end,1);
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_out[i],1);
   e=rmsd(2);
   if(TRES.logg>3)cerr<<"  "<<i<<" energy of "<<e<<"  "<<ent[i]<<endl;
   rms+=e/kep;
 }
 if(TRES.logg>3)cerr<<"average :"<<rms<<endl;
 delete [] ent;
}

int Loopy::checkomega() {
for(Res *r=segment->chn->res;r;r=r->next){
	Atm *a=r->isatmid(0);
	float d=a->gettorsionangle();
	if(d>500||d<-500) continue;
	if(d>0) d=180-d;
	else 	d=180+d;

	d=fabs(d);

	if(d>5) return 0;
}
return 1;
}

void Loopy::checkomega(FILE *fp) {
for(Res *r=segment->chn->res;r;r=r->next){
        Atm *a=r->isatmid(0);
        float e=a->gettorsionangle();
	float d=e;
        //if(d>500||d<-500) continue;
        if(d>0) d=180-d;
        else    d=180+d;
        d=fabs(d);
	if(fp) fprintf(fp,"%c%i %f %f\n",r->name,r->id,e,d);
}
}


