#include"source.h"

void Segen::initial() {
 randomtorsion=0;
 pdb=0;
 end=start=-1;
 segment=0;
 lat=new Lattice;
 step=0.17;
 cutoff=6.;
 last=0;
 next=0;
 direct=1;
 arbt=100;
 id=0;
 near=0;
 init=0;
 part=25;
 revs=1;
 pdbout=1;
 outlet=0.5;
 flat=1;
 cycle=1000;
 secd='-';
 cid='-';
 part0=50;
 discut=3.;
 
 //
 onlyenergy=0; //onlyenergy used to indicate when minimizing, energy is the only factor to accout for
	       //thorough distance geometry may exist.
 xyzsave=0;
 nomin=0;
 databaseonly=0;
 chiangle=0;
 fapr=0;
 numclose=50;
 bound=0;
 mutate=0;
 ranseed=18120;
 flex=0;
 flexrot=10;
 rotatm=0;
 onlysidechain=0;
 onlybound=0;
 smoothclash=0;
 strcpy(boundtag,"HP");
 //charge=0;
 disc=0;
 strcpy(dforce,"v");
 dace=0;
 randcoil=0;
 rmsdinsert=0;
}

Segen::Segen()
{
 initial();
}

Segen::Segen(Segen *ss)
{
 initial();
 end=ss->end;
 start=ss->start;
 step=ss->step; 
 cutoff=ss->cutoff;
 last=ss->last;
 direct=ss->direct;
 arbt=ss->arbt;
 id=ss->id+1;
 part=ss->part;
 revs=ss->revs;
 pdbout=ss->pdbout;
 outlet=ss->outlet;
 flat=ss->flat;
 cycle=ss->cycle;
 secd=ss->secd;
 cid=ss->cid;
 part0=ss->part0;
 discut=ss->discut;
 ranseed=ss->ranseed;
 nomin=ss->nomin;
 numclose=ss->numclose;
 //
 flex=ss->flex;
 onlyenergy=ss->onlyenergy;
 //charge=ss->charge;
 onlysidechain=ss->onlysidechain;
 onlybound=ss->onlybound;
 flexrot=ss->flexrot;
 fapr=ss->fapr;
 randcoil=ss->randcoil;
 databaseonly=ss->databaseonly;
 numclose=ss->numclose;
 smoothclash=ss->smoothclash;
 randomtorsion=ss->randomtorsion;
}

Segen::~Segen()
{
  if(segment) {delete segment;segment=0;}
  if(lat)     {delete lat;lat=0;}
  if(init) {Strhandler cc; cc.floatdel(init);init=0;}
  if(near) {Strhandler cc; cc.floatdel(near);near=0;}
  if(chiangle) {delete chiangle;chiangle=0;}
  if(bound)   delete bound;bound=0;
  if(rotatm) delete [] rotatm;rotatm=0;
  if(mutate) delete mutate;mutate=0;
  if(dace) delete dace;dace=0;
  if(disc) delete disc;disc=0;
  if(next)    {delete next;next=0;}
  if(xyzsave) {delete [] xyzsave;xyzsave=0;}
}

void Segen::clear() {
  if(segment) {delete segment;segment=0;}
  if(lat)     {delete lat;lat=0;}
  if(init) {Strhandler cc; cc.floatdel(init);init=0;}
  if(near) {Strhandler cc; cc.floatdel(near);near=0;}
  if(chiangle) {delete chiangle;chiangle=0;}
  if(bound)   delete bound;bound=0;
  if(rotatm) delete [] rotatm;rotatm=0;
  if(mutate) delete mutate;mutate=0;
  if(dace) delete dace;dace=0;
  if(disc) delete disc;disc=0;
  if(next)    {delete next;next=0;}
  if(xyzsave) {delete [] xyzsave;xyzsave=0;}
}

void Segen::ready()
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
  
  next=new Segen(this);
  next->id=this->id+1;
  next->pdb=pdb->next;
  next->last=this;
  next->next=0;

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

  //setup disc;
  if(disc) delete disc;disc=0;
  if(disc==0) disc=new Disc;
  disc->setupall(pdb,cutoff,1);

  if(next->disc) delete next->disc;next->disc=0;
  if(next->disc==0) next->disc=new Disc;
  next->disc->setupall(next->pdb,next->cutoff,1);
}

int Segen::create()
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
     if(TRES.logg)cerr<<"warning!!! this is not a loop prediction problem.\n"; 
     if(TRES.logg)cerr<<"it is a protein structure prediction"<<endl ;
     //exit(0); 
     direct=1;
     r=chn->res;
     start=chn->res->id0;
     end=chn->lastres()->id0;
     n=1;
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

void Segen::create(Res *s,int n)
{
  if(segment) delete segment;segment=0;
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
        if(r->temp) delete [] r->temp;r->temp=0;
	r->nemp=9;
	r->temp=new float[9];
	r1=(*s->chn)[r->id0];
	TRES.copy(r1->temp,r->temp,9);
  }
  //strcpy("unknown",segment->name);
  strcpy(segment->name,"unknown");
}
Pdb* Segen::mycreate(Res *s,int n)
{
  Pdb *segm;
  segm=new Pdb;
  segm->chn=new Chn;
  //segm->chn->pdb=segment;
  segm->chn->pdb=segm;
  segm->chn->create(s,n);
  segm->chn->id=s->chn->id;
  segm->chn->start=s->chn->start;
  segm->name=new char[100];
  segm->configure();
  Res *r1;
  for(Res *r=segm->chn->res;r;r=r->next) {
	r->id0+=s->id0;
	r->id+=s->id;
        if(r->temp) delete [] r->temp;
	r->nemp=9;
	r->temp=new float[9];
	r1=(*s->chn)[r->id0];
	TRES.copy(r1->temp,r->temp,9);
  }
  //strcpy("unknown",segment->name);
  strcpy(segm->name,"unknown");
  return segm;
}
void Segen::create(char *seqn)
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

void Segen::randmz(int rand)
{
  Res *r;
  Atm *a;
  int j,i,n;
  Rotate rot;

  //srandom(seed); //random initiializaton

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

void Segen::randmzHelix(int rand)
{
  Res *r;
  Atm *a,*a0;
  int j,i,n;
  Rotate rot;
  float rt;
  //srandom(seed); //random initiializaton

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
      if(direct==1) {
      if(a->tatm->id==1) ee=j-rt-60;
      else if(a->tatm->id==2) ee=j-rt-43;
      }
      else {
      if(a->tatm->id==1) ee=j+rt+60;
      else if(a->tatm->id==2) ee=j+rt+43;
      }
      if(direct==1) rot.rotate(a,ee,1);
      else          rot.rotate(a,ee,0);
    }
  }
  if(TRES.logg>3)segment->chn->write("s1");
  segment->chn->dihedral(stdout);
}

int Segen::randmzChiangle(int nu)
{
  Res *r;
  Atm *a,*a0;
  int i,n;
  Rotate rot;
  float rt;
  Chn *chn=(*pdb)[cid];
  if(segment==0||segment->chn==0) return 0;
  if(chn==0) return 0;
  Res *rr0=(*chn)[start];
  if(rr0==0) return 0;
  if(chiangle==0) return 0;
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
  if(last==0||last0==0) return 0;
  //return 1;
  //segment->chn->dihedral("s4"); 
  for(r=segment->chn->res;r;r=r->next)
  {
    if(r->id0>end) break;
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
	  if(TRES.logg)cerr<<"distance ..."<<TRES.distance(g.coo,r->isatmid(3)->xyz)<<" "<<TRES.distance(g.coo,w2)<<" "<<TRES.distance(g.coo,w3)<<" "<<rt<<endl;
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
      if(direct==1){
      ee=-rt+ee;
      }
      else {
	ee=rt-ee;
      }
      if(fabs(ee)>360) continue;
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
          if(TRES.logg)cerr<<" again distance ..."<<TRES.distance(g.coo,r->isatmid(3)->xyz)<<" "<<TRES.distance(g.coo,w2)<<" "<<TRES.distance(g.coo,w3)<<" "<<rt<<endl;
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
 
  if(TRES.logg>3) {
   	segment->chn->write("s3");
   	segment->chn->dihedral(stderr);
  	for(i=0;i<chia->number;i++){
		if(TRES.logg)cerr<<"angle: ";
		for(int j=0;j<chia->numa[i];j++)  if(TRES.logg)cerr<<" "<<chia->angle[i][j];
		if(TRES.logg)cerr<<endl;
  	}
  }
  
  return 1;
}

int Segen::randmzChianglefix(int ang1,int ang2)
{
  Res *r;
  Atm *a,*a0;
  Rotate rot;
  Res *last,*last0;
  float rt,ee;
  float *w1,*w2,*w3,*w4;

  Chn *chn=(*pdb)[cid];
  Chiangle *chia=chiangle;
  if(chia==0) return 0; 

  if(direct==1) {
    last=(*segment->chn)[start];
    last0=(*chn)[start];
  }
  else {
    last=(*segment->chn)[end];
    last0=(*chn)[end];
  }

  segment->chn->dihedral("s4"); 
  int nn=0;

  for(r=segment->chn->res;r;r=r->next)
  {	
    if(nn>=chia->number||chia->angle[nn]==0) continue;
    for(a=r->atm;a;a=a->next)
    {
      if(a->tatm->name[1]=='H') break;
      if(a->tatm->rotate==0) continue;
      if(a->tatm->id>=3) continue;
      if(r->name=='P'&&a->id==1) continue;
 
      if(direct==1&&r==segment->chn->res&&a->tatm->id==1) {
	  w1=last0->temp+3;
	  w2=r->atm->xyz;
	  w3=r->atm->next->xyz;
	  w4=r->atm->next->next->xyz;
          rt=TRES.dihedral(w1,w2,w3,w4);
      }
      else if(direct==1&&r->next==0&&a->tatm->id==2) {
	  if(databaseonly==0) continue;  //because the last is going to be superimposed
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
	  //cerr<<"distance ..."<<TRES.distance(g.coo,r->isatmid(3)->xyz)<<" "<<TRES.distance(g.coo,w2)<<" "<<TRES.distance(g.coo,w3)<<" "<<rt<<endl;
	  //continue;
      }
      else if(direct==0&&r==last&&a->tatm->id==2) {
	  a0=last0->atm->next->next;
          rt=a0->dihedral(0); 
      }
      else rt=a->dihedral(0);
      
      if(a->tatm->id==1) ee=chia->getangle(nn,0);
      else if(a->tatm->id==2) ee=chia->getangle(nn,1);
    
      if(a->tatm->id==1&&fabs(ee-rt)<ang1) continue;
      if(a->tatm->id==2&&fabs(ee-rt)<ang2) continue;

      if(direct==1) ee=-rt+ee;
      else          ee=rt-ee;
      if(direct==1) rot.rotate(a,ee,1);
      else          rot.rotate(a,ee,0);
    }
    nn++;
  } 
  segment->chn->dihedral("s3");
  return 1;
}


void Segen::randmzSheet(int rand)
{
  Res *r;
  Atm *a,*a0;
  int j,i,n;
  Rotate rot;
  float rt;
  //srandom(seed); //random initiializaton

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
      if(direct==1) {	
      if(a->tatm->id==1) ee=j-rt-100;
      else if(a->tatm->id==2) ee=j-rt+130;
      }
      else {
	if(a->tatm->id==1) ee=j+rt+100;
        else if(a->tatm->id==2) ee=j+rt-130;
      }
      if(direct==1) rot.rotate(a,ee,1);
      else          rot.rotate(a,ee,0);
      //rt=a->dihedral(0);

      //cerr<<a->name<<a->id<<" "<<rt<<endl;
    }
  }
}



float Segen::closeg(Segen *seg1,Segen *seg2)
{ 
  Segen *segen[2];
  Res *r,*last,*last0,*rr;
  Atm *aa[2],*aa0;
  float d,d1,dpre,cut;
  float xr[3],xr1[3],xr2[3];
  float **coeff=0,*unknown=0,lamba[12];
  int i,j,m,k,mitr,mseq,meg;
  int rec,rec1;
  Rotate rot; 
  coeff=new float*[12];
  for(i=0;i<12;i++) coeff[i]=0;
  unknown=0;
  last=0; 
  last0=0;
  mseq=0;

  segen[0]=seg1; segen[1]=seg2;
  if(segen[0]->end!=segen[1]->start) 
  {
    if(TRES.logg)cerr<<" the two segments not contigent"<<endl;
    if(coeff) delete [] coeff;coeff=0;
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
     if(coeff) delete [] coeff;coeff=0;
     return segen[1]->closeg(); 
  } 
  else if(segen[0]->segment&&segen[1]->segment==0) 
  { 
     if(coeff) delete [] coeff;coeff=0;
     return segen[0]->closeg(); 
  } 
  else if(segen[0]->segment==0&&segen[1]->segment==0) 
  { 
     if(TRES.logg)cerr<<" no existing of two segments for connecting"<<end; 
     if(coeff) delete [] coeff;coeff=0;
     return -1; 
  } 

  unknown=new float[mseq];
  for(i=0;i<12;i++) coeff[i]=new float[mseq];

  rec=0;rec1=0;dpre=1000;cut=1000;
  
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
	  {if(TRES.logg)cerr<<"warning in closeg(int) aa0=0"<<endl;exit(0);} 

          for(j=0;j<3;j++)xr[j]=aa[i]->xyz[j]-aa0->xyz[j];
          d=sqrt(xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2]);
          
	  if(d==0) {m++;continue;} 
          if(d<0.001) d=1.;

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
	  if(fabs(d)>360) {m++;continue;}
          if(meg==0) rot.rotate(aa[i],d,last->id0,1); 
          else       rot.rotate(aa[i],d,last0->id0,0);
          m++;
        }
      }
    }
    if(rec>100) break;
  }
  for(i=0;i<12;i++) if(coeff[i]){delete [] coeff[i]; coeff[i]=0;}
  if(coeff) delete [] coeff;coeff=0;
  if(unknown) {delete [] unknown; unknown=0;}
  //if(cut>2) cut=-2;
  //else if(cut>offset+outlet+0.1)
  if(cut<=0) return -2;
  else if(cut>outlet) return -2;
  return cut;
}
   

void Segen::writeseg(float **xyz_all,char *s,int nn)
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

void Segen::write(float **xyz_all,char *s,int nn)
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
  delete [] xyz;xyz=0;
}

float **Segen::predt0()
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

  xyz_out=new float*[revs*2];
  for(i=0;i<revs*2;i++) xyz_out[i]=0;
  xyz_out[0]=new float[nseq];
  chn->transfer(xyz_out[0],r0,end,0);
 
  if(end!=start) {
    //generate random initial conformation either based on ab-inito or database method
    //the total number of initial conformation would be arbt
    xyz_all=next->zipper(discut*(end-start+1),1);
  }
  else {
    xyz_all=new float*[2];
    xyz_all[0]=new float[nseq];
    chn->transfer(xyz_all[0],r0,end,0);
    xyz_all[1]=0;
  }
  if(end-start+1<=2) {
     strcpy(dforce,"v");
     ent=scap(xyz_all);if(ent) delete [] ent; ent=0;
     strcpy(dforce,"v");
     if(end!=start) ent=ensort(xyz_all);if(ent) delete [] ent; ent=0;
     strcpy(dforce,"shve");
     if(end!=start&&fapr==0) {ent=trim(xyz_out,arbt);if(ent) delete [] ent; ent=0; }
     strcpy(dforce,"hv");
     if(end!=start&&fapr==1) {ent=trim(xyz_out,arbt);if(ent) delete [] ent; ent=0; }
     kep=0;while(xyz_all[kep])kep++;
     for(i=0;i<revs&&i<kep;i++)
     {
       xyz_out[i]=xyz_all[i];
     }
     for(i=revs;i<kep;i++) delete [] xyz_all[i];
     delete [] xyz_all;
     xyz_all=copy(xyz_out);
     goto re1;
  }

  chn->transfer(xyzorg,r0,end,1);

  if(fapr==0) {
     next->cycle=200;next->outlet=0.3*(end-start);next->flat=300;
     ent=next->expand(xyz_all,0); if(ent) delete [] ent; ent=0;
  }

  next->cycle=200;next->outlet=0.3;next->flat=111;
  if(fapr==0) ent=next->expand(xyz_all,0); if(ent) delete [] ent; ent=0;
  strcpy(dforce,"v");
  ent=next->trim(xyz_all,arbt);if(ent) delete [] ent; ent=0;
  init=xyz_all;
  j=1;
  cycle=200;outlet=0.3;flat=111;
  
  kep=0;while(init&&init[kep])kep++;
  if(revs>kep) revs=kep;
  next->revs=revs;
  chn->transfer(xyzorg,r0,end,1); 
  for(i=0;i<revs;i++)
  {
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
          sprintf(namem,"%s_loopy.%i.pdb",pdb->name,i);
          pdb->write(namem);
      }   
  }

  chn->transfer(xyzorg,r0,end,1);  

  if(fapr==0) {
    strcpy(dforce,"v");
    ent=trim(xyz_out,arbt,1);if(ent) delete [] ent; ent=0;
    strcpy(dforce,"sv");
    ent=trim(xyz_out,arbt,1);if(ent) delete [] ent; ent=0;
    strcpy(dforce,"shv");
    ent=trim(xyz_out,arbt,1);if(ent) delete [] ent; ent=0;
    strcpy(dforce,"shvdcD");
    ent=trim(xyz_out,arbt,1);if(ent) delete [] ent; ent=0;
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

float **Segen::predt(int flg)
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
for(i=flg*kep/revs;i<(flg+1)*kep/revs;i++)
{
  xyz_all[j]=new float[nseq];
  TRES.copy(init[i],xyz_all[j],nseq);
  j++;
  if(j>=10&&flg>0) break;
}

//side-chain
strcpy(dforce,"v");
ent=scap(xyz_all);if(ent) delete [] ent; ent=0;

if(fapr==0) ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;

//filter on vw 
strcpy(dforce,"vw");
ent=trim(xyz_all,0);if(ent) delete [] ent; ent=0;

if(fapr==1) {
  chn->transfer(xyz,r0,end,1);
  if(xyz) {delete [] xyz; xyz=0;}
  return xyz_all;
}


//filter on next->vw
strcpy(dforce,"vw");
ent=next->trim(xyz_all,0); if(ent) delete [] ent; ent=0;

kep=0;while(xyz_all[kep])kep++;
if(kep>1.3*part) {ent=trim(xyz_all,part); if(ent) delete [] ent; ent=0;}

cop=1; 
 
for(i=0;i<2;i++)
{
  xyz_all=generate(xyz_all,cop,close);
  cerr<<"***cycle***: "<<i<<"  "<<j<<endl;
}

/*
for(i=0;i<2;i++)
{
  xyz_all=next->generate(xyz_all,cop,close);
  cerr<<"***cycle for backbone only***: "<<i<<" "<<j<<endl;
}

//side-chain
strcpy(dforce,"v");
ent=scap(xyz_all);if(ent) delete [] ent; ent=0;
ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;


if(part>10) part=part/2;
for(i=0;i<2;i++)
{
  xyz_all=generate(xyz_all,cop,close);
  cerr<<"***cycle***: "<<i<<"  "<<j<<endl;
}
*/

part=next->part;
strcpy(dforce,"vw");
ent=trim(xyz_all,arbt);if(ent) delete [] ent; ent=0;

chn->transfer(xyz,r0,end,1);
if(xyz) {delete [] xyz; xyz=0;}
return xyz_all;
}



float **Segen::predt0save()
{ 
  Res *r,*r0; 
  int nseq,k,i,j,kep;
  float **xyz_all,**xyz_out,*ent;
  Chn *chn=(*pdb)[cid];

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
  
  //create itself
  xyz_out=new float*[revs*2+100];
  for(i=0;i<revs*2+100;i++) xyz_out[i]=0;
  xyz_out[0]=new float[nseq];
  chn->transfer(xyz_out[0],r0,end,0);
  
  //create itself
  float *xyz0=chn->gettransfer(r0,end);
  
  //get sample
  xyz_all=next->zipper(30*(end-start),1);

  //delete
  strcpy(next->dforce,"v");
  ent=next->trim(xyz_all,next->part0);if(ent) delete [] ent; ent=0;

  //minimize
  next->cycle=20; next->numclose=20; next->step=0.17; next->flat=1; next->outlet=0.5;
  strcpy(next->dforce,"v");
  ent=next->minimizefix(xyz_all,1); if(ent) delete [] ent; ent=0;
  
  //delete
  strcpy(next->dforce,"v");
  ent=next->trim(xyz_all,2*part);if(ent) delete [] ent; ent=0;
  
  //sidechain
  strcpy(dforce,"v");
  ent=scap(xyz_all);if(ent) delete [] ent; ent=0;

  //minimize
  cycle=50; numclose=50;step=0.17; flat=111; outlet=0.5;
  strcpy(dforce,"v");
  ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;
  
  //delete
  strcpy(dforce,"v");
  ent=trim(xyz_all,part);if(ent) delete [] ent; ent=0;
   
  //check revs
  init=xyz_all;  
  kep=0;while(init&&init[kep])kep++;
  if(revs>kep) revs=kep;
  next->revs=revs;
  //end

  j=1;
  cycle=50; numclose=50;step=0.17; flat=111; outlet=0.5;
  
  for(i=0;i<revs;i++)
  {
    xyz_all=predt(i);
    
    xyz_out[j++]=xyz_all[0];
    kep=0;while(xyz_all[kep])kep++;
    for(k=1;k<kep;k++) delete [] xyz_all[k];
    delete [] xyz_all;
  }
  xyz_all=copy(xyz_out);
  strcpy(dforce,"v");
  ent=trim(xyz_out,arbt);if(ent) delete [] ent; ent=0;
  strcpy(dforce,"sv");
  ent=trim(xyz_out,arbt);if(ent) delete [] ent; ent=0;
  strcpy(dforce,"shv");
  ent=trim(xyz_out,arbt);if(ent) delete [] ent; ent=0;
  strcpy(dforce,"shve");
  ent=trim(xyz_out,arbt);if(ent) delete [] ent; ent=0;
  
  if(pdbout)
  {
      int totfile;
      char namem[100];
      totfile=0;while(xyz_all[totfile])totfile++;
      r0=(*chn)[start];
      for(i=1;i<totfile;i++)
      {
          chn->transfer(xyz_all[i],r0,end,1);
          if(totfile>1) sprintf(namem,"%s_loopy.%i.pdb",pdb->name,i);
	  else          sprintf(namem,"%s_loopy.pdb",pdb->name);
          pdb->write(namem);
      }   
  }

  //delete init
  kep=0;while(init&&init[kep])kep++;
  for(k=0;k<kep;k++) delete [] init[k];
  if(init) delete [] init;init=0;

  //delete xyz_out
  kep=0;while(xyz_out&&xyz_out[kep])kep++;
  for(k=0;k<kep;k++) delete [] xyz_out[k];
  if(xyz_out) delete [] xyz_out;xyz_out=0;
  
  //set original
  chn->transfer(xyz0,r0,end,1);
  if(xyz0) delete [] xyz0;xyz0=0;

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

float **Segen::predtsave(int flg)
{
Res *r,*r0;
float **xyz_all,*xyz,*ent;
int i,j,kep,nseq;
float close;
int cop;

//range
Chn *chn=(*pdb)[cid];
r0=(*chn)[start];
nseq=0;
for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;

//save original
xyz=new float[nseq];
chn->transfer(xyz,r0,end,0);

close=0.1*(end-start);

//choose needed colony

kep=0;while(init&&init[kep])kep++;
xyz_all=new float*[2*kep];
for(i=0;i<2*kep;i++)xyz_all[i]=0;
j=0;

for(i=flg*kep/revs;i<(flg+1)*kep/revs;i++)
{
  xyz_all[j]=new float[nseq];
  TRES.copy(init[i],xyz_all[j],nseq);
  j++;
}

if(flg) return xyz_all;

//loop fusion 
strcpy(dforce,"v");
cop=4; 
for(i=0;i<2;i++)
{
  xyz_all=generate(xyz_all,cop,close);
  strcpy(dforce,"v");
  ent=trim(xyz_all,part);if(ent) delete [] ent; ent=0;
  if(TRES.logg>3) cerr<<"***cycle***: "<<i<<"  "<<j<<endl;
}

//minimize
cycle=50; numclose=50;step=0.17; flat=111; outlet=0.5;
ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;

//delete and order
strcpy(dforce,"shvw");
ent=trim(xyz_all,part);if(ent) delete [] ent; ent=0;

//save to original
chn->transfer(xyz,r0,end,1);
if(xyz) {delete [] xyz; xyz=0;}
return xyz_all;
}



float **Segen::predt0old()
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

  xyz_out=new float*[revs*2];
  for(i=0;i<revs*2;i++) xyz_out[i]=0;
  
  if(end!=start) {
    //generate random initial conformation either based on ab-inito or database method
    //the total number of initial conformation would be arbt
    xyz_all=next->zipper(discut*(end-start+1),1);
  }
  else {
    xyz_all=new float*[2];
    xyz_all[0]=new float[nseq];
    chn->transfer(xyz_all[0],r0,end,0);
    xyz_all[1]=0;
  }
  if(end-start+1<=2) {
     strcpy(dforce,"v");
     ent=scap(xyz_all);if(ent) delete [] ent; ent=0;
     if(end!=start) ent=ensort(xyz_all);if(ent) delete [] ent; ent=0;
     if(end!=start&&fapr==0) {ent=trim(xyz_out,arbt);if(ent) delete [] ent; ent=0;}
     if(end!=start&&fapr==1) {ent=trim(xyz_out,arbt);if(ent) delete [] ent; ent=0;}
     kep=0;while(xyz_all[kep])kep++;
     for(i=0;i<revs&&i<kep;i++)
     {
       xyz_out[i]=xyz_all[i];
     }
     for(i=revs;i<kep;i++) delete [] xyz_all[i];
     delete [] xyz_all;xyz_all=0;
     xyz_all=copy(xyz_out);
     goto re1;
  }
  
  if(fapr==0) {
     next->cycle=200;next->outlet=0.3*(end-start);next->flat=300;
     strcpy(dforce,"v");
     ent=next->expand(xyz_all,0); if(ent) delete [] ent; ent=0;
  }

  next->cycle=200;next->outlet=0.3;next->flat=111;
  strcpy(dforce,"v");
  if(fapr==0) ent=next->expand(xyz_all,0); if(ent) delete [] ent; ent=0;
  strcpy(dforce,"v");
  ent=next->trim(xyz_all,arbt);if(ent) delete [] ent; ent=0;
  //ent=scap(xyz_all,"v");if(ent) delete [] ent; ent=0;
  strcpy(dforce,"v");
  ent=scap(xyz_all);if(ent) delete [] ent; ent=0;
  ent=minimize(xyz_all,1);if(ent) delete [] ent; ent=0;
  strcpy(dforce,"v");
  //ent=trim(xyz_all,"v",0);if(ent) delete [] ent; ent=0;
  ent=trim(xyz_all,0);if(ent) delete [] ent; ent=0;
  init=xyz_all;
  j=0;
  cycle=200;outlet=0.3;flat=111;
  
  for(i=0;i<revs;i++)
  {
    xyz_all=predtold(i);
    kep=0;while(xyz_all[kep])kep++;
    if(kep>0) xyz_out[j++]=xyz_all[0];
    for(k=1;k<kep;k++) delete [] xyz_all[k];
    delete [] xyz_all;xyz_all=0;
  }
  
  xyz_all=copy(xyz_out);

  re1:

  if(pdbout)
  {
      int totfile;
      char namem[100];
      totfile=0;while(xyz_all[totfile])totfile++;
      r0=(*chn)[start];
      for(i=0;i<totfile;i++)
      {
          chn->transfer(xyz_all[i],r0,end,1);
          sprintf(namem,"%s_looper_%i.pdb",pdb->name,i);
          pdb->write(namem);
      }   
  }

  if(fapr==0) {
    ent=trim(xyz_out,arbt,1);if(ent) delete [] ent; ent=0;
    ent=trim(xyz_out,arbt,1);if(ent) delete [] ent; ent=0;
    ent=trim(xyz_out,arbt,1);if(ent) delete [] ent; ent=0;
    ent=trim(xyz_out,arbt,1);if(ent) delete [] ent; ent=0;
  }


  kep=0;while(init&&init[kep])kep++;
  for(k=0;k<kep;k++) delete [] init[k];
  delete [] init;init=0;
  kep=0;while(xyz_out[kep])kep++;
  for(k=0;k<kep;k++) delete [] xyz_out[k];
  delete [] xyz_out; xyz_out=0;

  chn->transfer(xyzorg,r0,end,1);

  return xyz_all;
}

float **Segen::predtold(int flg)
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
for(i=flg*kep/revs;i<(flg+1)*kep/revs;i++)
{
  xyz_all[j]=new float[nseq];
  TRES.copy(init[i],xyz_all[j],nseq);
  j++;
  if(j>=10&&flg>0) break;
}

//side-chain

ent=scap(xyz_all);if(ent) delete [] ent; ent=0;

if(fapr==0) ent=minimize(xyz_all,1);if(ent) delete [] ent; ent=0;

//filter on vw 
ent=trim(xyz_all,0);if(ent) delete [] ent; ent=0;

if(fapr==1) {
  chn->transfer(xyz,r0,end,1);
  if(xyz) {delete [] xyz; xyz=0;}
  return xyz_all;
}


//filter on next->vw
ent=next->trim(xyz_all,0);if(ent) delete [] ent; ent=0;

kep=0;while(xyz_all[kep])kep++;
if(kep>1.3*part) {ent=trim(xyz_all,part); if(ent) delete [] ent; ent=0;}

cop=1; 
 
for(i=0;i<2;i++)
{
  xyz_all=generate(xyz_all,cop,close);
  if(TRES.logg)cerr<<"***cycle***: "<<i<<"  "<<j<<endl;
}

for(i=0;i<2;i++)
{
  xyz_all=next->generate(xyz_all,cop,close);
 if(TRES.logg) cerr<<"***cycle for backbone only***: "<<i<<" "<<j<<endl;
}

//side-chain
ent=scap(xyz_all);if(ent) delete [] ent; ent=0;
ent=minimize(xyz_all,1);if(ent) delete [] ent; ent=0;

if(part>10) part=part/2;
for(i=0;i<2;i++)
{
  xyz_all=generate(xyz_all,cop,close);
  if(TRES.logg)cerr<<"***cycle***: "<<i<<"  "<<j<<endl;
}
part=next->part;

ent=trim(xyz_all,arbt);if(ent) delete [] ent; ent=0;

chn->transfer(xyz,r0,end,1);
if(xyz) {delete [] xyz; xyz=0;}
return xyz_all;
}

int Segen::setend(int flg)
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

float Segen::closeg()
{
  Res *r0,*r,*r1,*last,*last0;
  Atm *aa[2],*aa0;
  float dn,dca,d;
  float xr[3],xr1[3],xr2[3],cut;
  float **coeff=0,*unknown=0,off,offset;
  float lamba[6];
  int i,j,k,m,mseq;
  int rec,rec1,*order=0,nuse;
  Qsort qsort;
  Rotate rot;
  Chn *chn=(*pdb)[cid];
  coeff=new float*[6];
  for(i=0;i<6;i++)coeff[i]=0;
  unknown=0;
//middle exist!    

  r0=(*segment->chn)[start];
  mseq=0;
  
  if(segment==0) {
    if(coeff) delete [] coeff;coeff=0;
    if(TRES.logg)cerr<<"no closing segment"<<endl;
    return -1;
  }

  if(direct==1) 
  {
    
    r1=(*chn)[end];
    if(r1->next==0) 
    {
        if(coeff) delete [] coeff;coeff=0;
	if(TRES.logg)cerr<<"no closing segment"<<endl;
	return -1;
    }
  }
  else if(direct==0) 
  {
    
    if(start==chn->res->id0)
    {
	if(coeff) delete [] coeff;coeff=0;
	if(TRES.logg)cerr<<"no closing segment"<<endl;
	return -1;
    }
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
    if(coeff) delete [] coeff;coeff=0;
    if(TRES.logg)cerr<<"no segment close process!"<<endl;
    return -1;
  }

  unknown=new float[mseq];
  order=new int[mseq];
  for(i=0;i<6;i++) coeff[i]=new float[mseq+1];

  rec=0;rec1=0;off=1000;cut=10000;
  //int nbig=0;

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
     
    if(TRES.logg>3)fprintf(stderr,"%3i %3i %c%d: %s:%8.3f  %s:%8.3f\n",rec,rec1,
            aa[1]->res->name,aa[1]->res->id,
            aa[1]->name,dn,aa[1]->next->name,dca);
    
    d=dn+dca;
    //cerr<<rec<<" "<<d<<endl;
    if(fabs(d-off)<0.01*d||d<outlet) {cut=d;break;} 
    if(fabs(off-d)<0.001) {cut=d; break;}
    //if(step<0.17&&d>off) {cut=d;break;}
    //if(d>off) step=0.085;
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
	if(d<0.001) d=1.;
	if(rotatm&&ifrotatmexist(aa[i])==0) {m++;continue;}
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
    if(TRES.logg>3)cerr<<rec1<<endl;
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
	if(fabs(d)>360) {m++;continue;}
        rot.rotate(aa[i],d,last->id0,direct);
        m++;
      }
    }

    if(rec>100) break;
  }

  for(i=0;i<6;i++) if(coeff[i]) {delete [] coeff[i];coeff[i]=0;}
  if(coeff) delete [] coeff;coeff=0;
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



float Segen::rmsd(int flg)
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

float Segen::rmsdmutate(int flg)
{
Res *r,*rr;
Atm *a,*aa;
int i,j;
float d;

Chn *chn=(*pdb)[cid];

int nlen=strlen(mutate->sqnto);

d=0;
j=0;
for(i=start;i<=end;i++)
{
  r=(*chn)[i];
  rr=(*segment->chn)[i];
  if(r==0||rr==0) continue;
  if(mutate&&rmsdinsert==0) {	
	int m;
	int mgo=0;
	int nsz=2;
	if(mutate->sharp>2) nsz=2;
	for(m=r->id0-nsz;m<=r->id0+nsz;m++) {
		Res *t=chn->isres0(m);
		if(t==0) continue;
		int n1=mutate->compare[t->id0];		
		int n2=mutate->owner->match[n1];
		if(n2!=-1) {
			mgo=1;
			break;
		}
	}
	if(mgo==0&&TRES.logg>3)  cerr<<"ignored in rmsdmtuate:"<<r->name<<r->id0<<endl;
	if(mgo==0) continue;	
  }
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

void Segen::printhelp()
{ 
printf("loopy -prm num -obj num-num -ini num -cid char -out num file.pdb\n");
printf("-prm   force parameters. default 1\n");
printf("-prm 1: CHARMM22 with all atom model\n");
printf("-prm 2: AMBER with all atom model\n");
printf("-prm 3: CHARMM with heavy atom model\n");
printf("-prm 4: AMBER with heavy atom model\n");
printf("-obj num-num: start and end id of loops to be predicted\n");
printf("-ini: number of initial conformations generated.default is 100\n");
printf("-cid: chain id.default is the first chain\n");
printf("-out: number of outputs. default is 1\n");
printf("file.pdb   pdb file\n");
}


void Segen::randmz()
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
  
  //srandom(seed+1); //skip random
 
  segp=segment->chn;

  r3=(*segp)[start];
  for(r=r3;r;r=r->next)
  {
    if(r->id0>end) break;
    
    //r1=TRES.rotamer->chn->isres(r->name,0);
    r1=TRES.findrotamer("backbone")->chn->isres(r->name,0);
    j=random()%r1->nummore;
    r2=r1->ismore(j);    
    if(r==r3)
    {
       if(r3->id0==segp->res->id0) r5=(*chn)[start-1];
       else r5=(*segp)[start-1];

       if(r5==0)
       {
          if(TRES.logg)cerr<<"strange in subroutine randmz..."<<endl;
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

    //r1=TRES.rotamer->chn->isres(r->name,0);
    r1=TRES.findrotamer("backbone")->chn->isres(r->name,0);
    j=random()%r1->nummore;
    r2=r1->ismore(j);
    if(r==last)
    {
       if(last->next) r5=last->next;
       else           r5=(*chn)[end+1];
       if(r5==0)
       {
           if(TRES.logg)cerr<<"strange in subroutine randmz..."<<endl;
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

float Segen::getenergy(Res *s,int to) {
	Res *t=(*s->chn)[to];
	if(t==0) t=s->chn->lastres();
	to=t->id0;
	pdb=s->chn->pdb;
	cid=s->chn->id;	
	start=s->id0;end=to;
	s->chn->header();
	readynewfix(); 
	float *xyz=s->chn->gettransfer(s,to);
	float d=ensort(xyz);	 
	s->chn->transfer(xyz,s,to,1);
	if(xyz) delete [] xyz;xyz=0;
	return d;
}
float Segen::clashold(Res *s,char *force) 
{ 
  Atm *a; 
  float e;
  e=0;
  for(a=s->atm;a;a=a->next)
  {
     if(a->flag<0) continue;
     e+=clashold(a,force);
  } 
  return e;
}

float Segen::clash(Res *s,int to)
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

float Segen::clash(Res *s,int from,int to)
{
  Atm *a;
  float e;
  e=0;
  for(a=s->atm;a;a=a->next)
  {
     if(a->tatm->name[1]=='H') break;
     if(a->tatm->id<from||a->tatm->id>to) continue;
     if(a->flag<0) continue;
     e+=clash(a);
  }
  return e;
}

 
float Segen::clash(Res *s) //original
{ 
  Atm *a; 
  float e;
  e=0;
  for(a=s->atm;a;a=a->next)
  {
     if(a->flag<0) continue;
     e+=clash(a);
  } 
  //e+=s->clashnoring(4,0);
  return e;
}

float Segen::clash(Atm *aa0)
{
  Atm *aa1,*all[50];
  float dr,dn,e,d,elon,xr[3];
  int i,j,m,n,isp;
  float ta,tb,tc,tn,tt;
  int mr; 
  float z[6];

  if(smoothclash) {
        ta=TRES.engcoeff->coeff[0];
        tb=TRES.engcoeff->coeff[1];
        tc=TRES.engcoeff->coeff[2];
        tt=TRES.engcoeff->coeff[3];
        tn=TRES.engcoeff->coeff[4];
  }
  else {
        ta=TRES.smooth[TRES.smt*4+0];
        tb=TRES.smooth[TRES.smt*4+1];
        tc=TRES.smooth[TRES.smt*4+2];
        tn=TRES.smooth[TRES.smt*4+3];
        tt=0;
  }
 
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

  int nte=0;
  int nget=0;
  nget=lat->nget;
  if(strchr(dforce,'v')) nte=1;
  for(m=0;nte&&m<nget;m++)
  {   
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
     if(aa0->tatm->hbond==5) d-=aa0->tatm->eng->radius;
     if(aa1->tatm->hbond==5) d-=aa1->tatm->eng->radius;
     if(d<=0) d=1;
     if(aa0->res->name=='C'&&aa1->res->name=='C')
     {
        if(aa0->tatm->id==5&&aa1->tatm->id==5) {d=2.0;mr=2;}
        else if(aa0->tatm->id==5&&aa1->tatm->id==4) d=3.;
        else if(aa0->tatm->id==4&&aa1->tatm->id==5) d=3.;
     }
     if(aa0->tatm->name[1]=='H'&&aa1->tatm->name[1]=='H') d=1.0;
     for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j])/d;
     dr=TRES.distsqr(xr);
     //ds=dr;
     if(dr<0.1) dr=0.1; 
     else if(dr>4) continue;
     if(aa0->res->name=='C'&&aa1->res->name=='C')
     {
        if(aa0->tatm->id==5&&aa1->tatm->id==5)e+=-10.*exp(-(dr-1)*(dr-1)); 
     }
     elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     elon=ta*elon*exp(-dr*tb);
     dr=1/(tt+dr);
     
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

  nte=0;
  if(strchr(dforce,'h')) nte=1;
  nget=lat->nget;  
 
  for(m=0;nte&&m<nget;m++)
  {  
     aa1=lat->obtain[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
     j=aa0->tatm->hbond*aa1->tatm->hbond;
     float ef=0;
     if(j==2||j==3||j==6||j==9)  {
	if(aa0->tatm->id<=4&&aa1->tatm->id<=4)      ef+=aa0->ishbond(aa1,1.5,50,-0.5);
        else if(aa0->tatm->id<=4&&aa1->tatm->id>=4) ef+=aa0->ishbond(aa1,1.5,50,-0.35);
        else if(aa0->tatm->id>=4&&aa1->tatm->id<=4) ef+=aa0->ishbond(aa1,1.5,50,-0.35);
        else    			            ef+=aa0->ishbond(aa1,1.5,50,-0.25);
     }     
     e+=ef; 
  }

  nte=0;
  nget=lat->nget;
  if(strchr(dforce,'u')) nte=1;
  for(m=0;nte&&m<nget;m++)
  {
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
     if(aa0->tatm->hbond==5) d-=aa0->tatm->eng->radius;
     if(aa1->tatm->hbond==5) d-=aa1->tatm->eng->radius;
     if(d<=0) d=1.0;
     if(aa0->res->name=='C'&&aa1->res->name=='C')
     {
        if(aa0->tatm->id==5&&aa1->tatm->id==5) {d=2.0;mr=2;}
        else if(aa0->tatm->id==5&&aa1->tatm->id==4) d=3.;
        else if(aa0->tatm->id==4&&aa1->tatm->id==5) d=3.;
     }
     if(aa0->tatm->name[1]=='H'&&aa1->tatm->name[1]=='H') d=1.0;
     for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j])/d;
     dr=TRES.distsqr(xr);
     //ds=dr;
     if(dr<0.1) dr=0.1; 
     else if(dr>4) continue;
     if(aa0->res->name=='C'&&aa1->res->name=='C')
     {
        if(aa0->tatm->id==5&&aa1->tatm->id==5)e+=-10.*exp(-(dr-1)*(dr-1)); 
     }
     elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     elon=ta*elon*exp(-dr*tb);
     dr=1/(tt+dr);
     
     dn=pow(dr,tn);  
     dn=elon*dn*(dn-tc);   
      
     dn=dn/10;
     if(aa0->tatm->id>=4&&aa1->tatm->id>=4&&dn<0) {
	if(aa0->tatm->name[1]=='H'&&aa0->tatm->bond[0]->id<4) {
		dn=dn;
	}	  
	else if(aa1->tatm->name[1]=='H'&&aa1->tatm->bond[0]->id<4) {
		dn=dn;
	}	  
	else if(strchr("AVLIFCMPW",aa0->res->name)&&strchr("AVLIFCMPW",aa1->res->name)) {
		dn=dn*10;
	} 
     }
     else if(dn>0) {
	dn=dn;
     }
     e+=dn;
  }

  Charge *ch0=disc->findcharge(aa0);
  Dipole *da0=disc->range[aa0->id0]; 
  
  nte=0;
  if(strchr(dforce,'D')) nte=1;

  float f0[3],f1[3],f2[3];
  float ff[2][3];

  if(da0) {
	for(i=0;i<3;i++) f0[i]=(da0->crg[0]->atm->xyz[i]+da0->crg[1]->atm->xyz[i])/2;
  }
  Dipole *dp0;
  for(m=0;da0&&nte&&m<nget;m++)  //dipole and dipole interaction
  {
     aa1=lat->obtain[m]->atm;     
     dp0=disc->range[aa1->id0];
     if(dp0==0) continue;  
     
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
 
     if(aa0->res==aa1->res) continue;
      
     for(i=0;i<3;i++) f1[i]=(dp0->crg[0]->atm->xyz[i]+dp0->crg[1]->atm->xyz[i])/2;      
     for(i=0;i<3;i++) f2[i]=f1[i]-f0[i];
     d=TRES.distance(f2);

     if(d<0.1) continue;
     else if(d<3.5) {	 
	for(j=0;j<3;j++)f2[j]=f2[j]*(3.5-d)/d;
	for(i=0;i<2;i++)  		
	for(j=0;j<3;j++) {
		ff[i][j]=dp0->crg[i]->atm->xyz[j]+f2[j];	
	}		 
     }
     else {
	for(i=0;i<2;i++) 
	for(j=0;j<3;j++) {
		ff[i][j]=dp0->crg[i]->atm->xyz[j];		
	}
     }

     int i0,j0;
      
     for(i0=0;i0<2;i0++) {
	Charge *c=da0->crg[i0];
	if(c==0) continue;
	for(j0=0;j0<2;j0++) {
		Charge *t=dp0->crg[j0];
		if(t==0) continue;
		for(i=0;i<3;i++) z[i]=c->atm->xyz[i]-ff[j0][i];
		float d=TRES.distance(z); 		 			
		float f=332*c->crg*t->crg/d/disc->eps;
		e+=f;
	}
     }     
  }

  nte=0;
  if(strchr(dforce,'d')) nte=1; //new energy function of dipole,charge interaction
  da0=disc->range[aa0->id0];
  int kep=0;while(disc->charge&&disc->charge[kep])kep++;
  for(i=0;da0&&i<3;i++) f0[i]=(da0->crg[0]->atm->xyz[i]+da0->crg[1]->atm->xyz[i])/2;
  for(m=0;da0&&nte&&m<kep;m++)  //dipole and charge interaction
  {
     aa1=disc->charge[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa0->res==aa1->res) continue;     
      
     for(i=0;i<3;i++) f1[i]=disc->charge[m]->atm->xyz[i];      
     for(i=0;i<3;i++) f2[i]=f1[i]-f0[i];
     d=TRES.distance(f2);

     if(d<0.1) continue;
     else if(d<3.5) {
	 for(i=0;i<3;i++) f2[i]=f2[i]*(3.5-d)/d; 
	 for(i=0;i<3;i++) f1[i]=disc->charge[m]->atm->xyz[i]+f2[i];	 
     }
     else {
	 for(i=0;i<3;i++) f1[i]=disc->charge[m]->atm->xyz[i];
     }
     

     int j0;
     for(j0=0;j0<2;j0++) {
	Charge *t=da0->crg[j0];
	if(t==0) continue;
	for(i=0;i<3;i++) z[i]=f1[i]-t->atm->xyz[i];
	float d=TRES.distance(z); 	 
	float f=332*disc->charge[m]->crg*t->crg/d/disc->eps;
	e+=f;	 
     }     
  }

  nte=0;
  if(strchr(dforce,'c')) nte=1; //new energy function of dipole,charge interaction
  
  if(nte) {	 
   	kep=0;while(disc->charge&&disc->charge[kep])kep++;
  }
  
  for(m=0;ch0&&nte&&m<kep;m++)  //charge and charge interaction
  {
     aa1=disc->charge[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
     if(aa0->res==aa1->res) continue;

     int i;
     for(i=0;i<3;i++) z[i]=ch0->atm->xyz[i]-disc->charge[m]->atm->xyz[i];
     float d=TRES.distance(z);     
     if(d<3.5) d=3.5;

     float f=332*ch0->crg*disc->charge[m]->crg/d/disc->eps;
     e+=f;	      	 
  }

  
  for(i=0;i<n;i++)lat->puton(all[i]);

  return e;  
}

float Segen::clashold(Atm *aa0,char *force)
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

float Segen::clashfix(Atm *aa0)
{
  Atm *aa1,*all[50];
  float dr,dn,e,d,elon,xr[3];
  int i,j,m,n,isp;
  float ta,tb,tc,tn,tt,ds,de;
  int mr; 
  
  if(smoothclash) {
        ta=TRES.engcoeff->coeff[0];
        tb=TRES.engcoeff->coeff[1];
        tc=TRES.engcoeff->coeff[2];
        tt=TRES.engcoeff->coeff[3];
        tn=TRES.engcoeff->coeff[4];
  }
  else {
        ta=TRES.smooth[TRES.smt*4+0];
        tb=TRES.smooth[TRES.smt*4+1];
        tc=TRES.smooth[TRES.smt*4+2];
        tn=TRES.smooth[TRES.smt*4+3];
        tt=0;
  }
 
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
     if(aa0->tatm->hbond==5) d-=aa0->tatm->eng->radius;
     if(aa1->tatm->hbond==5) d-=aa1->tatm->eng->radius;
     if(d<=0) d=1.0;
     if(aa0->res->name=='C'&&aa1->res->name=='C')
     {
        if(aa0->tatm->id==5&&aa1->tatm->id==5) {d=2.0;mr=2;}
        else if(aa0->tatm->id==5&&aa1->tatm->id==4) d=3.;
        else if(aa0->tatm->id==4&&aa1->tatm->id==5) d=3.;
     }
     if(aa0->tatm->name[1]=='H'&&aa1->tatm->name[1]=='H') d=1.0;
     for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j])/d;
     dr=TRES.distsqr(xr);
     ds=dr;
     if(dr<0.1) dr=0.1; 
     else if(dr>4) continue;
     if(aa0->res->name=='C'&&aa1->res->name=='C')
     {
        if(aa0->tatm->id==5&&aa1->tatm->id==5)e+=-10.*exp(-(dr-1)*(dr-1)); 
     }
     elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     elon=ta*elon*exp(-dr*tb);
     dr=1/(tt+dr);
     
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
    
     if(aa0->tatm->eng->charge!=0&&aa1->tatm->eng->charge!=0&&0)
     {
            de=1/(0.5+ds*d*d); 
            dn=aa0->tatm->eng->charge*aa1->tatm->eng->charge*de*332.3/20;	     
            e+=dn;
     }
  }
  
  if(strchr(dforce,'E')||strchr(dforce,'e')||strchr(dforce,'H')||\
     strchr(dforce,'h')||strchr(dforce,'i')) 
  {
    lat->resonly();
    for(m=0;m<lat->nget;m++)
    {
      if(aa0->tatm->eng->charge==0) break;
      if(strchr(dforce,'E')==0&&strchr(dforce,'e')==0&&strchr(dforce,'i')==0) 
      {
        if(aa0->tatm->name[1]!='O'&&aa0->tatm->name[1]!='N') break;
      }
      for(aa1=lat->obtain[m]->atm->res->atm;aa1;aa1=aa1->next)
      {
         if(aa1->tatm->eng->charge==0) continue;
         if(strchr(dforce,'E')==0&&strchr(dforce,'e')==0||strchr(dforce,'i')==0) 
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
         if(strchr(dforce,'E')||strchr(dforce,'e')||strchr(dforce,'i'))
         {
           e+=aa0->tatm->eng->charge*aa1->tatm->eng->charge/d/d*332.3/10;
         }
         if(strchr(dforce,'h')||strchr(dforce,'H'))
         {
           if(aa0->tatm->id<=4&&aa1->tatm->id<=4)      e+=aa0->ishbond(aa1,1.5,50,-0.5);
           else if(aa0->tatm->id<=4&&aa1->tatm->id>=4) e+=aa0->ishbond(aa1,1.5,50,-0.35);
           else if(aa0->tatm->id>=4&&aa1->tatm->id<=4) e+=aa0->ishbond(aa1,1.5,50,-0.35);
           else                                        e+=aa0->ishbond(aa1,1.5,50,-0.25);
         }
      }
    }
  }
 
  if(strchr(dforce,'e')||strchr(dforce,'h')||strchr(dforce,'i'))
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
      if(strchr(dforce,'e')==0) break; 
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
float*  Segen::scap(float **xyz_all)
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

 for(i=0;i<all;i++)
 {
   segment->chn->transfer(xyz_all[i],1);
   setend(1);
   d=scap();
   ent[i]=d;
   segment->chn->transfer(xyz_all[i],0);
   chn->transfer(xyz,r0,end,1);
   //cerr<<i<<" "<<"side-chain predicted energy: "<<rmsd(2)<<"  "<<d<<endl;
 }
 
 chn->transfer(xyz,r0,end,1);
 
 delete [] xyz;xyz=0;

 return ent;
}
float*  Segen::scapfix(float **xyz_all)
{
 Res *r,*r0;
 int nseq;
 float *xyz;
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
 mutate->setrotamerorder(0);
 for(i=0;i<all;i++)
 {
   Scap scprd;
   scprd.singletorsion=1;
   scprd.colonyline=1;
   scprd.colony=2;
   scprd.ncolony=1;
   scprd.nncolony=1;
   scprd.nummore=100;
   scprd.bmax=2;
   scprd.tormax=2;
   scprd.ring=1;	
   scprd.pdb=pdb;
   chn->transfer(xyz,r0,end,1);
   chn->transfer(xyz_all[i],r0,end,1);
   pdb->setflgr(-99999);
   pdb->chn->setflgr(r0,end,10000);
   pdb->setflgr('G',-99999);
   pdb->setflgr('A',-99999);
   pdb->setflgr('P',-99999);
   ent[i]=scprd.scpred(1);
   for(Chn *c=pdb->chn;c;c=c->next) c->setfreenextonmore();
   chn->transfer(xyz_all[i],r0,end,0);
   scprd.pdb=0;
 }
 pdb->setflgr(0);
 chn->transfer(xyz,r0,end,1);
 mutate->setrotamerorder(1);
 delete [] xyz;xyz=0;

 return ent;
}

void  Segen::hooksidechain() {

	Chn *chn=(*pdb)[cid];

	Res *r0=(*chn)[start];
	
	Res *r;
	
	for(r=r0;r;r=r->next) {
		if(r->id0>end) break;
		if(r->name=='G') continue;
		r->addsidechain();
	}
	pdb->configure();
	
	for(r=segment->chn->res;r;r=r->next) {
		if(r->id0>end) break;
		if(r->name=='G') continue;
		r->addsidechain();
	}
	segment->configure();
 
  	for(r=segment->chn->res;r;r=r->next) {
        	r->id0+=r0->id0;
        	r->id+=r0->id;
  	}

} 



float Segen::scap()
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
    //else        r1=TRES.rotamer->next->chn->isres(r->tres->name,0);
    //else        r1=TRES.findrotamer("sidechain")->chn->isres(r->tres->name,0);
    else        r1=TRES.findrotamername("side_small_rotamer")->chn->isres(r->tres->name,0);
    rot.hook(r1,r,4);
  }

  //preparing lattice
  lat->putoff();
  lat->flag=1;
  lat->grdsiz=2.;
  lat->radall=15.;
  lat->ready(pdb);
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

  crtmore();

  //puton lattice
  lat->putoff();
  for(r=r0;r;r=r->next)
  {
    if(r->id0>end) break;
    if(strchr("PGA",r->tres->name)) continue;   
    lat->puton(r,5,100);
  }

  //side-chain prediction

  d=scpred();

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

void Segen::crtmore()
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
  //for(r=TRES.rotamer->next->chn->res;r;r=r->next)
  for(r=TRES.findrotamer("sidechain")->chn->res;r;r=r->next)
  //for(r=TRES.findrotamername("side_small_rotamer")->chn->res;r;r=r->next)
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
    //else        r1=TRES.rotamer->next->chn->isres(r->tres->name,0);
    else        r1=TRES.findrotamer("sidechain")->chn->isres(r->tres->name,0);
    //else        r1=TRES.findrotamername("side_small_rotamer")->chn->isres(r->tres->name,0);
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
         e+=clash(a2);
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

  delete [] rm;rm=0;
}


float Segen::scpred()
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

    dmin=r->flag+clash(r,5,100);
    d=r->flag;
    r->transfer(xyz,0);
    n=0;
    for(r1=r->more;r1;r1=r1->more)
    {   
       n++;
       if(n>10) break; 
       r->transfer(r1);
       e=r1->flag+clash(r,5,100);
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

float Segen::minimize(float *xyz0,int flg)
{
  float **xyz_all;
  float *ent,d;
  xyz_all=new float*[2];
  xyz_all[0]=xyz0;
  xyz_all[1]=0;
  ent=minimize(xyz_all,1);
  d=ent[0];
  if(ent) delete [] ent; ent=0;
  delete [] xyz_all;  xyz_all=0;
  return d;
}

float Segen::minimizefix(float *xyz0,int flg)
{
  float **xyz_all;
  float *ent,d;
  xyz_all=new float*[2];
  xyz_all[0]=xyz0;
  xyz_all[1]=0;
  ent=minimizefix(xyz_all,1);
  d=ent[0];
  if(ent) delete [] ent; ent=0;
  delete [] xyz_all;  xyz_all=0;
  return d;
}


float* Segen::minimize(float **xyz_all,int flg)
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
  
  for(i=0;i<all;i++)
  {      
    segment->chn->transfer(xyz_all[i],1);                          
    ent[i]=segmin(flg); 
    segment->chn->transfer(xyz_all[i],0);          
  } 
  //rusage tt0;
  //getrusage(RUSAGE_SELF,&tt0);

  //cerr<<"\nminimize loos:"<<all<<" with time..."<<tt0.ru_stime.tv_usec-tt.ru_stime.tv_usec<<endl;
  cc.sort(ent,all,order);
  for(i=0;i<all;i++)
  {
    pot[i]=xyz_all[order[i]];
  }

  for(i=0;i<all;i++) xyz_all[i]=pot[i];

  delete [] pot;delete [] order;pot=0;order=0;

  chn->transfer(xyz,r0,end,1);
  delete [] xyz; xyz=0;

  return ent;
}




float* Segen::minimizefix(float **xyz_all,int flg)
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
  while(xyz_all&&xyz_all[all])all++;
  if(all==0) return 0;

  ent=new float[all];
  order=new int[all];
  pot=new float*[all];

  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;     
  
  xyz=new float[nseq];  

  chn->transfer(xyz,r0,end,0);

  //rusage tt;
  //getrusage(RUSAGE_SELF,&tt);
  
  for(i=0;i<all;i++)
  {      
    segment->chn->transfer(xyz_all[i],1);  
    chn->transfer(xyz_all[i],r0,end,1);  
    if(bound) {if(TRES.logg)cerr<<"....begin...bound....bound..."<<endl<<endl;if(TRES.logg)bound->printoutdistance();}
    if(TRES.logg>3)chn->write("s1");
    chn->transfer(xyz,r0,end,1);            
    ent[i]=boundmin(flg); 
    segment->chn->transfer(xyz_all[i],0);  
    chn->transfer(xyz_all[i],r0,end,1);  
    if(bound) {if(TRES.logg)cerr<<"....end...bound....bound..."<<endl<<endl;if(TRES.logg)bound->printoutdistance();}
    if(TRES.logg>3)chn->write("s2");
    chn->transfer(xyz,r0,end,1);         
  } 
  //rusage tt0;
  //getrusage(RUSAGE_SELF,&tt0);

  //cerr<<"\nminimize loos:"<<all<<" with time..."<<tt0.ru_stime.tv_usec-tt.ru_stime.tv_usec<<endl;
  cc.sort(ent,all,order);
  for(i=0;i<all;i++)
  {
    pot[i]=xyz_all[order[i]];
  }

  for(i=0;i<all;i++) xyz_all[i]=pot[i];

  delete [] pot;delete [] order;pot=0;order=0;

  chn->transfer(xyz,r0,end,1);
  delete [] xyz; xyz=0;

  return ent;
}

float* Segen::boundminimizefix(float **xyz_all,int flg)
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
  
  for(i=0;i<all;i++)
  {      
    segment->chn->transfer(xyz_all[i],1);                          
    ent[i]=boundmin(flg); 
    segment->chn->transfer(xyz_all[i],0);          
  } 
  //rusage tt0;
  //getrusage(RUSAGE_SELF,&tt0);

  //cerr<<"\nminimize loos:"<<all<<" with time..."<<tt0.ru_stime.tv_usec-tt.ru_stime.tv_usec<<endl;
  cc.sort(ent,all,order);
  for(i=0;i<all;i++)
  {
    pot[i]=xyz_all[order[i]];
  }

  for(i=0;i<all;i++) xyz_all[i]=pot[i];

  delete [] pot;delete [] order;pot=0;order=0;

  chn->transfer(xyz,r0,end,1);
  delete [] xyz; xyz=0;

  return ent;
}

float* Segen::minimizesidechain(float **xyz_all,int flg)
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
  
  for(i=0;i<all;i++)
  {      
    segment->chn->transfer(xyz_all[i],1);   
    if(TRES.logg>3)segment->chn->write("seg1");                       
    ent[i]=sidechainmin(flg); 
    if(TRES.logg>3)segment->chn->write("seg2");
    segment->chn->transfer(xyz_all[i],0);          
  } 
  //rusage tt0;
  //getrusage(RUSAGE_SELF,&tt0);

  //cerr<<"\nminimize loos:"<<all<<" with time..."<<tt0.ru_stime.tv_usec-tt.ru_stime.tv_usec<<endl;
  cc.sort(ent,all,order);
  for(i=0;i<all;i++)
  {
    pot[i]=xyz_all[order[i]];
  }

  for(i=0;i<all;i++) xyz_all[i]=pot[i];

  delete [] pot;delete [] order;pot=0;order=0;

  chn->transfer(xyz,r0,end,1);
  delete [] xyz; xyz=0;

  return ent;
}



float** Segen::zipper(float dis,int flg)
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
  //srandom(seed);
  nseq=0; 
  r0=(*chn)[start];
  for(r=r0;r&&r->id0<=end;r=r->next) 
  nseq+=r->tres->number*3;     
  
  xyz=new float[nseq];  
  xyz0=new float[nseq];
  int arbt0=arbt;
  if(chiangle) {
	arbt0=max(arbt0,chiangle->getsize());
  }
  xyz_all=new float*[arbt0+10];
  for(i=0;i<arbt0+10;i++) xyz_all[i]=0;

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
    //srandom(ranseed+n);
    randmz();
    
    //segment->chn->transfer(xyz0,0);
    /*
    if(checkomega()==0) {
	n--;continue; 
    }
    */
    //chn->transfer(xyz0,r0,end,1);
    //chn->dihedral(stdout);
    if(chiangle) {
	if(TRES.logg)cerr<<"rotate the dihedral to the specified in the database..."<<n-1<<"  "<<databaseonly<<endl;
	//checkomega(stdout);
	if(randmzChiangle(n-1)==0&&databaseonly==1) break;
	//checkomega(stdout);
        //FILE *fp=0;
	//segment->chn->dihedral(fp);
    }
    else if(secd=='h'&&end-start+1>3) randmzHelix(15);
    else if(secd=='e'&&end-start+1>3) randmzSheet(30);
    else if(chiangle==0&&databaseonly==1) {
	break;
    }
    //chn->dihedral(stdout);
    //test..start
    segment->chn->transfer(xyz0,0);
    chn->transfer(xyz0,r0,end,1);
    //segment->chn->write("se1");
    //chn->dihedral("s1");
    if(TRES.logg>3)chn->write("te1");
    //end
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
    //test
    if(TRES.logg>3)segment->chn->dihedral(stdout);    
    if(TRES.logg>3)pdb->chn->transfer(xyz_all[i],r0,end,1);
    if(TRES.logg>3) pdb->write("ss1");
    pdb->chn->transfer(xyz,r0,end,1);
    //end
    //char sf[100];
    //sprintf(sf,"_%i",i);
    //pdb->write(sf);
    i++;
    if(i>=arbt&&databaseonly==0) break;
  }
  if(TRES.logg)cerr<<"loop number :"<<n<<" "<<i<<endl;
  chn->transfer(xyz,r0,end,1);
  if(xyz) {delete [] xyz;xyz=0;}
  if(xyz0) {delete [] xyz0;xyz0=0;}
  return xyz_all;
}

float** Segen::cross(float **xyz_all,int cot,float cut)
{
  //loop cross
  Chn *cht;
  Res *r,*r0;  
  int all,i,nseq,j,m,m3,m4;
  float e,d,*xyz,**xyz_out;     
  Segen *seg1,*seg2;    
  Qsort cc;
  int *order,*used;
  float *dist;
  Chn *chn=(*pdb)[cid];

  if(end-start+1<=3) return xyz_all;
  
 //space
   
  seg1=new Segen(this); seg2=new Segen(this);
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
  xyz_out=new float*[all*(cot+10)];
  for(i=0;i<all*(cot+10);i++) xyz_out[i]=0;

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
  
  delete [] xyz_all;xyz_all=0;
  seg1->pdb=0;seg2->pdb=0;
  //seg1->rotamer=0;seg2->rotamer=0;
  //seg1->smooth=0; seg2->smooth=0;
  seg1->next=0;seg2->next=0;seg1->lat=0;seg2->lat=0;
  delete seg1;seg1=0;
  delete seg2;seg2=0;
  delete [] order;order=0;
  delete [] dist;dist=0;
  delete [] used;used=0;
  delete [] xyz;xyz=0;
  return xyz_out;
}

float Segen::ensort(float *xyz)
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
  ent=ensort(xyz_all);
  d=ent[0];
  if(ent) delete [] ent; ent=0;
  delete [] xyz_all;xyz_all=0;
  return d;
}

float* Segen::ensortold(float **xyz_all,char *flg)
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
  lat->grdsiz=2.;
  lat->radall=15.;
  lat->ready(pdb);
  for(Chn *cchn=pdb->chn;cchn;cchn=cchn->next)
  for(r=cchn->res;r;r=r->next)
  {
    if(r->id0>=start&&r->id0<=end&&cchn==chn) continue;
    lat->puton(r);
  }

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
       entv[i]+=clashold(r,flgv);
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
       ente[i]+=clashold(r,flgv);
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
           ente[i]+=clashold(r,flgs);
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

  if(strchr(flg,'w'))
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
  return ent;
}

float* Segen::ensort(float **xyz_all)
{
  Res *r,*r0;
  Atm *a;
  int i,j,nseq,all;
  float e,*xyz;
  float **pot,*ent,*entv;
  int  *order,*give;
  Qsort cc;
  
  if(xyz_all==0) return 0;
  Chn *chn=(*pdb)[cid];

  xyz=0;pot=0;order=0;ent=0;
  //space     
  all=0;
  while(xyz_all[all])all++;
  if(all==0) return 0;
  order =new int[all];
  r0=(*chn)[start];

  //double computing
 
  nseq=0; 
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;         
  xyz=new float[nseq];

  pot=new float*[all];
  for(i=0;i<all;i++) pot[i]=0; 
  
  ent=new float[all];
  entv=new float[all]; 
  for(i=0;i<all;i++) {ent[i]=0;entv[i]=0;}
  chn->transfer(xyz,r0,end,0);
     
  
  //test energy
  /*
  Disc disc0;
  disc0.pdb=pdb;
  disc0.ready();
  disc0.cutoff=8;
  disc0.setrange();
  disc0.setup();
  disc0.setdipoleascharge();
  disc0.setrange();
  strcpy(disc0.force,"udD");
  */
  //test end 

  //---------------------  
  i=0;
  for(r=chn->res;r;r=r->next)
  for(a=r->atm;a;a=a->next)i++;

  give=new int[i];
  //gives=new int[i];

  lat->putoff();
  lat->flag=1;
  lat->grdsiz=2.;
  lat->radall=15.;
  lat->ready(pdb);
  for(Chn *cchn=pdb->chn;cchn;cchn=cchn->next)
  for(r=cchn->res;r;r=r->next)
  {
    if(r->id0>=start&&r->id0<=end&&cchn==chn) continue;
    lat->puton(r);
  }
  
  //double computing
  r0=(*chn)[start];
  pdb->setflg(-1);
  chn->setflg(r0,end,1);
 
  if(r0->last) {
        r0->atm->flag=-1;
        r0->atm->next->flag=-1;
        a=(*r0)[" HN "];
        if(a)a->flag=-1;
  }
  r=(*chn)[end];
  if(r->next) {
        a=r->atm->next;
        a->flag=-1;
        a->next->flag=-1;
        a->next->next->flag=-1;
  }
   
  chn->setflg(give,chn->res,100000,0);
  //----------------------

 
  disc->cutoff=cutoff;
  disc->setrange(); 
  strcpy(disc->force,dforce);   
  chn->setflg(give,chn->res,100000,1);
  int nui=0;
  if(strchr(dforce,'D')||strchr(dforce,'d')||strchr(dforce,'c')) nui=1;
  for(i=0;i<all;i++) {	
  	chn->transfer(xyz_all[i],r0,end,1);  
	lat->putoff(r0,end);
        lat->puton(r0,end);	
	if(nui) { 	 
		disc->setup(r0,end,1);
		disc->setdipoleascharge(r0,end);
	}
	if(TRES.logg>3) r0->chn->write("b");	 
	entv[i]=disc->clash(lat,r0,end); 		        
  }
  
  //calculate surface
  if(strchr(dforce,'s'))
  {  
    float area,area0,area1,d;
    setoff(xyz_all,6.);
    chn->setflg(give,chn->res,1000000,0);
    for(i=0;i<all;i++)
    {
      chn->setflg(give,chn->res,1000000,1);
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
      if(strstr(dforce,"s1"))  d=0.025*(area0/4+area/2+area1);
      else                     d=0.05*(area0/4+area/2+area1);
      if(strchr(dforce,'s'))   entv[i]+=d;
      chn->transfer(xyz,r0,end,1);
      segment->chn->transfer(xyz_all[i],1);
     }
     chn->transfer(xyz,r0,end,1);        
  }
 
  pdb->setflg(0);
  
  for(i=0;i<all;i++) ent[i]=entv[i];
  
  
  if(strchr(dforce,'w'))
  {
    for(i=0;i<all;i++) ent[i]=0;
    float et=8.31*0.3;     
    float org=0,ne=0;
    float et0=et;
    et=1/et;
    for(i=0;i<all;i++)
    {
      chn->transfer(xyz_all[i],r0,end,1);
      org=entv[i];
      for(j=0;j<all;j++)
      {
        ne=entv[j];	 
        segment->chn->transfer(xyz_all[j],1);
        e=rmsd(2);
        //d=exp(-e*e*e/(end-start)/6);
        //d0+=d;
        //ent[i]+=-d*exp(-(entv[j])/(et));
	ent[i]+=exp(-ne*et-e*e*e/(end-start)/6);
      }
      if(ent[i]==0) {
	ent[i]=org;
      }
      else ent[i]=-et0*log(ent[i]);
    }  
  }
  
  
  if(strchr(dforce,'W'))
  {     
    float et=8.31*0.3;     
    float org=0,ne=0;
    float et0=et;
    et=1/et;
    float xx=-log(0.5);
    float avgrmsd=calcavgrmsd(xyz_all,20);
    if(TRES.logg)cerr<<"the average rmsd: "<<avgrmsd<<endl;
    if(avgrmsd==0) avgrmsd=1;
    avgrmsd=1/avgrmsd/avgrmsd; 
    for(i=0;i<all;i++) ent[i]=0;
    float emax=entv[0];
    for(i=0;i<all;i++) {
	if(entv[i]<emax) emax=entv[i];
    }
    
    for(i=0;i<all;i++)
    { 
      chn->transfer(xyz_all[i],r0,end,1);
      org=entv[i];
      int did=1;
      if(org-emax>20) did=0;
      for(j=0;did&&j<all;j++)
      {
        ne=entv[j];
	if(ne-org>20) continue;
        segment->chn->transfer(xyz_all[j],1);
        e=rmsd(2);       
	ent[i]+=exp(-ne*et-e*e*avgrmsd*xx);
      }
      if(ent[i]==0) {
	ent[i]=org;
      }
      else ent[i]=-et0*log(ent[i]);
    }  
    
  }

  if(strchr(dforce,'&'))
  {     
    float et=8.31*0.3;     
    float org=0,ne=0;
    float et0=et;
    et=1/et;
    float xx=-log(0.5);
    float avgrmsd=calcavgrmsd(xyz_all,20);
    if(TRES.logg)cerr<<"the average rmsd: "<<avgrmsd<<endl;
    if(avgrmsd==0) avgrmsd=1;
    avgrmsd=1/avgrmsd/avgrmsd; 
    for(i=0;i<all;i++) ent[i]=0;
    entv=calcenergycoeff(entv,all);
    float emax=entv[0];
    for(i=0;i<all;i++) {
	if(entv[i]<emax) emax=entv[i];
    }
    
    for(i=0;i<all;i++)
    { 
      chn->transfer(xyz_all[i],r0,end,1);
      org=entv[i];
      int did=1;
      if(org-emax>20) did=0;
      for(j=0;did&&j<all;j++)
      {
        ne=entv[j];
	if(ne-org>20) continue;
        segment->chn->transfer(xyz_all[j],1);
        e=rmsd(2);       
	ent[i]+=exp(-ne*et-e*e*avgrmsd*xx);
      }
      if(ent[i]==0) {
	ent[i]=org;
      }
      else ent[i]=-et0*log(ent[i]);
    }  
    
  }
  //first part of segment;

  for(i=0;i<all;i++)order[i]=i;

  if(!strstr(dforce,"?"))cc.sort(ent,all,order);

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
  delete [] order; order=0; 
  delete [] xyz; xyz=0;  
  delete [] pot;  pot=0;    
  delete [] entv; entv=0; 
  delete [] give; give=0; 
  
  return ent;
}


float Segen::calcavgrmsd(float **xyz_all,int num) {

    float *xyz;
    Res *r;
    int nseq=0;	
    int all=0;
    Chn *chn=(*pdb)[cid];
    while(xyz_all[all]) all++;

    Res *r0=(*chn)[start];    
    for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;

    xyz=new float[nseq];
    chn->transfer(xyz,r0,end,0);

    int nn=all/num;
    if(nn==0) nn=1;
    int total=0;
    float arm=0;
    int i,j;
    int *ip=new int[all+100];
    for(i=0;i<all;i++) ip[i]=i;
    for(i=0;i<all;i++) {
             int j=random();
             j=j%all;
             int k=ip[i];
             ip[i]=ip[j];
             ip[j]=k;
    }
    int ii,jj;	
    for(ii=0;ii<all;ii+=nn)
    {
      i=ip[ii];
      chn->transfer(xyz_all[i],r0,end,1);       
      for(jj=ii+1;jj<all;jj+=nn)
      {         
	j=ip[jj];
	if(i==j) continue;
        segment->chn->transfer(xyz_all[j],1);
        arm+=rmsd(2);       
	total++; 
      }
    }
    delete [] ip;ip=0;
    chn->transfer(xyz,r0,end,1);
    delete [] xyz;xyz=0;
    if(total==0) return 0;
    return arm/total;
}

float*Segen::calcenergycoeff(float *ent,int num) {

    	if(num<2) return ent;
 
    	int kep=num;
    	int *order=new int[kep];	
    	float *temp=new float[kep];

	int i;
	for(i=0;i<kep;i++) temp[i]=ent[i];
	Qsort cc;
	cc.sort(temp,kep,order);

	for(i=0;i<num;i++) {
    		int j=order[i];
		ent[j]=-8.31*0.3*log(num-i);
	}

	if(order) delete [] order;order=0;
	if(temp) delete [] temp;temp=0;
	return ent;
}


float Segen::noclose(float **xyz_all,float *xyz0) 
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
  delete [] xyz;xyz=0;
  return close;
}


int Segen::noclose(float **xyz_all,float close)
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


int Segen::nolarge(float **xyz_all,float close)
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
      if(d>close) kepp[j]=0;
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


float **Segen::around(float **xyz_all)
{
  float *ent;
  float **xyz_all0,**xyz_all00;
  int kep0,i,j,k,kep;
  Qsort cc;
  float **pot;
  int *order;

  if(TRES.logg)cerr<<"start expanding around..."<<endl;

  kep=0;while(xyz_all[kep])kep++;
  xyz_all0=new float*[kep*12];
  for(i=0;i<kep*12;i++) xyz_all0[i]=0;
  j=0;
  //srandom(10000);
  for(i=0;i<kep;i++)
  {
    xyz_all0[j++]=xyz_all[i];
    k=random()%100000;
    xyz_all00=expand(xyz_all[i],k,10,30);
    kep0=0;while(xyz_all00[kep0])kep0++;
    for(k=0;k<kep0;k++) xyz_all0[j++]=xyz_all00[k]; 
    delete [] xyz_all00;xyz_all00=0;
  }
  delete [] xyz_all;xyz_all=0;
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
    if(TRES.logg)cerr<<i<<" after expanding energy of v: "<<rmsd(2)<<" "<<ent[i]<<endl;
  }
  if(TRES.logg)cerr<<" the total number after refill: "<<kep<<" "<<d<<endl;
  delete [] order;delete [] pot; order=0;pot=0;
  return xyz_all;
}

float *Segen::expand(float **xyz_all,int flg)
{
  Res *r0,*r;
  int nseq;
  int i,kep;
  float *ent,e,*xyz;
  int kk;
  Qsort cc;
  int *order;
  float **pot;
 

  Chn *chn=(*pdb)[cid];
  
  strcpy(dforce,"v");

  r0=(*chn)[start];nseq=0;
  for(r=r0;r;r=r->next)if(r->id0<=end) nseq+=r->tres->number*3;
  xyz=new float[nseq];
  chn->transfer(xyz,r0,end,0);
 
  kep=0;
  while(xyz_all[kep])kep++;
  ent=ensort(xyz_all);

  kk=part0; 

  float rms,rms0;
  rms=0;
  for(i=0;i<kep;i++)
  {
    segment->chn->transfer(xyz_all[i],1);
    rms0=rmsd(2);
    rms+=rms0/kep;
    if(i<kk)
    {
      if(TRES.logg)cerr<<i<<" energy of "<<dforce<<": "<<rms0<<" "<<ent[i]<<" ";
      ent[i]=expand(xyz_all[i],flg);
      segment->chn->transfer(xyz_all[i],1);
      if(TRES.logg)cerr<<rmsd(2)<<" "<<ent[i]<<endl;
    }
    else    
    {
      if(TRES.logg)cerr<<i<<" energy of "<<dforce<<": "<<rms0<<" "<<ent[i]<<endl;
      delete [] xyz_all[i];xyz_all[i]=0;
    }
  }

  chn->transfer(xyz,r0,end,1);
  kep=0;
  while(xyz_all[kep])kep++;
  
  if(TRES.logg)cerr<<"the total number after "<<dforce<<": "<<kep<<"  "<<rms<<endl;
  
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
    if(TRES.logg)cerr<<i<<" energy of "<<dforce<<": "<<e<<" "<<ent[i]<<endl;
  }
  if(TRES.logg)cerr<<"the total number after "<<dforce<<": "<<kep<<"  "<<rms<<endl;
  delete [] pot;delete [] order;delete [] xyz;pot=0;order=0;xyz=0;
  return ent;
}

float Segen::expand(float *xyz0,int flg)
{
  Res *r0,*r;
  float d,*ent,**xyz_all;
  int kep,i,n;
  int nseq;

  Chn *chn=(*pdb)[cid];

  nseq=0;
  r0=(*chn)[start];
  for(r=r0;r;r=r->next)if(r->id0<=end) nseq+=r->tres->number*3;
 
  //srandom(100);
  //d=minimize(xyz0,1);
  d=minimizefix(xyz0,1);
  if(d<0||flg==0||fapr==1) return d;

  for(i=0;i<1;i++)
  {
    n=random()%100000;
    xyz_all=expand(xyz0,n,5,10);
    strcpy(dforce,"v");
    ent=ensort(xyz_all);if(ent) delete [] ent; ent=0;
    kep=0;
    while(xyz_all[kep])kep++;
    for(n=1;n<kep;n++) {delete [] xyz_all[n];xyz_all[n]=0;}
    ent=minimizefix(xyz_all,1);
    if(ent[0]<d) 
    {
      d=ent[0];
      TRES.copy(xyz_all[0],xyz0,nseq);
    } 
    delete [] xyz_all[0];
    delete [] xyz_all;xyz_all=0;
    if(ent) delete [] ent; ent=0;
    if(d<0) break;
  }
  return d;
}

float** Segen::expand(float *xyz0,int seed,int kep,int ang)
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
    randmz(ang);
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

float *Segen::peel(float **xyz_all,int got)
{
//recursively peel off
  float *ent;
  int i,kep;
 
  for(i=0;i<10;i++)
  {
    ent=trim(xyz_all,0); 
    kep=0;
    while(xyz_all[kep])kep++;
    if(kep<got+2) break;
  }
  return ent;
}

float **Segen::refill(float **xyz_all,int need)
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
  
  if(TRES.logg)cerr<<"start refilling..."<<endl;
  if(!strstr(dforce,"?")) {ent=trim(near,3*part);if(ent) delete [] ent; ent=0;}
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
  delete [] xyz_all;xyz_all=0;
  xyz_all=xyz_all0;
  kep=0;while(xyz_all[kep])kep++;
  d=0;
  for(i=0;i<kep;i++)
  {
    segment->chn->transfer(xyz_all[i],1);
    e=rmsd(2);
    d+=e/kep;
    if(TRES.logg)cerr<<i<<" after refill: "<<rmsd(2)<<endl;
  }
  if(TRES.logg)cerr<<" the total number after refill: "<<kep<<"  "<<d<<endl;
  return xyz_all;
}

float **Segen::fill(float **xyz_all,float cut)
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
  delete [] near;near=0;
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
  if(TRES.logg)cerr<<" the total number after fill: "<<kk<<" "<<cut<<"  "<<kep0<<"  "<<d<<endl;
  chn->transfer(xyz,r0,end,1);
  delete [] xyz;xyz=0;
  return xyz_all;
}

float **Segen::copy(float **xyz_all)
{
  int kep;
  kep=0;
  while(xyz_all[kep])kep++;
  return copy(xyz_all,kep);
}

float **Segen::copy(float **xyz_all,int all)
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

float* Segen::trim(float **xyz_all,float*ent,int haa)
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
  delete [] order;kept=0;pot=0;entt=0;order=0;
  return ent;
}

float* Segen::trim(float **xyz_all,float *ent,float cutf)
{
 int i,kk,kep;
 kep=0;
 while(xyz_all[kep])kep++;
 kk=0;
 for(i=0;i<kep;i++)if(ent[i]<cutf)kk++;
 ent=trim(xyz_all,ent,kk);
 return ent;
}

float* Segen::trimold(float **xyz_all,int haa)
{
 int kep,i;
 float e,rms,rms0,*ent,*entt;
 char dforce0[100];
 Qsort cc;
 int *order;
 Chn *chn=(*pdb)[cid];
 Res *r0=chn->isres0(start);
 float *xyz0=chn->gettransfer(r0,end);
 kep=0;
 while(xyz_all[kep])kep++;
 strcpy(dforce0,dforce); 
 i=strlen(dforce);
 dforce[i]='?';dforce[i+1]='\0';
   
 ent=ensortold(xyz_all,dforce);
 order=new int[kep];
 entt=new float[kep];
 for(i=0;i<kep;i++) {entt[i]=ent[i];order[i]=i;}
 if(!strstr(dforce0,"?"))cc.sort(entt,kep,order); 
 rms=0;
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[order[i]],1);
   e=rmsd(2);
   if(TRES.logg)cout<<"  "<<i<<" energy of "<<dforce0<<": "<<e<<"  "<<entt[i]<<endl;
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
 if(TRES.logg)cout<<" the total number after "<<dforce0<<" "<<kep<<"  "<<rms<<"  "<<rms0<<endl;
 delete [] entt;delete [] order;entt=0;order=0;
 chn->transfer(xyz0,r0,end,1);
 delete [] xyz0;xyz0=0;
 return ent;
}

float* Segen::trim(float **xyz_all,int haa)
{
 int kep,i;
 float e,rms,rms0,*ent,*entt;
 char dforce0[100];
 Qsort cc;
 int *order;
 Chn *chn=(*pdb)[cid];
 Res *r0=chn->isres0(start);
 float *xyz0=chn->gettransfer(r0,end);

 kep=0;
 while(xyz_all[kep])kep++;

 strcpy(dforce0,dforce); 
 i=strlen(dforce);
 dforce[i]='?';dforce[i+1]='\0';
   
 ent=ensort(xyz_all);
 order=new int[kep];
 entt=new float[kep];
 for(i=0;i<kep;i++) {entt[i]=ent[i];order[i]=i;}
 if(!strstr(dforce0,"?"))cc.sort(entt,kep,order); 
 rms=0;
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[order[i]],1);
   e=rmsd(2);
   if(TRES.logg>2)cout<<"  "<<i<<" energy of "<<dforce0<<": "<<e<<"  "<<entt[i]<<endl;
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
 if(TRES.logg>2)cout<<" the total number after "<<dforce0<<" "<<kep<<"  "<<rms<<"  "<<rms0<<endl;
 delete [] entt;delete [] order;entt=0;order=0;
 chn->transfer(xyz0,r0,end,1);
 delete [] xyz0;xyz0=0;
 return ent;
}

float* Segen::trim(float **xyz_all,int haa,int fe)
{
 int kep,i;
 float e,rms,rms0,*ent,*entt;
 char flg0[100];
 Qsort cc;
 int *order;

 kep=0;
 while(xyz_all[kep])kep++;
 strcpy(flg0,dforce);
 i=strlen(dforce);
 flg0[i]='?';flg0[i+1]='\0';
 ent=ensort(xyz_all);
 order=new int[kep];
 entt=new float[kep];
 for(i=0;i<kep;i++) {entt[i]=ent[i];order[i]=i;}
 if(!strstr(dforce,"?"))cc.sort(entt,kep,order);
 rms=0;
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[order[i]],1);
   e=rmsd(2);
   if(fe==0) {
   if(TRES.logg)cerr<<"  "<<i<<" energy of "<<dforce<<": "<<e<<"  "<<entt[i]<<endl;
   }
   else {
   if(TRES.logg)cerr<<"  "<<i<<" energy of "<<dforce<<": "<<e<<"  "<<entt[i]<<endl; 
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
 if(TRES.logg)cerr<<" the total number after "<<dforce<<" "<<kep<<"  "<<rms<<"  "<<rms0<<endl;
 delete [] entt;delete [] order;entt=0;order=0;
 return ent;
}


void Segen::setoff(float **xyz_all,float ct)
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


float **Segen::generateold(float **xyz_all,int cop,float close)
{

   Chn *chn=(*pdb)[cid];
   Res *r0=chn->isres0(start);

   //save original
   float *xyz0=chn->gettransfer(r0,end);
   
   //fusion
   
   xyz_all=cross(xyz_all,cop,close);
    
   //save back
   chn->transfer(xyz0,r0,end,1);
   if(xyz0) delete [] xyz0;xyz0=0;

   //
   return xyz_all;
}
 
float **Segen::generate(float **xyz_all,int cop,float close)
{
   float *ent;
   int i,kep,kep0;
   kep0=0;while(xyz_all[kep0])kep0++;
   xyz_all=cross(xyz_all,cop,close);
   kep=0;while(xyz_all[kep])kep++;

   for(i=kep0;i<kep;i++) minimizefix(xyz_all[i],1);

   strcpy(dforce,"vw");
   if(kep>1.3*part) {ent=trim(xyz_all,0); if(ent) delete [] ent; ent=0;}
   kep=noclose(xyz_all,close);
   cout<<" the total number after noclose of "<<close<<" :"<<kep<<endl;
   for(i=(int)(1.3*part);i<kep;i++) {delete [] xyz_all[i];xyz_all[i]=0;}
   return xyz_all;
}

int **Segen::allconf(int nall)
{
  Res *r,*r1;
  int **outid;
  int k,n,m,i,j,jj;
 
  n=1; 
  for(r=segment->chn->res;r;r=r->next)
  {
    //r1=TRES.rotamer->chn->isres(r->name,0);
    r1=TRES.findrotamer("backbone")->chn->isres(r->name,0);
    n=n*r1->nummore;
  }
  if(TRES.logg)cerr<<" the total number of conformation is :"<<n<<endl;
  if(nall>n) nall=n;
  n=end-start+1;
  outid=new int*[nall+10];
  for(i=0;i<nall+10;i++) outid[i]=0;
  for(i=0;i<nall;i++) outid[i]=new int[n];
  for(i=0;i<nall;i++) for(j=0;j<n;j++) outid[i][j]=0;

  k=0;n=1;
  for(r=segment->chn->res;r;r=r->next)
  {    
    //r1=TRES.rotamer->chn->isres(r->name,0);
     r1=TRES.findrotamer("backbone")->chn->isres(r->name,0);
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

float Segen::segmin(int fast)
{
  Res *r,*r0,*r1,*last,*last0;
  Atm *aa2,*aa0,*aa1,*all[50];
  float dn,dca,d,dmin;
  float xr[3],xr1[3],xr2[3];
  float **coeff=0,*unknown=0,offset;
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
  coeff=new float*[6];
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
    dcut=-1;
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
	if(d<0.001) d=1;
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
          if(di>10000) di=10000;
	  else if(di<-10000) di=-10000;
          if(dca>10000) dca=10000;
	  else if(dca<-10000) dca=-10000;
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
	  //if(di>10000) di=10000;
	  //else if(di<-10000) di=-10000;
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

    if(fast==0&&rec2>6&&qut) break;
    if(fast==1&&rec2>1&&qut) break; 
    if(fast==2&&rec2>1&&qut) break; 
    if(qut&&fast==1&&fabs(exp1st[mseq]-dmin)<0.1&&dcut<outlet) break;
    if(qut&&fast==2&&fabs(exp1st[mseq]-dmin)<0.1&&dcut<outlet) break;

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
      if(dmin<0&&fast==2) break;
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

    if(TRES.logg>3)cerr<<rec<<"energy:  "<<only<<"  "<<exp1st[mseq]<<"  "<<dcut<<"  "<<d<<endl;


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
	if(fabs(d)>360) {m++;continue;}
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
  for(i=0;i<6;i++)   if(coeff[i]) { delete [] coeff[i];coeff[i]=0;}
  if(coeff) delete [] coeff;coeff=0;
  delete [] unknown;unknown=0;
  delete [] order;order=0;
  delete [] exp2st;exp2st=0;
  delete [] exp1st;exp1st=0;
  delete [] xyz;xyz=0;
  delete [] xyz0;xyz0=0;
  for(i=0;i<nseq;i++) delete [] effco[i];
  delete [] effco;effco=0;
  for(i=0;i<nseq/3;i++) delete [] affect[i];
  delete [] affect;affect=0;
  delete [] xyz1;xyz1=0; delete [] mcsc;mcsc=0;
  delete [] dist;dist=0;
  if(qut==0) dmin=pow(10,25);
  return dmin;
}

float Segen::boundmin(int fast)
{
  Res *r,*r0,*r1,*last,*last0;
  Atm *aa2,*aa0,*aa1,*all[50];
  float dn,dca,d,dmin;
  float xr[3],xr1[3],xr2[3];
  float **coeff,*unknown,offset;
  float **effco,*exp1st,*xyz1,*xyz,*xyz0,lamba[6],*exp2st;
  int i,j,k,m,mseq,nseq,*mcsc,mdel; 
  int rec,rec1,rec2,*order,nuse,endcon;
  int **affect,n;
  Qsort qsort;
  Rotate rot;
  float ta,tb,tc,tn,ds,de;
  float dcut,di; 
  int isp;
  float *dist;
  int qut;
  int mr;
  int only,onlysum;
  Chn *chn=(*pdb)[cid];
  Chn *cchn;
  float engprev=0;
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
  coeff=new float*[6];
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
  
  float tt=0;
  if(smoothclash) {
	ta=TRES.engcoeff->coeff[0];
	tb=TRES.engcoeff->coeff[1];
	tc=TRES.engcoeff->coeff[2];
	tt=TRES.engcoeff->coeff[3];
	tn=TRES.engcoeff->coeff[4];
  }
  else {
  	ta=TRES.smooth[TRES.smt*4+0];
  	tb=TRES.smooth[TRES.smt*4+1];
  	tc=TRES.smooth[TRES.smt*4+2];
  	tn=TRES.smooth[TRES.smt*4+3];
	tt=0;
  }
 
  int minget=0; //number of structure get, modified 3-15-02
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
    dcut=-1;
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
      if(dcut>outlet) only=1;
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
	if(d<0.001) d=1;
        if(onlysidechain==1&&aa2->tatm->id<=3) {m++;continue;}
  	if(rotatm&&ifrotatmexist(aa2)==0) {m++;continue;}
        for(j=0;j<3;j++)xr[j]=xr[j]/d; //unit vector of axis

        // create equation satifying end constraint
        rec1=0;
        if(aa2->tatm->id>3) goto re300;
        for(k=0;k<3;k++)
        {
          if(onlysidechain==1) break;
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
          if(onlysidechain&&mcsc[m]==1) break;
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
    if(endcon&&onlysidechain==0) 
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
    if(dace&&only==0) {
	dace->update(xyz1);
    } 
    for(i=0;i<3;i++) xr1[i]=0;
    for(r=r0;r;r=r->next)
    {
      if(only) break;
      if(r->id0>end) break; 
      for(aa0=r->atm;aa0;aa0=aa0->next)
      {
        if(aa0->flag<0) continue;
	//if(onlybound==1) break;
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

	if(onlybound==1) {
		lat->nget=0;
	}	
	else {
        	lat->getcell(aa0,cutoff);
        	lat->cutoff(aa0,cutoff);
	}
	Charge *ch0=disc->findcharge(aa0);
	 
 	int nte=0;
	if(strchr(dforce,'d')||strchr(dforce,'D')) nte=1;
	if(nte) {
		int io=pdb->maxatmid0()+100;
		if(disc->change) delete [] disc->change; disc->change=0;
		disc->change=new float[io];
		int i; for(i=0;i<io;i++) disc->change[i]=0;
	}
	
	for(m=0;nte&&m<lat->nget;m++) {				
		aa1=lat->obtain[m]->atm;
		if(aa1==0) continue;		
		Dipole *t=disc->range[aa1->id0];
		if(t==0) continue;
		if(t->crg==0||t->crg[0]==0||t->crg[1]==0) {
			cerr<<"error in dipole's charge"<<endl;
			exit(0);
		}
		lat->addifnotexist(t->crg[0]->atm);
		lat->addifnotexist(t->crg[1]->atm);
		disc->change[t->crg[0]->atm->id0]+=t->crg[0]->crg;
        	disc->change[t->crg[1]->atm->id0]+=t->crg[1]->crg;
	}

	int ntr=0;
	if(strchr(dforce,'u')) ntr++;
	if(strchr(dforce,'v')) ntr++;
	if(strchr(dforce,'d')) ntr++;
	if(strchr(dforce,'D')) ntr++;
	if(TRES.logg>4) lat->printoutLattice();
        for(m=0;ntr&&m<lat->nget;m++)
        {
          aa1=lat->obtain[m]->atm;
          if(aa0->flag<0&&aa1->flag<0) continue;
          if(aa0->flag>=0&&aa1->flag>=0){if(aa1->flag<aa0->flag) continue;}
	  if(aa0==aa1||aa0->isnear(aa1)) continue;
	 
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
  	  ds=dn;
          dn=1/(tt+dn);
          xr2[1]=pow(dn,tn);  
          xr2[2]=xr2[1]*xr2[1];
          dca=tn*xr2[1]*dn*(tc-2*xr2[1])
             +tb*xr2[1]*(tc-xr2[1]);
          dca=dca*xr2[0];
          di=xr2[0]*xr2[1]*(xr2[1]-tc);
          if(di>10000) di=10000;
	  else if(di<-10000) di=-10000;
          if(dca>10000) dca=10000;
	  else if(dca<-10000) dca=-10000;
          isp=aa0->tatm->ispolar*aa1->tatm->ispolar;
          int nte=0;
	  if(strchr(dforce,'v')) nte=1;
          if(nte&&di<0)
          {
             if(mr==1)       {di=1.50*di;dca=1.50*dca;}
             else if(mr==2)  {di=5*di;dca=5*dca;}
             else if(isp==1) {di=1.50*di;dca=1.50*dca;}
             else if(isp==2) {di=1.25*di;dca=1.25*dca;}
             else if(isp==3) {di=0.75*di;dca=0.75*dca;}	    
          }
          else if(nte)
          {
             if(mr==1)       {di=0.50*di;dca=0.50*dca;}
             else if(mr==2)  {di=0.20*di;dca=0.20*dca;}
             else if(isp==1) {di=0.50*di;dca=0.50*dca;}
             else if(isp==2) {di=0.75*di;dca=0.75*dca;}
             else if(isp==3) {di=1.25*di;dca=1.25*dca;}	     
          }
	 
	  if(nte) {
		if(dace==0) {
			exp1st[mseq]+=di;
             		for(j=0;j<3;j++) xr1[j]=2*xr[j]*dca/d; 
		}
		else {
			exp1st[mseq]+=di*dace->coeff;
			for(j=0;j<3;j++) xr1[j]=2*xr[j]*dca/d*dace->weight;
		}
	  }

	  nte=0;
	  if(strchr(dforce,'u')) {		
		nte=1;
	  }
	  int hp=0;
	  if(aa0->tatm->id>=4&&aa1->tatm->id>=4) {
		if(aa0->tatm->name[1]=='H'&&aa0->tatm->bond[0]->id<4) {
			hp=0;
		}	  
		else if(aa1->tatm->name[1]=='H'&&aa1->tatm->bond[0]->id<4) {
			hp=0;
		}	  
		else if(strchr("AVLIFCMPW",aa0->res->name)&&strchr("AVLIFCMPW",aa1->res->name)) {
			hp=1;
		} 
     	  }
	  if(di>0) {
  		if(hp==1) {di=di/10; dca=dca/10;}
          }
 	  else {
		if(hp==0) {di=di/10; dca=dca/10;}
	  }
     	  if(nte) {
		if(dace==0) {
			exp1st[mseq]+=di;
             		for(j=0;j<3;j++) xr1[j]=2*xr[j]*dca/d; 
		}
		else {
			exp1st[mseq]+=di*dace->coeff;
			for(j=0;j<3;j++) xr1[j]=2*xr[j]*dca/d*dace->weight;
		}
	  } 
          //exp1st[mseq]+=di; 
	  //cerr<<m<<" "<<aa0->res->name<<" "<<aa0->name<<" "<<aa0->xyz[1]<<" "<<aa1->res->name<<" "<<aa1->name<<" "<<aa1->xyz[1]<<" "<<di<<endl;	
          //for(j=0;j<3;j++) xr1[j]=2*xr[j]*dca/d; 
          
          nte=0;
	  if(strchr(dforce,'D')) nte=1;
          if(nte&&disc->change[aa0->id0]!=0&&disc->change[aa1->id0]!=0&&aa0->res!=aa1->res)
          {	    
	    de=sqrt((ds+tt)*d*d);     
	    if(de<3.5) de=3.5;     
            dca=disc->change[aa0->id0]*disc->change[aa1->id0]*332*0.25/disc->eps/de;
            di=-0.5*dca/de/de;
	    if(di>10000) di=10000;
	    else if(di<-10000) di=-10000;
	    if(dca>10000) dca=10000;
	    else if(dca<-10000) dca=-10000;
	    if(dace==0) { 	   
            	exp1st[mseq]+=dca;
            	for(j=0;j<3;j++) xr1[j]+=2*xr[j]*d*di;
	    }
	    else {
		exp1st[mseq]+=dca*dace->coeff;
            	for(j=0;j<3;j++) xr1[j]+=2*xr[j]*d*di*dace->weight;
	    }
          }            
          	  
	  nte=0;
	  if(strchr(dforce,'d')) nte=1;
          if(ch0&&nte&&disc->change[aa1->id0]!=0&&aa0->res!=aa1->res)
          {	   
	    de=sqrt((ds+tt)*d*d);     
	    if(de<3.5) de=3.5;     
            dca=ch0->crg*disc->change[aa1->id0]*332*0.25/disc->eps/de;
            di=-0.5*dca/de/de; 	   
	    if(di>10000) di=10000;
	    else if(di<-10000) di=-10000;
	    if(dca>10000) dca=10000;
	    else if(dca<-10000) dca=-10000;
	    if(dace==0) {
            	exp1st[mseq]+=dca;
            	for(j=0;j<3;j++) xr1[j]+=2*xr[j]*d*di; 
            }
	    else {
		exp1st[mseq]+=dca*dace->coeff;
            	for(j=0;j<3;j++) xr1[j]+=2*xr[j]*d*di*dace->weight;
	    }
          }          
	  nte=0;

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

        nte=0;
	if(strchr(dforce,'c')) nte=1;
	int kep=0;while(disc->charge[kep])kep++;
	for(m=0;nte&&ch0&&m<kep;m++) {
		aa1=disc->charge[m]->atm;
		d=aa0->tatm->eng->radius+aa1->tatm->eng->radius;
		for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j])/d;
          	ds=TRES.distsqr(xr);
		
		de=sqrt((ds+tt)*d*d);     
	  	if(de<3.5) de=3.5;     
            	dca=disc->charge[m]->crg*ch0->crg*332/disc->eps/de;
            	di=-0.5*dca/de/de; 	
		if(di>10000) di=10000;
	    	else if(di<-10000) di=-10000;
	    	if(dca>10000) dca=10000;
	    	else if(dca<-10000) dca=-10000;   
		if(dace==0) {
            		exp1st[mseq]+=dca;
            		for(j=0;j<3;j++) xr1[j]=2*xr[j]*d*di;
	        }
		else {
			exp1st[mseq]+=dca*dace->coeff;
            		for(j=0;j<3;j++) xr1[j]=2*xr[j]*d*di*dace->weight;
		}
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

	if(dace&&dace->more) {
		Dace *s;		 
		for(s=dace->more;s;s=s->more) {
			if(aa0->tatm->id>4) break;
			Atm *aa1=s->atoms[aa0->id0];
			if(aa1==0) continue;
			if(aa1->tatm->id>4) continue;
			if(aa1->res->id0!=aa0->res->id0) {
				if(TRES.logg) cerr<<"the problem detected..."<<endl;
				if(TRES.logg) exit(0);
				continue;
			}
			d=aa0->tatm->eng->radius+aa1->tatm->eng->radius;
			for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j])/d;			
			for(j=0;j<3;j++) xr1[j]=2*xr[j]*d*s->weight;
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
	}
	//bound constraint
	float a1;
	
        if(bound) {
		
		int bt=0;
		for(bt=bound->tag[2*aa0->id0];bt<=bound->tag[2*aa0->id0+1];bt++) {
			if(bt==-1) break;
			if(aa0->flag<0) break;
			if(bound->predt[aa0->id0]==0) break;
			aa1=bound->atoms[bound->atmpair[2*bt+1]];
			//if(aa1==0) continue;
			if(aa1->res->chn==aa0->res->chn) {
				if(aa0->flag<0&&aa1->flag<0) continue;
          			if(aa0->flag>=0&&aa1->flag>=0){if(aa1->flag<aa0->flag) continue;}				
			}
			for(j=0;j<3;j++) xr[j]=(aa0->xyz[j]-aa1->xyz[j]);
          		dn=TRES.distsqr(xr);
			float d=sqrt(dn);
			float mid=bound->dist[3*bt];
			float low=bound->dist[3*bt+1];
			float high=bound->dist[3*bt+2];
			float weight=bound->weight[bt];

			if(onlybound==0||1) {
				if(d==0) d=0.1;	
					
				if(d>=low&&d<=high) {				 
					if(d<mid) {
						a1=mid-d;
						if(a1>10) a1=10;
						dca=-weight*exp(-a1)/8;
						di=-weight*exp(-a1)/8;					 
					}
					else if(d>mid) {
						a1=d-mid;
						if(a1>10) a1=10;
						dca=weight*exp(-a1)/8;
						di=-weight*exp(-a1)/8;					 
					}
					else {
						di=-weight*0.125;dca=weight*0;
					}
				 									 
				}
				else if(d<low) {
					a1=low-d;
					if(a1>10) a1=10;
					dca=-weight*exp(2*a1);
					di=weight*exp(2*a1);
				 
				}
				else if(d>high) {
					a1=d-high;
					if(a1>10) a1=10;
					dca=weight*exp(2*a1);
					di=weight*exp(2*a1);				 
				}		 	 
			}
			else {
				if(d>mid) {
					a1=d*d-mid*mid;
					dca=weight*1;
					di=weight*a1;
				}
				else {
					a1=mid*mid-d*d;
					dca=-weight*1;
					di=weight*a1;
				}
				//dca*=100000;
				//di*=100000;
			}
			if(di>10000) di=10000;
	    		else if(di<-10000) di=-10000;
	    		if(dca>10000) dca=10000;
	    		else if(dca<-10000) dca=-10000;
			for(j=0;j<3;j++) xr1[j]=2*xr[j]*dca;	
			if(onlyenergy==0) exp1st[mseq]+=di; 
			for(n=0;n<mseq;n++)
          		{				 
             			i=affect[aa0->flag][n];
             			if(i==-1) continue;
             			k=order[i];
             			for(j=0;j<3;j++)    xr2[j]=effco[aa0->flag*3+j][k];
             			exp1st[i]+=xr1[0]*xr2[0]+xr1[1]*xr2[1]+xr1[2]*xr2[2];
          		}
			
          		for(n=0;n<mseq;n++)
          		{
				if(aa1->res->chn!=aa0->res->chn) break;
				if(bound->predt[aa1->id0]==0) break;    
				if(aa1->flag<0) break;         			
             			i=affect[aa1->flag][n];
             			if(i==-1) break;
             			k=order[i];
             			for(j=0;j<3;j++) xr2[j]=-effco[aa1->flag*3+j][k];
             			exp1st[i]+=xr1[0]*xr2[0]+xr1[1]*xr2[1]+xr1[2]*xr2[2];
          		}
		}		
		    	
  	}
	//end

        for(j=0;j<mdel;j++) lat->puton(all[j]);
      }
    }


   if(dace&&dace->more) {	 
	exp1st[mseq]=dace->energy;	
   }
    //exit(0);
    if(fast==0&&rec2>8&&qut) break;
    if(fast==1&&rec2>1&&qut) break;
    if(fast==4&&rec2>3&&qut) break;
    if(qut&&fast==1&&fabs(exp1st[mseq]-dmin)<0.1&&dcut<outlet) break;
    if(qut&&fast==4&&fabs(exp1st[mseq]-dmin)<0.1&&dcut<outlet) break;

    if(fast==0&&rec2>8&&qut&&only==0) break;
    if((fast==1||fast==3)&&rec2>3&&qut&&only==0) break; 
    if((fast==4||fast==3)&&rec2>3&&qut&&only==0) break; 
    if(fast==2&&rec2>3&&qut&&only==0) break;
    if(qut&&(fast==3||fast==1||fast==4)&&fabs(exp1st[mseq]-dmin)<0.1&&fabs(exp1st[mseq]-dmin)<0.05*fabs(exp1st[mseq])&&dcut<outlet&&only==0) break;
    if(qut&&fast==2&&fabs(exp1st[mseq]-dmin)<0.1&&fabs(exp1st[mseq]-dmin)<0.05*fabs(exp1st[mseq])&&dcut<outlet&&only==0) break;
    
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
	 float dxe=(50+random()%50)*0.01;
	 if(randomtorsion==0||only==1)  dxe=1;
         if(d>step)d=step*dxe;
         else if(d<-step) d=-step*dxe;
         for(i=0;i<nuse;i++) lamba[i]=lamba[i]+coeff[i][j]*d;
         if(fabs(d-unknown[j])>0.01*fabs(d))m=1;
         unknown[j]=d;
       }
       if(m==0) break;
    }
    
    for(i=0;i<nuse;i++) 
    {
       d=lamba[i];
       float dxe=(50+random()%50)*0.01;
       if(randomtorsion==0||only==1)  dxe=1;

       if(d>step)d=step*dxe;
       else if(d<-step) d=-step*dxe;
       unknown[i]=d;
    }
    
    d=0; //calculating decreased energy 
    for(i=0;i<mseq;i++)d+=exp1st[i]*unknown[i];

    if(rec==0&&only==0&&dcut<outlet) {
      segment->chn->transfer(xyz,0);dmin=exp1st[mseq];rec2=0;qut=1;minget++;
      if(dmin<0&&fast==2) break;
    }
    else if(dmin>exp1st[mseq]&&dcut<outlet&&only==0) //store lower energy conformation
    {
      segment->chn->transfer(xyz,0);dmin=exp1st[mseq];rec2=0;qut=1;minget++;
      if(dmin<0&&fast==2) break;
    }
    else if(engprev>exp1st[mseq]&&dcut<outlet&&only==0&&fast==3) {
      rec2=0;
    }
    else if(dmin<exp1st[mseq]&&dcut<outlet&&only==0)
    {
      rec2++;
    }
    if(dcut<outlet&&only==0) engprev=exp1st[mseq];
    /*
    if(dcut<outlet&&minget==0&&only==0) { //to make sure the first time is assigned
	dmin=exp1st[mseq];qut=1;minget++;
    }
    */

    for(i=0;i<mseq;i++)
    {
       k=order[i];
       exp2st[k]=unknown[i];
    }
    for(i=0;i<mseq;i++) unknown[i]=exp2st[i];

    if(TRES.logg>3) cerr<<rec<<"energy:  "<<only<<"  "<<exp1st[mseq]<<"  "<<dcut<<"  "<<d<<endl;


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
	if(fabs(d)>360) {m++;continue;}
	
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
  for(i=0;i<6;i++) if(coeff[i]) {delete [] coeff[i];coeff[i]=0;}
  delete [] coeff;coeff=0;
  delete [] unknown;unknown=0;
  delete [] order;order=0;
  delete [] exp2st;exp2st=0;
  delete [] exp1st;exp1st=0;
  delete [] xyz;xyz=0;
  delete [] xyz0;xyz0=0;
  for(i=0;i<nseq;i++) delete [] effco[i];
  delete [] effco;effco=0;
  for(i=0;i<nseq/3;i++) delete [] affect[i];
  delete [] affect;affect=0;
  delete [] xyz1;xyz1=0;delete [] mcsc;mcsc=0;
  delete [] dist;dist=0;
  if(qut==0) dmin=pow(10,25);
  return dmin;
}



float Segen::sidechainmin(int fast)
{
  Res *r,*r0,*r1,*last,*last0;
  Atm *aa2,*aa0,*aa1,*all[50];
  float dn,dca,d,dmin;
  float xr[3],xr1[3],xr2[3];
  float **coeff,*unknown,offset;
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
  coeff=new float*[6];
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

  int sidechainonly=1;
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
	if(d<0.001) d=1;
        for(j=0;j<3;j++)xr[j]=xr[j]/d; //unit vector of axis

        // create equation satifying end constraint
        rec1=0;
        if(aa2->tatm->id>3) goto re300;
        for(k=0;k<3;k++)
        {
          if(sidechainonly==1) break;
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
	  if(sidechainonly&&mcsc[m]==1) break;
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
    if(endcon&&sidechainonly==0) 
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
	  if(di>10000) di=10000;
	  else if(di<-10000) di=-10000;
	  if(dca>10000) dca=10000;
	  else if(dca<-10000) dca=-10000;
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

    if(fast==0&&rec2>6&&qut) break;
    if(fast==1&&rec2>1&&qut) break; 
    if(fast==2&&rec2>1&&qut) break;
    if(qut&&fast==1&&fabs(exp1st[mseq]-dmin)<0.1&&dcut<outlet) break;
    if(qut&&fast==2&&fabs(exp1st[mseq]-dmin)<0.1&&dcut<outlet) break; 
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

    if(TRES.logg>3)cerr<<rec<<"energy:  "<<only<<"  "<<exp1st[mseq]<<"  "<<dcut<<"  "<<d<<endl;


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
	if(fabs(d)>360) {m++;continue;}
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
  for(i=0;i<6;i++) if(coeff[i]) {delete [] coeff[i];coeff[i]=0;}
  delete [] coeff;coeff=0;
  delete [] unknown;unknown=0;
  delete [] order;order=0;
  delete [] exp2st;exp2st=0;
  delete [] exp1st;exp1st=0;
  delete [] xyz;xyz=0;
  delete [] xyz0;xyz0=0;
  for(i=0;i<nseq;i++) delete [] effco[i];
  delete [] effco;effco=0;
  for(i=0;i<nseq/3;i++) delete [] affect[i];
  delete [] affect;affect=0;
  delete [] xyz1;xyz1=0;delete [] mcsc;mcsc=0;
  delete [] dist;dist=0;
  if(qut==0) dmin=pow(10,25);
  return dmin;
}



void Segen::testclosure()
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

  xyz_all=zipper(30,0);
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
  if(TRES.logg)cerr<<"**********initial**************"<<endl; 
  ent=ensort(xyz_out);//ent=ensort(xyz_out,"v?");
  kep=0;while(xyz_out[kep])kep++;
  rms=0;
  chn->transfer(xyz,r0,end,1);
  for(i=0;i<kep;i++)
  {
   segment->chn->transfer(xyz_out[i],1);
   e=rmsd(2);
   if(TRES.logg)cerr<<"  "<<i<<" energy of "<<e<<"  "<<ent[i]<<endl;
   rms+=e/kep;
  }
  if(TRES.logg)cerr<<"average :"<<rms<<endl;
  if(ent) delete [] ent; ent=0;

  if(TRES.logg)cerr<<"**********random treawk****** **"<<endl;
  ent=ensort(xyz_all);//ent=ensort(xyz_all,"v?");
  
  rms=0;
  chn->transfer(xyz,r0,end,1);
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[i],1);
   e=rmsd(2);
   if(TRES.logg)cerr<<"  "<<i<<" energy of "<<e<<"  "<<ent[i]<<endl;
   rms+=e/kep;
 }
 if(TRES.logg)cerr<<"average :"<<rms<<endl;
 if(ent) delete [] ent; ent=0;


 for(i=0;i<kep;i++)
 {
    chn->transfer(xyz,r0,end,1);   
    segment->chn->transfer(xyz_all[i],1);
    cycle=200;outlet=0.3;flat=300;
    segmin(1);
    segment->chn->transfer(xyz_all[i],0);
 }
  if(TRES.logg)cerr<<"**********after minimization********"<<endl;
  ent=ensort(xyz_all);//ent=ensort(xyz_all,"v?");//
  rms=0;
  chn->transfer(xyz,r0,end,1);
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_all[i],1);
   e=rmsd(2);
   if(TRES.logg)cerr<<"  "<<i<<" energy of "<<e<<"  "<<ent[i]<<endl;
   rms+=e/kep;
 }
 if(TRES.logg)cerr<<"average :"<<rms<<endl;
 if(ent) delete [] ent; ent=0;


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
  if(TRES.logg)cerr<<"**********channeled treak*****"<<endl;
  ent=ensort(xyz_out);//ent=ensort(xyz_out,"v?");
  rms=0;
  chn->transfer(xyz,r0,end,1);
 for(i=0;i<kep;i++)
 {
   segment->chn->transfer(xyz_out[i],1);
   e=rmsd(2);
   if(TRES.logg)cerr<<"  "<<i<<" energy of "<<e<<"  "<<ent[i]<<endl;
   rms+=e/kep;
 }
 if(TRES.logg)cerr<<"average :"<<rms<<endl;
 if(ent) delete [] ent; ent=0;
}

int Segen::checkomega() {
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

void Segen::checkomega(FILE *fp) {
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

 

int Segen::linkviaboth(int first,int last,int mid) {
	Rotate rot;
	start=first;
	end=last;
	direct=1; 
	int n=getdirection(); 
	Chn *chn; 
	Res *r0,*r,*t0;//*t;  
	chn=(*pdb)[cid];
	if(n==0||n==1) {
  		
  		r0=(*chn)[start];
  		t0=(*chn)[end];
  		r=(*chn)[start-1];
  		//t=(*chn)[end+1];
		
		Res *rr;
		if(r==0) {
			if(TRES.logg>3)chn->write("s0");
			for(rr=t0;rr;rr=rr->last) {
				if(rr->last==0) continue;
				rot.link(rr->last,rr,rr->last->id0,0);
			}
			if(TRES.logg>3)chn->write("s1");				
		}
		else {
			for(rr=r0;rr;rr=rr->next) {
				if(rr->next==0) continue;
				rot.link(rr,rr->next,rr->next->id0,1);
			}
		}
		return 1;
	}

	Segen seg1,seg2; 

	seg1.pdb=pdb;seg1.cid=cid;
	seg2.pdb=pdb;seg2.cid=cid;
	
	//int mid0;
	if(mid-start+1>end-mid+1) {
		seg1.start=start;
		seg1.end=mid;
		seg2.start=mid;
		seg2.end=end;
	}
	else {
		seg1.start=start;
		seg1.end=mid+1;
		seg2.start=mid+1;
		seg2.end=end;
	}
	
	seg1.create();
	seg2.create();
	seg1.direct=1;
	seg2.direct=0;
	
	 
	if(TRES.logg>3)seg1.segment->write("seg1");
	if(TRES.logg>3)seg2.segment->write("seg2");
	if(mid-start+1>end-mid+1) {
		Res *r0,*r1;
		r0=seg2.segment->chn->res;
		r1=seg2.segment->chn->res->next; 
		rot.link(r0,r1,r0->id0,0);
	}
	else {
		Res *r0,*r1;
		r0=seg1.segment->chn->isres0(mid);
		r1=seg1.segment->chn->isres0(mid+1); 
		rot.link(r0,r1,r1->id0,1);  
	}
	if(TRES.logg>3)seg1.segment->write("seg1a");
	if(TRES.logg>3)seg2.segment->write("seg2a");
	float d=closeg(&seg1,&seg2);	
	
	if(fabs(d-2)>0.001) {
		
		
		int nseq,m3,m4;
		nseq=0;m3=0;
		r0=(*chn)[start];
  		nseq=0;m3=0;
  		for(r=r0;r;r=r->next)
  		{
    			if(r->id0<=end) nseq+=r->tres->number*3;
    			if(r->id0<seg2.start) m3+=r->tres->number*3;
  		}
		m4=m3+seg2.segment->chn->res->tres->number*3;
		float *xyz_out=new float[nseq];
		Chn *cht=seg1.segment->chn; 
		r=cht->res;
		seg1.segment->chn->transfer(xyz_out,r,seg1.end,0);
		cht=seg2.segment->chn; r=cht->res->next;
     		cht->transfer(xyz_out+m4,r,seg2.end,0);
		seg1.pdb=0;seg2.pdb=0;
		r0=(*chn)[start];
		pdb->chn->transfer(xyz_out,r0,end,1);
		if(TRES.logg>3)pdb->write("s2");
		chn->header(r0,end+1);
		delete [] xyz_out;xyz_out=0;
		return 1;
	}
	else return 0;
} 

int Segen::getdirection()
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
     cerr<<"warning!!! this is not a loop prediction problem.no stems\n"; 
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
 
  return n;
}
 

void Segen::updatexyzsave(int *stem) {

	if(xyzsave) delete [] xyzsave;xyzsave=0;

	if(stem==0) return;
 
	Strhandler cc;
	Chn *chn=(*pdb)[cid];	
	int nlen=strlen(mutate->sqnto);
 
	int n1=	stem[0];
	int n2= stem[1];
	
	int start0,end0;

	if(n1==-1) start0=chn->res->id0;
	else if(n1==nlen) start0=chn->lastres()->id0;
	else	start0=mutate->match[n1];
	
	if(n2==-1) end0=chn->res->id0;
	else if(n2==nlen) end0=chn->lastres()->id0;	
	else	end0=mutate->match[n2];
		
	if(start0==-1||end0==-1) return;
	
	Res *r0=(*chn)[start0];
	if(r0==0) return;

	Res *r1=(*chn)[end0];
	if(r1==0) return;

	xyzsave=chn->gettransfer(r0,end0);
	return;
}

float *Segen::myfixsegment(int *stem) {

	if(stem==0) return 0;	
	Res *r;	
	Strhandler cc;
	Chn *chn=(*pdb)[cid];	
	int nlen=strlen(mutate->sqnto);
	
 
	int n1=	stem[0];
	int n2= stem[1];

	if(n1==-1) start=chn->res->id0;
	else if(n1==nlen) start=chn->lastres()->id0;
	else	start=mutate->match[n1];
	
	if(n2==-1) end=chn->res->id0;
	else if(n2==nlen) end=chn->lastres()->id0;	
	else	end=mutate->match[n2];
		
	if(start==-1||end==-1) return 0;
	
	int end0=end;
	int i;
	
	Res *r0=(*chn)[start];
  	int nseq=0;
  	for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
 	float *xyzout=new float[nseq];
	chn->transfer(xyzout,r0,end0,0);	

	readynewfix();
	
	float **out=predtfix();
	chn->transfer(xyzout,r0,end0,1);	

	r0=(*chn)[start];
	nseq=0;
  	for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
	delete [] xyzout;xyzout=0;
	xyzout=new float[nseq];
	chn->transfer(xyzout,r0,end,0);
	 
	if(out==0) {
		chn->transfer(xyzout,r0,end,1);
		delete [] xyzout;xyzout=0;
		return 0;
	}
	
	chn->transfer(out[0],r0,end,1);		 
	if(TRES.logg>1)chn->write("temp.pdb");	
	chn->transfer(xyzout,r0,end,1);
	if(TRES.logg>1)chn->write("semp.pdb");	

	for(i=0;i<nseq;i++) xyzout[i]=out[0][i];
	
	out=cc.floatdel(out);	
	
 	return xyzout;
} 

float **Segen::newfixsegment(int *stem) {

	if(stem==0) return 0;	
	Res *r;	
	Strhandler cc;
	Chn *chn=(*pdb)[cid];	
	int nlen=strlen(mutate->sqnto);
	
 
	int n1=	stem[0];
	int n2= stem[1];

	if(n1==-1) start=chn->res->id0;
	else if(n1==nlen) start=chn->lastres()->id0;
	else	start=mutate->match[n1];
	
	if(n2==-1) end=chn->res->id0;
	else if(n2==nlen) end=chn->lastres()->id0;	
	else	end=mutate->match[n2];
		
	if(start==-1||end==-1) return 0;
	
	int end0=end;
 
	
	Res *r0=(*chn)[start];
  	int nseq=0;
  	for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
 	float *xyzout=new float[nseq];
	chn->transfer(xyzout,r0,end0,0);	
	if(TRES.logg>3) chn->write("temp1.pdb");
	readynewfix();
	
	float **out=predtfix();
	chn->transfer(xyzout,r0,end0,1);	

	r0=(*chn)[start];
	nseq=0;
  	for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
	delete [] xyzout;xyzout=0;
	xyzout=new float[nseq];
	chn->transfer(xyzout,r0,end,0);
	 
	if(out==0) {
		chn->transfer(xyzout,r0,end,1);
		delete [] xyzout;xyzout=0;
		return 0;
	}
	
	chn->transfer(out[0],r0,end,1);		 
	if(TRES.logg>3)chn->write("temp.pdb");	
	chn->transfer(xyzout,r0,end,1);
	if(TRES.logg>3)chn->write("semp.pdb");	

	//for(i=0;i<nseq;i++) xyzout[i]=out[0][i];
	if(xyzout) delete [] xyzout; xyzout=0;
	//out=cc.floatdel(out);	
	
 	return out;
} 

int Segen::fixsegment(int *stem) {

	if(stem==0) return 0;	
	 
 	
	Strhandler cc;
	Chn *chn=(*pdb)[cid];	
	int nlen=strlen(mutate->sqnto);
	
	int nec=0;
	int n1=	stem[0];
	int n2= stem[1];

	if(n1==-1) start=chn->res->id0;
	else if(n1==nlen) start=chn->lastres()->id0;
	else	start=mutate->match[n1];
	
	if(n2==-1) end=chn->res->id0;
	else if(n2==nlen) end=chn->lastres()->id0;	
	else	end=mutate->match[n2];
		
	if(start==-1||end==-1) return 0;
	
	float **out=gofix();
  
	if(out) {
		nec++;
		Res *r=(*chn)[start];
		chn->transfer(out[0],r,end,1);	
		out=cc.floatdel(out);
 		r=(*chn)[start];
		chn->header(r,end+1);
	}
	return nec; 
}

int Segen::resort(SegBed *segbed) {

 	int i;
	
	Chn *chn=(*pdb)[cid];	
	 	
	start=segbed->first;
	end=segbed->last;
		
	if(start==-1||end==-1) return 0;
	
	readynewfix();
	Res *r0=chn->isres0(start);
	float *xyzorg=chn->gettransfer(r0,end);

	int kep=segbed->num;//int kep=0;while(segbed->xyzout[kep]) kep++;
	float **out=new float*[kep+1];
	
	int ii=0;
	for(i=0;i<kep;i++) {
		if(segbed->xyzout[i]==0) continue;
		out[ii++]=segbed->xyzout[i];
		segbed->xyzout[i]=0;
	}
	out[ii]=0; 
	segbed->num=0;
	float *ent;
	int nseq;
	Res *r;
	//int insert=mutate->calcinsert();
	if(0&&ii>1) {
		chn->transfer(xyzorg,r0,end,1);
		strcpy(boundtag,"@HPL");
        	out=avgminimize(out);
		if(bound) delete bound;bound=0;
        	r0=(*chn)[start];
        	nseq=0;
        	for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
        	delete [] xyzorg;xyzorg=0;
        	xyzorg = new float[nseq];
        	chn->transfer(xyzorg,r0,end,0);
	}

	int same=1; 
	for(i=0;i<segbed->num;i++) {
		if(segbed->start[i]!=start||segbed->end[i]!=end) {
			same=0;
			break;
		}	
	}
	
	//since segbed->num==0, always on
	if(0&&same) {
		float avgrmsd=calcavgrmsd(out,20);
		int kep=0;while(out[kep]) kep++;
		out=cross(out,kep,avgrmsd/2);		
	}
	//strcpy(dforce,"uDd@");
	float cu=cutoff;
	cutoff=10.;
	strcpy(dforce,"udD");
	ent=ensort(out);//ent=ensort(out,"vW"); 
	cutoff=cu;
	r0=chn->isres0(start);
	for(i=0;i<kep;i++) {	
		int n1=segbed->start[i];
		int n2=segbed->end[i];	 
		char line[1000];
		sprintf(line,"last.%i.pdb",i);
		chn->transfer(out[i],r0,end,1);	
		if(TRES.logg)cerr<<i<<":"<<n1<<" "<<n2<<" "<<" "<<ent[i]<<endl;
		if(TRES.logg>3)chn->write(line);		
	}
	if(ent) delete [] ent;ent=0;
	chn->transfer(out[0],r0,end,1);	
	chn->header(r0,end+1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
	return 1;
}
int Segen::resortmean(SegBed *segbed) {

 	int i;
	
	Chn *chn=(*pdb)[cid];	
	 	
	start=segbed->first;
	end=segbed->last;
		
	if(start==-1||end==-1) return 0;
	
	readynewfix();
	Res *r0=chn->isres0(start);
	float *xyzorg=chn->gettransfer(r0,end);

	int kep=segbed->num;//int kep=0;while(segbed->xyzout[kep]) kep++;
	float **out=new float*[kep+1];
	
	int ii=0;
	for(i=0;i<kep;i++) {
		if(segbed->xyzout[i]==0) continue;
		out[ii++]=segbed->xyzout[i];
		segbed->xyzout[i]=0;
	}
	out[ii]=0; 
	segbed->num=0;
	float *ent;
	int nseq;
	Res *r;
	//int insert=mutate->calcinsert();
	/*
	//buildhbond
	chn->transfer(xyzorg,r0,end,1);
	chn->clearhbond();
        chn->header();
        chn->buildhbond();
        chn->setdsspstr();
        chn->buildssbond();
        chn->setthreestatesec();
	//endbuildhbone
	*/
	int ne;
	cerr<<endl;
	cerr<<"minimize with constraints from the average structure..."<<endl;
	cerr<<endl;
	for(ne=0;ne<3;ne++) {
		int n=mutate->refine;
		randomtorsion=1;
		mutate->refine=1;
		chn->transfer(xyzorg,r0,end,1);
		strcpy(boundtag,"@P*");
        	out=avgminimize(out);
		if(bound) delete bound;bound=0;
        	r0=(*chn)[start];
        	nseq=0;
        	for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
        	delete [] xyzorg;xyzorg=0;
        	xyzorg = new float[nseq];
        	chn->transfer(xyzorg,r0,end,0);
		mutate->refine=n;
		randomtorsion=0;
	}
	//
	//chn->clearhbond();
	//

	int same=1; 
	for(i=0;i<segbed->num;i++) {
		if(segbed->start[i]!=start||segbed->end[i]!=end) {
			same=0;
			break;
		}	
	}
	
	//since segbed->num==0, always on
	if(0&&same) {
		float avgrmsd=calcavgrmsd(out,20);
		int kep=0;while(out[kep]) kep++;
		out=cross(out,kep,avgrmsd/2);		
	}
	//strcpy(dforce,"uDd@");
	float cu=cutoff;
	cutoff=10.;
	strcpy(dforce,"udDw");
	ent=ensort(out);//ent=ensort(out,"vW"); 
	cutoff=cu;
	r0=chn->isres0(start);
	for(i=0;i<kep;i++) {	
		int n1=segbed->start[i];
		int n2=segbed->end[i];	 
		char line[1000];
		sprintf(line,"last.%i.pdb",i);
		chn->transfer(out[i],r0,end,1);	
		if(TRES.logg)cerr<<i<<":"<<n1<<" "<<n2<<" "<<" "<<ent[i]<<endl;
		if(TRES.logg>3)chn->write(line);		
	}
	if(ent) delete [] ent;ent=0;
	chn->transfer(out[0],r0,end,1);	
	chn->header(r0,end+1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
	return 1;
}

float **Segen::gofix() {
  
	readyfix();
	float **tt=predtfix();		
	return tt;
}

float **Segen::gofix(int site) {
  
	readyfix();
	
	int n1=mutate->compare[start];
	int n2=mutate->compare[end];

	int i;
	int p1=0,p2=0;
	for(i=n1+1;i<n2;i++) {
		if(mutate->owner->seqngap[i]!='-') p1++;
		if(mutate->sqnto[i]!='-') p2++;
	}
	
	if(p1==p2) {
		if(TRES.logg)cerr<<"bad alignment...."<<endl;
		return 0;
	}
	float **tt=predtfix();		
	return tt;
}

void Segen::readyfix()
{
  Chn *chn;

  if(cid=='1')  cid=pdb->chn->id;

  //create another one segment class
	
  chn=(*pdb)[cid];
 
  //create segments
  Res *r=(*chn)[start];
  	
  chn->transform(r,end,3);
  pdb->configure();
  if(TRES.logg>3)pdb->write("seg1");
  
  //create segments
  create();
   

  Atm *a;
  for(r=chn->res;r;r=r->next) {
	if(r->id0<start-1||r->id0>end+1) continue;
	for(a=r->atm;a;a=a->next) a->allnear(3,0);
  }
  
  //segmentfix();
  //setup disc;
  if(disc) delete disc;disc=0;
  if(disc==0) disc=new Disc;
  disc->setupall(pdb,cutoff,1);
}

void Segen::readynewfix()
{
  Chn *chn;

  if(cid=='1')  cid=pdb->chn->id;

  //create another one segment class
	
  chn=(*pdb)[cid];
 
  //create segments
  Res *r=(*chn)[start];
 
  //create segments
  create();
   
  Atm *a;
  for(r=chn->res;r;r=r->next) {
	if(r->id0<start-1||r->id0>end+1) continue;
	for(a=r->atm;a;a=a->next) a->allnear(3,0);
  }

  //setup disc;
  if(disc) delete disc;disc=0;
  if(disc==0) disc=new Disc;
  disc->setupall(pdb,cutoff,1);
}


void Segen::setallnear() {

	Chn *chn =(*pdb)[cid];
	Atm *a;Res *r;
  	for(r=chn->res;r;r=r->next) {
		//if(r->id0<start-1||r->id0>end+1) continue;
		if(r->id0<start||r->id0>end) continue;
		for(a=r->atm;a;a=a->next) a->allnear(3,0);
  	}
  
} 

void Segen::segmentfix(){
	
 
	Res *r0=segment->chn->isres0(start);

	Res *r;
	Rotate rot;

	int h=0;
	int e=0;
	int c=0;
	
	if(TRES.logg>3)segment->write("seg1");
	for(r=r0;r;r=r->next) {
		int n=mutate->compare[r->id0];
		if(mutate->dssp[n]=='h') h++;
		if(mutate->dssp[n]=='e') e++;
		if(mutate->dssp[n]=='-') c++;
	} 
	
	if(e<=1&&direct==1) { //the segment is only alpha and coil
		Res *t=0;
		Res *x=0;
		for(r=r0;r;r=r->next) {		
			if(r->id0>end) break;							
			int n=mutate->findlongestreliable(r);		
			rot.link(t,r,n,1);
			x=segment->chn->isres0(n);			
			t=x;
			r=x;			
		}
		if(TRES.logg>3)segment->write("seg2");		
	}
 	else if(e<=1&&direct==0) {
		Res *t=0;
		Res *x=0;
		int i=0;
		for(i=end;i>=start;i--) {	
			r=segment->chn->isres0(i);										
			int n=mutate->findminuslongestreliable(r);		
			rot.link(r,t,n,0);
			//rot.link(t,r,n,0);
			x=segment->chn->isres0(n);			
			t=x;
			i=n;			
		}
		if(TRES.logg>3)segment->write("seg2");
	}
}


void Segen::setbound(int flg) {

	if(mutate==0) return;
	StrFmt *parent=mutate->owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;
	
	int i;
	int m=0;
	Chn *chn=(*pdb)[cid];

	Pdb *s;
	Chn *chn0=pdb->lastchain();	
	for(s=segment->next;s;s=s->next) {
		Chn *c=s->chn;
		if(c==0) continue;
		chn0->next=c;
		chn0=chn0->next;
	}
	chn0->next=mutate->owner->pdb->chn;
	chn0=pdb->lastchain();		
	if(mutate->mold) chn0->next=mutate->mold->mpdb->chn;  
	pdb->configureatmid0();
	bound=new Bound;
	bound->model=pdb;
	bound->flag|=TRES.constant->pdbundeletable;
	bound->easyready();
	
	int nlen=strlen(mutate->sqnto);
	float dd=0; 

	if(strstr(boundtag,"lllll")) dd=5*0.5;
	else if(strstr(boundtag,"llll")) dd=4*0.5;
	else if(strstr(boundtag,"lll")) dd=3*0.5;
	else if(strstr(boundtag,"ll")) dd=2*0.5;
	else if(strstr(boundtag,"l")) dd=1*0.5;
 
	float varn=dd;

	dd=0;
	if(strstr(boundtag,"bbbbb")) dd=5*0.5;
	else if(strstr(boundtag,"bbbb")) dd=4*0.5;
	else if(strstr(boundtag,"bbb")) dd=3*0.5;
	else if(strstr(boundtag,"bb")) dd=2*0.5;
	else if(strstr(boundtag,"b")) dd=1*0.5;
  
	float bvarn=dd;	

	Res *r0=(*chn)[start];
	Res *r,*r1;
	Atm *a,*b;
	
	//check for hydrogen bond constraint among pdb of owner
	DistPopular *hdist=TRES.popbin->gethbond();

	int nte=0;
	if(strchr(boundtag,'H')) nte=1;	
	else nte=0;
	
	//# means backbone hydrogen bond only
	for(r=r0;r&&nte;r=r->next) {
		
		if(r->id0>end) break;
		i=mutate->compare[r->id0];		 
		int j=mutate->owner->match[i];
		if(j==-1) continue;
		r1=mutate->owner->resn[j];
		if(r1==0) continue;
		HBondList *hlist;
		for(hlist=r1->hbond;hlist;hlist=hlist->next) {
			if(hlist->donor->res==hlist->acceptor->res) continue; //same residue
			a=0;b=0;
			if(hlist->donor->res==r1) {
				a=hlist->donor;
				b=hlist->acceptor;
			}
			else if(hlist->acceptor->res==r1) {
				a=hlist->acceptor;
				b=hlist->donor;
			}
			else {
				a=0;b=0;
			}
			if(a==0||b==0) continue;
			if(strstr(boundtag,"#")&&(a->tatm->id>3||b->tatm->id>3)) continue; 			

			Res *t1=b->res;
			
			Res *t=0;
			
			if(abs(t1->id0-r1->id0)<5) { //local hydrogen bond
				int nt1=t1->id0-r1->id0;
				int nt2=r->id0+nt1;
				t=pdb->chn->isres0(nt2);		
			}
			else {  		     //distant hydrogen bond
 				int nt1=mutate->owner->findpost(t1);			 			 
				if(nt1==-1) continue; 
				int nt2=mutate->match[nt1];
				if(nt2!=-1) t=mutate->resn[nt2];
			}
			
			if(t==0) continue;

			Atm *a1,*b1;			
			 			 
			a1=r->isatm(a->name);
			b1=t->isatm(b->name);											 			
			if(a1==0||b1==0) continue;
			if(a1->tatm->id!=a->tatm->id) continue;
			if(b1->tatm->id!=b->tatm->id) continue;		

			DistPopular *p;
			float x;
			if(t1->sec==r1->sec&&t1->sec!='-'&&a->tatm->id<=3&&b->tatm->id<=3) {
				p=hdist->getDistPopular(a1->res->tres,a1->name,b1->name,b1->res->tres->name);
				//p=hdist->getDistPopular(a1->name,b1->name,0);				
				//if(!p) p=hdist->getDistPopular(b1->name,a1->name,0);
				if(!p) p=hdist->getDistPopular(b1->res->tres,b1->name,a1->name,a1->res->tres->name);
				if(p) { 
					x=p->findmostpopulardistance();
					bound->addbounds(a1->id0,b1->id0,x,p->dist[0],p->dist[1]+bvarn,1);
				}	
				//
				if(a1->bond[0]==0||b1->bond[0]==0) continue;
				if(a->bond[0]==0||b->bond[0]==0) continue;
				//
				p=0;				
				//p=hdist->getDistPopular(a1->bond[0]->name,b1->name,0);
				p=hdist->getDistPopular(a1->res->tres,a1->bond[0]->name,b1->name,b1->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b1->name,a1->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b1->res->tres,b1->name,a1->bond[0]->name,a1->res->tres->name);
				if(p&&a1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					bound->addbounds(a1->bond[0]->id0,b1->id0,x,p->dist[0],p->dist[1]+bvarn,1);
				}
				//
				//p=hdist->getDistPopular(a1->name,b1->bond[0]->name,0);
				p=hdist->getDistPopular(a1->res->tres,a1->name,b1->bond[0]->name,b1->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b1->bond[0]->name,a1->name,0);
				if(!p) p=hdist->getDistPopular(b1->res->tres,b1->bond[0]->name,a1->name,a1->res->tres->name);
				if(p&&b1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					bound->addbounds(a1->id0,b1->bond[0]->id0,x,p->dist[0],p->dist[1]+bvarn,1);	
				}
				//
				//p=hdist->getDistPopular(a1->bond[0]->name,b1->bond[0]->name,0);
				p=hdist->getDistPopular(a1->res->tres,a1->bond[0]->name,b1->bond[0]->name,b1->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b1->bond[0]->name,a1->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b1->res->tres,b1->bond[0]->name,a1->bond[0]->name,a1->res->tres->name);
				if(p&&a1->tatm->id!=0&&b1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					bound->addbounds(a1->bond[0]->id0,b1->bond[0]->id0,x,p->dist[0],p->dist[1]+bvarn,1);
				}			
			} 	
			else {				
				//p=hdist->getDistPopular(a->name,b->name,0);
				p=hdist->getDistPopular(a->res->tres,a->name,b->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->name,0);
			 	if(!p) p=hdist->getDistPopular(b->res->tres,b->name,a->name,a->res->tres->name);
				if(p) {
					x=p->findmostpopulardistance();					
					float y=TRES.distance(a1,b1);
					if(strchr(boundtag,'@')) y=p->dist[1];
					float z=max(p->dist[1],y)+varn;	
					//z=100;					
					bound->addbounds(a1->id0,b1->id0,x,p->dist[0],z,1);
				}
				//
				if(a1->bond[0]==0||b1->bond[0]==0) continue;
				if(a->bond[0]==0||b->bond[0]==0) continue;
				//
				//p=hdist->getDistPopular(a->bond[0]->name,b->name,0);
				p=hdist->getDistPopular(a->res->tres,a->bond[0]->name,b->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->name,a->bond[0]->name,a->res->tres->name);
				if(p&&a1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					float y=TRES.distance(a1->bond[0],b1);
					if(strchr(boundtag,'@')) y=p->dist[1];
					float z=max(p->dist[1],y)+varn;	
					//z=100;
					bound->addbounds(a1->bond[0]->id0,b1->id0,x,p->dist[0],z,1);
				}
				//
				//p=hdist->getDistPopular(a1->name,b1->bond[0]->name,0);
				p=hdist->getDistPopular(a1->res->tres,a1->name,b1->bond[0]->name,b1->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->bond[0]->name,a->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->bond[0]->name,a->name,a->res->tres->name); 
				if(p&&b1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					float y=TRES.distance(a1,b1->bond[0]);
					if(strchr(boundtag,'@')) y=p->dist[1];
					float z=max(p->dist[1],y)+varn;	
					//z=100;
					bound->addbounds(a1->id0,b1->bond[0]->id0,x,p->dist[0],z,1);	
				}
				//
				//p=hdist->getDistPopular(a->bond[0]->name,b->bond[0]->name,0);
				p=hdist->getDistPopular(a->res->tres,a->bond[0]->name,b->bond[0]->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->bond[0]->name,a->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->bond[0]->name,a->bond[0]->name,a->res->tres->name); 
				if(p&&a1->tatm->id!=0&&b1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					float y=TRES.distance(a1->bond[0],b1->bond[0]);
					if(strchr(boundtag,'@')) y=p->dist[1];
					float z=max(p->dist[1],y)+varn;	
					//z=100;
					bound->addbounds(a1->bond[0]->id0,b1->bond[0]->id0,x,p->dist[0],z,1);
				}
			} 			
			m++;	
		}					
	}
	//end
	
	//hydrogen bond in mold template
	nte=0;
	if(strchr(boundtag,'h')&&mutate->mold) nte=1;	
	else nte=0;
	
	for(r=r0;r&&nte;r=r->next) {
		if(r->id0>end) break;
		r1=mutate->mold->findresid(r->id0);
		if(r1==0) continue;
		HBondList *hlist;
		for(hlist=r1->hbond;hlist;hlist=hlist->next) {
			if(hlist->donor->res==hlist->acceptor->res) continue; //same residue
			a=0;b=0;
			if(hlist->donor->res==r1) {
				a=hlist->donor;
				b=hlist->acceptor;
			}
			else if(hlist->acceptor->res==r1) {
				a=hlist->acceptor;
				b=hlist->donor;
			}
			else {
				a=0;b=0;
			}
			if(a==0||b==0) continue;
			if(strstr(boundtag,"#")&&(a->tatm->id>3||b->tatm->id>3)) continue; 
			//if(a->tatm->id>3||b->tatm->id>3) continue; 			

			Res *t1=b->res;
			
			Res *t=0;
			
			if(abs(t1->id0-r1->id0)<5) { //local hydrogen bond
				int nt1=mutate->mold->findresid(t1);
				t=pdb->chn->isres0(nt1);
			}
			else {  		     //distant hydrogen bond
				int nt1=mutate->mold->findresid(t1);
				t=pdb->chn->isres0(nt1); 				
			}
			
			if(t==0) continue;

			Atm *a1,*b1;			
			 
			a1=r->isatm(a->name);
			b1=t->isatm(b->name);											 			
			if(a1==0||b1==0) continue;
			if(a1->tatm->id!=a->tatm->id) continue;
			if(b1->tatm->id!=b->tatm->id) continue; 			
			
			DistPopular *p;
			float x,y,z;
			if(t1->sec==r1->sec&&t1->sec!='-'&&a->tatm->id<=3&&b->tatm->id<=3&&t1->name==t->name) {
				//p=hdist->getDistPopular(a1->name,b1->name,0);
				p=hdist->getDistPopular(a1->res->tres,a1->name,b1->name,b1->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b1->name,a1->name,0);
				if(!p) p=hdist->getDistPopular(b1->res->tres,b1->name,a1->name,a1->res->tres->name);
				if(p) { 										
					x=p->findmostpopulardistance();
					y=p->dist[0];
					z=p->dist[1];
					if(abs(t1->id0-r1->id0)>=5) {
						float e=TRES.distance(a1->xyz,b1->xyz);	
						if(strchr(boundtag,'@'))  e=p->dist[1];	
						e=max(e,z)+varn;					 
						bound->addbounds(a1->id0,b1->id0,x,y,e,1);
					}
					else bound->addbounds(a1->id0,b1->id0,x,y,z,1);
				}	
				//
				if(a1->bond[0]==0||b1->bond[0]==0) continue;
				if(a->bond[0]==0||b->bond[0]==0) continue;
				//
				p=0;				
				//p=hdist->getDistPopular(a1->bond[0]->name,b1->name,0);		
				p=hdist->getDistPopular(a1->res->tres,a1->bond[0]->name,b1->name,b1->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b1->name,a1->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b1->res->tres,b1->name,a1->bond[0]->name,a1->res->tres->name);
				if(p&&a1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					y=p->dist[0];
					z=p->dist[1];
					if(abs(t1->id0-r1->id0)>5) {
						float e=TRES.distance(a1->bond[0]->xyz,b1->xyz);
						if(strchr(boundtag,'@')) e=p->dist[1];   	
						e=max(e,z)+varn;							 
						bound->addbounds(a1->bond[0]->id0,b1->id0,x,y,e,1);
					}
					else bound->addbounds(a1->bond[0]->id0,b1->id0,x,y,z,1);
				}
				//
				//p=hdist->getDistPopular(a1->name,b1->bond[0]->name,0);
				p=hdist->getDistPopular(a1->res->tres,a1->name,b1->bond[0]->name,b1->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b1->bond[0]->name,a1->name,0);
				if(!p) p=hdist->getDistPopular(b1->res->tres,b1->bond[0]->name,a1->name,a1->res->tres->name);
				if(p&&b1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					y=p->dist[0];
					z=p->dist[1];
					if(abs(t1->id0-r1->id0)>5) {
						float e=TRES.distance(a1->xyz,b1->bond[0]->xyz);
						if(strchr(boundtag,'@')) e=p->dist[1]; 
						e=max(e,z)+varn;								 
						bound->addbounds(a1->id0,b1->bond[0]->id0,x,y,e,1);
					}
					else bound->addbounds(a1->id0,b1->bond[0]->id0,x,y,z,1);	
				}
				//
				//p=hdist->getDistPopular(a1->bond[0]->name,b1->bond[0]->name,0);
				p=hdist->getDistPopular(a1->res->tres,a1->bond[0]->name,b1->bond[0]->name,b1->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b1->bond[0]->name,a1->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b1->res->tres,b1->bond[0]->name,a1->bond[0]->name,a1->res->tres->name);
				if(p&&a1->tatm->id!=0&&b1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					y=p->dist[0];
					z=p->dist[1];
					if(abs(t1->id0-r1->id0)>5) {
						float e=TRES.distance(a1->bond[0]->xyz,b1->bond[0]->xyz);
						if(strchr(boundtag,'@')) e=p->dist[1]; 
						e=max(e,z)+varn;								 
						bound->addbounds(a1->bond[0]->id0,b1->bond[0]->id0,x,y,e,1);
					}
					else bound->addbounds(a1->bond[0]->id0,b1->bond[0]->id0,x,y,z,1);
				}			
			} 	
			else {				
				//p=hdist->getDistPopular(a->name,b->name,0);
				p=hdist->getDistPopular(a->res->tres,a->name,b->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->name,a->name,a->res->tres->name);
				if(p) {
					x=p->findmostpopulardistance();
					y=p->dist[0]; 
					z=p->dist[1];
					float e=TRES.distance(a1,b1);
					if(strchr(boundtag,'@')) e=p->dist[1]; 
					e=max(e,z)+varn;						
					bound->addbounds(a1->id0,b1->id0,x,y,e,1);
				}
				//
				if(a1->bond[0]==0||b1->bond[0]==0) continue;
				if(a->bond[0]==0||b->bond[0]==0) continue;
				//
				//p=hdist->getDistPopular(a->bond[0]->name,b->name,0);
				p=hdist->getDistPopular(a->res->tres,a->bond[0]->name,b->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->name,a->bond[0]->name,a->res->tres->name);
				if(p&&a1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					y=p->dist[0]; 	
					z=p->dist[1];
					float e=TRES.distance(a1->bond[0],b1);
					if(strchr(boundtag,'@')) e=p->dist[1]; 
					e=max(e,z)+varn;					
					bound->addbounds(a1->bond[0]->id0,b1->id0,x,y,e,1);
				}
				//
				//p=hdist->getDistPopular(a->name,b->bond[0]->name,0);
				p=hdist->getDistPopular(a->res->tres,a->name,b->bond[0]->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->bond[0]->name,a->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->bond[0]->name,a->name,a->res->tres->name); 
				if(p&&b1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					y=p->dist[0]; 	
					z=p->dist[1];
					float e=TRES.distance(b1->bond[0],a1);
					if(strchr(boundtag,'@')) e=p->dist[1]; 
					e=max(e,z)+varn;		
					bound->addbounds(a1->id0,b1->bond[0]->id0,x,y,e,1);	
				}
				//
				//p=hdist->getDistPopular(a->bond[0]->name,b->bond[0]->name,0);
				p=hdist->getDistPopular(a->res->tres,a->bond[0]->name,b->bond[0]->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->bond[0]->name,a->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->bond[0]->name,a->bond[0]->name,a->res->tres->name); 
				if(p&&a1->tatm->id!=0&&b1->tatm->id!=0) {
					x=p->findmostpopulardistance();
					y=p->dist[0]; 	
					z=p->dist[1];
					float e=TRES.distance(b1->bond[0],a1->bond[0]);
					if(strchr(boundtag,'@')) e=p->dist[1]; 
					e=max(e,z)+varn;	
					bound->addbounds(a1->bond[0]->id0,b1->bond[0]->id0,x,y,e,1);
				}
			} 			
			m++;	
		}					
	}
	//end



	//check hydrogen bond within itself
	nte=0;
	if(strchr(boundtag,'R')) nte=1;	
	else nte=0;
	
	for(r=r0;r&&nte;r=r->next) {
		
		if(r->id0>end) break;
		 	 
		HBondList *hlist;
		for(hlist=r->hbond;hlist;hlist=hlist->next) {
			if(hlist->donor->res==hlist->acceptor->res) continue; //same residue
 
			a=0;b=0;
			if(hlist->donor->res==r) {
				a=hlist->donor;
				b=hlist->acceptor;
			}
			else if(hlist->acceptor->res==r) {
				a=hlist->acceptor;
				b=hlist->donor;
			}
			else {
				a=0;b=0;
			}
			if(a==0||b==0) continue;
			if(strstr(boundtag,"#")&&(a->tatm->id>3||b->tatm->id>3)) continue; 
			//if(a->tatm->id>3||b->tatm->id>3) continue; 			

			Res *t=b->res;
 
			if(t==0) continue;
 
			DistPopular *p;
			float x;
			if(t->sec==r->sec&&t->sec!='-'&&a->tatm->id<=3&&b->tatm->id<=3) {
				//p=hdist->getDistPopular(a->name,b->name,0);
				p=hdist->getDistPopular(a->res->tres,a->name,b->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->name,a->name,a->res->tres->name);
				if(p) { 
					x=p->findmostpopulardistance();
					bound->addbounds(a->id0,b->id0,x,p->dist[0],p->dist[1]+bvarn,1);
				}	
				//
				if(a->bond[0]==0||b->bond[0]==0) continue;
				 
				//
				p=0;				
				//p=hdist->getDistPopular(a->bond[0]->name,b->name,0);
				p=hdist->getDistPopular(a->res->tres,a->bond[0]->name,b->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->name,a->bond[0]->name,a->res->tres->name);
				if(p&&a->tatm->id!=0) {
					x=p->findmostpopulardistance();
					bound->addbounds(a->bond[0]->id0,b->id0,x,p->dist[0],p->dist[1]+bvarn,1);
				}
				//
				//p=hdist->getDistPopular(a->name,b->bond[0]->name,0);
				p=hdist->getDistPopular(a->res->tres,a->name,b->bond[0]->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->bond[0]->name,a->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->bond[0]->name,a->name,a->res->tres->name);
				if(p&&b->tatm->id!=0) {
					x=p->findmostpopulardistance();
					bound->addbounds(a->id0,b->bond[0]->id0,x,p->dist[0],p->dist[1]+bvarn,1);	
				}
				//
				//p=hdist->getDistPopular(a->bond[0]->name,b->bond[0]->name,0);
				p=hdist->getDistPopular(a->res->tres,a->bond[0]->name,b->bond[0]->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->bond[0]->name,a->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->bond[0]->name,a->bond[0]->name,a->res->tres->name);
				if(p&&a->tatm->id!=0&&b->tatm->id!=0) {
					x=p->findmostpopulardistance();
					bound->addbounds(a->bond[0]->id0,b->bond[0]->id0,x,p->dist[0],p->dist[1]+bvarn,1);
				}			
			} 	
			else {				
				//p=hdist->getDistPopular(a->name,b->name,0);
				p=hdist->getDistPopular(a->res->tres,a->name,b->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->name,a->name,a->res->tres->name); 
				if(p) {
					x=p->findmostpopulardistance();
					float y=TRES.distance(a,b);
					if(strchr(boundtag,'@')) y=p->dist[1]; 
					float z=max(p->dist[1],y)+varn;
					bound->addbounds(a->id0,b->id0,x,p->dist[0],z,1);
				}
				//
				if(a->bond[0]==0||b->bond[0]==0) continue;
				//
				//p=hdist->getDistPopular(a->bond[0]->name,b->name,0);
				p=hdist->getDistPopular(a->res->tres,a->bond[0]->name,b->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->name,a->bond[0]->name,a->res->tres->name);
				if(p&&a->tatm->id!=0) {
					x=p->findmostpopulardistance();
					float y=TRES.distance(a->bond[0],b);
					if(strchr(boundtag,'@')) y=p->dist[1]; 
					y=max(p->dist[1],y)+varn;
					bound->addbounds(a->bond[0]->id0,b->id0,x,p->dist[0],y,1);
				}
				//
				//p=hdist->getDistPopular(a->name,b->bond[0]->name,0);
				p=hdist->getDistPopular(a->res->tres,a->name,b->bond[0]->name,b->res->tres->name);
				//if(!p) p=hdist->getDistPopular(b->bond[0]->name,a->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->bond[0]->name,a->name,a->res->tres->name);
				if(p&&b->tatm->id!=0) {
					x=p->findmostpopulardistance();
					float y=TRES.distance(a,b->bond[0]);
					if(strchr(boundtag,'@')) y=p->dist[1]; 
					y=max(p->dist[1],y)+varn;
					bound->addbounds(a->id0,b->bond[0]->id0,x,p->dist[0],y,1);	
				}
				//
				//p=hdist->getDistPopular(a->bond[0]->name,b->bond[0]->name,0);
				p=hdist->getDistPopular(a->bond[0]->res->tres,a->bond[0]->name,b->bond[0]->name,b->res->name);
				//if(!p) p=hdist->getDistPopular(b->name,a->bond[0]->name,0);
				if(!p) p=hdist->getDistPopular(b->res->tres,b->bond[0]->name,a->bond[0]->name,a->res->name);
				if(p&&a->tatm->id!=0&&b->tatm->id!=0) {
					x=p->findmostpopulardistance();
					float y=TRES.distance(a->bond[0],b->bond[0]);
					if(strchr(boundtag,'@')) y=p->dist[1]; 
					y=max(p->dist[1],y)+varn;
					bound->addbounds(a->bond[0]->id0,b->bond[0]->id0,x,p->dist[0],y,1);
				}
			} 			
			m++;	
		}					
	}
	//end

	//distance constraint in itself
	nte=0;
	if(strchr(boundtag,'U')) nte=1;	 
	else nte=0;
	
	for(r=r0;r&&nte;r=r->next) {		
		if(r->id0>end) break;
		 
		Res *t;
		for(t=r->next;t;t=t->next) {
			if(t->id0>end) break;
			
			if(t->id0-r->id0<2) continue; 
			
			for(i=0;i<3;i++) {
				a=r->isatmid(i);
				b=t->isatmid(i);
				float x=TRES.distance(a,b);
				bound->addbounds(a->id0,b->id0,x,x-varn,x+varn,1);
			}
		}		
	}
	//end

	//check for position constraint which comes from segment

	if(strchr(boundtag,'P')) nte=1;	
	else nte=0; 

	for(s=segment->next;s&&nte;s=s->next) {
	 	Chn *c=s->chn;
	    	for(r1=c->res;r1;r1=r1->next) {
			if(r1->id0>end) break;
			r=chn->isres0(r1->id0);
			if(r==0) continue;
 
			for(a=r1->atm;a;a=a->next) {
				if(a->tatm->id!=1&&flg==1) continue;
				else if(a->tatm->id>2&&flg==2) continue;
				else if(a->tatm->id>3&&flg==3) continue;			
				Atm *b=r->isatmid(a->tatm->id);
				if(b==0) continue;			 
				float dd=0.5; 
				if(strstr(boundtag,"LLLLL")) dd=5;
				else if(strstr(boundtag,"LLLL")) dd=4;
				else if(strstr(boundtag,"LLL")) dd=3;
				else if(strstr(boundtag,"LL")) dd=2;
				else if(strstr(boundtag,"L")) dd=1;
			
				bound->addbounds(a->id0,b->id0,0,0,dd,1);
				bound->setpredt(a->id0,0);	
				m++;	
			}
		}					
	}
	//end
 
	if(m==0) {
		delete bound;bound=0;
	}
	else {
		bound->adjust();
		bound->reorder();
		//bound->clearredundancy();
		bound->settag();
		bound->setbuffertag();
	}
	if(bound&&TRES.logg>3) bound->printoutdistance();
	
	
	chn0=chn;
	while(chn0) {
		Chn *s=chn0->next;
		if(s==0) {
			chn0->next=0;
			chn0=s;
		}
		else if(s->pdb==chn0->pdb) {
			chn0=s;
		}
		else {
			chn0->next=0;
			chn0=s;
		}
	}
}
 
void Segen::setbound(int flg,int flg0) {

	if(mutate==0) return;
	StrFmt *parent=mutate->owner->getparentStrFmt();
        StrFmt *se = parent->getsequenceStrFmt();
        if(se==0) return;
 
	int i;
	int m=0;
	Chn *chn=(*pdb)[cid];
		
	chn->next=mutate->pdbcopy->chn;
	chn->next->next=mutate->owner->pdb->chn;
	pdb->configureatmid0();
	bound=new Bound;
	bound->model=pdb;
	bound->flag|=TRES.constant->pdbundeletable;
	bound->easyready();
	
	int nlen=strlen(mutate->sqnto);

	
	Res *r0=(*chn)[start];
	Res *r;
	Atm *a;
	for(r=r0;r;r=r->next) {
		if(r->id0>end) break;
		i=mutate->compare[r->id0];
		 
		
		if(mutate->owner->seqngap[i]=='-'||se->seqngap[i]=='-') continue;
		if(mutate->rely[i]==0) continue;		

		int j=mutate->owner->match[i];
		Res *r1;
		
		if(mutate->refine==1) {
			if(mutate->owner->seqngap[i]!='-'&&se->seqngap[i]!='-') {
				r1=mutate->owner->resn[j];
			}
			else if(mutate->owner->seqngap[i]=='-'&&se->seqngap[i]!='-') {
				r1=mutate->pdbcopy->chn->isres0(r->id0);
			}
		}
		
		
		if(r1==0) continue;
		for(a=r->atm;a;a=a->next) {
			if(a->tatm->id!=1&&flg==1) continue;
			else if(a->tatm->id>2&&flg==2) continue;
			else if(a->tatm->id>3&&flg==3) continue;			
			Atm *b=r1->isatmid(a->tatm->id);
			if(b==0) continue;			 
			float dd=mutate->boundscore[i];			
			bound->addbounds(a->id0,b->id0,0,0,dd,parent->sitescore[i]);
			bound->setpredt(b->id0,0);	
			m++;	
		}					
	}
	if(m==0) {
		delete bound;bound=0;
	}
	else {
		bound->adjust();
		bound->reorder();
		//bound->clearredundancy();
		bound->settag();
	}
	if(bound) bound->printoutdistance();
	chn->next->next=0;chn->next=0;
}
float **Segen::zipperallpossible(){
  
	 if(mutate->refine==1) {
		if(rotatm) delete [] rotatm; rotatm=0;
	 	if(bound)  delete bound;bound=0;
		setbound(2);
		return zipperboth(discut*(end-start+1)*100);
	 }
	 else if(mutate->refine==2) {
		if(rotatm) delete [] rotatm; rotatm=0;
	 	if(bound)  delete bound;bound=0;
		setbound(2);
		return zipperfix(discut*(end-start+1)*100,1);
	 }
	 else if(mutate->refine==3) {
		if(rotatm) delete [] rotatm; rotatm=0;
	 	if(bound)  delete bound;bound=0;
		setbound(2);
		return zipperboth(discut*(end-start+1)*100);
	 }
	 else if(mutate->refine==4) {
		if(rotatm) delete [] rotatm; rotatm=0;
	 	if(bound)  delete bound;bound=0;
		return zipperboth(discut*(end-start+1)*100);
	 }
	 else {
		if(rotatm) delete [] rotatm; rotatm=0;
	 	if(bound)  delete bound;bound=0;
		return zipperboth(discut*(end-start+1)*100);
	 }
	 
}

float **Segen::zipperit(){
  	 
	return zipperboth(discut*(end-start+1)*100);	 	 
}
float **Segen::zipperit(float d){
  	 
	return zipperboth(d);	 	 
}
float **Segen::addup(float **x,float **y) {

	if(x==0) return y;
	if(y==0) return x;
	int n1=0;while(x[n1])n1++;
	int n2=0;while(y[n2])n2++;

	float **z=new float*[n1+n2+10];

	int i=0;
	for(i=0;i<n1;i++) z[i]=x[i];
	int t=i;
	for(i=0;i<n2;i++) z[i+t]=y[i];
	z[n2+n1]=0;
	delete [] x;delete [] y;x=0;y=0;
	return z;
}

float **Segen::addup(float **x,float *y) {
	
	if(y==0) return x;
	float **yy=new float*[2];
	yy[0]=y;
	yy[1]=0;
	return addup(x,yy);
}

float **Segen::indavgminimize(float **xyz_all) {

	Chn *chn=(*pdb)[cid];
	Res *r0=(*chn)[start];	
	float *xyzorg=chn->gettransfer(r0,end);
	float **xyz_temp=new float*[2];
	int kep=0;while(xyz_all[kep])kep++;
	
	int i=0;
	for(i=0;i<kep;i++) {
		xyz_temp[0]=xyz_all[i];
		xyz_temp[1]=0;
		chn->transfer(xyzorg,r0,end,1);
		xyz_temp=avgminimize(xyz_temp);
		xyz_all[i]=xyz_temp[0];
	}
	chn->transfer(xyzorg,r0,end,1);
	if(xyzorg) delete [] xyzorg;xyzorg=0;
	if(xyz_temp) delete [] xyz_temp;xyz_temp=0;
	return xyz_all;
} 

float *Segen::avgminrmsd(float **xyz_all,float *avg) { 

	Chn *chn=(*pdb)[cid];
	Res *r0=(*chn)[start];
	
	float *xyz0=chn->gettransfer(r0,end);
	int kep=0;while(xyz_all[kep])kep++;
	int i;
	float a=1000000;
	int  n=0;
	chn->transfer(avg,r0,end,1);
	for(i=0;i<kep;i++) {
		segment->chn->transfer(xyz_all[i],1);
		float d=rmsd(2);
		if(d<a) {
			a=d;
			n=i;
		}
	}
	if(a>10000) {
		chn->transfer(xyz0,r0,end,1);
		if(xyz0) delete [] xyz0; xyz0=0;
		return avg;
	}
	chn->transfer(xyz_all[n],r0,end,1);
	chn->transfer(avg,r0,end,0);
	chn->transfer(xyz0,r0,end,1);
	if(xyz0) delete [] xyz0; xyz0=0;
	if(TRES.logg) cerr<<"the minimal rmsd is:"<<a<<endl;
	return avg;
}
float **Segen::avgminimize(float **xyz_all) {
 
	Chn *chn; 
  	chn=(*pdb)[cid];
 	Res *r1,*r2;
 
	r1=chn->isres0(start);r2=chn->isres0(end);

	if(strchr(boundtag,'M')) {
		
		int nw=0;
		if(r1&&r1->last&&r1->last->last) nw=1;
	        int nf=0;	//take care of too long segments		 
		for(r1=r1->last;nw&&r1->last;r1=r1->last) {
			nf++;
			if(nf>5) {
			break;
			}
			else {
			if(mutate->rely[r1->id0]!=2) continue;
			Res *re=mutate->getres(r1);
			if(re==0) continue;
			float f=r1->directrmsdanyway(re,0,3);
			if(f<1) break;			
			}
		}
		if(nw==0) r1=chn->isres0(start);
	 	nw=0;
		if(r2&&r2->next&&r2->next->next) nw=1;
		nf=0;       //take care of too long segments
		for(r2=r2->next;nw&&r2->next;r2=r2->next) {
			nf++;
                        if(nf>5) {
                        break;
                        }
                        else {
			if(mutate->rely[r2->id0]!=2) continue;
			Res *re=mutate->getres(r2);
			if(re==0) continue;
			float f=r2->directrmsdanyway(re,0,3);
			if(f<1) break;			
			}
		}
 		if(nw==0) r2=chn->isres0(end);
		if((r2->id0-r1->id0+1)-(end-start+1)<4) {
			if(r2->id0-end<2) {
				r2=chn->isres0(end+2);
				if(r2==0) r2=pdb->chn->lastres();
			}
			if(start-r1->id0<2) {		 		  	
				r1=chn->isres0(start-2);
				if(r1==0) r1=pdb->chn->res;
			}			
		} 	

		if(r1==0) r1=chn->isres0(start);
		if(r2==0) r2=chn->isres0(end);
		int k1=mutate->compare[r1->id0];
		int k2=mutate->compare[r2->id0];
		int kk[2];
		kk[0]=k1;kk[1]=k2;
		if(TRES.logg) mutate->printsegment(kk);
	}
	//end
	
	float *xyzorg=chn->gettransfer(r1,r2->id0);	 	
	xyz_all=getnewxyz(xyz_all,r1,r2->id0);
 		
	start=r1->id0;
	end=r2->id0;	
 	
	readynewfix();
	segment->next=mycreate(r1,r2->id0);	
	float *avg=getaveragepost(xyz_all);
	//avg=avgminrmsd(xyz_all,avg); //find average rmsd close to the average rmsd
	if(mutate->refine<10) { //fix
		segment->next->transfer(avg,1);	 
		chn->transfer(avg,r1,r2->id0,1);
	}	
	else if(mutate->refine>=10&&mutate->refine<=90) { //assemble
		segment->next->transfer(xyzorg,1);	 
		chn->transfer(xyzorg,r1,r2->id0,1);
	}
	else if(mutate->refine>=100&&mutate->refine<=900) { //loopy
		segment->next->transfer(avg,1);	 
		chn->transfer(avg,r1,r2->id0,1);
	}
	else if(mutate->refine>=1000&&mutate->refine<=9000) { //refinement
		segment->next->transfer(mutate->xyzself,1);	 
		chn->transfer(mutate->xyzself,r1,r2->id0,1);
	}
	else {
		cerr<<"refinement not defined"<<endl;
		exit(0);
	} 
	if(TRES.logg>3)pdb->chn->write("avg.pdb");
	chn->transfer(xyzorg,r1,r2->id0,1);
	 
 	delete [] avg;avg=0;

	int nte=0;
	if(strchr(boundtag,'*')) nte=1;
	if(bound==0&&nte==0) setbound(2); 
	else if(bound==0&&nte==1) setbound(1);

	float *ent;
	if(0&&mutate->fapr){
		ent=ensort(xyz_all);if(ent) delete [] ent; ent=0;
		int kep=0;while(xyz_all[kep])kep++;
		clearxyz(xyz_all,max(1,kep/2));
	}
	if(mutate->refine>=10&&mutate->refine<=90&&TRES.logg>3){
		chn->transfer(xyz_all[0],r1,r2->id0,1);
		if(bound)bound->printoutdistance();
		chn->transfer(xyzorg,r1,r2->id0,1);
		if(bound)bound->printoutdistance();
	}
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;
	
	if(bound) delete bound;bound=0;
	chn->transfer(xyzorg,r1,r2->id0,1);
	delete [] xyzorg;xyzorg=0;
	return xyz_all;
}
 

float **Segen::allavgminimize(float **xyz_all) {
 
	Chn *chn; 
  	chn=(*pdb)[cid];
 	Res *r1,*r2;
 
	
	/* 
	//former
	n1=start-2;
	n2=end+2;
	
	r1=chn->isres0(n1);
	r2=chn->isres0(n2);
	if(r1==0) r1=pdb->chn->res;
	if(r2==0) r2=pdb->chn->lastres();
	//end
	*/ 
	 
	//test
	r1=chn->isres0(start);r2=chn->isres0(end);
	if(strchr(boundtag,'M')) {
		for(r1=r1;r1->last;r1=r1->last) {
			if(mutate->rely[r1->id0]==2) break;
		}
	 
		for(r2=r2;r2->next;r2=r2->next) {
			if(mutate->rely[r2->id0]==2) break;
		}

		if((r2->id0-r1->id0+1)-(end-start+1)<4) {
			if(r2->id0-end<2) {
				r2=chn->isres0(end+2);
				if(r2==0) r2=pdb->chn->lastres();
			}
			if(start-r1->id0<2) {		 		  	
				r1=chn->isres0(start-2);
				if(r1==0) r1=pdb->chn->res;
			}			
		} 	
	}
	 
	//end
	
	float *xyzorg=pdb->chn->gettransfer(r1,r2->id0);	 	
	xyz_all=getnewxyz(xyz_all,r1,r2->id0);
 		
	start=r1->id0;
	end=r2->id0;	
 	
	readynewfix();
	Pdb *sg=segment;
	int kep=0;while(xyz_all[kep])kep++;
	int i;
	for(i=0;i<kep;i++) {
		sg->next=mycreate(r1,r2->id0);		
		sg=sg->next;
		sg->transfer(xyz_all[i],1);
	}
	//segment->next=mycreate(r1,r2->id0);	
	float *avg=getaveragepost(xyz_all);
	segment->next->transfer(avg,1);
	 
	pdb->chn->transfer(avg,r1,r2->id0,1);
	if(TRES.logg>3)pdb->chn->write("avg.pdb");
	pdb->chn->transfer(xyzorg,r1,r2->id0,1);
	 
 	delete [] avg;avg=0;
	if(bound==0) setbound(2); 
	float *ent;
	//ent=ensort(xyz_all,"v");if(ent) delete [] ent; ent=0;
	//int kep;while(xyz_all[kep])kep++;
	//clearxyz(xyz_all,max(1,kep/2));
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;
	
	if(bound) delete bound;bound=0;
	pdb->chn->transfer(xyzorg,r1,r2->id0,1);
	delete [] xyzorg;xyzorg=0;
	return xyz_all;
}


void Segen::clearxyz(float **xyz_all,int n) {

	int kep=0;while(xyz_all[kep])kep++;
	for(int i=n;i<kep;i++) {
		delete [] xyz_all[i];
		xyz_all[i]=0;
	}
}

float **Segen::create2d(int num) {

	Chn *chn;
  	chn=(*pdb)[cid];
	Res *r0=(*chn)[start];
	Res *r,*t;
	float **xyz_all=0;
	int arbt0=arbt;
	//arbt=num;
	int databaseonly0=databaseonly;
	databaseonly=1;
	
	int tote=0;
 
	int rd=random();
	//helix
	int tt=1;
	int nj=end-start+1;
	 
	if(chiangle) delete chiangle;chiangle=0;
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;	
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r;t;t=t->next) {
				if(t->id0>end) break;
				if(t->id0+tt>end) break;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%15;
				if(fo%2==1) fo=-fo;		
				if(t->name!='P') ang[k][0]=-60+fo;
				fo=random()%15;
				if(fo%2==1) fo=-fo;
				ang[k][1]=-43+fo;			 
				gott++;got++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		tt++;
	}
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	//if(chiangle) xyz_all=zipperfix(rd,100,1);
		

	//sheet
	tt=1;
	nj=end-start+1;
	//if(chiangle) delete chiangle;chiangle=0;
	tote=0;
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r;t;t=t->next) {
				if(t->id0>end) break;
				if(t->id0+tt>end) break;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%30;
				if(fo%2==1) fo=-fo;			
				if(t->name!='P') ang[k][0]=-110+fo;
				fo=random()%30;
				if(fo%2==1) fo=-fo;
				ang[k][1]=130+fo;			 
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
			
		}
		if(got==0) break;
		tt++;
	}
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	//float **tf=0;
	//if(chiangle) tf=zipperfix(rd,100,1);
	//xyz_all=addup(xyz_all,tf);
	//itself

	tt=1;
	nj=end-start+1;
	//if(chiangle) delete chiangle;chiangle=0;
	tote=0;
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			if(tote>arbt) break;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r0;t;t=t->next) {
				if(t->id0>end) break;
				//if(t->id0+tt>end) break;
				int k=t->id0-start;	
				if(k>=nj) break;	
				Atm *a=(*chn)[t->id0]->isatmid(1);
				int fo=random()%10;
				if(fo%2==1) fo=-fo;	
				if(a) {
					ang[k][0]=a->chi+fo;	
				}
				a=(*chn)[t->id0]->isatmid(2);
				fo=random()%30;
				if(fo%2==1) fo=-fo;
				if(a) {
					ang[k][1]=a->chi+fo;
				}
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		if(tote>arbt) break;
		tt++;
	}
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	//if(chiangle) tf=zipperfix(rd,100,1);
	 
	if(chiangle) 
	{
		arbt=chiangle->getsize();
		arbt=max(3,arbt);
		databaseonly=1;
		int ranseed0=ranseed;
		ranseed=random();
		xyz_all=zipperfix(100,1);
		ranseed=ranseed0;
	}
	else xyz_all=0;
	
	arbt=arbt0;
	databaseonly=databaseonly0;
	if(chiangle) delete chiangle;chiangle=0;
	return xyz_all;
}

float **Segen::create3d(int num) {

	Chn *chn;
  	chn=(*pdb)[cid];
	Res *r0=(*chn)[start];
	Res *r,*t;
	float **xyz_all=0;
	int arbt0=arbt;
	//arbt=num;
	int databaseonly0=databaseonly;
	databaseonly=1;
	float *xyzorg=chn->gettransfer(r0,end);
	int tote=0;
 
	int rd=random();
	//helix
	int tt=1;
	int nj=end-start+1;
	chn->transfer(xyzorg,r0,end,1); 
	if(chiangle) delete chiangle;chiangle=0;
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;	
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r;t;t=t->next) {
				if(t->id0>end) break;
				if(t->id0+tt>end) break;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%15;
				if(fo%2==1) fo=-fo;		
				if(t->name!='P') ang[k][0]=-60+fo;
				fo=random()%15;
				if(fo%2==1) fo=-fo;
				ang[k][1]=-43+fo;			 
				gott++;got++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		tt++;
	}
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
 	chn->transfer(xyzorg,r0,end,1);
	//sheet
	tt=1;
	nj=end-start+1;
	chn->transfer(xyzorg,r0,end,1); 
	tote=0;
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r;t;t=t->next) {
				if(t->id0>end) break;
				if(t->id0+tt>end) break;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%30;
				if(fo%2==1) fo=-fo;			
				if(t->name!='P') ang[k][0]=-110+fo;
				fo=random()%30;
				if(fo%2==1) fo=-fo;
				ang[k][1]=130+fo;			 
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
			
		}
		if(got==0) break;
		tt++;
	}
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	
	chn->transfer(xyzorg,r0,end,1);
	tt=1;
	nj=end-start+1;
	//if(chiangle) delete chiangle;chiangle=0;
	tote=0;
	chn->transfer(xyzorg,r0,end,1);
	chn->dihedral(r0,end);
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			if(tote>arbt) break;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r0;t;t=t->next) {
				if(t->id0>end) break;				
				int k=t->id0-start;	
				if(k>=nj) break;	
				Atm *a=(*chn)[t->id0]->isatmid(1);
				int fo=random()%10;
				if(tote>arbt/2) fo=random()%5;
				if(fo%2==1) fo=-fo;	
				if(a) {
					ang[k][0]=a->chi+fo;	
				}
				a=(*chn)[t->id0]->isatmid(2);
				fo=random()%10;
				if(tote>arbt/2) fo=random()%5;
				if(fo%2==1) fo=-fo;
				if(a) {
					ang[k][1]=a->chi+fo;
				}
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		if(tote>arbt) break;
		tt++;
	}
	chn->transfer(xyzorg,r0,end,1);
	//chn->dihedral(r0,end);
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	//if(chiangle) tf=zipperfix(rd,100,1);
	 
	if(chiangle) 
	{
		arbt=chiangle->getsize();
		arbt=max(3,arbt);
		databaseonly=1;
		int ranseed0=ranseed;
		ranseed=random();
		xyz_all=zipperfix(100,1);
		ranseed=ranseed0;
	}
	else xyz_all=0;
	chn->transfer(xyzorg,r0,end,1);
	arbt=arbt0;
	databaseonly=databaseonly0;
	if(chiangle) delete chiangle;chiangle=0;
	return xyz_all;
}

float **Segen::create6d(int num) {
	Chn *chn;
  	chn=(*pdb)[cid];
	Res *r0=(*chn)[start];
	Res *r,*t;
	float **xyz_all=0;
	int arbt0=arbt;
	//arbt=num;
	int databaseonly0=databaseonly;
	databaseonly=1;
	float *xyzorg=chn->gettransfer(r0,end);

	int tote=0;
 
	int rd=random();
	//helix
	int tt=1;
	int nj=end-start+1;
	 
	if(chiangle) delete chiangle;chiangle=0;
	chn->transfer(xyzorg,r0,end,1);
	while(0) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;	
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r0;t;t=t->next) {
				if(t->id0>end) break;
				//if(t->id0+tt>end) break;
				if(mutate->rely[t->id0]==2) continue;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%15;
				if(fo%2==1) fo=-fo;	
				int rot=ifrotatmexist(r->isatmid(1));	
				if(t->name!='P'&&rot==0) ang[k][0]=-60+fo;
				fo=random()%15;
				if(fo%2==1) fo=-fo;
				rot=ifrotatmexist(r->isatmid(2));	
				if(rot==0) ang[k][1]=-43+fo;			 
				gott++;got++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		if(tote>arbt) break;
		tt++;
	}
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
 	chn->transfer(xyzorg,r0,end,1);
	//sheet
	tt=1;
	nj=end-start+1;
	chn->transfer(xyzorg,r0,end,1); 
	tote=0;
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r0;t;t=t->next) {
				if(t->id0>end) break;
				//if(t->id0+tt>end) break;
				if(mutate->rely[t->id0]==2) continue;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%15;
				if(fo%2==1) fo=-fo;			
				if(t->name!='P') ang[k][0]=-110+fo;
				fo=random()%15;
				if(fo%2==1) fo=-fo;
				ang[k][1]=130+fo;			 
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
			
		}
		if(got==0) break;
		if(tote>arbt) break;
		tt++;
	}
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	chn->transfer(xyzorg,r0,end,1);

	tt=1;
	nj=end-start+1;
	//if(chiangle) delete chiangle;chiangle=0;
	chn->transfer(xyzorg,r0,end,1);
	chn->dihedral(r0,end);
	tote=0;
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			if(tote>arbt) break;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r0;t;t=t->next) {
				if(t->id0>end) break;				
				int k=t->id0-start;	
				if(k>=nj) break;	
				Atm *a=(*chn)[t->id0]->isatmid(1);
				int fo=random()%5;
				//if(tote>arbt/2) fo=random()%3;
				if(fo%2==1) fo=-fo;	
				if(a) {
					ang[k][0]=a->chi+fo;	
				}
				a=(*chn)[t->id0]->isatmid(2);
				fo=random()%5;
				//if(tote>arbt/2) fo=random()%3;
				if(fo%2==1) fo=-fo;
				
				if(a) {
					ang[k][1]=a->chi+fo;
				}
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		if(tote>arbt) break;
		tt++;
	}
	chn->transfer(xyzorg,r0,end,1); 
	rd=random();
	//chn->dihedral(r0,end);
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	//if(chiangle) tf=zipperfix(rd,100,1);
	 
	if(chiangle) 
	{
		arbt=chiangle->getsize();
		arbt=max(3,arbt);
		databaseonly=1;
		int ranseed0=ranseed;
		ranseed=random();
		xyz_all=zipperfix(100,1);
		ranseed=ranseed0;
	}
	else xyz_all=0;
	chn->transfer(xyzorg,r0,end,1);
	arbt=arbt0;
	databaseonly=databaseonly0;
	if(xyzorg) delete [] xyzorg;
	if(chiangle) delete chiangle;chiangle=0;
	return xyz_all;
}
float **Segen::create5d(int num) {
	Chn *chn;
  	chn=(*pdb)[cid];
	Res *r0=(*chn)[start];
	Res *r,*t;
	float **xyz_all=0;
	int arbt0=arbt;
	//arbt=num;
	int databaseonly0=databaseonly;
	databaseonly=1;
	float *xyzorg=chn->gettransfer(r0,end);	
	
	int tote=0;
 
	int rd=random();
	//helix
	int tt=1;
	int nj=end-start+1;
	 
	if(chiangle) delete chiangle;chiangle=0;
	chn->transfer(xyzorg,r0,end,1);
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;	
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r0;t;t=t->next) {
				if(t->id0>end) break;
				//if(t->id0+tt>end) break;
				if(mutate->rely[t->id0]==2) continue;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%15;
				if(fo%2==1) fo=-fo;	
				int rot=ifrotatmexist(r->isatmid(1));	
				if(t->name!='P'&&rot==0) ang[k][0]=-60+fo;
				fo=random()%15;
				if(fo%2==1) fo=-fo;
				rot=ifrotatmexist(r->isatmid(2));	
				if(rot==0) ang[k][1]=-43+fo;			 
				gott++;got++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		if(tote>arbt) break;
		tt++;
	}
	chn->transfer(xyzorg,r0,end,1);
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
 
	//sheet
	tt=1;
	nj=end-start+1;
	 
	tote=0;
	chn->transfer(xyzorg,r0,end,1);
	while(0) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r;t;t=t->next) {
				if(t->id0>end) break;
				if(t->id0+tt>end) break;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%30;
				if(fo%2==1) fo=-fo;			
				if(t->name!='P') ang[k][0]=-110+fo;
				fo=random()%30;
				if(fo%2==1) fo=-fo;
				ang[k][1]=130+fo;			 
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
			
		}
		if(got==0) break;
		
		tt++;
	}
	chn->transfer(xyzorg,r0,end,1);
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	

	tt=1;
	nj=end-start+1;
	//if(chiangle) delete chiangle;chiangle=0;
	tote=0;
	chn->transfer(xyzorg,r0,end,1);
	chn->dihedral(r0,end);
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			if(tote>arbt) break;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r0;t;t=t->next) {
				if(t->id0>end) break;				
				int k=t->id0-start;	
				if(k>=nj) break;	
				Atm *a=(*chn)[t->id0]->isatmid(1);
				int fo=random()%5;
				//if(tote>arbt/2) fo=random()%3;
				if(fo%2==1) fo=-fo;	
				if(a) {
					ang[k][0]=a->chi+fo;	
				}
				a=(*chn)[t->id0]->isatmid(2);
				fo=random()%5;
				//if(tote>arbt/2) fo=random()%3;
				if(fo%2==1) fo=-fo;
				if(a) {
					ang[k][1]=a->chi+fo;
				}
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		if(tote>arbt) break;
		tt++;
	}
	
	chn->transfer(xyzorg,r0,end,1);
	//chn->dihedral(r0,end);
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	//if(chiangle) tf=zipperfix(rd,100,1);
	 
	if(chiangle) 
	{
		arbt=chiangle->getsize();
		arbt=max(3,arbt);
		databaseonly=1;
		int ranseed0=ranseed;
		ranseed=random();
		xyz_all=zipperfix(20000,1);
		ranseed=ranseed0;
	}
	else xyz_all=0;
	chn->transfer(xyzorg,r0,end,1);
	arbt=arbt0;
	databaseonly=databaseonly0;
	if(chiangle) delete chiangle;chiangle=0;
	if(xyzorg) delete [] xyzorg;xyzorg=0;
	return xyz_all;
}

float **Segen::create4d(int num) {

	Chn *chn;
  	chn=(*pdb)[cid];
	Res *r0=(*chn)[start];
	Res *r,*t;
	float **xyz_all=0;
	int arbt0=arbt;
	//arbt=num;
	int databaseonly0=databaseonly;
	databaseonly=1;
	
	int tote=0;
 
	int rd=random();
	//helix
	int tt=1;
	int nj=end-start+1;
	 
	if(chiangle) delete chiangle;chiangle=0;
	while(0) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;	
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r;t;t=t->next) {
				if(t->id0>end) break;
				if(t->id0+tt>end) break;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%15;
				if(fo%2==1) fo=-fo;		
				if(t->name!='P') ang[k][0]=-60+fo;
				fo=random()%15;
				if(fo%2==1) fo=-fo;
				ang[k][1]=-43+fo;			 
				gott++;got++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		tt++;
	}
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
 
	//sheet
	tt=1;
	nj=end-start+1;
	 
	tote=0;
	while(0) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r;t;t=t->next) {
				if(t->id0>end) break;
				if(t->id0+tt>end) break;
				int k=t->id0-start;	
				if(k>=nj) break;	
				int fo=random()%30;
				if(fo%2==1) fo=-fo;			
				if(t->name!='P') ang[k][0]=-110+fo;
				fo=random()%30;
				if(fo%2==1) fo=-fo;
				ang[k][1]=130+fo;			 
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
			
		}
		if(got==0) break;
		tt++;
	}
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	

	tt=1;
	nj=end-start+1;
	//if(chiangle) delete chiangle;chiangle=0;
	tote=0;
	while(1) {
		int got=0;
		for(r=r0;r;r=r->next) {
			if(r->id0>end) break;
			
			int j;
			tote++;
			if(tote>arbt) break;
			Chiangle *chi=new Chiangle;
			float **ang=new float*[nj+10];
			for(j=0;j<nj+10;j++) ang[j]=0;				
			for(j=0;j<nj;j++) {
				ang[j]=new float[10];
				int k;
				for(k=0;k<10;k++) ang[j][k]=999;
				ang[j][0]=999;
				ang[j][1]=999;
			}	
			int *numa=new int[nj+10];
			for(j=0;j<nj+10;j++) numa[j]=2;
		 	chi->angle=ang;
			chi->numa=numa;
			chi->number=nj;
			int gott=0;
			for(t=r0;t;t=t->next) {
				if(t->id0>end) break;				
				int k=t->id0-start;	
				if(k>=nj) break;	
				Atm *a=(*chn)[t->id0]->isatmid(1);
				int fo=random()%4;
				if(fo%2==1) fo=-fo;	
				if(a) {
					ang[k][0]=a->chi+fo;	
				}
				a=(*chn)[t->id0]->isatmid(2);
				fo=random()%4;
				if(fo%2==1) fo=-fo;
				if(a) {
					ang[k][1]=a->chi+fo;
				}
				got++;gott++;
			}
			if(gott<num) {
				delete chi;chi=0;
			}
			else if(gott&&chiangle) {
				chiangle->add(chi);chi=0;
			}
			else if(gott&&chiangle==0) {
				chiangle=chi;chi=0;
			}
			else {
				delete chi;chi=0;	
			}
		}
		if(got==0) break;
		if(tote>arbt) break;
		tt++;
	}
	rd=random();
	if(TRES.logg) cerr<<"the total number of predefined conformation:" <<chiangle->getsize()<<endl;
	//if(chiangle) tf=zipperfix(rd,100,1);
	 
	if(chiangle) 
	{
		arbt=chiangle->getsize();
		arbt=max(3,arbt);
		databaseonly=1;
		int ranseed0=ranseed;
		ranseed=random();
		xyz_all=zipperfix(100,1);
		ranseed=ranseed0;
	}
	else xyz_all=0;
	
	arbt=arbt0;
	databaseonly=databaseonly0;
	if(chiangle) delete chiangle;chiangle=0;
	return xyz_all;
}

float **Segen::cutlargermsd(float **xyz_all,float rcut) {

	if(xyz_all==0) return 0;

	Chn *chn=(*pdb)[cid];
	Res *r0=chn->isres0(start);
	float *xyz0=chn->gettransfer(r0,end);
	int kep=0;while(xyz_all[kep])kep++;	
	int i;
	for(i=0;i<kep;i++) {
		segment->chn->transfer(xyz_all[i],1);
		float d=rmsdmutate(2);
		if(d>rcut) {
			delete [] xyz_all[i]; xyz_all[i]=0;
		}
		if(TRES.logg>3) cerr<<i<<" "<<"the rmsd from the original: "<<d<<endl;
	}

	int m=0;
	for(i=0;i<kep;i++) {
		if(xyz_all[i]==0) continue;
		xyz_all[m++]=xyz_all[i];
	}
	xyz_all[m]=0;
	chn->transfer(xyz0,r0,end,1);
	if(xyz0) delete [] xyz0;xyz0=0;
	return xyz_all;
}

float **Segen::cutlargermsd(float **xyz_all,float rcut,int ntot) {

	if(xyz_all==0) return 0;

	Chn *chn=(*pdb)[cid];
	Res *r0=chn->isres0(start);
	float *xyz0=chn->gettransfer(r0,end);
	int kep=0;while(xyz_all[kep])kep++;	
	float *temp=new float[kep+100];
	int   *order=new int[kep+100];
	
	int i;
	for(i=0;i<kep;i++) {
		segment->chn->transfer(xyz_all[i],1);
		float d=rmsdmutate(2);
		temp[i]=d;		
	}

	Qsort cc;
	cc.sort(temp,kep,order);
	
	for(i=ntot;i<kep;i++) {
		int j=order[i];
		if(xyz_all[j]) {
			delete [] xyz_all[j];
			xyz_all[j]=0;
		}  
		if(TRES.logg>3) cerr<<i<<" "<<ntot<<" "<<"the rmsd from the original: "<<temp[i]<<endl;
	}

	int m=0;
	for(i=0;i<kep;i++) {
		
		if(xyz_all[i]==0) continue;
		xyz_all[m++]=xyz_all[i];
	}
	xyz_all[m]=0;
	if(temp) delete [] temp;temp=0;
	if(order) delete [] order;order=0;
	chn->transfer(xyz0,r0,end,1);
	if(xyz0) delete [] xyz0;xyz0=0;
	return xyz_all;
}

int Segen::iscompositesegment() {
 
	int n1=mutate->compare[start];
	int n2=mutate->compare[end];
	StrFmt *cp=mutate->owner->getStrFmt("composite");
	StrFmt *se=mutate->owner->getStrFmt("sequence");
	StrFmt *parent=mutate->owner->getparentStrFmt();
	if(parent==0||parent->iscom==0) return 0;
	if(cp==0||se==0) {		 
		return 0;
	}
 
	int i;
	for(i=n1;i<=n2;i++) {
		if(se->seqngap[i]=='-') continue;
		if(cp->seqngap[i]!=se->zipcode) return 1;
	} 
	return 0;
}

void Segen::setboundtag(int f) {
 
	if(f==0)  {
		strcpy(boundtag," ");
		return;
	}
	else if(f==1) {
		strcpy(boundtag,"#H@"); //only restraints from backbone and standard
		return;
	}
	else if(f==2) { //composite case
		int n1=mutate->compare[start];
		int n2=mutate->compare[end];
		StrFmt *cp=mutate->owner->getStrFmt("composite");
		StrFmt *se=mutate->owner->getStrFmt("sequence");
		if(cp==0||se==0) {
			strcpy(boundtag,"#H@");	
			return;
		}
		int n=0;
		int i;
		for(i=n1;i<=n2;i++) {
			if(cp->seqngap[i]!=se->zipcode) n++;
		} 
		if(n==0) strcpy(boundtag,"#H@");
		else 	 strcpy(boundtag,"#H@");
		return;
	}
	else {
		strcpy(boundtag,"#H@");	
		return;
	}

}

float **Segen::predtfix()
{ 
  Res *r,*r0;//,*rr0; 
  int nseq,i,kep;
  float **xyz_all,*ent;
  Chn *chn; 
  chn=(*pdb)[cid];

  r0=(*chn)[start];
  //rr0=(*chn)[end];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
  
  float *xyzorg = new float[nseq];
  chn->transfer(xyzorg,r0,end,0);
   
  step=0.17;numclose=15;cycle=50;flat=1;
  //srandom(ranseed);
  if(mutate->refine==0) {//artifical assembly
	randcoil=0;
	int n=random();	 
	xyz_all=zipperit();
  	if(xyz_all==0) goto re200; 
  }
  else if(mutate->refine==100) {//loopy   
	int n1=mutate->calcinsert(r0,end);
	int n2=mutate->calcdelete(r0,end); 
	int arbt0=arbt;
	if(n1==0&&n2==0) {
		arbt=max(1,arbt0/3);
	}	
	else {
		arbt=max(1,arbt0/2);
	}     
	//start	 
	randcoil=2;	 
        xyz_all=zipperit();
  	if(xyz_all==0) {arbt=arbt0; goto re200;} 
	//end
 	
	if(n1==0&&n2==0) {
		arbt=max(1,2*arbt0/3);
	}
	else {
		arbt=max(1,arbt0/2);
	}
	chn->transfer(xyzorg,r0,end,1);         	 
	float **t= create3d(3);
	xyz_all=addup(xyz_all,t);
	t=0;
	arbt=arbt0;
  }  
  else if(mutate->refine==200) {//loopy
	int ranseed0=ranseed;
        int arbt0=arbt;
	Res *rend=(*chn)[end];
        if(direct==1&&rend->next) {
		arbt=arbt/4;
        }
	else arbt=arbt/3;

	//start
	randcoil=2;
	ranseed=random();
        xyz_all=zipperit();
  	if(xyz_all==0) goto re200; 
	//end
	 
	//start
	ranseed=random();
	float **t=zipperfix(100,1);
	xyz_all=addup(xyz_all,t);	
	//end

	//start
	ranseed=random();	
	randcoil=0;
	t=zipperit();
	xyz_all=addup(xyz_all,t);	
	//end
 
	if(direct==1&&rend->next) {
		randcoil=2;
		direct=0;ranseed=random();
		float **t=zipperfix(100,1);
		xyz_all=addup(xyz_all,t);
		direct=1;
	} 
 
	chn->transfer(xyzorg,r0,end,1);         	 
	randcoil=2;
	t= create2d(3);
	xyz_all=addup(xyz_all,t);
	arbt=arbt0;
	ranseed=ranseed0;
  } 
  else if(mutate->refine==1000) {//helix refinement
	
	//int ranseed0=ranseed;
        int arbt0=arbt;
	//Res *rend=(*chn)[end];
        
	//start
	arbt=arbt0;
	randcoil=0;
	ranseed=random();
	arbt=max(1,arbt0/3);
        xyz_all=zipperfix(100,1);
  	if(xyz_all==0) goto re200; 
	//end
 
	randcoil=0;
	arbt=max(1,2*arbt0/3);
	chn->transfer(xyzorg,r0,end,1);         	 
	segment->chn->transfer(xyzorg,1);
	float **t= create5d(3);
	xyz_all=addup(xyz_all,t);	 
	t=0;
	
	arbt=arbt0;	 
  }
  else if(mutate->refine==2000) {//helix refinement
	
	//int ranseed0=ranseed;
        int arbt0=arbt;
	//Res *rend=(*chn)[end];
        
	//start
	arbt=arbt0;
	randcoil=0;
	ranseed=random();
	arbt=max(1,arbt0/3);
        xyz_all=zipperit();
  	if(xyz_all==0) goto re200; 
	//end
 
	randcoil=0;
	arbt=max(1,2*arbt0/3);
	chn->transfer(xyzorg,r0,end,1);         	 
	segment->chn->transfer(xyzorg,1);
	float **t= create6d(3);
	xyz_all=addup(xyz_all,t);	 
	t=0;
	
	arbt=arbt0;	  
       
  }
  else if(mutate->refine==2002) {//sheet refinement

	int ranseed0=ranseed;
        int arbt0=arbt;
	Res *rend=(*chn)[end];
        
	//start
	arbt=arbt0;
	randcoil=0;
	ranseed=random();
        xyz_all=zipperit();
  	if(xyz_all==0) goto re200; 
	//end
  
	arbt=arbt0;
	ranseed=ranseed0;
  }
  else if(mutate->refine==10) {//composite
	randcoil=0;
	xyz_all=zipperfix(100,1);  	 
	Res *rend=(*chn)[end];
	if(direct==1&&rend->next) {
		direct=0;
		//n=random(); 
		float **t=zipperfix(100,1);
		xyz_all=addup(xyz_all,t);
		direct=1;
	} 
	
	chn->transfer(xyzorg,r0,end,1); 
	if(xyz_all==0) goto re200;
  }  
  else if(mutate->refine==20) {//composite
	int n=random();
	xyz_all=zipperit();
  	if(xyz_all==0) goto re200; 
  }

  if(xyz_all==0) goto re200; 
  kep=0;while(xyz_all[kep])kep++;
  if(kep==0) goto re200;

  

  //assemble
  int insert;insert=mutate->calcinsert();
  if(TRES.logg>3) {
	float *xyz0=chn->gettransfer(r0,end);
	strcpy(dforce,"v");  
	ent=ensort(xyz_all); 
  	kep=0;while(xyz_all[kep])kep++;
  	for(i=0;i<kep;i++) {
     		chn->transfer(xyz_all[i],r0,end,1);	
     		char nn[100];
     		if(TRES.logg)cerr<<"segment energy: "<<i<<" "<<ent[i]<<endl;
     		sprintf(nn,"semp0.%i.pdb",i);
     		if(TRES.logg)chn->write(nn);
  	}
	chn->transfer(xyz0,r0,end,1);
	if(xyz0) delete [] xyz0;xyz0=0;
  	if(ent) delete [] ent;ent=0;
  }
  
 
  if(insert>=3&&mutate->refine>=0&&mutate->refine<10&&mutate->fapr!=5) { //builder
	strcpy(boundtag,"@HPML");	 //"HPMR"
	//H: hydrogen bond
	//P: segment
	//M: moving start and end
        //R: hydrogen in itself
	//h: hydrogen in composite
	//U: position in itself
	//#: only backbone hydrogen
	//@: standard constraints
	//*: setbound(1);
	step=0.17;numclose=25;cycle=30;flat=1;
	strcpy(dforce,"udD");
  	chn->transfer(xyzorg,r0,end,1); 
  	xyz_all=avgminimize(xyz_all);			
  	if(bound) delete bound;bound=0;
	r0=(*chn)[start];
	//rr0=(*chn)[end];
        nseq=0;
        for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
        delete [] xyzorg;xyzorg=0;
        xyzorg = new float[nseq];
        chn->transfer(xyzorg,r0,end,0);	
  }
  
  chn->transfer(xyzorg,r0,end,1);

  //remove high energy case
  
  if(mutate->refine>=0&&mutate->refine<10) {
  	strcpy(dforce,"udD");
  	ent=ensort(xyz_all);
  	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
  	for(i=max(1,kep/2);i<kep;i++) {  
     		if(ent[i]>10&&ent[i]-ent[0]>10) {
     			delete [] xyz_all[i];xyz_all[i]=0;
     		}		
		else if(i>5) {
			delete [] xyz_all[i];xyz_all[i]=0;
		}		
  	}
	if(ent) delete [] ent;ent=0;
	
  	int ij=0;
  	for(i=0;i<kep;i++) {
		if(xyz_all[i]==0) continue;
		xyz_all[ij++]=xyz_all[i];
  	}
  	xyz_all[ij]=0;
  	chn->transfer(xyzorg,r0,end,1);    

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;	
  }

  //composite
  //remove bad candidates first
  if(mutate->refine>=10&&mutate->refine<=90) {
 
	//remove high energy candidates
	chn->transfer(xyzorg,r0,end,1); 
	strcpy(dforce,"v");  
	ent=ensort(xyz_all); if(ent) delete [] ent;ent=0;
  	kep=0;while(xyz_all[kep])kep++;	
	int ii=min(30,kep/2+1);
	for(i=ii;i<kep;i++) {
		delete [] xyz_all[i];xyz_all[i]=0;
	}
	chn->transfer(xyzorg,r0,end,1); 

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
 
	//remove too close candidates
	//float avgrmsd=calcavgrmsd(xyz_all,20);
	int pn=part;part=0;
	noclose(xyz_all,0.25);
	part=pn;
	kep=0;while(xyz_all[kep])kep++;
  	chn->transfer(xyzorg,r0,end,1); 

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;

	//remove too many candidates
	ii=max(20,kep/2+1);
	for(i=ii;i<kep;i++) {
		delete [] xyz_all[i];xyz_all[i]=0;
	}
	chn->transfer(xyzorg,r0,end,1); 
	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
  }

  //constrained minimization on loops
  if(mutate->refine>=10&&mutate->refine<=90) {
 
	step=0.17;numclose=50;cycle=50;flat=111; 
	//minimize again
	if(bound) delete bound;bound=0;  
	chn->transfer(xyzorg,r0,end,1); 
	strcpy(dforce,"udD"); 
	strcpy(boundtag,"PLL");
	segment->next=mycreate(r0,end);
	segment->next->transfer(xyzsave,1);
	setbound(2);
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;	
	chn->transfer(xyzorg,r0,end,1); 
	if(bound) delete bound;bound=0; 
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;	
	//remove so that the number is no more 10
	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
  	for(i=max(1,kep/2);i<kep;i++) {  
     		if(i>5) {
			delete [] xyz_all[i];xyz_all[i]=0;
		}		
  	}		    
  }

  
  //loop refinement
  //remove bad candidates first
  if(mutate->refine>=100&&mutate->refine<=900) {

	//remove those larger rmsd
	r0=(*chn)[start]; 
	chn->transfer(mutate->xyzself,r0,end,1);
	xyz_all=cutlargermsd(xyz_all,mutate->rmsd); 
	chn->transfer(xyzorg,r0,end,1);
	
	kep=0;while(xyz_all&&xyz_all[kep])kep++;	
	if(kep==0)  goto re200;
 
	//remove high energy candidates
	chn->transfer(xyzorg,r0,end,1); 
	strcpy(dforce,"v");  
	ent=ensort(xyz_all); if(ent) delete [] ent;ent=0;
  	kep=0;while(xyz_all[kep])kep++;	
	int ii=min(30,kep/2+1);
	for(i=ii;i<kep;i++) {
		delete [] xyz_all[i];xyz_all[i]=0;
	}
	chn->transfer(xyzorg,r0,end,1); 

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
 
	//remove too close candidates
	//float avgrmsd=calcavgrmsd(xyz_all,20);
	int pn=part;part=0;
	noclose(xyz_all,0.5);
	part=pn;
	kep=0;while(xyz_all[kep])kep++;
  	chn->transfer(xyzorg,r0,end,1); 

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;

	//remove too many candidates
	ii=max(20,kep/2+1);
	for(i=ii;i<kep;i++) {
		delete [] xyz_all[i];xyz_all[i]=0;
	}
	chn->transfer(xyzorg,r0,end,1); 
	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
  }

  //add original back to list.
  if(mutate->refine>=100&&mutate->refine<=900&&xyz_all) {
	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
	float **tmp=new float*[kep+100];
	int ie;
	for(ie=0;ie<kep+100;ie++) tmp[ie]=0;
	for(ie=0;ie<kep;ie++) {
		tmp[ie]=xyz_all[ie];xyz_all[ie]=0;
	}
	chn->transfer(xyzsave,r0,end,1); 
	tmp[kep]= chn->gettransfer(r0,end);
	tmp[kep+1]=0;	 
	delete [] xyz_all;
	xyz_all=tmp;
  }


  //constrained minimization on loops
  if(mutate->refine>=100&&mutate->refine<=900) {
	step=0.17;numclose=50;cycle=20;flat=1; 	 

	//minimize for hbonds


	//only restraints from backbone and standard
	int iscs=iscompositesegment();
	if(mutate->restraint==1) {
		strcpy(boundtag,"#H@");
	}
	else if(mutate->restraint==2&&iscs==0) {
		strcpy(boundtag,"#H@");
	}
	else if(mutate->restraint==2&&iscs==1) {
		strcpy(boundtag,"#R@");
		//modified, take care of bus error in composite refinement. 07/11/2002
		//chn->transfer(xyzorg,r0,end,1); 
		//mutate->buildhbond();
		//end
	}
	else {
		strcpy(boundtag," ");
	}
	  
	chn->transfer(xyzorg,r0,end,1); 
	strcpy(dforce,"udD"); 
	setbound(2);
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;
	if(bound) delete bound;bound=0; 
	chn->transfer(xyzorg,r0,end,1); 	
	
	//remove larger rmsd
  	r0=(*chn)[start]; 	 
	chn->transfer(mutate->xyzself,r0,end,1);
	xyz_all=cutlargermsd(xyz_all,mutate->rmsd); 
	chn->transfer(xyzorg,r0,end,1);	 

	kep=0;while(xyz_all&&xyz_all[kep])kep++;	
	if(kep==0) goto re200;
 
	//cross	
	int ic=0;

	//repeat!!!
	for(ic=0;ic<2;ic++) {//repeat	
	//repeat!!!

	float avgrmsd=calcavgrmsd(xyz_all,20); 
	xyz_all=cross(xyz_all,1,avgrmsd/2);
	chn->transfer(xyzorg,r0,end,1);	 
	 
	//ensort
	strcpy(dforce,"udD"); 
	ent=ensort(xyz_all);
	chn->transfer(xyzorg,r0,end,1);	 

	//remove larger rmsd
  	r0=(*chn)[start]; 	 
	chn->transfer(mutate->xyzself,r0,end,1);
	xyz_all=cutlargermsd(xyz_all,mutate->rmsd); 
	chn->transfer(xyzorg,r0,end,1);

 	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;

	//remove too close candidates
	//float avgrmsd=calcavgrmsd(xyz_all,20);
	int pn=part;part=0;
	noclose(xyz_all,0.5);
	part=pn;
	kep=0;while(xyz_all[kep])kep++;
  	chn->transfer(xyzorg,r0,end,1); 

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;


	//remove higher energy
	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
  	for(i=max(1,kep/2);i<kep;i++) {  
     		if(ent[i]>10&&ent[i]-ent[0]>10) {
     			delete [] xyz_all[i];xyz_all[i]=0;
     		}		
		else if(i>20) {
			delete [] xyz_all[i];xyz_all[i]=0;
		}		
  	}
	if(ent) delete [] ent;ent=0;
	
  	int ij=0;
  	for(i=0;i<kep;i++) {
		if(xyz_all[i]==0) continue;
		xyz_all[ij++]=xyz_all[i];
  	}
  	xyz_all[ij]=0;
  	chn->transfer(xyzorg,r0,end,1);    

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;

	//end repeat!!!		
	}
	//end repeat!!!	

	//remove so that the number is no more 10
	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
  	for(i=max(1,kep/2);i<kep;i++) {  
     		if(i>5) {
			delete [] xyz_all[i];xyz_all[i]=0;
		}		
  	}
	
	step=0.17;numclose=50;cycle=50;flat=111; 
	//minimize again
	if(bound) delete bound;bound=0;  
	chn->transfer(xyzorg,r0,end,1); 
	strcpy(dforce,"udD"); 
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;	
	chn->transfer(xyzorg,r0,end,1); 

		    
  }

  //secondary minimization
  //remove bad candidates first
  if(mutate->refine>=1000&&mutate->refine<=9000) {
	step=0.17;numclose=50;cycle=20;flat=111;
	//remove those larger rmsd
	r0=(*chn)[start]; 
	chn->transfer(mutate->xyzself,r0,end,1);
	xyz_all=cutlargermsd(xyz_all,mutate->rmsd); 
	chn->transfer(xyzorg,r0,end,1);
	
	kep=0;while(xyz_all&&xyz_all[kep])kep++;	
	if(kep==0)  goto re200;

	//remove high energy candidates
	chn->transfer(xyzorg,r0,end,1); 
	strcpy(dforce,"v");  
	ent=ensort(xyz_all); if(ent) delete [] ent;ent=0;
  	kep=0;while(xyz_all[kep])kep++;	
	int ii=min(30,kep/2+1);
	for(i=ii;i<kep;i++) {
		delete [] xyz_all[i];xyz_all[i]=0;
	}
	chn->transfer(xyzorg,r0,end,1); 

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
 
	//remove too close candidates
	//float avgrmsd=calcavgrmsd(xyz_all,20);
	int pn=part;part=0;
	noclose(xyz_all,0.5);
	part=pn;
	kep=0;while(xyz_all[kep])kep++;
  	chn->transfer(xyzorg,r0,end,1); 

	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;

	//remove too many candidates
	ii=max(20,kep/2+1);
	for(i=ii;i<kep;i++) {
		delete [] xyz_all[i];xyz_all[i]=0;
	}
	chn->transfer(xyzorg,r0,end,1); 
	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
  }

  //add original back to list.
  if(mutate->refine>=1000&&mutate->refine<=9000&&xyz_all) {
	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
	float **tmp=new float*[kep+100];
	int ie;
	for(ie=0;ie<kep+100;ie++) tmp[ie]=0;
	for(ie=0;ie<kep;ie++) {
		tmp[ie]=xyz_all[ie];xyz_all[ie]=0;
	}
	chn->transfer(xyzsave,r0,end,1); 
	tmp[kep]= chn->gettransfer(r0,end);
	tmp[kep+1]=0;	 
	delete [] xyz_all;
	xyz_all=tmp;
  }


  //unconstrained minimization on secondary 
  if(mutate->refine>=1000&&mutate->refine<=9000) {
	step=0.17;numclose=50;cycle=20;flat=111; 	 

	//minimize for hbonds
	//only restraints from backbone and standard
	int iscs=iscompositesegment();
	if(mutate->restraint==1) {
		strcpy(boundtag,"#H@");
	}
	else if(mutate->restraint==2&&iscs==0) {
		strcpy(boundtag,"#H@");
	}
	else if(mutate->restraint==2&&iscs==1) {
		strcpy(boundtag,"#R@");
		//modified, take care of bus error in composite refinement. 07/11/2002
		//chn->transfer(xyzorg,r0,end,1); 
		//mutate->buildhbond();
		//end
	}
	else {
		strcpy(boundtag," ");
	}
	chn->transfer(xyzorg,r0,end,1); 
	strcpy(dforce,"udD"); 
	setbound(2);
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;
	if(bound) delete bound;bound=0; 
	chn->transfer(xyzorg,r0,end,1); 	
	
	//remove larger rmsd
  	r0=(*chn)[start]; 	 
	chn->transfer(mutate->xyzself,r0,end,1);
	xyz_all=cutlargermsd(xyz_all,mutate->rmsd); 
	chn->transfer(xyzorg,r0,end,1);	 

	kep=0;while(xyz_all&&xyz_all[kep])kep++;	
	if(kep==0) goto re200;
 
	//remove so that the number is no more 10
	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
  	for(i=max(1,kep/2);i<kep;i++) {  
     		if(i>5) {
			delete [] xyz_all[i];xyz_all[i]=0;
		}		
  	}
	
	step=0.17;numclose=50;cycle=50;flat=111; 
	//minimize again
	if(bound) delete bound;bound=0;  
	chn->transfer(xyzorg,r0,end,1); 
	strcpy(dforce,"udD"); 
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;	
	chn->transfer(xyzorg,r0,end,1); 		    
  }

  
  //debug
  if(TRES.logg>3) {
	chn->transfer(xyzorg,r0,end,1);
  	strcpy(dforce,"udD");
  	ent=ensort(xyz_all);//ent=ensortfix(xyz_all,"ucdD");

  	kep=0;while(xyz_all[kep])kep++;
  	for(i=0;i<kep;i++) {
     		float *xyz0=chn->gettransfer(r0,end);
     		chn->transfer(xyz_all[i],r0,end,1);	
     		char nn[100];
     		if(TRES.logg)cerr<<"segment energy: "<<i<<" "<<ent[i]<<endl;
     		sprintf(nn,"semp.%i.pdb",i);
     		if(TRES.logg)chn->write(nn);
     		chn->transfer(xyz0,r0,end,1);
     		delete [] xyz0;xyz0=0;
  	}  
  	if(ent) delete [] ent;ent=0;
	chn->transfer(xyzorg,r0,end,1);
  }
  
  

  
  hooksidechain();
  pdb->configure();

  if(xyzsave) { 
  	delete [] xyzorg;
	chn->transfer(xyzsave,r0,end,1);
	xyzorg=chn->gettransfer(r0,end); 
  }
  setallnear();
 
  strcpy(dforce,"u");
   
  kep=0; while(xyz_all&&xyz_all[kep]) kep++;
  if(kep>100) {
	int ip;
	for(ip=100;ip<kep;ip++) {
		if(xyz_all[ip]) delete [] xyz_all[ip];
		xyz_all[ip]=0;
	}
  }
  kep=0; while(xyz_all&&xyz_all[kep]) kep++;
  if(kep==0) {
	if(xyz_all) delete [] xyz_all;xyz_all=0;
	goto re200; 
  }

  
  ent=scap(xyz_all); if(ent) delete [] ent; ent=0;
  
   //minimize again only on sidechain
  if(mutate->refine>=0&&mutate->refine<10) {
	step=0.17/2;numclose=50;cycle=20;flat=111;
	onlysidechain=1;
	strcpy(dforce,"udD");
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;	
	onlysidechain=0;
	chn->transfer(xyzorg,r0,end,1);	

 	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
  }

   //minimize again only on sidechain
  if(mutate->refine>=10&&mutate->refine<90) {
	step=0.17/2;numclose=50;cycle=20;flat=111;
	onlysidechain=1;
	strcpy(dforce,"udD");
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;	
	onlysidechain=0;
	chn->transfer(xyzorg,r0,end,1);	

 	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
  }

  //add original back to list.
  if(mutate->refine>=100&&mutate->refine<=900&&xyz_all) {
	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
	float **tmp=new float*[kep+100];
	int ie;
	for(ie=0;ie<kep+100;ie++) tmp[ie]=0;
	for(ie=0;ie<kep;ie++) {
		tmp[ie]=xyz_all[ie];xyz_all[ie]=0;
	}
	chn->transfer(xyzsave,r0,end,1); 
	tmp[kep]= chn->gettransfer(r0,end);
	tmp[kep+1]=0;	 
	delete [] xyz_all;
	xyz_all=tmp;
  }

  //minimize again only on sidechain
  if(mutate->refine>=100&&mutate->refine<=900) {
	step=0.17/2;numclose=50;cycle=20;flat=111;
	onlysidechain=1;
	strcpy(dforce,"udD");
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;	
	onlysidechain=0;
	chn->transfer(xyzorg,r0,end,1); 
	
	//remove larger rmsd
  	r0=(*chn)[start]; 	 
	chn->transfer(mutate->xyzself,r0,end,1);
	xyz_all=cutlargermsd(xyz_all,mutate->rmsd); 
	chn->transfer(xyzorg,r0,end,1);

 	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
  }
 
   //add original back to list.
  if(mutate->refine>=1000&&mutate->refine<=9000&&xyz_all) {
	kep=0; while(xyz_all&&xyz_all[kep]) kep++;
	float **tmp=new float*[kep+100];
	int ie;
	for(ie=0;ie<kep+100;ie++) tmp[ie]=0;
	for(ie=0;ie<kep;ie++) {
		tmp[ie]=xyz_all[ie];xyz_all[ie]=0;
	}
	chn->transfer(xyzsave,r0,end,1); 
	tmp[kep]= chn->gettransfer(r0,end);
	tmp[kep+1]=0;	 
	delete [] xyz_all;
	xyz_all=tmp;
  }

  //minimize again only on sidechain
  if(mutate->refine>=1000&&mutate->refine<=9000) {
	step=0.17/2;numclose=50;cycle=20;flat=111;
	onlysidechain=1;
	strcpy(dforce,"udD");
	ent=minimizefix(xyz_all,1);if(ent) delete [] ent; ent=0;	
	onlysidechain=0;
	chn->transfer(xyzorg,r0,end,1); 
	
	//remove larger rmsd
  	r0=(*chn)[start]; 	 
	chn->transfer(mutate->xyzself,r0,end,1);
	xyz_all=cutlargermsd(xyz_all,mutate->rmsd); 
	chn->transfer(xyzorg,r0,end,1);

 	kep=0;while(xyz_all[kep])kep++;	
	if(kep==0) goto re200;
  }
 

  //check exist!
  kep=0; while(xyz_all&&xyz_all[kep]) kep++;
  if(kep==0) {
	if(xyz_all) delete [] xyz_all;xyz_all=0;
	goto re200; 
  }
  
  //cluster and sort!
  float cu;cu=cutoff;
  cutoff=10.0;
  strcpy(dforce,"udDW");
  ent=ensort(xyz_all); 
  cutoff=cu;

  kep=0;while(xyz_all[kep])kep++;
 
  //write debug
  for(i=0;i<kep&&TRES.logg>3;i++) {
     float *xyz0=chn->gettransfer(r0,end);
     chn->transfer(xyz_all[i],r0,end,1);	
     char nn[100];
     if(TRES.logg)cerr<<"sidechain segment energy: "<<i<<" "<<ent[i]<<endl;
     sprintf(nn,"temp.%i.pdb",i);
     if(TRES.logg)chn->write(nn);
     chn->transfer(xyz0,r0,end,1);
     if(xyz0) delete [] xyz0;xyz0=0;
  }
  //writeout temporary file
  if(0&&mutate->refine>=10&&mutate->refine<=90&&mutate->out>1&&mutate->code) {

	chn->transfer(xyzorg,r0,end,1); 	 
	int pn=part;part=0;
	noclose(xyz_all,0.5);
	part=pn;
	kep=0;while(xyz_all[kep])kep++;
	if(kep==0) goto re200;

  	chn->transfer(xyzorg,r0,end,1); 	 
	kep=0;while(xyz_all[kep])kep++;	
        segment->chn->transfer(xyzorg,1);
	 
	for(i=0;i<kep;i++) {
		chn->transfer(xyz_all[i],r0,end,1);
		float d=rmsd(2);
		if(d<0.1) continue;
		char nn[100];		
		sprintf(nn,"%s_refine_composite.%i.pdb",mutate->code,mutate->refineid++);    				
		char c=pdb->chn->id;
        	StrFmt *parent=mutate->owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	pdb->chn->id=se->cid;
        	if(pdb->chn->id=='0') pdb->chn->id=' ';
     		pdb->write(nn);		
		pdb->chn->id=c;
 	} 	 
  }
 
  //writeout temporary file
  if(mutate->refine>=100&&mutate->refine<=900&&mutate->out>1&&mutate->code) {

	chn->transfer(xyzorg,r0,end,1); 	 
	int pn=part;part=0;
	noclose(xyz_all,0.5);
	part=pn;
	kep=0;while(xyz_all[kep])kep++;
	if(kep==0) goto re200;

  	chn->transfer(xyzorg,r0,end,1); 	 
	kep=0;while(xyz_all[kep])kep++;	
        segment->chn->transfer(xyzorg,1);
	 
	for(i=0;i<kep;i++) {
		chn->transfer(xyz_all[i],r0,end,1);
		float d=rmsd(2);
		if(d<0.1) continue;
		char nn[100];		
		sprintf(nn,"%s_loop_refine_tmp.%i.pdb",mutate->code,mutate->refineid++);    				
		cerr<<"write down intermediate structures in refinement of loop regions: "<<nn<<endl;
		char c=pdb->chn->id;
        	StrFmt *parent=mutate->owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	pdb->chn->id=se->cid;
        	if(pdb->chn->id=='0') pdb->chn->id=' ';
     		pdb->write(nn);		
		pdb->chn->id=c;
 	} 	 
  }
 
  if(mutate->refine>=1000&&mutate->refine<=9000&&mutate->out>1&&mutate->code) {

	chn->transfer(xyzorg,r0,end,1); 	 
	int pn=part;part=0;
	noclose(xyz_all,0.5);
	part=pn;
	kep=0;while(xyz_all[kep])kep++;
	if(kep==0) goto re200;

  	chn->transfer(xyzorg,r0,end,1); 	 
	kep=0;while(xyz_all[kep])kep++;	
        segment->chn->transfer(xyzorg,1);
	 
	for(i=0;i<kep;i++) {
		chn->transfer(xyz_all[i],r0,end,1);
		float d=rmsd(2);
		if(d<0.1) continue;
		char nn[100];	
		sprintf(nn,"%s_sec_refine_tmp.%i.pdb",mutate->code,mutate->refineid++);   
		cerr<<"write down intermediate structures in refinement of secondary regions: "<<nn<<endl;
		char c=pdb->chn->id;
        	StrFmt *parent=mutate->owner->getparentStrFmt();
        	StrFmt *se = parent->getsequenceStrFmt();
        	pdb->chn->id=se->cid;
        	if(pdb->chn->id=='0') pdb->chn->id=' ';
     		pdb->write(nn);	
		pdb->chn->id=c;	
 	} 	 
  }

  if(ent) delete [] ent;ent=0;
  
  re200:
  chn->transfer(xyzorg,r0,end,1);
  if(xyzorg) delete [] xyzorg;xyzorg=0;
  kep=0;while(xyz_all&&xyz_all[kep])kep++;	
  if(kep==0&&xyz_all) {
	delete [] xyz_all; 
	xyz_all=0;
  }
  return xyz_all;
}

float **Segen::getnewxyz(float **xyz_all,Res *ss,int nn){

   int i;
   Chn *chn;
   chn=(*pdb)[cid];

   int nseq=0; 
   Res *r; 
   for(r=ss;r&&r->id0<=nn;r=r->next) nseq+=r->tres->number*3;

   int kep=0; while(xyz_all[kep])kep++;
   float **out=new float*[kep+10];
   for(i=0;i<kep+10;i++) {
	out[i]=0;
   }    

   Res *r0=chn->isres0(start);
   for(i=0;i<kep;i++) {
	chn->transfer(xyz_all[i],r0,end,1);
	if(TRES.logg>3)chn->write("s");
	out[i]=new float[nseq];
	chn->transfer(out[i],ss,nn,0);
   }
   
   Strhandler cc;
   xyz_all=cc.floatdel(xyz_all) ;
   return out;
}

void Segen::readyminfix(Res *s,int n) {
  
  //create segments
  if(s==0) s=pdb->chn->res;
  Chn *chn=(*pdb)[cid];
 
  start=s->id0;  
  
  Res *r=pdb->chn->isres0(n);
  if(r==0) r=pdb->chn->lastres();

  end=r->id0;
  Atm *aa;
  for(r=s;r;r=r->next) {
	if(r->id0>n) break;
	for(aa=r->atm;aa;aa=aa->next) {
		aa->allnear(3,0);
	}
  }
  
  create();

  //
  //setup disc;
  if(disc) delete disc;disc=0;
  if(disc==0) disc=new Disc;
  disc->setupall(pdb,cutoff,1);
}

float **Segen::minfix(float **xyz_all) {

  Res *r,*r0; 
  int nseq,i,kep;
 
  Chn *chn; 

  chn=(*pdb)[cid];

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
  
  float *xyzorg = new float[nseq];
  chn->transfer(xyzorg,r0,end,0);
   
  step=0.17;numclose=15;cycle=20;
  

  kep=0; while(xyz_all[kep]) kep++;
 
  for(i=0;i<kep;i++) {
	chn->transfer(xyz_all[i],r0,end,1);
        
  } 

  if(xyzorg) delete [] xyzorg;xyzorg=0;

  return xyz_all;
}

float Segen::predtmin()
{ 
  Res *r,*r0; 
  int nseq;
  float ent;
  Chn *chn; 
 
  chn=(*pdb)[cid];

  r0=(*chn)[start];
  nseq=0;
  for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;
  
  float *xyzorg = new float[nseq];
  chn->transfer(xyzorg,r0,end,0);  

  step=0.17/2;numclose=55;cycle=200;flat=111;
  
  //write orginal structure
  char nn[100];
  if(TRES.logg>3) sprintf(nn,"semp.%i.pdb",0);
  if(TRES.logg>3)chn->write(nn);
  
  
  if(mutate->refine>=-10&&mutate->refine<0&&onlysidechain==0) {
	segment->next=mycreate(r0,end);
	segment->next->transfer(mutate->xyzself,1);
	strcpy(boundtag,"#@HPL");
	if(mutate->restraint==0)  strcpy(boundtag," ");
	setbound(2);
  }
  else if(mutate->refine>=-900&&mutate->refine<=-100&&onlysidechain==0) {
	if(rotatm) delete [] rotatm;rotatm=0;
	if(bound) delete bound;bound=0;
  }
  else if(mutate->refine>=-9000&&mutate->refine<=-1000&&onlysidechain==0) {
	if(rotatm) delete [] rotatm;rotatm=0;
	if(bound) delete bound;bound=0;
  }
  else if(mutate->refine>=10&&mutate->refine<=90&&onlysidechain==0) {
	segment->next=mycreate(r0,end);
	strcpy(boundtag,"#@HPLL");
	if(mutate->restraint==0)  strcpy(boundtag," ");
	setbound(2);
  }
  else if(mutate->refine>=100&&mutate->refine<=900&&onlysidechain==0) {
	segment->next=mycreate(r0,end);
	segment->next->transfer(mutate->xyzself,1);
	//only restraints from backbone and standard
	int iscs=iscompositesegment();
	if(mutate->restraint==1) {
		strcpy(boundtag,"#H@PLL");
	}
	else if(mutate->restraint==2&&iscs==0) {
		strcpy(boundtag,"#H@PLL");
	}
	else if(mutate->restraint==2&&iscs==1) {
		strcpy(boundtag,"#R@PLL");
		//modified, take care of bus error in composite refinement. 07/11/2002
		//chn->transfer(xyzorg,r0,end,1); 
		//mutate->buildhbond();
		//end
	}
	else {
		strcpy(boundtag," ");
	}
	setbound(2); 
  }
  else if(mutate->refine>=1000&&mutate->refine<=9000&&onlysidechain==0) {
	segment->next=mycreate(r0,end);
	segment->next->transfer(mutate->xyzself,1);
	int iscs=iscompositesegment();
	if(mutate->restraint==1) {
		strcpy(boundtag,"#H@PLL");
	}
	else if(mutate->restraint==2&&iscs==0) {
		strcpy(boundtag,"#H@PLL");
	}
	else if(mutate->restraint==2&&iscs==1) {
		strcpy(boundtag,"#R@PLL");
		//modified, take care of bus error in composite refinement. 07/11/2002
		//chn->transfer(xyzorg,r0,end,1); 
		//mutate->buildhbond();
		//end
	}
	else {
		strcpy(boundtag," ");
	}
	setbound(2); 
  }   
  else {
        if(rotatm) delete [] rotatm;rotatm=0;
	if(bound) delete bound;bound=0;
  }
   
  if(bound&&TRES.logg>3) mutate->printoutdistance();
  
  strcpy(dforce,"udD");
  ent=minimizefix(xyzorg,1);
  chn->transfer(xyzorg,r0,end,1);  
  if(bound&&TRES.logg>3) mutate->printoutdistance();
  if(bound) delete bound;bound=0;
  if(TRES.logg>3)sprintf(nn,"temp.%i.pdb",0);
  if(TRES.logg>3)chn->write(nn);
  
  if(xyzorg) delete [] xyzorg;xyzorg=0;
  return ent;
}


int Segen::findmidresidue(int a,int b) {

	Chn *chn=(*pdb)[cid];

	Res *r0=(*chn)[start];
	Res *t0=(*chn)[end];

        Res *r=getmidresidue(r0,t0);

        if(r) return r->id0;
        else  return -1;
}

int Segen::findmidresidue(Res *a,Res *b) {

        Res *r=getmidresidue(a,b);

        if(r) return r->id0;
        else  return -1;
}

Res *Segen::getmidresidue(Res *r1,Res *r2) {

        if(r1->last==0||r2->next==0) return 0;

        Chn *c=r1->chn;
        Res *r;
        for(r=r1;r;r=r->next) {
                if(r->id0>=r2->id0) break;
                if(c->islinked(r,r->next)==0) return r;
        }

        int n=(r1->id0+r2->id0)/2;

        for(r=r1;r;r=r->next) {
                if(r->id0>r2->id0) break;
                if(r->id0==n) return r;
        }

        return 0;
}

float** Segen::zipperboth(float dis)
{
  //generate loop conformations
  Segen seg1,seg2;
  Res  *r1,*r2;
  Res  *r0,*r;  
  Atm  *a,*a1;
  int ne,i,n,nseq,m3,m4,m5;
  float d;
  float *xyz,*xyz0,**xyz_all;
  Chn *chn=(*pdb)[cid];
  Chn *cht;

  StrFmt *parent=mutate->owner->getparentStrFmt();
 
  int mid=findmidresidue(start,end);

  if(mid==-1) return zipperfix(dis,1);
  //srandom(seed);
  seg1.pdb=pdb;seg1.cid=cid;
  seg2.pdb=pdb;seg2.cid=cid;
  
  int nj1=mutate->compare[start];
  int nj2=mutate->compare[end];
  if(mid-start+1>end-mid+1) {
     	seg1.start=start;
   	seg1.end=mid;
    	seg2.start=mid;
     	seg2.end=end;
  }
  else if(mid-start+1==end-mid+1){
	if(parent->sitescore[nj1]>=parent->sitescore[nj2]) {
		seg1.start=start;
      		seg1.end=mid+1;
      		seg2.start=mid+1;
    		seg2.end=end;
	}
	else {
		seg1.start=start;
   		seg1.end=mid;
    		seg2.start=mid;
     		seg2.end=end;
	}
  }
  else {
    	seg1.start=start;
      	seg1.end=mid+1;
      	seg2.start=mid+1;
    	seg2.end=end;
  }
   
  
  seg1.create(); 
  seg2.create();
  seg1.segment->chn->setlastresorder();
  seg2.segment->chn->setlastresorder();
  seg1.direct=1;  
  seg2.direct=0;  
 
  r0=(*chn)[start];
  m3=0;m4=0;
  for(r=r0;r;r=r->next) {                         
     if(r->id0<seg2.start) m3+=r->tres->number*3;
     if(r->id0<seg2.start) m4+=9;
  }
  m5=m3+seg2.segment->chn->res->tres->number*3;
  xyz=0;xyz_all=0;

  //coordinate buffer

  nseq=0; 
  
  for(r=r0;r&&r->id0<=end;r=r->next) 
  nseq+=r->tres->number*3;     
  
  xyz=new float[nseq];  
  xyz0=new float[nseq*3];
  int arbt0=arbt;
  if(chiangle) {
	arbt0=max(arbt0,chiangle->getsize());
  }
  xyz_all=new float*[arbt0+10];
  for(i=0;i<arbt0+10;i++) xyz_all[i]=0;

  chn->transfer(xyz,r0,end,0);
  chn->transfertemp(xyz0,r0,end,0); 

  r1=(*seg1.segment->chn)[seg1.end];
  r2=(*seg2.segment->chn)[seg2.start]; 

  //decides the segment to be kept

  //find possible hits

  i=0;n=0;
   
  ne=0;
  
  for(; ;)
  {
    if(flex==0&&i==0&&n>20) {
	if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	delete [] xyz_all;
        xyz_all=0;
	break;
    }
    else if(flex==1&&i==0&&n>10) {
	if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	delete [] xyz_all;
        xyz_all=0;
	break;
    }
    else if(flex==1&&i&&n/i>10) {
        if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	for(int j=0;j<i;j++) if(xyz_all[j]) delete [] xyz_all[j];
	delete [] xyz_all;
        xyz_all=0;
	break;
    }
    else if(flex==0&&i&&n/i>10) {
        if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	for(int j=0;j<i;j++) if(xyz_all[j]) delete [] xyz_all[j];
	delete [] xyz_all;
	xyz_all=0;
	break;
    }
    else if(flex==2&&i==0&&n>20) {
	if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	delete [] xyz_all;
        xyz_all=0;
	break;
    }
    else if(flex==2&&i&&n/i>5) {
        if(TRES.logg)cerr<<"hard to connect the segment between:"<<start<<" "<<end<<endl;
	break;
    }
    else if(flex==3&&i==0&&n>5) {
	if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	delete [] xyz_all;
        xyz_all=0;
	break;
    }
    else if(flex==3&&i&&n/i>5) {
        if(TRES.logg)cerr<<"hard to connect the segment between:"<<start<<" "<<end<<endl;
	
	break;
    }
    n++;ne++;
    seg1.segment->chn->transfer(xyz,seg1.segment->chn->res,seg1.end,1);    
    seg1.segment->chn->transfertemp(xyz0,seg1.segment->chn->res,seg1.end,1);
    seg2.segment->chn->transfer(xyz+m3,seg2.segment->chn->res,seg2.end,1);
    seg2.segment->chn->transfertemp(xyz0+m4,seg2.segment->chn->res,seg2.end,1);
    if(TRES.logg>3) segment->write("seg1");
    //segment->chn->transfer(xyz,t0,end,1);
    //segment->chn->transfertemp(xyz0,t0,end,1);
    //segment->write("seg2");
    //int nrd=random();srandom(nrd+1);
    //srandom(ranseed+ne+1);
    randmzfix(&seg1,&seg2);
    if(TRES.logg>3) segment->write("seg3");
    if(r1->id0==128&&0) {
	if(TRES.logg>3)seg1.segment->chn->write("aa");
	if(TRES.logg>3)seg2.segment->chn->write("bb");
        d=TRES.distance(r1->last->atm->next->next->xyz,r1->atm->xyz);
        if(TRES.logg)cerr<<"the distance..."<<d<<endl;
    }
    /*
    if(chiangle&&chiangle->next==0&&arbt>1){
	randmzChianglefix(15,15); 
    }
    */
    if(chiangle) {
	//cerr<<"rotate the dihedral to the specified in the database..."<<n-1<<"  "<<databaseonly<<endl;
	//checkomega(stdout);
	if(randmzChiangle(n-1)==0&&databaseonly==1) break;
	//checkomega(stdout);
    }
    else if(secd=='h'&&end-start+1>3) randmzHelix(15);
    else if(secd=='e'&&end-start+1>3) randmzSheet(30);
    else if(chiangle==0&&databaseonly==1) {
	break;
    }
    chn->transfer(xyz,r0,end,1);
  
    d=TRES.distance(r1->atm->xyz,r2->atm->xyz)+
      TRES.distance(r1->atm->next->xyz,r2->atm->next->xyz);
    
    if(d>dis) continue;

    if(r1->last->id0==38) {
	d=TRES.distance(r1->last->atm->next->next->xyz,r1->atm->xyz);
        if(TRES.logg)cerr<<"the distance..."<<d<<endl;
    }

    d=closeg(&seg1,&seg2);  
     if(TRES.logg>3) segment->write("seg4");	
    if(r1->last->id0==38) {
        d=TRES.distance(r1->last->atm->next->next->xyz,r1->atm->xyz);
        if(TRES.logg)cerr<<"the new distance..."<<d<<endl;
    }

    if(d==-2) continue;  

    if(r1&&r2&&end-start+1<=3) {
   	a=(*r1)[2];a1=(*r2)[2];
        d=0;
	if(a&&a1) d+=TRES.distance(a,a1);
	a=(*r1)[3];a1=(*r2)[3];
 	if(a&&a1) d+=TRES.distance(a,a1);
	if(d>1.5) continue;	
    }
	
    chn->transfer(xyz,r0,end,1);
    if(outlet<0.5)    {if(setend(1)==0) continue;}
    xyz_all[i]=new float[nseq];
    cht=seg1.segment->chn;r=cht->res;
    cht->transfer(xyz_all[i],r,seg1.end,0);
    cht=seg2.segment->chn; r=cht->res->next;
    seg2.segment->chn->transfer(xyz_all[i]+m5,r,seg2.end,0);
    segment->chn->transfer(xyz_all[i],segment->chn->res,end,1);
    if(TRES.logg>3)segment->chn->write("sss");
    i++;
    if(i>=arbt) break;    
  }
  if(TRES.logg)cerr<<"loop number :"<<n<<" "<<i<<endl;
  chn->transfer(xyz,r0,end,1);
  chn->transfertemp(xyz0,r0,end,1);
  if(xyz) {delete [] xyz;xyz=0;}
  if(xyz0) {delete [] xyz0;xyz0=0;}
  return xyz_all;
}


float** Segen::zipperfix(float dis,int flg)
{
  //generate loop conformations
  Res  *r1,*r2;
  Res  *r0,*t0,*r;  
  Atm  *a,*a1;
  int ne,i,n,nseq;
  float d;
  float *xyz,*xyz0,**xyz_all;
  Chn *chn=(*pdb)[cid];

  xyz=0;xyz_all=0;
  //srandom(seed);
  //coordinate buffer

  nseq=0; 
  r0=(*chn)[start];t0=(*segment->chn)[start];
  for(r=r0;r&&r->id0<=end;r=r->next) 
  nseq+=r->tres->number*3;     
  
  xyz=new float[nseq];  
  xyz0=new float[nseq*3];
  int arbt0=arbt;
  if(chiangle) {
	arbt0=max(arbt0,chiangle->getsize());
  }
  xyz_all=new float*[arbt0+10];
  for(i=0;i<arbt0+10;i++) xyz_all[i]=0;

  chn->transfer(xyz,r0,end,0);
  chn->transfertemp(xyz0,r0,end,0); 
  r1=(*chn)[end];r2=(*segment->chn)[end]; 

  //decides the segment to be kept

  //find possible hits

  i=0;n=0;
   
  ne=0;

  
  for(; ;)
  {
    if(flex==0&&i==0&&n>20) {
	if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	delete [] xyz_all;
        xyz_all=0;
	break;
    }
    else if(flex==1&&i==0&&n>10) {
	if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	delete [] xyz_all;
        xyz_all=0;
	break;
    }
    else if(flex==1&&i&&n/i>10) {
        if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	for(int j=0;j<i;j++) if(xyz_all[j]) delete [] xyz_all[j];
	delete [] xyz_all;
        xyz_all=0;
	break;
    }
    else if(flex==0&&i&&n/i>10) {
        if(TRES.logg)cerr<<"could not connect the segment between:"<<start<<" "<<end<<endl;
	for(int j=0;j<i;j++) if(xyz_all[j]) delete [] xyz_all[j];
	delete [] xyz_all;
	xyz_all=0;
	break;
    }

    n++;ne++;
    if(TRES.logg>3)segment->write("seg1");
    segment->chn->transfer(xyz,t0,end,1);
    segment->chn->transfertemp(xyz0,t0,end,1);
    if(TRES.logg>3){segment->write("seg2");}//segment->chn->dihedral(stderr);}
    //int nrd=random();
    //srandom(ranseed+ne+1);
    if(chiangle==0||databaseonly!=1)  randmzfix();
    if(TRES.logg>3){segment->write("seg3");}//segment->chn->dihedral(stderr);}
 
    if(chiangle) {
	//cerr<<"rotate the dihedral to the specified in the database..."<<n-1<<"  "<<databaseonly<<endl;
	//checkomega(stdout);
	if(randmzChiangle(n-1)==0&&databaseonly==1) break;
	//checkomega(stdout);
    }
    else if(secd=='h'&&end-start+1>3) randmzHelix(15);
    else if(secd=='e'&&end-start+1>3) randmzSheet(30);
    else if(chiangle==0&&databaseonly==1) {
	break;
    }
    chn->transfer(xyz,r0,end,1);
  
    d=TRES.distance(r1->atm->xyz,r2->atm->xyz)+
      TRES.distance(r1->atm->next->xyz,r2->atm->next->xyz);
   
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
    
    i++;
    //if(i>=arbt&&databaseonly==0) break;
    if(i>=arbt) break;
  }
  if(TRES.logg)cerr<<"loop number :"<<n<<" "<<i<<endl;
  chn->transfer(xyz,r0,end,1);
  if(xyz) {delete [] xyz;xyz=0;}
  if(xyz0) {delete [] xyz0;xyz0=0;}
  return xyz_all;
}

void Segen::randmzfix(Segen *seg1,Segen *seg2)
{
  Res *r,*r1,*r2,*r3;//,*r5;
  Res *last;
  //Atm *a;
  int j;
  Chn *segp;
  Rotate rot;
  Chn *rotamer=TRES.findrotamer("backbone")->chn;
  //float xr[3],xr1[3],dn;//,dca;
  Chn *chn=(*pdb)[cid];
  int e;
   

  StrFmt *parent=mutate->owner->getparentStrFmt();
  StrFmt *se = parent->getsequenceStrFmt();
  if(se==0) return;

  //for the first segment

  //if(direct==0) goto re200;
  
  //srandom(seed+1); //skip random
 
  //segp=segment->chn;

  segp=seg1->segment->chn; 
  r3=(*segp)[seg1->start];
  int n=0;

  for(r=r3;r;r=r->next)
  {
    if(r->id0>seg1->end) break;
    n=mutate->compare[r->id0];
    
    if(r==r3) {  
	 
	if(mutate->rely[r->id0]>0) {
		r3=r; continue;
	}
	else  {
		e=random()%5;
                if(e%2==0) e=-e;		 
                rot.rotate(r->atm->next->next,e,r->id0,1);				
	} 
    }  
    
    else if(randcoil==2) {
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;	
	r2=r1->ismore(j); 
 
	r->transfer(r2);
        r->copytemp(r2); 
        rot.link(r3,r,r->id0,1);    
    } 
    else if(r->name=='P') {
	rot.link(r3,r,r->id0,1); 	
	e=random()%10; 
	if(e%2==0) e=-e;
	rot.rotate(r->atm->next->next,e,r->id0,1);	
    } 
    else if(randcoil==-1&&mutate->rely[r->id0]==0) {
	r1=rotamer->isres(r->name,0);
	Res *rt=r3; 
	if(rt==0) rt=r;
	rt=chn->isres0(rt->id0);
	float a=rt->atm->next->chi;
	float b=rt->atm->next->next->chi;
	
	r2=r1->findbestmore(a,b);
	r->transfer(r2);
        r->copytemp(r2);

	rot.link(r3,r,r->id0,1);    		
	e=random()%15; 
	if(e%2==0) e=-e;
	if(r->name!='P') rot.rotate(r->atm->next,e,r->id0,1);
	e=random()%15; 
	if(e%2==0) e=-e;
	rot.rotate(r->atm->next->next,e,r->id0,1);	    
    }
    else if(r->next==0) {
     
    	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;	
	r2=r1->ismore(j); 
	 
	r->transfer(r2);
        r->copytemp(r2); 
        rot.link(r3,r,r->id0,1);    
    } 
    else if(randcoil==1)
    {    
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;	
	r2=r1->ismore(j); 
	r->transfer(r2);
        r->copytemp(r2);
        rot.link(r3,r,r->id0,1);       	 
    } 
    else if(mutate->rely[r->id0]>0) { 
	rot.link(r3,r,r->id0,1);	
    } 
    else if(mutate->dssp[n]!='-') {
	//int n=mutate->findlongestreliable(r);
       	rot.link(r3,r,r->id0,1);       
	if(mutate->rely[r->id0]==1) { 	 
		e=random()%5; 
		if(e%2==0) e=-e;
		if(r->name!='P') rot.rotate(r->atm->next,e,r->id0,1);
		e=random()%5; 
		if(e%2==0) e=-e;
		rot.rotate(r->atm->next->next,e,r->id0,1); 		 	
	}
	else if(mutate->dssp[n]=='h') {

		r1=rotamer->isres(r->name,0);
		r2=r1->findbestmore(-60,-45);
		r->transfer(r2);
        	r->copytemp(r2);

		rot.link(r3,r,r->id0,1);    		
		e=random()%15; 
		if(e%2==0) e=-e;
		if(r->name!='P') rot.rotate(r->atm->next,e,r->id0,1);
		e=random()%15; 
		if(e%2==0) e=-e;
		rot.rotate(r->atm->next->next,e,r->id0,1);	    
	}
	else if(mutate->dssp[n]=='e') {

		r1=rotamer->isres(r->name,0);
		r2=r1->findbestmore(-110,130);
		r->transfer(r2);
        	r->copytemp(r2);
		 
		rot.link(r3,r,r->id0,1);    		
		e=random()%30; 
		if(e%2==0) e=-e;
		if(r->name!='P') rot.rotate(r->atm->next,e,r->id0,1);
		e=random()%30; 
		if(e%2==0) e=-e;
		rot.rotate(r->atm->next->next,e,r->id0,1);			
	} 
    } 
    else  
    {    
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;	
	r2=r1->ismore(j); 
	r->transfer(r2);
        r->copytemp(r2);
        rot.link(r3,r,r->id0,1);       	 
    }
     
    r3=r;
  }
  
  //for the next segment;

  segp=seg2->segment->chn;
  if(segp->res==0) return;

  last=(*segp)[end];
  r3=last;
  
  //for(m=seg2->end;m>=seg2->start;m--)
  for(r=last;r;r=r->last)
  {
    n=mutate->compare[r->id0];
   
    
    if(r==r3) {	
	
	if(mutate->rely[r->id0]>0) {
		r3=r;
		continue;
	}
	else {
		e=random()%5;
                if(e%2==0) e=-e;		
                rot.rotate(r->atm->next,e,r->id0,0);		
	}  
    } 
   
    else if(randcoil==2) {
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;		
	r2=r1->ismore(j); 		 	 
	r->transfer(r2);
        r->copytemp(r2);	 
	rot.link(r,r3,r->id0,0);     
    } 
    else if(r->name=='P') {
	rot.link(r,r3,r->id0,0); 	 
        e=random()%10;
        if(e%2==0) e=-e;
        rot.rotate(r->atm->next->next,e,r->id0,0);        
    } 
    else if(randcoil==-1&&mutate->rely[r->id0]==0) {
	r1=rotamer->isres(r->name,0);
	Res *rt=r3; 
	if(rt==0) rt=r;
	rt=chn->isres0(rt->id0);
	float a=rt->atm->next->chi;
	float b=rt->atm->next->next->chi;
	r2=r1->findbestmore(a,b);
	r->transfer(r2);
        r->copytemp(r2);
	rot.link(r,r3,r->id0,0);  		
	e=random()%15;
        if(e%2==0) e=-e;
        rot.rotate(r->atm->next->next,e,r->id0,0);
	e=random()%15;
        if(e%2==0) e=-e;
        if(r->name!='P') rot.rotate(r->atm->next,e,r->id0,0);    
    }
    else if(r->last==0) {	        
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;		
	r2=r1->ismore(j); 		 	 
	r->transfer(r2);
        r->copytemp(r2);	 
	rot.link(r,r3,r->id0,0);     
    } 
    else if(randcoil==1)
    {      
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;	
	r2=r1->ismore(j); 
	r->transfer(r2);
        r->copytemp(r2);
	rot.link(r,r3,r->id0,0);        
    } 
    else if(mutate->rely[r->id0]>0) {
	rot.link(r,r3,r->id0,0);
	
    } 
    else if(mutate->dssp[n]!='-') {
	rot.link(r,r3,r->id0,0);
	if(mutate->rely[r->id0]==1) {		 
                e=random()%5;
                if(e%2==0) e=-e;
                rot.rotate(r->atm->next->next,e,r->id0,0);
		e=random()%5;
                if(e%2==0) e=-e;
                if(r->name!='P') rot.rotate(r->atm->next,e,r->id0,0);
		 
        }
	else if(mutate->dssp[n]=='h') {
		r1=rotamer->isres(r->name,0);
		r2=r1->findbestmore(-60,-45);
		r->transfer(r2);
        	r->copytemp(r2);
		rot.link(r,r3,r->id0,0);  
		
		e=random()%15;
                if(e%2==0) e=-e;
                rot.rotate(r->atm->next->next,e,r->id0,0);
		e=random()%15;
                if(e%2==0) e=-e;
                if(r->name!='P') rot.rotate(r->atm->next,e,r->id0,0);
    
	}
	else if(mutate->dssp[n]=='e') {

		r1=rotamer->isres(r->name,0);
		r2=r1->findbestmore(-110,130);
		r->transfer(r2);
        	r->copytemp(r2);
		rot.link(r,r3,r->id0,0);  
		
		e=random()%30;
                if(e%2==0) e=-e;
                rot.rotate(r->atm->next->next,e,r->id0,0);
		e=random()%30;
                if(e%2==0) e=-e;
                if(r->name!='P') rot.rotate(r->atm->next,e,r->id0,0);
	}
    } 
    else  
    {      
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;	
	r2=r1->ismore(j); 
	r->transfer(r2);
        r->copytemp(r2);
	rot.link(r,r3,r->id0,0);        
    }
    r3=r; 
  }
}

 
 
void Segen::randmzfix()
{
  Res *r,*r1,*r2,*r3; 
  Res *last;
  int j,m;
  Chn *segp;
  Rotate rot;
  Chn *rotamer=TRES.findrotamer("backbone")->chn;
  Chn *chn=(*pdb)[cid];
  int e;

  Res *rb=(*chn)[start-1];
  Res *re=(*chn)[end+1];
  //for the first segment
  if(direct==0) goto re200;
  

  //srandom(seed+1); //skip random
 
  segp=segment->chn;

  r3=(*segp)[start];
  //int n=0;

  for(r=r3;r;r=r->next)
  {
    if(r->id0>end) break;
  
    if(r==r3) {  	 
	if(mutate->rely[r->id0]==2) {
		r3=r; continue;
	}
	else {				 
		e=random()%(flexrot/2); 
		if(e%2==0) e=-e;
		rot.rotate(r->atm->next->next,e,r->id0,1);		    		
	} 
    }    
    else if(randcoil==2) {
        r1=rotamer->isres(r->name,0);
        j=random()%r1->nummore;
        r2=r1->ismore(j);

        r->transfer(r2);
        r->copytemp(r2);
        rot.link(r3,r,r->id0,1);
    }         
    else if(r->name=='P'&&mutate&&mutate->rely[r->id0]==2) {
	rot.link(r3,r,r->id0,1); 
    }
    else if(r->name=='P'&&mutate&&mutate->rely[r->id0]==1) {
        rot.link(r3,r,r->id0,1);
        e=random()%flexrot;
        if(e%2==0) e=-e;
        rot.rotate(r->atm->next->next,e,r->id0,1);       
    }
    else if(r->next==0&&re) {
	r1=rotamer->isres(r->name,0);
        j=random()%r1->nummore;
        r2=r1->ismore(j);
        r->transfer(r2);
        r->copytemp(r2);
        rot.link(r3,r,r->id0,1);
    }
    else if(randcoil==1) {
        r1=rotamer->isres(r->name,0);
        j=random()%r1->nummore;
        r2=r1->ismore(j);
        r->transfer(r2);
        r->copytemp(r2);
        rot.link(r3,r,r->id0,1);
    } 
    else if(mutate->rely[r->id0]==2) {	 
       	rot.link(r3,r,r->id0,1);       	 
    }   
    else if(mutate->rely[r->id0]==1) {
       	rot.link(r3,r,r->id0,1);       	
	e=random()%flexrot; 
	if(e%2==0) e=-e;
	rot.rotate(r->atm->next,e,r->id0,1);	
	//
	e=random()%flexrot; 
	if(e%2==0) e=-e;
	rot.rotate(r->atm->next->next,e,r->id0,1);	 
    }   
    else  
    {    
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;	
	r2=r1->ismore(j); 
	r->transfer(r2);
        r->copytemp(r2);
        rot.link(r3,r,r->id0,1);       	
    }
    r3=r;   
  }


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
	 
    //if(r->id0>end) break;
	
    if(r==r3) {  
	 
	if(mutate->rely[r->id0]==2) {
		r3=r; continue;
	}
	else {		 		 
		e=random()%(flexrot/2);
                if(e%2==0) e=-e;
		rot.rotate(r->atm->next,e,r->id0,0);              		
	} 
    }
    else if(randcoil==2) {
        r1=rotamer->isres(r->name,0);
        j=random()%r1->nummore;
        r2=r1->ismore(j);
        r->transfer(r2);
        r->copytemp(r2);
        rot.link(r,r3,r->id0,0);
    }   
    else if(r->name=='P'&&mutate&&mutate->rely[r->id0]==2) {
	rot.link(r,r3,r->id0,0);
    }
    else if(r->name=='P'&&mutate&&mutate->rely[r->id0]==1) {
        rot.link(r,r3,r->id0,0);
        e=random()%flexrot;
        if(e%2==0) e=-e;
        rot.rotate(r->atm->next,e,r->id0,0);
    }    
    else if(r->id0==start&&rb) {
	r1=rotamer->isres(r->name,0);
        j=random()%r1->nummore;
        r2=r1->ismore(j);
        r->transfer(r2);
        r->copytemp(r2);
        rot.link(r,r3,r->id0,0);
    }
    else if(randcoil==1)
    {
        r1=rotamer->isres(r->name,0);
        j=random()%r1->nummore;
        r2=r1->ismore(j);
        r->transfer(r2);
        r->copytemp(r2);
        rot.link(r,r3,r->id0,0);
    }
    else if(mutate->rely[r->id0]==2) {
	rot.link(r,r3,r->id0,0);       
    }
    else if(mutate->rely[r->id0]==1) {
	rot.link(r,r3,r->id0,0);
        e=random()%flexrot;
        if(e%2==0) e=-e;	
        rot.rotate(r->atm->next->next,e,r->id0,0);
	//
	e=random()%flexrot;
        if(e%2==0) e=-e;	
        rot.rotate(r->atm->next,e,r->id0,0);
    }
    else  
    {      
	r1=rotamer->isres(r->name,0);
	j=random()%r1->nummore;	
	r2=r1->ismore(j); 
	r->transfer(r2);
        r->copytemp(r2);
	rot.link(r,r3,r->id0,0);    
    }
    r3=r;  
  }
}
 

float *Segen::getaveragepost(float **xyz_all) {

	Res *r0,*r;
 
	Chn *chn=(*pdb)[cid];
	r0=(*chn)[start];

	int nseq=0; 
  	for(r=r0;r&&r->id0<=end;r=r->next) nseq+=r->tres->number*3;         

  	float *xyz=new float[nseq];
	 

	int all=0;
  	while(xyz_all[all])all++;
 
	int n,i;
	for(i=0;i<nseq;i++) xyz[i]=0;

	for(n=0;n<all;n++) {
		for(i=0;i<nseq;i++) xyz[i]+=xyz_all[n][i]/all;		
	}
	 
	return xyz;
}

 
 
int Segen::ifrotatmexist(Atm *aa2) {

	if(rotatm==0) return 0;
	int n=0;
	while(rotatm[n]) {
		if(rotatm[n]==aa2) return 1;
		if(rotatm[n]->res->id0==aa2->res->id0&&rotatm[n]->tatm->id==aa2->tatm->id) return 1;
		n++;
	}
	return 0;
}
