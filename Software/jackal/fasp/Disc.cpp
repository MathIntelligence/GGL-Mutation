#include"source.h"
Disc::Disc()
{
 dipole=0;
 charge=0;
 pdb=0;
 strcpy(force,"udD");
 cutoff=6.0;
 range=0;
 change=0;
 int i=0;
 for(i=0;i<200;i++) energy[i]=0;
 eps=20.;
 smoothclash=0;
}

Disc::~Disc()
{
 int i,kep;
 if(dipole) {
	kep=0;
	while(dipole[kep])kep++;
	for(i=0;i<kep;i++) delete dipole[i];
        delete [] dipole;
 }
 if(charge) {
        kep=0;
        while(charge[kep])kep++;        
        for(i=0;i<kep;i++) delete charge[i];
        delete [] charge;
 }
 if(range) delete [] range;
 if(change) delete [] change;
}

void Disc::clear() {
 int i,kep;
 if(dipole) {
	kep=0;
	while(dipole[kep])kep++;
	for(i=0;i<kep;i++) delete dipole[i];
        delete [] dipole;dipole=0;
 }
 if(charge) {
        kep=0;
        while(charge[kep])kep++;        
        for(i=0;i<kep;i++) delete charge[i];
        delete [] charge;charge=0;
 }
 if(range) delete [] range;range=0;
}

void Disc::setupall(Pdb *s,float off,int iscrg) {
	pdb=s;
	ready();
	cutoff=off;
	setrange();
	setup(iscrg);
	setdipoleascharge();
	setrange();
}

void Disc::write(FILE *fp) {

	int n=0;while(charge&&charge[n])n++;

  	int i;
	for(i=0;i<n;i++) {
		fprintf(fp,"%c %5i %s\n",charge[i]->atm->res->name,charge[i]->atm->res->id0,charge[i]->atm->name); 
	}
	n=0;while(dipole&&dipole[n])n++;

 	
	for(i=0;i<n;i++) {
		fprintf(fp,"%c %5i %s: ",dipole[i]->atm->res->name,dipole[i]->atm->res->id0,dipole[i]->atm->name); 
		fprintf(fp,"%s %s ",dipole[i]->crg[0]->atm->name,dipole[i]->crg[1]->atm->name);
		float d=TRES.distance(dipole[i]->crg[0]->atm->xyz,dipole[i]->crg[0]->xyz);
		float d1=TRES.distance(dipole[i]->crg[1]->atm->xyz,dipole[i]->crg[1]->xyz);
		fprintf(fp,"%f %f\n",d,d1);
		//fprintf(fp,"%f %f %f\n",dipole[i]->crg[1]->xyz[0],dipole[i]->crg[1]->xyz[1],dipole[i]->crg[1]->xyz[2]);
	}
}


float Disc::clash(Lattice *lat, Res *rr,int n) {
 
	Res *r;
	Atm *a;
	float e=0;

	int i=0;
 	for(i=0;i<200;i++) energy[i]=0;
	for(r=rr;r;r=r->next) {	
		if(r->id0>n) break;
		float ee=0;
		for(a=r->atm;a;a=a->next) {
			if(a->flag<0) continue;
			ee+=clash(lat,a);	
		}
		if(TRES.logg>3&&strchr("AVLIFCMPWY",r->name)) 
		if(TRES.logg>3)cerr<<r->name<<r->id0<<" "<<ee<<endl;
  		e+=ee;
		//cerr<<r->name<<r->id0<<e<<" "<<energy['u']<<endl;
	}
	return e;
}

float Disc::clash(Res *rr,int n) {
 
	Lattice lat;
  	lat.putoff();
	lat.flag=1;
	lat.grdsiz=2;
	lat.radall=15;
	lat.ready(pdb);
	lat.puton(pdb);
	
	pdb->setflg(-1);
	rr->chn->setflg(rr,n,1);

	int i=0;
 	for(i=0;i<200;i++) energy[i]=0;

	return clash(&lat,rr,n);
}

Dipole *Disc::finddipole(Atm *a) {

	int kep=0;while(dipole[kep])kep++;
	int i;
	for(i=0;i<kep;i++) {
		if(dipole[i]->atm==a) return dipole[i];
	}
	return 0;
}

Charge *Disc::findcharge(Atm *a) {

	int kep=0;while(charge[kep])kep++;
	int i;
	for(i=0;i<kep;i++) {
		if(charge[i]->atm==a) return charge[i];
	}
	return 0;
}

float Disc::clash(Lattice *lat,Atm *aa0) {
 
  Atm *aa1,*all[50];
  float dr,dn,e,d,elon,xr[3];
  int i,j,m,n,isp;
  float ta,tb,tc,tn,tt;
  int mr; 
  int nget=0;
  float z[3];

  e=0;
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
  
  int nte=0;
  if(strchr(force,'v')) nte=1;
   
  if(nte&&nget==0) { 
  	lat->getcell(aa0,cutoff);
  	lat->cutoff(aa0,cutoff);
	nget=lat->nget;
  }
 
  if(TRES.logg>4) lat->printoutLattice();
  for(m=0;nte&&m<nget;m++)
  {
     aa1=lat->obtain[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
     mr=0;float ef=0;
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
        if(aa0->tatm->id==5&&aa1->tatm->id==5)ef+=-10.*exp(-(dr-1)*(dr-1)); 
     }
     elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     elon=ta*elon*exp(-dr*tb);
     dr=1/(tt+dr);
     
     dn=pow(dr,tn);  
     dn=elon*dn*(dn-tc);   
     
     isp=aa0->tatm->ispolar*aa1->tatm->ispolar;
     
     if(dn<0)
     {
      if(mr==1)       ef+= dn*1.5;
      else if(mr==2)  ef+= dn*5;
      else if(isp==1) ef+= dn*1.50;  
      else if(isp==2) ef+= dn*1.25; 
      else if(isp==3) ef+= dn*0.75; 
      else            ef+= dn;
     }
     else
     {
      if(mr==1)       ef+= dn*0.50;
      else if(mr==2)  ef+= dn*0.20;
      else if(isp==1) ef+= dn*0.5; 
      else if(isp==2) ef+= dn*0.75;
      else if(isp==3) ef+= dn*1.25;
      else            ef+= dn;
     }

     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

     if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	ef=ef/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	ef=ef/5;
     }
    	
     energy['v']+=ef;
     e+=ef; 
     //if(dn>1000||dn<-1000) cerr<<"dn..."<<dn<<" "<<aa0->name<<" "<<aa0->res->id<<" : "<<aa1->name<<" "<<aa1->res->id<<endl;
  }
  
  nte=0;
  if(strchr(force,'h')) nte=1;
   
  if(nte&&nget==0) { 
  	lat->getcell(aa0,cutoff);
  	lat->cutoff(aa0,cutoff);
	nget=lat->nget;
  }
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
 
     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

     if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	ef=ef/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	ef=ef/5;
     }

     energy['h']+=ef;
     e+=ef; 
  }


  nte=0;
  if(strchr(force,'u')) nte=1;
   
  if(nte&&nget==0) { 
  	lat->getcell(aa0,cutoff);
  	lat->cutoff(aa0,cutoff);
	nget=lat->nget;
  }
 
  for(m=0;nte&&m<nget;m++)
  {
     aa1=lat->obtain[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
     mr=0;float ef=0;
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
        if(aa0->tatm->id==5&&aa1->tatm->id==5)ef+=-10.*exp(-(dr-1)*(dr-1)); 
     }
     elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     elon=ta*elon*exp(-dr*tb);
     dr=1/(tt+dr);
     
     dn=pow(dr,tn);  
     dn=elon*dn*(dn-tc);   
     dn=dn;

     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

     if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	dn=dn/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	dn=dn/5;
     }
    
     //

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
     if(dn>0) {
	if(hp==1) dn=dn/10;
     }
     else {
	if(hp==1) dn=dn;
	else dn=dn/10;
     }

     ef+=dn;
     energy['u']+=ef;
     e+=ef;
  }


  nte=0;
  if(strchr(force,'V')) nte=1;
   
  if(nte&&nget==0) { 
  	lat->getcell(aa0,cutoff);
  	lat->cutoff(aa0,cutoff);
	nget=lat->nget;
  }
 
  for(m=0;nte&&m<nget;m++)
  {
     aa1=lat->obtain[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
     mr=0;float ef=0;
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
        if(aa0->tatm->id==5&&aa1->tatm->id==5)ef+=-10.*exp(-(dr-1)*(dr-1)); 
     }
     elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     elon=ta*elon*exp(-dr*tb);
     dr=1/(tt+dr);
     
     dn=pow(dr,tn);  
     dn=elon*dn*(dn-tc);   

     ef+=dn;

     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

     if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	ef=ef/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	ef=ef/5;
     }

     e+=ef;
     energy['V']+=ef;
     //if(dn>1000||dn<-1000) cerr<<"dn..."<<dn<<" "<<aa0->name<<" "<<aa0->res->id<<" : "<<aa1->name<<" "<<aa1->res->id<<endl;
  }
  
  nte=0;
  if(strchr(force,'D')) nte=1; //new energy function of dipole,dipole interaction
  
  Dipole *da0=0;
  float f0[3],f1[3],f2[3];
  float ff[2][3];

  int kep=0;
  
  if(nte) {
	//da0=finddipole(aa0);
	da0=range[aa0->id0];
   	//kep=0;while(dipole&&dipole[kep])kep++;	
  }

  if(nte&&nget==0&&da0) { 
  	lat->getcell(aa0,cutoff);
  	lat->cutoff(aa0,cutoff);
	nget=lat->nget;
	for(i=0;i<3;i++) f0[i]=(da0->crg[0]->xyz[i]+da0->crg[1]->xyz[i])/2;
  }

  Dipole *dp0;
 
  for(m=0;da0&&nte&&m<nget;m++)  //dipole and dipole interaction
  {
     aa1=lat->obtain[m]->atm;     
     dp0=range[aa1->id0];
     if(dp0==0) continue;  
     
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
 
     if(aa0->res==aa1->res) continue;
      
     for(i=0;i<3;i++) f1[i]=(dp0->crg[0]->xyz[i]+dp0->crg[1]->xyz[i])/2;      
     for(i=0;i<3;i++) f2[i]=f1[i]-f0[i];
     d=TRES.distance(f2);

     if(d<0.1) continue;
     else if(d<3.5) {	 
	for(j=0;j<3;j++)f2[j]=f2[j]*(3.5-d)/d;
	for(i=0;i<2;i++)  		
	for(j=0;j<3;j++) {
		ff[i][j]=dp0->crg[i]->xyz[j]+f2[j];	
	}		 
     }
     else {
	for(i=0;i<2;i++) 
	for(j=0;j<3;j++) {
		ff[i][j]=dp0->crg[i]->xyz[j];		
	}
     }

     int i0,j0;
     float ef=0; 
     for(i0=0;i0<2;i0++) {
	Charge *c=da0->crg[i0];
	if(c==0) continue;
	for(j0=0;j0<2;j0++) {
		Charge *t=dp0->crg[j0];
		if(t==0) continue;
		for(i=0;i<3;i++) z[i]=c->xyz[i]-ff[j0][i];
		float d=TRES.distance(z); 		 			
		float f=332*c->crg*t->crg/d/eps;
		ef+=f;
		
	}
     }     

     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

    if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	ef=ef/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	ef=ef/5;
     }

     e+=ef;
     energy['D']+=ef;   
  }
   
  nte=0;
  if(strchr(force,'d')) nte=1; //new energy function of dipole,charge interaction
  //da0=finddipole(aa0);
  da0=range[aa0->id0];
  kep=0;while(charge&&charge[kep])kep++;
  for(i=0;da0&&i<3;i++) f0[i]=(da0->crg[0]->xyz[i]+da0->crg[1]->xyz[i])/2;
  for(m=0;da0&&nte&&m<kep;m++)  //dipole and charge interaction
  {
     aa1=charge[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa0->res==aa1->res) continue;     
      
     for(i=0;i<3;i++) f1[i]=charge[m]->xyz[i];      
     for(i=0;i<3;i++) f2[i]=f1[i]-f0[i];
     d=TRES.distance(f2);

     if(d<0.1) continue;
     else if(d<3.5) {
	 for(i=0;i<3;i++) f2[i]=f2[i]*(3.5-d)/d; 
	 for(i=0;i<3;i++) f1[i]=charge[m]->xyz[i]+f2[i];	 
     }
     else {
	 for(i=0;i<3;i++) f1[i]=charge[m]->xyz[i];
     }
     

     int j0;
     float ef=0;
     for(j0=0;j0<2;j0++) {
	Charge *t=da0->crg[j0];
	if(t==0) continue;
	for(i=0;i<3;i++) z[i]=f1[i]-t->xyz[i];
	float d=TRES.distance(z); 	 
	float f=332*charge[m]->crg*t->crg/d/eps;
	ef+=f;	 
     }     

     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

    if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	ef=ef/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	ef=ef/5;
     }

     e+=ef;
     energy['d']+=ef;
  }

  nte=0;
  if(strchr(force,'c')) nte=1; //new energy function of dipole,charge interaction
  Charge *ch0=0;
  if(nte) {
	ch0=findcharge(aa0);
   	kep=0;while(charge&&charge[kep])kep++;
  }
  
  for(m=0;ch0&&nte&&m<kep;m++)  //charge and charge interaction
  {
     aa1=charge[m]->atm;
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
     if(aa0->res==aa1->res) continue;

     int i;
     for(i=0;i<3;i++) z[i]=ch0->xyz[i]-charge[m]->xyz[i];
     float d=TRES.distance(z);     
     if(d<3.5) d=3.5;

     float f=332*ch0->crg*charge[m]->crg/d/eps;

     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

     if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	f=f/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	f=f/5;
     }

     e+=f;	
     energy['c']+=f;
     if(TRES.logg>3)cerr<<aa0->res->name<<aa0->res->id0<<" "<<aa0->name<<" "<<aa1->res->name<<aa1->res->id0<<" "<<aa1->name<<" ";
     if(TRES.logg>3)cerr<<f<<":"<<d<<endl;	 
  }

      
  nte=0;
  if(strchr(force,'p')) nte=1; //hydrophobic contact energy
  if(nte&&nget==0) {
	lat->getcell(aa0,cutoff);
  	lat->cutoff(aa0,cutoff);
	nget=lat->nget;
  }
 
  for(m=0;nte&&m<nget;m++) //hydrophobic contact
  {
     if(aa0->tatm->id<4) break;     
     if(aa0->tatm->name[1]=='H'&&aa0->tatm->bond[0]->id<4) break;
     if(strchr("AVLIFCMPW",aa0->res->name)==0)break;

     aa1=lat->obtain[m]->atm;
     //if(aa0->res==aa1->res) continue;
     if(aa1->tatm->id<4) continue;     
     if(strchr("AVLIFCMPW",aa1->res->name)==0)continue;
     if(aa1->tatm->name[1]=='H'&&aa1->tatm->bond[0]->id<4) continue;

     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
       
 
     float d=TRES.distance(aa0,aa1);
     float ef=0;
     if(d<aa0->tatm->eng->radius+aa1->tatm->eng->radius+2.8) {
	if(TRES.logg>3) cerr<<aa0->res->name<<aa0->res->id0<<" "<<aa0->name<<" "<<aa1->res->name<<aa1->res->id0<<" "<<aa1->name<<" "<<d<<" "<<aa0->tatm->eng->radius+aa1->tatm->eng->radius+2.8<<endl;
	ef+=-1;
        
     }

     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

    if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	ef=ef/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	ef=ef/5;
     }

     e+=ef;
     energy['p']+=ef;
  } 

  nte=0;
  if(strchr(force,'P')) nte=1; //hydrophobic contact energy
  if(nte&&nget==0) {
	lat->getcell(aa0,cutoff);
  	lat->cutoff(aa0,cutoff);
	nget=lat->nget;
  }
 
  for(m=0;nte&&m<nget;m++) //hydrophobic contact
  {
     if(aa0->tatm->id<4) break;
     if(aa0->tatm->name[1]=='H'&&aa0->tatm->bond[0]->id<4) break; 
     if(strchr("AVLIFCMPW",aa0->res->name)==0)break;
  
     aa1=lat->obtain[m]->atm;

     if(aa1->tatm->id<4) continue;
     if(aa1->tatm->name[1]=='H'&&aa1->tatm->bond[0]->id<4) continue;
     if(strchr("AVLIFCMPW",aa1->res->name)==0)continue;
 
     if(aa1->flag<0&&aa0->flag<0) continue;
     if(aa1->flag>=0&&aa0->flag>=0)  //kick out doubling computing.
     {
       if(aa1->id0<aa0->id0) continue;
     }
  
     mr=0;float ef=0;
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
        if(aa0->tatm->id==5&&aa1->tatm->id==5)ef+=-10.*exp(-(dr-1)*(dr-1)); 
     }
     elon=sqrt(aa0->tatm->eng->epslon*aa1->tatm->eng->epslon);
     elon=ta*elon*exp(-dr*tb);
     dr=1/(tt+dr);
     
     dn=pow(dr,tn);  
     dn=elon*dn*(dn-tc);  
     ef+=dn; 


     //set backbone flag
     int bk1,bk2;

     if(aa0->tatm->name[1]=='H') bk1=0;
     else if(aa0->tatm->id>4) bk1=0;
     else bk1=1;

     if(aa1->tatm->name[1]=='H') bk2=0;
     else if(aa1->tatm->id>4) bk2=0;
     else bk2=1;

     if(strchr(force,'@')&&(bk1==0&&bk2==0)) {
	ef=ef/10;
     }
     else if(strchr(force,'@')&&(bk1==0||bk2==0)) {
	ef=ef/5;
     }

     energy['P']+=ef;
     e+=ef;     
  } 
   
  for(i=0;i<n;i++)lat->puton(all[i]);

  return e;  
 
}

void Disc::setrange() {
  
  if(range) delete [] range;range=0;
  int nn=pdb->maxatmid0()+100;
  range=new Dipole*[nn];
  
 
  int i;
  for(i=0;i<nn;i++) range[i]=0;

  int kep=0;while(dipole[kep])kep++;

  for(i=0;i<kep;i++) {
	Atm *a=dipole[i]->atm;
        range[a->id0]=dipole[i];	
  }
}


void Disc::ready() {

clear();

if(pdb==0) return;
pdb->setlastresorder();
int n=pdb->manyatm();

dipole=new Dipole*[2*n];
charge=new Charge*[2*n];

int i;
for(i=0;i<2*n;i++) {
dipole[i]=0;charge[i]=0; 
}
 
Chn *c;
Res *r;
Atm *a;

//for(c=pdb->chn;c;c=c->next) c->header(); 

for(c=pdb->chn;c;c=c->next) 
for(r=c->res;r;r=r->next)
for(a=r->atm;a;a=a->next) {
	a->allnear(3,0);
}
}

void Disc::setup(int iscrg) {
	Chn *c;
	for(c=pdb->chn;c;c=c->next) 
		setup(c->res,c->lastres()->id0+1000,iscrg);
}

void Disc::setdipoleascharge(Dipole *c) {
	
		Res *s=c->atm->res;
		if((s->name=='Y'||s->name=='T'||s->name=='S')&&c->atm->tatm->id>4) {
			if(c->crg[0]==0) c->crg[0]=new Charge;
			if(c->crg[1]==0) c->crg[1]=new Charge;
			c->crg[0]->crg=-1;  
			c->crg[1]->crg=1;  
			
			int j;
			Atm *b=0;
			for(j=0;j<4;j++) {
				if(c->atm->bond[j]==0) continue;
				if(c->atm->bond[j]->name[1]=='H') {b=c->atm->bond[j];break;}
			}	
			if(b==0) {
				c->crg[0]->atm=c->atm->bond[0];
				c->crg[1]->atm=c->atm;
			}
			else {
				c->crg[0]->atm=c->atm;
				c->crg[1]->atm=b;
			}	
			for(j=0;j<3;j++) c->crg[0]->xyz[j]=c->atm->xyz[j]-0.25*c->pot[j]; 
			for(j=0;j<3;j++) c->crg[1]->xyz[j]=c->atm->xyz[j]+0.25*c->pot[j];
				
		}				
		else if(c->atm->tatm->name[1]=='N') {
			if(c->crg[0]==0) c->crg[0]=new Charge;
			if(c->crg[1]==0) c->crg[1]=new Charge;
			c->crg[0]->crg=-1;  
			c->crg[1]->crg=1;  
			
			int j;
			Atm *b=0;
			for(j=0;j<4;j++) {
				if(c->atm->bond[j]==0) continue;
				if(c->atm->bond[j]->name[1]=='H') {
					b=c->atm->bond[j]; 
					break;
				}
			}	
			if(b==0) {
				if(c->atm->tatm->id==0) c->crg[0]->atm=c->atm->next;
				else c->crg[0]->atm=c->atm->bond[0];
				c->crg[1]->atm=c->atm;
			}
			else {
				c->crg[0]->atm=c->atm;
				c->crg[1]->atm=b;
			}	
			for(j=0;j<3;j++) c->crg[0]->xyz[j]=c->atm->xyz[j]+0.5*c->pot[j]; 
			for(j=0;j<3;j++) c->crg[1]->xyz[j]=c->atm->xyz[j]+1.0*c->pot[j];			 
		}
		else if(c->atm->tatm->name[1]=='O') {
			if(c->crg[0]==0) c->crg[0]=new Charge;
			if(c->crg[1]==0) c->crg[1]=new Charge;
			c->crg[0]->crg=-1; c->crg[0]->atm=c->atm;
			c->crg[1]->crg=1;  c->crg[1]->atm=c->atm->bond[0];
			int j;
			for(j=0;j<3;j++) c->crg[0]->xyz[j]=c->atm->xyz[j]+0.0*c->pot[j];
			for(j=0;j<3;j++) c->crg[1]->xyz[j]=c->atm->xyz[j]+0.5*c->pot[j];
		}	
		else if(c->atm->tatm->name[1]=='S') {
			if(c->crg[0]==0) c->crg[0]=new Charge;
			if(c->crg[1]==0) c->crg[1]=new Charge;
			c->crg[0]->crg=-10; c->crg[0]->atm=c->atm;
			c->crg[1]->crg=10;  c->crg[1]->atm=c->atm->bond[0];
			int j;
			for(j=0;j<3;j++) c->crg[0]->xyz[j]=c->crg[0]->atm->xyz[j];
			for(j=0;j<3;j++) c->crg[1]->xyz[j]=c->crg[1]->atm->xyz[j];
		}
}

void Disc::setdipoleascharge() {

	int kep=0;while(dipole[kep])kep++;

	int i;

	for(i=0;i<kep;i++) {
		Dipole *c=dipole[i];
		setdipoleascharge(c);		
	}

}

void Disc::setdipoleascharge(Res *rr,int nn) {

	 
	Res *r;
	Atm *a;
	 
	
	for(r=rr;r;r=r->next) {
		if(r->id0>nn) break;
		for(a=r->atm;a;a=a->next) {
			Dipole *c=range[a->id0];
			if(c==0) continue;
			setdipoleascharge(c);		
		}
	}

}

void  Disc::setup(Res *rr,int nn,int iscrg) {
 
Res *r,*s;
Atm *a,*b,*a1,*a2;
int crg=0;
int dip=0;
float t;
float xyz[6];
 
Dipole *dp0;
Charge *ch0;

crg=0;while(charge&&charge[crg]) crg++;
dip=0;while(dipole&&dipole[dip]) dip++;

int i;
for(r=rr;r&&r->id0<=nn;r=r->next)
for(a=r->atm;a;a=a->next) { 
   if(iscrg==0&&a->name[1]!='N'&&a->name[1]!='O'&&a->name[1]=='S') continue;
   if(iscrg==1&&a->name[1]!='N'&&a->name[1]!='O'&&a->name[1]=='S'&&strchr("DER",r->name)==0) continue;
   if(r->name=='P'&&a->name[1]=='N') continue;
   
   s=r; 
   if(a->tatm->id==0&&s->last) {
	b=a->ischildatm(" HN ");
	if(b==0) {			
		Atm *t1=s->isatmid(1); 	
		Atm *t2=s->last->isatmid(2);
		if(t1==0||t2==0) continue;
		b=new Atm;
		for(i=0;i<3;i++)b->xyz[i]=(t1->xyz[i]+t2->xyz[i])/2;
		a1=b;a2=a; 							
	}
	else {
		a1=a;a2=b;
	}  
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.01) {
		if(b->res==0) delete b;b=0;
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}
	
	for(i=3;i<6;i++) {
		if(b->res==0) xyz[i]=a->xyz[i-3]+xyz[i-3]; 
		else	      xyz[i]=b->xyz[i-3]; 
	}
	if(b->res==0) delete b;b=0;

	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];			 
   }	
   else if(iscrg==0&&a->tatm->id==0&&s->last==0) {

	b=a->next;

	if(b==0||b->tatm->id!=1) continue;
	 	 
	a1=b;a2=a;
	   
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.01) {		
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}
	
	for(i=3;i<6;i++)  xyz[i]=a->xyz[i-3]+xyz[i-3]; 
 
	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];			 
   }					
   else if(a->tatm->id==3&&(s->next||iscrg==0)) {
	b=a->bond[0];
	if(b==0) continue;
	a1=a;a2=b;	 
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.1) {
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}	
	for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]; 
	 
	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];		 	 				 
   } 
   else if((s->name=='Q'||s->name=='N')&&a->name[1]=='N'&&a->tatm->id>4) {	
	b=a->bond[0]; 
	if(b==0) continue;
	a1=b;a2=a;
	   
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.01) {
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}
	
	for(i=3;i<6;i++)  xyz[i]=a->xyz[i-3]+xyz[i-3]; 
 
	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];		
	continue;
   }
   else if((s->name=='N'||s->name=='Q')&&a->name[1]=='O'&&a->tatm->id>4) {
	b=a->bond[0];  	
	if(b==0) continue;
	a1=a;a2=b;

	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.1) {
		continue;
	}
	for(i=0;i<3;i++) xyz[i]=xyz[i]/t;   
	for(i=3;i<6;i++) xyz[i]=a->xyz[i-3];
	
	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];		
   } 
   else if(s->name=='R'&&strcmp(a->name," NE ")==0) {
	continue; 
   }
   else if(iscrg==0&&(s->name=='R'||s->name=='K')&&a->name[1]=='N'&&a->tatm->id>4) {
	b=a->bond[0];
	if(b==0) continue;
	a1=b;a2=a;
	   
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.01) {
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}
	
	for(i=3;i<6;i++)  xyz[i]=a->xyz[i-3]+xyz[i-3]; 
 
	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];
	continue;
   }   
   else if(iscrg==0&&(s->name=='D'||s->name=='E')&&a->name[1]=='O'&&a->tatm->id>4) {
	b=a->bond[0];  
	a1=a;a2=b;
	if(a2==0)continue;

	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.1) {
		if(b->res==0) delete b;b=0;
		continue;
	}
	for(i=0;i<3;i++) xyz[i]=xyz[i]/t;   
	for(i=3;i<6;i++) xyz[i]=a->xyz[i-3];

	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];		
   }
   else if(s->name=='H'&&strcmp(a->name," ND1")==0) {
	b=a->ischildatm(" HD1");
	if(b==0) {			
		Atm *t1=a->res->isatm(" CG "); 	
		Atm *t2=a->res->isatm(" CE1");
		if(t1==0||t2==0) continue;
		b=new Atm;
		for(i=0;i<3;i++)b->xyz[i]=(t1->xyz[i]+t2->xyz[i])/2;
		a1=b;a2=a; 							
	}
	else {
		a1=a;a2=b;
	} 
	
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.01) {
		if(b->res==0) delete b;b=0;
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}
	if(b->res==0) delete b;b=0;
	for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2; 

	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];		
   }
   else if(s->name=='H'&&strcmp(a->name," NE2")==0) {
        continue;
   }
   else if(s->name=='W'&&a->name[1]=='N'&&a->tatm->id>4) {
	b=a->ischildatm(" HE1");
	if(b==0) {			
		Atm *t1=a->res->isatm(" CD1"); 	
		Atm *t2=a->res->isatm(" CE2");
		if(t1==0||t2==0) continue;
		b=new Atm;
		for(i=0;i<3;i++)b->xyz[i]=(t1->xyz[i]+t2->xyz[i])/2;
		a1=b;a2=a; 							
	}
	else {
		a1=a;a2=b;
	} 
	
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.01) {
		if(b->res==0) delete b;b=0;
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}
	for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3];
	if(b->res==0) delete b;b=0;
	//dp0=finddipole(a);
	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];		
   } 					
   else if((s->name=='T'||s->name=='S'||s->name=='Y')&&a->name[1]=='O'&&a->tatm->id>4) {	
	continue; 	 
	b=a->bond[0];
	if(b==0) continue;
	a1=b;a2=a;	 
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.1) {
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}	
	for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]; 

	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];		 		
   }
   else if(r->name=='C'&&a->name[1]=='S') {
	continue;
	b=a->bond[0];
	if(b==0) continue;
	a1=b;a2=a;	 
	for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
	t=TRES.distance(xyz);
	if(t<0.1) {
		continue;
	}
	else {
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
	}	
	for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]; 

	dp0=range[a->id0];
	if(dp0==0) {
		dipole[dip]=new Dipole;
		dipole[dip]->db=0.5;	
		dipole[dip]->atm=a;
		dp0=dipole[dip];
		range[a->id0]=dp0;
		dip++;	
	}
	for(i=0;i<3;i++) dp0->xyz[i]=xyz[i+3];
	for(i=0;i<3;i++) dp0->pot[i]=xyz[i];	
   }	
   else if(iscrg&&r->last==0&&a->tatm->id==0) {	 
	//continue;	
	ch0=findcharge(a);
        if(ch0==0) { 
		charge[crg]=new Charge;
        	charge[crg]->crg=1;
		charge[crg]->atm=a;
		ch0=charge[crg];
		crg++;
	}
        for(i=0;i<3;i++)ch0->xyz[i]=a->xyz[i];				
   }
   else if(iscrg&&r->next==0&&a->tatm->id==3) {
	//continue;	
        b=r->isatmid(2); //c
        if(b==0)continue;
        a1=s->isatmid(1); //ca
        if(a1==0) continue;
	ch0=findcharge(a);
        if(ch0==0) { 
		charge[crg]=new Charge;
        	charge[crg]->crg=-1;
		charge[crg]->atm=a;
		ch0=charge[crg];
		crg++;
	}
        for(i=0;i<3;i++)ch0->xyz[i]=b->xyz[i]+(b->xyz[i]-a1->xyz[i])/2;	 	
   }
   else if(iscrg&&s->name=='D'&&strcmp(a->name," CG ")==0) {
	Atm *a1=a->ischildatm(" OD1");
	Atm *a2=a->ischildatm(" OD2");
	if(a1==0||a2==0) continue;
	ch0=findcharge(a);
	if(ch0==0) { 
		charge[crg]=new Charge;
        	charge[crg]->crg=-1;
		charge[crg]->atm=a;
		ch0=charge[crg];
		crg++;
	} 
	for(i=0;i<3;i++) ch0->xyz[i]=(a->xyz[i]+a2->xyz[i]+a1->xyz[i])/3; 	 	 
   }
   else if(iscrg&&s->name=='E'&&strcmp(a->name," CD ")==0) {
	Atm *a1=a->ischildatm(" OE1");
	Atm *a2=a->ischildatm(" OE2");
	if(a1==0||a2==0) continue;
	ch0=findcharge(a);
	if(ch0==0) {
		charge[crg]=new Charge;
       	 	charge[crg]->crg=-1;
		charge[crg]->atm=a;
		ch0=charge[crg];
		crg++; 
	}	 
	for(i=0;i<3;i++) ch0->xyz[i]=(a->xyz[i]+a2->xyz[i]+a1->xyz[i])/3; 	 	
   } 
   else if(iscrg&&s->name=='K'&&strcmp(a->name," NZ ")==0) {
	ch0=findcharge(a);
	if(ch0==0) {
		charge[crg]=new Charge;
        	charge[crg]->crg=1;
		charge[crg]->atm=a; 
		ch0=charge[crg];
		crg++;
	}
	for(i=0;i<3;i++) ch0->xyz[i]=a->xyz[i]; 
   }
   else if(iscrg&&s->name=='R'&&strcmp(a->name," CZ ")==0) {
	ch0=findcharge(a);
	if(ch0==0) {
		charge[crg]=new Charge;
        	charge[crg]->crg=1;
		charge[crg]->atm=a; 
		ch0=charge[crg];
		crg++;
	}
	for(i=0;i<3;i++) ch0->xyz[i]=a->xyz[i]; 
   }						   
}
if(TRES.logg>3)cerr<<"the total number of charges and dipoles:"<<crg<<" "<<dip<<endl;
}


void Disc::calcenglist(FILE *fpp,char *file,char *org){

	Strhandler te;
        FILE *fp;

	Pdb porg;

	porg.read(org,'1');

        fp=fopen(file,"r");

        if(fp==0) {
                cerr<<"could not open file:"<<file<<endl;
                return;
        }

        char line0[1000];
	
        while(fgets(line0,1000,fp)!=NULL) {
	
		char *line=strdup(line0);
		line=te.clearendchar(line,"\r\t\n ");
		char *s=strchr(line,' ');
		if(s) *s='\0';
               

                Pdb pdb;

                te.lowercase(line);

                pdb.read(line,'1');
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->addhatoms();
		pdb.configure();
		if(TRES.logg>3) pdb.write("b");
		Disc disc;
		disc.cutoff=10;
		disc.pdb=&pdb;
		disc.ready();
		disc.setrange();
		disc.setup(pdb.chn->res,1000,1);
		disc.setdipoleascharge();
		disc.write(stdout);
		strcpy(disc.force,force);
		float d=disc.clash(pdb.chn->res,10000);
		Stralg algn;
		algn.flag=1;
		float e=algn.superimpose(pdb.chn,porg.chn,0);
		fprintf(fpp,"%s %f %f",line,e,d);
		int nn=strlen(disc.force);
		int i;
		for(i=0;i<nn;i++) {			
			fprintf(fpp," %c %f",disc.force[i],disc.energy[disc.force[i]]);
		}
		
		fprintf(fpp,"\n");
		if(line) delete [] line;line=0;
	}

}
