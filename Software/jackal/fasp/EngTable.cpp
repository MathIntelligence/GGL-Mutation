#include"source.h"
EngTable::EngTable()
{
	initial();
}

void EngTable::initial() {
	id=-1;
	size=0;
	engtable=0;
	numsim=-1;
	next=0;
        more=0;
        table=0;
        expand=1;
        wall=40;
        resolution=100;
        distance=15;
        name=0;
        energy=0;
	radius=0;
	charge=0;
	epslon=0;
}
EngTable::EngTable(EngTable *s)
{
	initial();     
        expand=s->expand;
        wall=s->wall;
        resolution=s->resolution;
        distance=s->distance;
	numsim=s->numsim;	
}


EngTable::~EngTable()
{
	if(engtable) delete [] engtable;
	if(more) delete more;more=0;
	if(table) delete table;table=0; 
 	if(next) delete next;next=0;
	if(name) delete name;name=0;
	if(energy) delete energy;energy=0;
}




void EngTable::createwalltable(float dist) {


	float r;
	EngTable *e=0;	
	for(r=dist;r>0;r-=0.2){
		if(e==0) {
                        e=this;
                }
                else {
                        e->more=new EngTable(this);
                        e=e->more;
                }
		e->distance=r;
		e->resolution=r/100;		 		 
		e->createwalltable(); 	
		//exit(0);	
	}
}
 
void EngTable::createwalltable() {

	//create table for each pair atoms interaction
 
	int leng=(int)((distance+distance)/resolution+1);
	table=new float[leng];
	energy=new float[leng];

	int i;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
	
	//treatment of radius to favor hydrogen bond and ss bond
	FILE *fp=0,*fp1=0;
	if(distance<0.5&&distance>0.3) {
 	   fp=fopen("s2","w");fp1=fopen("s3","w");
	}
	float r0;
	for(i=0;i<leng;i++) {

		r0=resolution*(i+0.5);
		if(r0<0) {
			table[i]=2*wall*r0;
                        energy[i]=wall*r0*r0;
		}
		else if(r0>distance) {
			
			table[i]=2*wall*(r0-distance);
			energy[i]=wall*(r0-distance)*(r0-distance);
		} 
		if(fp) fprintf(fp,"%f %f\n",r0,table[i]);
		if(fp1) fprintf(fp1,"%f %f\n",r0,energy[i]);
	}
	size=leng;

	if(fp) fclose(fp);
	if(fp1) fclose(fp1);
}

void EngTable::createdefwalltable(float dist) {


	 
	EngTable *e=0;	
	float r;
	for(r=dist;r>0;r-=0.2){
		if(e==0) {
                        e=this;
                }
                else {
                        e->more=new EngTable(this);
                        e=e->more;
                }
		e->distance=r;
		e->resolution=r/100;		 		 
		e->createdefwalltable(); 	
		//exit(0);	
	}
}
 
void EngTable::createdefwalltable() {

	//create table for each pair atoms interaction
 
	int leng=(int)((distance+distance)/resolution+1);
	table=new float[leng];
	energy=new float[leng];

	int i,j;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
	
	
	//treatment of radius to favor hydrogen bond and ss bond
	FILE *fp=0;
	FILE *fp1=0;
	if(distance>0.3&&distance<0.5) {
 		fp=fopen("s2","w");
		fp1=fopen("s3","w");
	}
	float r0,r,d,e0,f,g,r1;
	for(i=-leng/2;i<leng;i++) {

		r0=resolution*(i+0.5);

		float nom=0;

		float denom=0;
		float flat=0;
		int nu=(int)((distance)/(0.5*resolution)+1);
	 
		for(j=-4*nu;j<0;j++) {
                        r=0.5*resolution*(j+0.5);
                        d=r*r*wall;			
                        f=d*exp(-expand*(r-r0)*(r-r0));
                        e0=f;
			denom+=e0;
                        g=e0*(-2*expand)*(r0-r); 
			flat+=exp(-expand*(r-r0)*(r-r0));		 
                        nom+=g;
                } 
				
		int note=(int)(distance/(0.5*resolution)+0.5);
		for(j=note;j<note+4*nu;j++) {			
			r=0.5*resolution*(j+0.5)-distance;
                        d=r*r*wall;			
			r1=0.5*resolution*(j+0.5);
                        f=d*exp(-expand*(r1-r0)*(r1-r0));
                        e0=f; 
			denom+=e0;
                        g=e0*(-2*expand)*(r0-r1); 
			flat+=exp(-expand*(r1-r0)*(r1-r0));
                        nom+=g;
                }  
		if(i>=0) {
			table[i]=nom/flat;
			energy[i]=denom/flat;
		}
		if(fp) {
			fprintf(fp,"%f %f\n",r0,nom/flat);
			fprintf(fp1,"%f %f\n",r0,denom/flat);
		}
	}	
	if(fp) {
		fclose(fp);
		fclose(fp1);
	}
	size=leng;

}

void EngTable::createcolonywalltable(float dist) {


	EngTable *e=0;	
	float r;
	for(r=dist;r>0;r-=0.2){
		if(e==0) {
                        e=this;
                }
                else {
                        e->more=new EngTable(this);
                        e=e->more;
                }
		e->distance=r;
		//e->distance=1.2;
		e->resolution=r/100;		 		 
		e->createcolonywalltable(); 	
		cerr<<r<<" "<<e->resolution<<endl;
		//exit(0);
	}
}
 
void EngTable::createcolonywalltable() {

	//create table for each pair atoms interaction
 
	int leng=(int)((distance+distance)/resolution+1);
	table=new float[leng];
	energy=new float[leng];

	int i,j;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
	
	
	//treatment of radius to favor hydrogen bond and ss bond
	
	FILE *fp=0,*fp1=0;
	if(distance<0.5&&distance>0.3) {
 		fp=fopen("s2","w");
		fp1=fopen("s3","w");
	}
	float r0,r,e0,f,g,r1;
	float gc=1/(8.31*0.298);
	for(i=-leng/2;i<leng;i++) {

		r0=resolution*(i+0.5);

		float nom=0;

		float denom=0;
		float flat=0;
		int nu=(int)(4*(distance)/(0.5*resolution)+1);
		//cerr<<nu<<" "<<i<<endl; 
		for(j=-nu;j<0;j++) {
                        r=0.5*resolution*(j+0.5);                       		
                        f=exp(-gc*r*r*wall-expand*(r-r0)*(r-r0));
                        e0=f;
			denom+=e0;
                        g=e0*(-2*expand)*(r0-r); 
			flat+=exp(-expand*(r-r0)*(r-r0));		 
                        nom+=g;
                } 
				
		
		int note=(int)(distance/(0.5*resolution)+0.5);
		for(j=0;j<note;j++) {
			r=0.5*resolution*(j+0.5);
                        f=exp(-expand*(r-r0)*(r-r0));
                        e0=f;
                        denom+=e0;
                        g=e0*(-2*expand)*(r0-r);
                        flat+=exp(-expand*(r-r0)*(r-r0));
                        nom+=g;
		}
		for(j=note;j<note+nu;j++) {			
			r=0.5*resolution*(j+0.5)-distance;                  			
			r1=0.5*resolution*(j+0.5);
                        f=exp(-gc*r*r*wall-expand*(r1-r0)*(r1-r0));
                        e0=f; 
			denom+=e0;
                        g=e0*(-2*expand)*(r0-r1); 
			flat+=exp(-expand*(r1-r0)*(r1-r0));
                        nom+=g;
                }  
		if(i>=0){
		table[i]=-1/gc*nom/denom;
		energy[i]=-1/gc*log(denom)/flat;
		}
		if(fp)fprintf(fp,"%f %f %f\n",r0,-1/gc*nom/denom,distance);
		if(fp1)fprintf(fp1,"%f %f %f\n",r0,-1/gc*log(denom)/flat,distance);
	}
	size=leng;
	if(fp) fclose(fp);
	if(fp1) fclose(fp1);
}

void EngTable::createcolonytable() { 
 
	EngTable *e=0;
	SimAtm *a,*b;
	for(a=TRES.simatm;a;a=a->next)
	for(b=TRES.simatm;b;b=b->next){
		if(e==0) {
			e=this;
		}
		else {
			e->more=new EngTable(this);
			e=e->more;
		}
		e->createcolonytable(a->id,b->id);
		//return;
	}
} 

void EngTable::createcolonytable(int a,int b) {

	//create table for each pair atoms interaction

	
	SimAtm *sa=TRES.simatm->find(a);
	SimAtm *sb=TRES.simatm->find(b);
	 
	int leng=(int)(distance/resolution+1);
	table=new float[leng];
	energy=new float[leng];

	int i,j;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
	
	//treatment of radius to favor hydrogen bond and ss bond
	float dr=sa->radius+sb->radius;
	float eps=sqrt(sa->epslon*sb->epslon);
	int ndr=0;
	if(sa->name=='N'||sa->name=='O') ndr++;
	if(sb->name=='N'||sb->name=='O') ndr++;	
	if(ndr==2)  dr=3.0;
	else if(sa->name=='S'&&sb->name=='S') dr=2.0;
	else if(sa->name=='S'&&sb->name=='C') dr=3.0;
	else if(sa->name=='C'&&sb->name=='S') dr=3.0;
	else if(sa->name=='H'&&sb->name=='H') dr=1.0;
	else if(sa->name=='H') dr-=sa->radius/2;
	else if(sb->name=='H') dr-=sb->radius/2;
 	
	float crg=sa->charge*sb->charge*332/20;
	float ec=1/(8.31*0.298);
	float r0,r,e,d,e0,f,g;
	
	
	FILE *fp,*fp0;
	fp=fopen("s","w");
	fp0=fopen("s1","w");//eps=1;
	
	for(i=0;i<leng;i++) {
		r0=resolution*(i+0.5);	
		float nom=0;
		float denom=0;
		int nu=(int)((distance+distance)/(resolution*0.5));
		float flat=0;
		//float zero=0;
		for(j=0;j<nu;j++) {
			r=0.5*resolution*(j+0.5);
			d=dr/r;
			d=d*d;
			d=d*d*d; 
			e=-(eps*(d*d-2*d)+crg/r)*ec;
			f=e-expand*(r-r0)*(r-r0);
			//f=e-expand*(r*r-2*r*r0);
			//if(j==0) zero=f;
			e0=exp(f);
			//e0=exp(f-zero);
			
			denom+=e0; 			
			g=e0*(-2*expand)*(r0-r);
			nom+=g;	
			flat+=exp(-expand*(r-r0)*(r-r0));		
			//flat+=exp(-expand*(r*r-2*r*r0)); 
		}		
		table[i]=-nom/denom/ec;
		energy[i]=-log(denom/flat)/ec;
		fprintf(fp,"%f %f\n",r0,table[i]);
		fprintf(fp0,"%f %f %f\n",r0,-log(denom/flat)/ec,denom/flat);
		//cout<<r0<<" "<<table[note][i]<<" "<<denom<<" "<<endl;		
	}
	size=leng;
	if(fp) {
		fclose(fp);
		fclose(fp0);
	}
}


void EngTable::createdeftable() {//simple diffusion table
 
	EngTable *e=0;

	SimAtm *a,*b;
	for(a=TRES.simatm;a;a=a->next)
	for(b=TRES.simatm;b;b=b->next){
		if(e==0) {
			e=this;
		}
		else {
			e->more=new EngTable(this);
			e=e->more;
		}
		e->createdeftable(a->id,b->id);
		// exit(0);
	}
} 

void EngTable::createdeftable(int a,int b) {

	//create table for each pair atoms interaction

	SimAtm *sa=TRES.simatm->find(a);
	SimAtm *sb=TRES.simatm->find(b);
	
	
	int leng=(int)(distance/resolution+1);
	table=new float[leng];
	energy=new float[leng];

	int i,j;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
	 
	
	//treatment of radius to favor hydrogen bond and ss bond
	float dr=sa->radius+sb->radius;
	float eps=sqrt(sa->epslon*sb->epslon);
	int ndr=0;
	if(sa->name=='N'||sa->name=='O') ndr++;
	if(sb->name=='N'||sb->name=='O') ndr++;	
	if(ndr==2)  dr=3.0;
	else if(sa->name=='S'&&sb->name=='S') dr=2.0;
	else if(sa->name=='S'&&sb->name=='C') dr=3.0;
	else if(sa->name=='C'&&sb->name=='S') dr=3.0;
	else if(sa->name=='H'&&sb->name=='H') dr=1.0;
	else if(sa->name=='H') dr-=sa->radius/2;
	else if(sb->name=='H') dr-=sb->radius/2;
 	
	float crg=sa->charge*sb->charge*332/20;
	//float ec=1/(8.31*0.298);
	float r0,r,e,d,e0,f,g;
	
	//FILE *fp,*fp0;
	//fp=fopen("s","w");
	//fp0=fopen("s1","w");//eps=1;
	for(i=0;i<leng;i++) {
		r0=resolution*(i+0.5);
		
		float nom=0;
		float denom=0;
		int nu=(int)((distance+distance)/(resolution*0.5));
		float flat=0;
		for(j=0;j<nu;j++) {
			r=0.5*resolution*(j+0.5);
			d=dr/r;
			//if(d>2) d=2;
			d=d*d;
			d=d*d*d; 
			e=eps*(d*d-2*d)+crg/r;
			
			//e=-(eps*(d*d-2*d)+crg/r/r)*ec;
			//f=e-expand*(r-r0)*(r-r0);
			f=-expand*(r-r0)*(r-r0);
			e0=e*exp(f);
			denom+=e0; 			
			g=e0*(-2*expand)*(r0-r);
			nom+=g;	
			flat+=exp(-expand*(r-r0)*(r-r0));		
		}
		
		table[i]=nom/flat;
		energy[i]=denom/flat;
		//fprintf(fp,"%f %f\n",r0,table[note][i]);
		//fprintf(fp0,"%f %f\n",r0,denom/flat);
		//cout<<r0<<" "<<table[note][i]<<" "<<denom<<" "<<endl;		
	}
	size=leng;
}


void EngTable::createvdwtable() {//simple diffusion table

	 

	EngTable *e=0;

	SimAtm *a,*b;
	for(a=TRES.simatm;a;a=a->next)
	for(b=TRES.simatm;b;b=b->next){
		if(e==0) {
                        e=this;
                }
                else {
                        e->more=new EngTable(this);
                        e=e->more;
                }

		e->createvdwtable(a->id,b->id);
		//exit(0);
	}
} 

void EngTable::createvdwtable(int a,int b) {

	//create table for each pair atoms interaction

	SimAtm *sa=TRES.simatm->find(a);
	SimAtm *sb=TRES.simatm->find(b);
	int leng=(int)(distance/resolution+1);

	table=new float[leng];
	energy=new float[leng];

	int i;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
 
	//treatment of radius to favor hydrogen bond and ss bond
	float dr=sa->radius+sb->radius;
	float eps=sqrt(sa->epslon*sb->epslon);
	int ndr=0;
	if(sa->name=='N'||sa->name=='O') ndr++;
	if(sb->name=='N'||sb->name=='O') ndr++;	
	if(ndr==2)  dr=3.0;
	else if(sa->name=='S'&&sb->name=='S') dr=2.0;
	else if(sa->name=='S'&&sb->name=='C') dr=3.0;
	else if(sa->name=='C'&&sb->name=='S') dr=3.0;
	else if(sa->name=='H'&&sb->name=='H') dr=1.0;
	else if(sa->name=='H') dr-=sa->radius/2;
	else if(sb->name=='H') dr-=sb->radius/2;
 	
	float crg=sa->charge*sb->charge*332/20;
	//float ec=1/(8.31*0.298);
	float r,e,d,f;
	
	/*
	FILE *fp,*fp0;
	fp=fopen("s","w");
	fp0=fopen("s1","w");//eps=1;
	*/
	for(i=0;i<leng;i++) {
		r=resolution*(i+0.5);				 
		d=dr/r;		
		d=d*d;
		d=d*d*d; 
		e=eps*12*(-d*d+d)/r-crg/r;
 		f=eps*(d*d-2*d)+crg/r;
		table[i]=e;
		energy[i]=f;
		//fprintf(fp,"%f %f\n",r,table[i]);
		//fprintf(fp0,"%f %f\n",r,f);
		//cout<<r0<<" "<<table[note][i]<<" "<<denom<<" "<<endl;		
	}
	size=leng;
}



void EngTable::createsoftvdwtable() {//simple diffusion table

	 
	EngTable *e=0;
	SimAtm *a,*b;
	for(a=TRES.simatm;a;a=a->next)
	for(b=TRES.simatm;b;b=b->next){
		if(e==0) {
                        e=this;
                }
                else {
                        e->more=new EngTable(this);
                        e=e->more;
                }

		e->createsoftvdwtable(a->id,b->id);
		//return;
	}
} 

void EngTable::createsoftvdwtable(int a,int b) {

	//create table for each pair atoms interaction

	SimAtm *sa=TRES.simatm->find(a);
	SimAtm *sb=TRES.simatm->find(b);
	//int n=TRES.simatm->getnumber();
	//int note=a*n+b;
	int leng=(int)(distance/resolution+1);
	table=new float[leng];
	energy=new float[leng];

	int i;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
	 
	
	//treatment of radius to favor hydrogen bond and ss bond
	float dr=sa->radius+sb->radius;
	float eps=sqrt(sa->epslon*sb->epslon);
	int ndr=0;
	if(sa->name=='N'||sa->name=='O') ndr++;
	if(sb->name=='N'||sb->name=='O') ndr++;	
	if(ndr==2)  dr=3.0;
	else if(sa->name=='S'&&sb->name=='S') dr=2.0;
	else if(sa->name=='S'&&sb->name=='C') dr=3.0;
	else if(sa->name=='C'&&sb->name=='S') dr=3.0;
	else if(sa->name=='H'&&sb->name=='H') dr=1.0;
	else if(sa->name=='H') dr-=sa->radius/2;
	else if(sb->name=='H') dr-=sb->radius/2;
 	
	float crg=sa->charge*sb->charge*332/20;
	//float ec=1/(8.31*0.298);
	float r,e,d,f,g;

	/*	
	FILE *fp,*fp0;
	fp=fopen("s","w");
	fp0=fopen("s1","w");//eps=1;
	*/
	
	float ta,tb,tc,tn;
	i=4*TRES.smt;
  	ta=TRES.smooth[i+0];
  	tb=TRES.smooth[i+1];
  	tc=TRES.smooth[i+2];
  	tn=TRES.smooth[i+3];

	for(i=0;i<leng;i++) {
		r=resolution*(i+0.5);				 
		d=dr/r;	
		g=pow(d,tn);
		e=eps*ta*exp(-tb/d/d)*(g*g-tc*g)*(-2*tb/dr/dr*r);
		e+=eps*ta*exp(-tb*d*d)*(-2*tn*g*g/r+tn*tc*g/r);
		e+=-crg/r;
		f=eps*ta*exp(-tb/d/d)*(g*g-tc*g)+crg/r;	
		table[i]=e;
		energy[i]=f;
		//fprintf(fp,"%f %f\n",r,table[i]);
		//fprintf(fp0,"%f %f\n",r,energy[i]);
		//cout<<r0<<" "<<table[note][i]<<" "<<denom<<" "<<endl;		
	}
	size=leng;
	//fclose(fp);
	//fclose(fp0);
}

void EngTable::createdefsoftvdwtable() {//simple diffusion table

	 
	EngTable *e=0;
	SimAtm *a,*b;
	for(a=TRES.simatm;a;a=a->next)
	for(b=TRES.simatm;b;b=b->next){
		if(e==0) {
                        e=this;
                }
                else {
                        e->more=new EngTable(this);
                        e=e->more;
                }

		e->createdefsoftvdwtable(a->id,b->id);
		//exit(0);
	}
} 

void EngTable::createdefsoftvdwtable(int a,int b) {

	//create table for each pair atoms interaction

	SimAtm *sa=TRES.simatm->find(a);
	SimAtm *sb=TRES.simatm->find(b);
	//int n=TRES.simatm->getnumber();
	//int note=a*n+b;
	int leng=(int)(distance/resolution+1);
	table=new float[leng];
	energy=new float[leng];

	int i,j;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
	 
	
	//treatment of radius to favor hydrogen bond and ss bond
	float dr=sa->radius+sb->radius;
	float eps=sqrt(sa->epslon*sb->epslon);
	int ndr=0;
	if(sa->name=='N'||sa->name=='O') ndr++;
	if(sb->name=='N'||sb->name=='O') ndr++;	
	if(ndr==2)  dr=3.0;
	else if(sa->name=='S'&&sb->name=='S') dr=2.0;
	else if(sa->name=='S'&&sb->name=='C') dr=3.0;
	else if(sa->name=='C'&&sb->name=='S') dr=3.0;
	else if(sa->name=='H'&&sb->name=='H') dr=1.0;
	else if(sa->name=='H') dr-=sa->radius/2;
	else if(sb->name=='H') dr-=sb->radius/2;
 	
	float crg=sa->charge*sb->charge*332/20;
	//float ec=1/(8.31*0.298);
	float r0,r,d,e0,g;
	
	//FILE *fp,*fp0;
	//fp=fopen("s","w");
	//fp0=fopen("s1","w");//eps=1;
	
	float ta,tb,tc,tn;
	i=4*TRES.smt;
  	ta=TRES.smooth[i+0];
  	tb=TRES.smooth[i+1];
  	tc=TRES.smooth[i+2];
  	tn=TRES.smooth[i+3];

	for(i=0;i<leng;i++) {

		r0=resolution*(i+0.5);
		
		float nom=0;
		float denom=0;
		int nu=(int)((distance+distance)/(resolution*0.5));
		float flat=0;
		for(j=0;j<nu;j++) {
			r=0.5*resolution*(i+0.5);				 
			d=dr/r;	
			g=pow(d,tn);			
			e0=eps*ta*exp(-tb/d/d)*(g*g-tc*g)+crg/r;	
			denom+=e0; 
			g=e0*(-2*expand)*(r0-r);
			nom+=g;	
			flat+=exp(-expand*(r-r0)*(r-r0));	
		}
		table[i]=nom/flat;
		energy[i]=denom/flat;
		//fprintf(fp,"%f %f\n",r0,table[note][i]);
		//fprintf(fp0,"%f %f\n",r0,denom/flat);
		//cout<<r0<<" "<<table[note][i]<<" "<<denom<<" "<<endl;				
	}
	size=leng;
}

void EngTable::createcolonysoftvdwtable() {//simple diffusion table

	 
	EngTable *e=0;
	SimAtm *a,*b;
	for(a=TRES.simatm;a;a=a->next)
	for(b=TRES.simatm;b;b=b->next){
		if(e==0) {
                        e=this;
                }
                else {
                        e->more=new EngTable(this);
                        e=e->more;
                }

		e->createcolonysoftvdwtable(a->id,b->id);
		//exit(0);
	}
} 

void EngTable::createcolonysoftvdwtable(int a,int b) {

	//create table for each pair atoms interaction

	SimAtm *sa=TRES.simatm->find(a);
	SimAtm *sb=TRES.simatm->find(b);
	//int n=TRES.simatm->getnumber();
	//int note=a*n+b;
	int leng=(int)(distance/resolution+1);
	table=new float[leng];
	energy=new float[leng];

	int i,j;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
	 
	
	//treatment of radius to favor hydrogen bond and ss bond
	float dr=sa->radius+sb->radius;
	float eps=sqrt(sa->epslon*sb->epslon);
	int ndr=0;
	if(sa->name=='N'||sa->name=='O') ndr++;
	if(sb->name=='N'||sb->name=='O') ndr++;	
	if(ndr==2)  dr=3.0;
	else if(sa->name=='S'&&sb->name=='S') dr=2.0;
	else if(sa->name=='S'&&sb->name=='C') dr=3.0;
	else if(sa->name=='C'&&sb->name=='S') dr=3.0;
	else if(sa->name=='H'&&sb->name=='H') dr=1.0;
	else if(sa->name=='H') dr-=sa->radius/2;
	else if(sb->name=='H') dr-=sb->radius/2;
 	
	float crg=sa->charge*sb->charge*332/20;
	float ec=1/(8.31*0.298);
	float r0,r,e,d,e0,f,g;
	
	//FILE *fp,*fp0;
	//fp=fopen("s","w");
	//fp0=fopen("s1","w");//eps=1;
	
	float ta,tb,tc,tn;
	i=4*TRES.smt;
  	ta=TRES.smooth[i+0];
  	tb=TRES.smooth[i+1];
  	tc=TRES.smooth[i+2];
  	tn=TRES.smooth[i+3];

	for(i=0;i<leng;i++) {

		r0=resolution*(i+0.5);
		
		float nom=0;
		float denom=0;
		int nu=(int)((distance+distance)/(resolution*0.5));
		float flat=0;
		for(j=0;j<nu;j++) {

			r=0.5*resolution*(j+0.5);
			d=dr/r;
			//d=d*d;
			//d=d*d*d;
			g=pow(d,tn);			
			e0=eps*ta*exp(-tb/d/d)*(g*g-tc*g)+crg/r;	 
			e=-e0*ec;
			f=e-expand*(r-r0)*(r-r0);
			e0=exp(f);
			denom+=e0; 			
			g=e0*(-2*expand)*(r0-r);
			nom+=g;	
			flat+=exp(-expand*(r-r0)*(r-r0));				 
		}
		table[i]=nom;
		energy[i]=denom/flat;
		//fprintf(fp,"%f %f\n",r0,table[note][i]);
		//fprintf(fp0,"%f %f\n",r0,denom/flat);
		//cout<<r0<<" "<<table[note][i]<<" "<<denom<<" "<<endl;				
	}
	size=leng;

}


void EngTable::chargecurve(FILE *fp,float f){


	for(int i=0;i<1000;i++) {

		float d=(i+0.5)*0.01;
		float e=f/d;
		float p=exp(-e/8.31/0.3);
		fprintf(fp,"%f %f %f\n",d,p,e);
	}
}

void EngTable::vdwcurve(FILE *fp,float f){


        for(int i=0;i<1000;i++) {

                float d=(i+0.5)*0.01;
		float f1=f/d;
                float e=pow(f1,12)-2*pow(f1,6);
                float p=exp(-e/8.31/0.3);
                fprintf(fp,"%f %f %f\n",d,p,e);
        }
}
 

EngTable **EngTable::getallengtable() {

	EngTable **out=0;

	int n=getnumber();

	out=new EngTable*[n+1];
	for(int i=0;i<n+1;i++) out[i]=0;

	int k=0; 
	for(EngTable *e=this;e;e=e->more) {
		out[e->id]=e;
		if(e->id!=k) {
			cerr<<"warning..."<<endl;
		}
		k++;
	}	
	return out;
}

int EngTable::getnumber() {
	int n=0;
	for(EngTable *e=this;e;e=e->more) n++;
	return n;
}
void EngTable::setchargezero() {

	SimAtm *a;
        for(a=TRES.simatm;a;a=a->next) a->charge=0;

}

float EngTable::getvalue(int a,int b,float d) {

	if(numsim==-1) numsim=TRES.simatm->getnumber();
	if(engtable==0) engtable=getallengtable();

	int n=a*numsim+b;
	EngTable *aa=engtable[n];
	if(aa==0) return 0;

	n=(int)(d/resolution+0.5);
	if(n>=size) n=size-1;
	return aa->table[n];
}

float EngTable::getenergy(int a,int b,float d) {
	
	if(numsim==-1) numsim=TRES.simatm->getnumber();
        if(engtable==0) engtable=getallengtable();
	
	int n=a*numsim+b;
        EngTable *aa=engtable[n];
        if(aa==0) return 0;

        n=(int)(d/resolution+0.5);
	if(n>=size) n=size-1;
        return aa->energy[n];
}


float EngTable::getvalue(float dis,float d) {

	float dis0=dis;

	if(engtable==0) engtable=getallengtable();
	
	if(dis>distance) dis=distance;
	int n=(int)((distance-dis)/0.2+0.5);
	
	EngTable *a=engtable[n];

	float d0=d;	
	if(d<0) d=a->distance-d;
	n=(int)(d/a->resolution+0.5);
	
	if(n>=200) {
		n=199;
		//if(d0<0) return d0*2*wall;
		//else     return (d0-a->distance)*2*wall;
	}
	
	if(d0<0) return -a->table[n];
	else return a->table[n];
}
	

float EngTable::getenergy(float dis,float d) {

	float dis0=dis;

        if(engtable==0) engtable=getallengtable();

        if(dis>distance) dis=distance;
        int n=(int)((distance-dis)/0.2+0.5);

        EngTable *a=engtable[n];

        float d0=d;
        if(d<0) d=a->distance-d;
        n=(int)(d/a->resolution+0.5);

        if(n>=200) {
                n=199;
                //if(d0<0) return d0*d0*wall;
                //else     return (d0-a->distance)*(d0-a->distance)*wall;
        }
        if(d0<0) return a->energy[n];
        else return a->energy[n];
}

void EngTable::setvdwtable() {//simple diffusion table
 
	EngTable *e=0;

	if(numsim==-1) numsim=TRES.simatm->getnumber();

	SimAtm *a,*b;

	for(a=TRES.simatm;a;a=a->next)
	for(b=TRES.simatm;b;b=b->next){
		if(e==0) {
                        e=this;
                }
                else {
                        e->more=new EngTable(this);
                        e=e->more;
                }
		e->id=a->id*numsim+b->id;
		e->setvdwtable(a->id,b->id);	
	}
} 

void EngTable::setvdwtable(int a,int b) {

	//create table for each pair atoms interaction

	SimAtm *sa=TRES.simatm->find(a);
	SimAtm *sb=TRES.simatm->find(b);
	int leng=(int)(distance*resolution+1);

	table=new float[leng];
	energy=new float[leng];

	int i;

	for(i=0;i<leng;i++) {
		table[i]=0;
		energy[i]=0;
	}
 
	//treatment of radius to favor hydrogen bond and ss bond
	float dr=sa->radius+sb->radius;
	float eps=sqrt(sa->epslon*sb->epslon);
	int ndr=0;
	if(sa->name=='N'||sa->name=='O') ndr++;
	if(sb->name=='N'||sb->name=='O') ndr++;	
	if(ndr==2)  dr=3.0;
	else if(sa->name=='S'&&sb->name=='S') dr=2.0;
	else if(sa->name=='S'&&sb->name=='C') dr=3.0;
	else if(sa->name=='C'&&sb->name=='S') dr=3.0;
	else if(sa->name=='H'&&sb->name=='H') dr=1.0;
	else if(sa->name=='H') dr-=sa->radius/2;
	else if(sb->name=='H') dr-=sb->radius/2; 	
	float crg=sa->charge*sb->charge*332/4;
	
	radius=dr;
	charge=crg;
	epslon=eps;

	//float ec=1/(8.31*0.298);
	float r,e,d,f;
	
	float reso=1./resolution;	
	//FILE *fp,*fp0;
	//fp=fopen("s","w");
	//fp0=fopen("s1","w");//eps=1;
	
	for(i=0;i<leng;i++) {
		//r=(i+0.5)*reso;	
		if(i==0) r=(i+0.5)*reso;
		else 	 r=i*reso;
		energy[i]=calcvdwvalue(radius,charge,epslon,r);	
		table[i]=calcvdwtable(radius,charge,epslon,r);
		/*
		e=TRES.engcoeff->getvdw1storder(r/radius);
		d=TRES.engcoeff->getcrg1storder(r);
		table[i]=eps*e/(radius*r*radius)+charge*d/(radius*r);			 
		e=TRES.engcoeff->getvdwvalue(r/radius);
		d=TRES.engcoeff->getcrgvalue(r);
		energy[i]=epslon*e+charge*d;
		*/
		//fprintf(fp,"%f %f %f\n",r*radius,energy[i],table[i]);
		//fprintf(fp0,"%f %f\n",r*radius,f);
		//cout<<r0<<" "<<table[note][i]<<" "<<denom<<" "<<endl;		
	}
	//fclose(fp);//fclose(fp0);
	size=leng;
}

float EngTable::calcvdwvalue(float radius,float charge,float epslon,float r){

	float e=TRES.engcoeff->getvdwvalue(r/radius);
	float d=TRES.engcoeff->getcrgvalue(r);
	float t=epslon*e+charge*d;
	return t;
}

float EngTable::calcvdwtable(float radius,float charge,float epslon,float r){

        float e=TRES.engcoeff->getvdw1storder(r/radius);
        float d=TRES.engcoeff->getcrg1storder(r);
        float t=epslon*e/(radius*r)+charge*d/r;
        return t;
}

float EngTable::gettablevalue(int a,int b,float d) {

	if(numsim==-1) numsim=TRES.simatm->getnumber();
	if(engtable==0) engtable=getallengtable();

	int n=a*numsim+b;
	EngTable *aa=engtable[n];
	if(aa==0) return 0;

	float e,x;
        e=d*resolution;
        n=(int)e;
        x=e-n;

	if(n>=size-1) n=size-2;
	return aa->table[n]*(1-x)+aa->table[n+1]*x;
}

float EngTable::getenergyvalue(int a,int b,float d) {
	
	int n=a*numsim+b;
        EngTable *aa=engtable[n];
        if(aa==0) return 0;

	float e,x;
        e=d*resolution;
	n=(int)e;
	x=e-n;

	if(n>=size-1) n=size-2;
        return aa->energy[n]*(1-x)+aa->energy[n+1]*x;
}

void EngTable::preparecall() {

	if(numsim==-1) numsim=TRES.simatm->getnumber();
	if(engtable==0) engtable=getallengtable();
}
