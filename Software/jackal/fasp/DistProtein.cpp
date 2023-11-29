#include"source.h"

DistProtein::DistProtein()
{
tres=0;
distdb=0;
int i;
for(i=0;i<200;i++) tresnum[i]=0;
}

DistProtein::~DistProtein()
{
if(distdb) delete distdb;
}

void DistProtein::reordernearnext() {

	Qsort cc;
	DistDatabase *target[50000],*target0[50000],*target1[50000];
	int order[50000];
	float temp[50000];
        DistDatabase *s;//,*s0;

	int i=0;
        for(i=0;i<50000;i++) target[i]=0;

        for(s=distdb;s;s=s->next) {
		if(s->tres) target[s->tres->name]=s;
		else target[0]=s;
        }

	s=0;
	for(i=0;i<50000;i++) {
		if(target[i]==0) continue;
		int j=0;
		for(j=0;j<50000;j++) target0[j]=0;
		if(target[i]==0||target[i]->aim==0) continue;
		j=0;
		for(s=target[i];s;s=s->more) {
			target0[j]=s;
			
			if(s->resn==0) temp[j]=-10000;
			else if(s->resn==-1) temp[j]=0;
			else if(s->resn==1) temp[j]=10000;
			else  temp[j]=10000*s->resn;
			Tatm *aa;
			if(s->tres) aa=s->tres->isatm(s->aim);
			else 	    aa=TRES.findanyatmwithname(s->aim);
			Tatm *bb;
			bb=TRES.findanyatmwithname(s->des);
			if(aa) temp[j]+=aa->id*50;
			if(bb) temp[j]+=bb->id;
			j++;
		}

		cc.sort(temp,j,order);
	
		int ii;	
		for(ii=0;ii<j;ii++) {
			int nn=order[ii];
			target1[ii]=target0[nn];
			target1[ii]->more=0;
			target1[ii]->next=0;
		}

		target[i]=target1[0];

		s=target[i];
		for(ii=1;ii<j;ii++) {
			s->more=target1[ii];
			s=s->more;
		}
		s->more=0;
	}
	s=0;
	for(i=0;i<50000;i++) {
                if(target[i]==0) continue;
		if(s==0) { 
			distdb=target[i];
			s=distdb;
		}
		else {
			s->next=target[i];
			s=s->next;
		}
	}
	s->next=0;
}

void DistProtein::buildnearnextdatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		pdb.read(line,c);

		cerr<<line<<" "<<c<<endl;
		Res *r,*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			for(rr=pdb.chn->res;rr;rr=rr->next) {
				if(abs(rr->id-r->id)>=2) continue;

				if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				  
				for(Atm *a=r->atm;a;a=a->next) 
				for(Atm *b=rr->atm;b;b=b->next) {
					
					if(a==b) continue;
					if(rr!=r&&b->tatm->id>4) continue;
					if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					if(distdb==0) distdb = new DistDatabase(r->tres);
					int n=rr->id-r->id;
					distdb->addbounds(r->tres,a->name,b->name,n,d);
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			}
		}
	}

}


void DistProtein::buildnearnextfulldatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		pdb.read(line,c);

		cerr<<line<<" "<<c<<endl;
		Res *r,*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			for(rr=pdb.chn->res;rr;rr=rr->next) {
				if(abs(rr->id-r->id)>=2) continue;

				if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				  
				for(Atm *a=r->atm;a;a=a->next) 
				for(Atm *b=rr->atm;b;b=b->next) {
					
					if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					if(distdb==0) distdb = new DistDatabase(r->tres);
					int n=rr->id-r->id;
					distdb->addbounds(r->tres,a->name,b->name,n,d);
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			}
		}
	}

}

void DistProtein::buildnearnextfulldatabaseh(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		line[4]='_';
		line[5]='m';
		line[6]='o';
		line[7]='d';	
		line[8]='\0';
		te.lowercase(line);

		pdb.read(line,c);

		cerr<<line<<" "<<c<<endl;
		Res *r,*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			for(rr=pdb.chn->res;rr;rr=rr->next) {
				if(abs(rr->id-r->id)>=2) continue;

				if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				  
				for(Atm *a=r->atm;a;a=a->next) 
				for(Atm *b=rr->atm;b;b=b->next) {
					
					if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					if(distdb==0) distdb = new DistDatabase(r->tres);
					int n=rr->id-r->id;
					if(n==1) continue;
					if(n==-1&&a->tatm->name[1]=='H')continue;
					if(n==-1&&b->tatm->name[1]=='H') continue;
					if(n==-1&&b->tatm->id>4) continue;
					distdb->addbounds(r->tres,a->name,b->name,n,d);
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			}
		}
	}

}
void DistProtein::buildfarawaydatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		pdb.read(line,c);

		cerr<<line<<" "<<c<<endl;
		Res *r,*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			for(rr=pdb.chn->res;rr;rr=rr->next) {
				//if(abs(rr->id-r->id)>=2) continue;
				if(r==rr) continue;
				
				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				int f=pdb.chn->isnearnext(r,rr);	
				if(f==1) continue;
				for(Atm *a=r->atm;a;a=a->next) 
				for(Atm *b=rr->atm;b;b=b->next) {
					
					if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					/*if(f) {
						if(a->tatm->id<=4&&b->tatm->id<=4) continue;	
					}*/
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					if(distdb==0) distdb = new DistDatabase(r->tres);
					int n=0;//n=rr->id-r->id;
					distdb->addbounds(r->tres,a->name,b->name,n,d);
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			}
		}
	}

}

void DistProtein::buildfarawaytruedatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		pdb.read(line,c);

		cerr<<line<<" "<<c<<endl;
		Res *r,*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			for(rr=pdb.chn->res;rr;rr=rr->next) {
				if(abs(rr->id-r->id)<=1) continue;
				if(r==rr) continue;
				
				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				int f=pdb.chn->isnearnext(r,rr);	
				if(f==1) continue;
				
				for(Atm *a=r->atm;a;a=a->next) 
				for(Atm *b=rr->atm;b;b=b->next) {
					
					if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					/*if(f) {
						if(a->tatm->id<=4&&b->tatm->id<=4) continue;	
					}*/
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					if(distdb==0) distdb = new DistDatabase(r->tres);
					int n=0;//n=rr->id-r->id;
					distdb->addbounds(r->tres,a->name,b->name,n,d);
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			}
		}
	}

}

void DistProtein::buildallfarawaytruedatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		pdb.read(line,c);

		cerr<<line<<" "<<c<<endl;
		Res *r,*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			for(rr=r->next;rr;rr=rr->next) {
				if(abs(rr->id-r->id)<=1) continue;
				if(r==rr) continue;
				
				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				int f=pdb.chn->isnearnext(r,rr);	
				if(f==1) continue;
				
				for(Atm *a=r->atm;a;a=a->next) 
				for(Atm *b=rr->atm;b;b=b->next) {
					
					if(a==b) continue;
					if(a->id0>=b->id0) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					/*if(f) {
						if(a->tatm->id<=4&&b->tatm->id<=4) continue;	
					}*/
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					if(distdb==0) distdb = new DistDatabase(0);
					int n=0;//n=rr->id-r->id;
					distdb->addbounds(0,a->name,b->name,n,d);
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			}
		}
	}

}


void DistProtein::buildhbonddatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();
		pdb.chn->buildssbond();
		cerr<<line<<" "<<c<<endl;
		Res *r;//*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			//for(rr=pdb.chn->res;rr;rr=rr->next) {
				//if(abs(rr->id-r->id)<=1) continue;
				//if(r==rr) continue;
				
				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				//int f=pdb.chn->isnearnext(r,rr);	
				//if(f==1) continue;
				
				HBondList *h;
				for(h=r->hbond;h;h=h->next) {
				//for(Atm *a=r->atm;a;a=a->next) 
				//for(Atm *b=rr->atm;b;b=b->next) {
					Atm *a,*b;
					//if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					/*if(f) {
						if(a->tatm->id<=4&&b->tatm->id<=4) continue;	
					}*/

					

					if(h->donor->res==r) {a=h->donor;b=h->acceptor;}
					else 	{b=h->donor;a=h->acceptor;}
	
					if(a==0||b==0) continue;
					if(fabs(a->res->id-b->res->id)<=1) continue;
					//if(a->res->sec!=b->res->sec||a->res->sec=='-') continue;
					Atm *a0,*b0;
				
					if(a->tatm->id==0) a0=a->bond[1];
					else 	a0=a->bond[0];

					if(b->tatm->id==0) b0=b->bond[1];       
                                        else    b0=b->bond[0];  
					
					if(a0==0||b0==0) continue;	
					if(a->id0>=b->id0) continue;
					if(distdb==0) distdb = new DistDatabase(r->tres);
					int n=0;//n=rr->id-r->id;
					
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a->tatm->name,b->tatm->name,n,d);

					
					d=TRES.distance(a0->xyz,b->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a0->tatm->name,b->tatm->name,n,d);

					d=TRES.distance(a->xyz,b0->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a->tatm->name,b0->tatm->name,n,d);

					d=TRES.distance(a0->xyz,b0->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a0->tatm->name,b0->tatm->name,n,d);

					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			//}
		}
	}

}


void DistProtein::buildloopdatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;
		
		te.lowercase(line);
		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();		
		cerr<<line<<" "<<c<<endl;

		Res *r,*r0,*r1; 
 
		r=pdb.chn->res;

		
		
		r0=0;r1=0;
		while(r) {
			 
			if(r->sec!='-') {				
				r0=r;
				r=r->next;
				continue; 
			}
			if(r0==0||r0->sec=='-') {
				
				r0=r;
				r=r->next;
				continue;
			}
			 
			Res *r2=0;
			for(r1=r;r1;r1=r1->next) {
				if(r1->sec!='-') break;	
				r2=r1;		 
			}
			if(r1==0) break;
			int n=r2->id-r->id+1;
			Atm *a1=r0->atm->next;
			Atm *a2=r1->atm->next;
			if(a1==0) a1=r0->atm;
			if(a2==0) a2=r1->atm;
			if(a1==0||a2==0) {
				r0=r1;
				r=r1->next;
				continue;
			}
			float d=TRES.distance(a1,a2);
			if(distdb==0) distdb = new DistDatabase(0);
			int dd=(int)d;
			char name[5];
			sprintf(name,"%4i",dd);name[4]='\0';
			distdb->addbounds(0,name,"   ",0,n*1.0);
			r0=r1;
			r=r1->next;
			continue;
		}
	}

}

void DistProtein::buildnewloopdatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;
		
		te.lowercase(line);
		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();		
		cerr<<line<<" "<<c<<endl;

		Res *r,*r0,*r1; 
 
		r=pdb.chn->res;

		
		
		r0=0;r1=0;
		while(r) {
			 
			if(r->sec!='-') {				
				r0=r;
				r=r->next;
				continue; 
			}
			if(r0==0||r0->sec=='-') {
				
				r0=r;
				r=r->next;
				continue;
			}
			 
			Res *r2=0;
			for(r1=r;r1;r1=r1->next) {
				if(r1->sec!='-') break;	
				r2=r1;		 
			}
			if(r1==0) break;
			int n=r2->id-r->id+1;
			Atm *a1=r0->atm->next;
			Atm *a2=r1->atm->next;
			if(a1==0) a1=r0->atm;
			if(a2==0) a2=r1->atm;
			if(a1==0||a2==0) {
				r0=r1;
				r=r1->next;
				continue;
			}
			float d=TRES.distance(a1,a2);
			if(distdb==0) distdb = new DistDatabase(0);
			//int dd=d;
			char name[5];
			sprintf(name,"%4i",n);name[4]='\0';
			char name0[5];
			if(a1->res->sec=='e'&&a2->res->sec=='e') {
				strcpy(name0,"ee  ");
			} 
			else if(a1->res->sec=='e'&&a2->res->sec=='h') {
				strcpy(name0,"he  ");
			} 
			else if(a1->res->sec=='h'&&a2->res->sec=='h') {
				strcpy(name0,"hh  ");
			} 
			else if(a1->res->sec=='h'&&a2->res->sec=='e') {
				strcpy(name0,"he  ");
			} 
			else {
				cerr<<"warning...no match of secondary structure\n";
			}
			distdb->addbounds(0,name,name0,0,d);
			r0=r1;
			r=r1->next;
			continue;
		}
	}

}
void DistProtein::buildthenewloopdatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;
		
		te.lowercase(line);
		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();		
		cerr<<line<<" "<<c<<endl;

		Res *r,*r0,*r1; 
 
		r=pdb.chn->res;

		
		
		r0=0;r1=0;
		while(r) {
			 
			if(r->sec!='-') {				
				r0=r;
				r=r->next;
				continue; 
			}
			if(r0==0||r0->sec=='-') {
				
				r0=r;
				r=r->next;
				continue;
			}
			 
			Res *r2=0;
			for(r1=r;r1;r1=r1->next) {
				if(r1->sec!='-') break;	
				r2=r1;		 
			}
			if(r1==0) break;
			int n=r2->id-r->id+1;
			Atm *a1=r0->atm->next;
			Atm *a2=r1->atm->next;
			if(a1==0) a1=r0->atm;
			if(a2==0) a2=r1->atm;
			if(a1==0||a2==0) {
				r0=r1;
				r=r1->next;
				continue;
			}
			float d=TRES.distance(a1,a2);
			if(distdb==0) distdb = new DistDatabase(0);
			//int dd=d;
			char name[5];
			sprintf(name,"%4i",n);name[4]='\0';
			char name0[5];
			if(a1->res->sec=='e'&&a2->res->sec=='e') {
				strcpy(name0,"ee  ");
			} 
			else if(a1->res->sec=='e'&&a2->res->sec=='h') {
				strcpy(name0,"eh  ");
			} 
			else if(a1->res->sec=='h'&&a2->res->sec=='h') {
				strcpy(name0,"hh  ");
			} 
			else if(a1->res->sec=='h'&&a2->res->sec=='e') {
				strcpy(name0,"he  ");
			} 
			else {
				cerr<<"warning...no match of secondary structure\n";
			}
			distdb->addbounds(0,name,name0,0,d);
			r0=r1;
			r=r1->next;
			continue;
		}
	}

}
void DistProtein::buildallnewloopdatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;
		
		te.lowercase(line);
		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();		
		cerr<<line<<" "<<c<<endl;

		Res *r,*r0,*r1; 
 
		r=pdb.chn->res;

		
		
		r0=0;r1=0;
		while(r) {
			 
			if(r->sec!='-') {				
				r0=r;
				r=r->next;
				continue; 
			}
			if(r0==0||r0->sec=='-') {
				
				r0=r;
				r=r->next;
				continue;
			}
			 
			Res *r2=0;
			for(r1=r;r1;r1=r1->next) {
				if(r1->sec!='-') break;	
				r2=r1;		 
			}
			if(r1==0) break;
			int n=r2->id-r->id+1;
			Atm *a1=r0->atm->next;
			Atm *a2=r1->atm->next;
			if(a1==0) a1=r0->atm;
			if(a2==0) a2=r1->atm;
			if(a1==0||a2==0) {
				r0=r1;
				r=r1->next;
				continue;
			}
			float d=TRES.distance(a1,a2);
			if(distdb==0) distdb = new DistDatabase(0);
			//int dd=d;
			char name[5];
			sprintf(name,"%4i",n);name[4]='\0';
			char name0[5]; 
			if(a1->res->sec=='e'&&a2->res->sec=='e') {
				strcpy(name0,"0   ");
			} 
			else if(a1->res->sec=='e'&&a2->res->sec=='h') {
				strcpy(name0,"0   ");
			} 
			else if(a1->res->sec=='h'&&a2->res->sec=='h') {
				strcpy(name0,"0   ");
			} 
			else if(a1->res->sec=='h'&&a2->res->sec=='e') {
				strcpy(name0,"0   ");
			} 
			else {
				cerr<<"warning...no match of secondary structure\n";
			}
			distdb->addbounds(0,name,name0,0,d);
			r0=r1;
			r=r1->next;
			continue;
		}
	}

}


void DistProtein::buildallhbonddatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();
		pdb.chn->buildssbond();
		cerr<<line<<" "<<c<<endl;
		Res *r;//*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			//for(rr=pdb.chn->res;rr;rr=rr->next) {
				//if(abs(rr->id-r->id)<=1) continue;
				//if(r==rr) continue;
				
				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				//int f=pdb.chn->isnearnext(r,rr);	
				//if(f==1) continue;
				
				HBondList *h;
				for(h=r->hbond;h;h=h->next) {
				//for(Atm *a=r->atm;a;a=a->next) 
				//for(Atm *b=rr->atm;b;b=b->next) {
					Atm *a,*b;
					//if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					/*if(f) {
						if(a->tatm->id<=4&&b->tatm->id<=4) continue;	
					}*/

					

					if(h->donor->res==r) {a=h->donor;b=h->acceptor;}
					else 	{b=h->donor;a=h->acceptor;}
	
					if(a==0||b==0) continue;
					if(fabs(a->res->id-b->res->id)<=1) continue;
					//if(a->res->sec!=b->res->sec||a->res->sec=='-') continue;
					Atm *a0,*b0;
				
					if(a->tatm->id==0) a0=a->bond[1];
					else 	a0=a->bond[0];

					if(b->tatm->id==0) b0=b->bond[1];       
                                        else    b0=b->bond[0];  
					
					if(a0==0||b0==0) continue;	
					if(a->id0>=b->id0) continue;
					if(distdb==0) distdb = new DistDatabase(0);
					int n=0;//n=rr->id-r->id;
					
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					distdb->addbounds(0,a->tatm->name,b->tatm->name,n,d);

					
					d=TRES.distance(a0->xyz,b->xyz);
					if(d==0) continue;
					distdb->addbounds(0,a0->tatm->name,b->tatm->name,n,d);

					d=TRES.distance(a->xyz,b0->xyz);
					if(d==0) continue;
					distdb->addbounds(0,a->tatm->name,b0->tatm->name,n,d);

					d=TRES.distance(a0->xyz,b0->xyz);
					if(d==0) continue;
					distdb->addbounds(0,a0->tatm->name,b0->tatm->name,n,d);

					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			//}
		}
	}

}

void DistProtein::indbuildallhbonddatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();
		pdb.chn->buildssbond();
		cerr<<line<<" "<<c<<endl;
		Res *r;//*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			tresnum[r->name]++;
			if(tres&&r->tres!=tres) continue;
			//for(rr=pdb.chn->res;rr;rr=rr->next) {
				//if(abs(rr->id-r->id)<=1) continue;
				//if(r==rr) continue;
				
				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				//int f=pdb.chn->isnearnext(r,rr);	
				//if(f==1) continue;
				
				HBondList *h;
				for(h=r->hbond;h;h=h->next) {
				//for(Atm *a=r->atm;a;a=a->next) 
				//for(Atm *b=rr->atm;b;b=b->next) {
					Atm *a,*b;
					//if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					/*if(f) {
						if(a->tatm->id<=4&&b->tatm->id<=4) continue;	
					}*/

					

					if(h->donor->res==r) {a=h->donor;b=h->acceptor;}
					else 	{b=h->donor;a=h->acceptor;}
	
					if(a==0||b==0) continue;
					if(fabs(a->res->id-b->res->id)<=1) continue;
					//if(a->res->sec!=b->res->sec||a->res->sec=='-') continue;
					Atm *a0,*b0;
				
					if(a->tatm->id==0) a0=a->bond[1];
					else 	a0=a->bond[0];

					if(b->tatm->id==0) b0=b->bond[1];       
                                        else    b0=b->bond[0];  
					
					if(a0==0||b0==0) continue;	
					if(a->id0>=b->id0) continue;
					if(distdb==0) distdb = new DistDatabase(0);
					int n=b->tatm->tres->name;//n=rr->id-r->id;
					
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a->tatm->name,b->tatm->name,b->tatm->tres->name,d);

					
					d=TRES.distance(a0->xyz,b->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a0->tatm->name,b->tatm->name,b->tatm->tres->name,d);

					d=TRES.distance(a->xyz,b0->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a->tatm->name,b0->tatm->name,b->tatm->tres->name,d);

					d=TRES.distance(a0->xyz,b0->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a0->tatm->name,b0->tatm->name,b->tatm->tres->name,d);

					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<b->tatm->tres->name<<" "<<d<<endl;
				}
			//}
		}
	}

}

void DistProtein::indONbuildallhbonddatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();
		pdb.chn->buildssbond();
		cerr<<line<<" "<<c<<endl;
		Res *r;//*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			tresnum[r->name]++;
			if(tres&&r->tres!=tres) continue;
			//for(rr=pdb.chn->res;rr;rr=rr->next) {
				//if(abs(rr->id-r->id)<=1) continue;
				//if(r==rr) continue;
				
				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				//int f=pdb.chn->isnearnext(r,rr);	
				//if(f==1) continue;
				
				HBondList *h;
				for(h=r->hbond;h;h=h->next) {
				//for(Atm *a=r->atm;a;a=a->next) 
				//for(Atm *b=rr->atm;b;b=b->next) {
					Atm *a,*b;
					//if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					/*if(f) {
						if(a->tatm->id<=4&&b->tatm->id<=4) continue;	
					}*/

					

					if(h->donor->res==r) {a=h->donor;b=h->acceptor;}
					else 	{b=h->donor;a=h->acceptor;}
	
					if(a==0||b==0) continue;
					if(fabs(a->res->id-b->res->id)<=1) continue;
					//if(a->res->sec!=b->res->sec||a->res->sec=='-') continue;
					Atm *a0,*b0;
				
					if(a->tatm->id==0) a0=a->bond[1];
					else 	a0=a->bond[0];

					if(b->tatm->id==0) b0=b->bond[1];       
                                        else    b0=b->bond[0];  
					
					if(a0==0||b0==0) continue;	
					if(a->id0>=b->id0) continue;
					if(distdb==0) distdb = new DistDatabase(0);
					int n=b->tatm->tres->name;//n=rr->id-r->id;
					
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a->tatm->name,b->tatm->name,b->tatm->tres->name,d);

					/*
					d=TRES.distance(a0->xyz,b->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a0->tatm->name,b->tatm->name,b->tatm->tres->name,d);

					d=TRES.distance(a->xyz,b0->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a->tatm->name,b0->tatm->name,b->tatm->tres->name,d);

					d=TRES.distance(a0->xyz,b0->xyz);
					if(d==0) continue;
					distdb->addbounds(r->tres,a0->tatm->name,b0->tatm->name,b->tatm->tres->name,d);
					*/
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<b->tatm->tres->name<<" "<<d<<endl;
				}
			//}
		}
	}

}

void DistProtein::buildneighboratomdatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];
	
	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';
		line[4]='\0';

		Pdb pdb;
		
		te.lowercase(line);
		
		pdb.read(line,c);
		if(pdb.chn==0) continue;
		cerr<<line<<" "<<c<<endl;

		if(pdb.chn==0) continue;	
		LatSurfv latsurfv;
		latsurfv.ready(&pdb,2,1.4);
		latsurfv.calcarea(&pdb);		

		//puton lattice
		Lattice lat;
		lat.putoff();
		lat.flag=1; 
		lat.grdsiz=2.;
		lat.radall=15.;
		lat.ready(&pdb);
		lat.puton(&pdb);
		if(pdb.manychain()>1) {	
			cerr<<"too many chains"<<endl;	
			continue;
		}
 
		Res *r; 
		Atm *a;					
		for(r=pdb.chn->res;r;r=r->next) 
		for(a=r->atm;a;a=a->next) {
			lat.getcell(a,8);
  	 		int nget;
			int m;
			nget=lat.nget;
			int e=0;
			for(m=0;m<nget;m++) {
				Atm *b=lat.obtain[m]->atm;
				float d=TRES.distance(a,b);
				float r=a->tatm->eng->radius+b->tatm->eng->radius;
				if(d<r+0.0) e++;
			}
			
			e=e/5*5;
			if(e>=100) {
				cerr<<"the percentage >100"<<endl;
				continue;
			}

			float x=(a->tatm->eng->radius+1.4)*(1.4+a->tatm->eng->radius)*12.56; 
			float y=a->area/x*100;
			 
			char ss[10];
			sprintf(ss,"%4i",e);
						 
			if(distdb==0) distdb = new DistDatabase(0);
			distdb->addbounds(r->tres,a->tatm->name,ss,'0',y);			
		}		
	}

}

void DistProtein::testhbonddatabase(char *file) {

	Strhandler te;
	/*
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char *line[1000];

	//while(fgets(line,1000,fp)!=NULL) {
	*/
	char *line=file;
	line=file;
 	if(1){
		char c=line[4];

		//if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		
		pdb.read(line,'0');
		if(pdb.chn==0) return;
		pdb.chn->header();
		pdb.chn->buildhbond();		
		pdb.chn->setdsspstr();
		pdb.chn->setthreestatesec();
		pdb.chn->buildssbond();
		pdb.chn->addhatoms();
		pdb.chn->write("a");
		cerr<<line<<" "<<c<<endl;
		Res *r;//*rr;
		pdb.setlastresorder();
		if(pdb.chn==0) return;

		for(r=pdb.chn->res;r;r=r->next) {
			tresnum[r->name]++;
			if(tres&&r->tres!=tres) continue;
			//for(rr=pdb.chn->res;rr;rr=rr->next) {
				//if(abs(rr->id-r->id)<=1) continue;
				//if(r==rr) continue;
				
				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				//int f=pdb.chn->isnearnext(r,rr);	
				//if(f==1) continue;
				
				HBondList *h;
				for(h=r->hbond;h;h=h->next) {
					if(h->donor->tatm->eng->dipole==0||h->acceptor->tatm->eng->dipole==0) continue;
					float dd=TRES.distance(h->donor,h->acceptor);
					if(dd>3.5||dd<2.5) continue;
					if(h->donor->tatm->eng->dipole==0||h->acceptor->tatm->eng->dipole==0) continue;
				  	float *x=TRES.getdipolevector(h->donor,1);
					float *y=TRES.getdipolevector(h->acceptor,1);
					if(x==0) {
						cerr<<h->donor->res->name<<" "<<h->donor->name<<" does not exist!"<<endl;
						continue;
					}
					if(y==0) {
						cerr<<h->acceptor->res->name<<" "<<h->acceptor->name<<" does not exist!"<<endl;
						continue;
					}
					float z=x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; 
					 
					cerr<<h->donor->res->name<<h->donor->res->id+pdb.chn->start<<" "<<h->donor->name;
					cerr<<h->acceptor->res->name<<h->acceptor->res->id+pdb.chn->start<<" "<<h->acceptor->name<<" :"<<acos(z)/3.14*180<<endl;
					float e1=TRES.distance(x);
					float e2=TRES.distance(x);
					cerr<<h->donor->res->name<<" "<<h->donor->name;
					cerr<<h->acceptor->res->name<<" "<<h->acceptor->name<<" :"<<e1<<" "<<e2<<endl;
				}
			//}
		}
	}

}
void DistProtein::writetresnum(FILE *fp) {

	int i;

	for(i=0;i<200;i++) {
		if(tresnum[i]==0) continue;
		fprintf(fp,"%c %6i\n",char(i),tresnum[i]);
	}
}


void DistProtein::buildssyesdatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		pdb.read(line,c);

		cerr<<line<<" "<<c<<endl;
		Res *r,*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			if(r->tres->name!='C') continue;
			for(rr=r;rr;rr=rr->next) {
				if(rr->tres->name!='C') continue;
				if(rr==r) continue;				
				if(abs(rr->id-r->id)<=1) continue;

				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				if(pdb.chn->isssbond(r,rr)==0) continue;  
				
				for(Atm *a=r->atm;a;a=a->next) 
				for(Atm *b=rr->atm;b;b=b->next) {
					
					if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					if(distdb==0) distdb = new DistDatabase(r->tres);
					//int n=rr->id-r->id;
					int n=0;
					distdb->addbounds(r->tres,a->name,b->name,n,d);
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			}
		}
	}

}

void DistProtein::buildssnodatabase(char *file) {

	Strhandler te;
	FILE *fp;

	fp=fopen(file,"r");

	if(fp==0) {
		cerr<<"could not open file:"<<file<<endl;
		return;
	}

	char line[1000];

	while(fgets(line,1000,fp)!=NULL) {

		char c=line[4];

		if(c=='0') c='1';

		line[4]='\0';

		Pdb pdb;

		
		te.lowercase(line);

		pdb.read(line,c);

		cerr<<line<<" "<<c<<endl;
		Res *r,*rr;
		
		if(pdb.chn==0) continue;

		for(r=pdb.chn->res;r;r=r->next) {
			if(tres&&r->tres!=tres) continue;
			if(r->tres->name!='C') continue;
			for(rr=r;rr;rr=rr->next) {
				if(rr->tres->name!='C') continue;
				if(rr==r) continue;				
				if(abs(rr->id-r->id)<=1) continue;

				//if(!(rr==r||pdb.chn->isnearnext(r,rr))) continue;
				if(pdb.chn->isssbond(r,rr)==1) continue;  
				
				for(Atm *a=r->atm;a;a=a->next) 
				for(Atm *b=rr->atm;b;b=b->next) {
					
					if(a==b) continue;
					//if(rr!=r&&b->tatm->id>4) continue;
					//if(rr==r&&b->tatm->id<=a->tatm->id) continue;
					
					float d=TRES.distance(a->xyz,b->xyz);
					if(d==0) continue;
					if(distdb==0) distdb = new DistDatabase(r->tres);
					//int n=rr->id-r->id;
					int n=0;
					distdb->addbounds(r->tres,a->name,b->name,n,d);
					//cerr<<r->id<<" "<<a->id<<" "<<rr->id<<" "<<b->id<<" "<<n<<" "<<d<<endl;
				}
			}
		}
	}

}

