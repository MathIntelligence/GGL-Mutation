#include"source.h"

DistPopular::DistPopular()
{
	dist[0]=100000; dist[1]=-10000;
	popular=0;
	aim=0;
	des=0;
	resn=0;
	tres=0;
	more=0;
	next=0;
	size=0;
	num=0;
	file=0;
	bond=-1;
}

DistPopular::~DistPopular()
{
	if(popular) delete [] popular;
	if(aim) delete [] aim;
	if(des)  delete [] des;
	if(more) delete more;
	if(next) delete next;
	if(file) delete [] file;
}

DistPopular::DistPopular(Tres *tr)
{
        dist[0]=100000; dist[1]=-10000;
        popular=0;
        aim=0;
        des=0;
        resn=0;
        tres=tr;
        more=0;
        next=0;
	size=0;
	num=0;
	file=0;
	bond=-1;
}

void DistPopular::write(FILE *fp) {

	DistPopular *target=0;

	for(DistPopular *target0=this;target0;target0=target0->next) {

	fprintf(fp,"!start residue internal:%c\n",target0->tres->name);

	for(target=target0;target;target=target->more) {
		int m=0;
		int nnn=target->num;
		fprintf(fp,"%c %s %s %2i %4.1f %4.1f %7i",target->tres->name,target->aim,target->des,target->resn,target->dist[0],target->dist[1],nnn);
		for(m=0;m<10;m++) {
			fprintf(fp," %4.1f",target->popular[m]);
		}
		fprintf(fp,"\n");


	}

	fprintf(fp,"!end residue internal:%c\n",target0->tres->name);
	}
}

float DistPopular::findmax() {

	DistPopular *e=0,*target=0;
	float d=0;
        for(DistPopular *target0=this;target0;target0=target0->next) {
        	for(target=target0;target;target=target->more) {
			if(target->dist[1]-target->dist[0]>d) {
				d=target->dist[1]-target->dist[0];
				e=target;
			}
        	}
        }
 	if(TRES.logg>3) {
		if(e&&e->tres) cerr<<"the maximum distance is: "<<e->tres->name<<" "<<e->aim<<" "<<e->des<<" "<<e->resn<<" "<<d<<endl;
		else  cerr<<"the maximum distance is: "<<d<<endl;
	}
	return d;
}

void DistPopular::writelarge(FILE *fp) {

        DistPopular *target=0;

        for(DistPopular *target0=this;target0;target0=target0->next) {

        fprintf(fp,"!start residue internal:%c\n",target0->tres->name);

        for(target=target0;target;target=target->more) {
		if(target->size==0) continue;
                int m=0;
                int nnn=target->num;
                fprintf(fp,"%c %s %s %2i %4.1f %4.1f %7i",target->tres->name,target->aim,target->des,target->resn,target->dist[0],target->dist[1],nnn);
                for(m=0;m<10;m++) {
                        fprintf(fp," %4.1f",target->popular[m]);
                }
                fprintf(fp,"\n");


        }

        fprintf(fp,"!end residue internal:%c\n",target0->tres->name);
        }
}


DistPopular::DistPopular(Tres *tr,char *a,char *b,int n)
{
        dist[0]=100000; dist[1]=-10000;
        popular=0;
        aim=strdup(a);
        des=strdup(b);
        resn=n;
        tres=tr;
        more=0;
        next=0;
	size=0;
	num=0;
	file=0;
}

DistPopular *DistPopular::getDistPopular(Tres *tr,char *a,char *b,int n) {

	
	DistPopular *target0=getDistPopular(tr);

	if(target0==0) return 0;

	return target0->getDistPopular(a,b,n);
}

DistPopular *DistPopular::getDistPopular(Tres *tr) {

       for(DistPopular *target0=this;target0;target0=target0->next) {
		
		if(target0->tres==tr) return target0;
	} 

	return 0;
}


DistPopular *DistPopular::getDistPopular(char *a,char *b,int n) {

	
	for(DistPopular *target0=this;target0;target0=target0->more) {
		if(target0->resn!=n) continue;
		if(target0->aim[1]!=a[1]||target0->des[1]!=b[1]) continue;
		if(a&&strcmp(target0->aim,a)!=0) continue;
		if(b&&strcmp(target0->des,b)!=0) continue;
		return target0;
        }
	return 0; 
}


void DistPopular::read(char *filnam) {

	Strhandler cc;
	char *mlib,line[200];
	FILE *fp;
	int i;
	Tres *tr;
	char nm1[5],nm2[5];
	int nid;
	float dis1,dis2;
	int num0; 
	float pop[10];

	file=strdup(filnam);

	mlib=RCS["library"];

	fp=0;i=0;
	sprintf(line,"%s/%s",mlib,filnam);
	i=strlen(mlib);

	while((fp=fopen(line,"r"))==NULL) {	

		mlib=mlib+i+1;
		i=strlen(mlib);
		if(*mlib=='\0') {
			cerr<<"warning:could not open: "<<filnam<<endl;return;
			cerr<<"check the jackal.dir file first in your current directory then in the environment path:JACKALDIR"<<endl;
        		cerr<<"something wrong in specifying the library"<<endl;
        		cerr<<"make sure the library points to the abosulte path"<<endl;
			exit(0);
		}
		sprintf(line,"%s/%s",mlib,filnam);
	}
	cerr<<endl;
	cerr<<"reading file:"<<line<<endl<<endl;
	while(fgets(line,200,fp)!=NULL) {
		
		cc.clearendemptyspace(line);
		
		if(line[0]=='!') continue;

		tr=TRES[line[0]];

		if(tr==0&&line[0]!='0') continue;

		strncpy(nm1,line+2,4);
		strncpy(nm2,line+7,4);
		nm1[4]='\0';
		nm2[4]='\0';	

		if((line[12]>='A'&&line[12]<='Z')||(line[12]>='a'&&line[12]<='z')) {
			nid=line[12];
			sscanf(line+13,"%f %f %d %f %f %f %f %f %f %f %f %f %f",&dis1,&dis2,&num0,\
                       	pop+0,pop+1,pop+2,pop+3,pop+4,pop+5,pop+6,pop+7,pop+8,pop+9);
		}
		else {
			sscanf(line+12,"%d %f %f %d %f %f %f %f %f %f %f %f %f %f",&nid,&dis1,&dis2,&num0,\
                       	pop+0,pop+1,pop+2,pop+3,pop+4,pop+5,pop+6,pop+7,pop+8,pop+9);
		}
		if(tres==0) tres=tr;
		if(TRES.findanyatmwithname(nm1)==0||TRES.findanyatmwithname(nm2)==0) continue;
		addbounds(tr,nm1,nm2,nid,num0,dis1,dis2,pop,10);
	}
	findmax();
}

void DistPopular::addbounds(Tres *tr,char *a, char *b,int n,int num0,float d0,float d1,float *pop,int siz) {

	//find target with the same tr
	DistPopular *target=0;
	for(target=this;target;target=target->next) {
 		if(target->tres==tr) break;
	}

	//if no such target
	if(target==0) {
		for(target=this;target->next;target=target->next);
		if(target==0) {cerr<<"target is zero"<<endl;return;}
		target->next=new DistPopular(tr,a,b,n);
		target=target->next;
	}

	if(target==0) {
		cerr<<"target should not be zero..."<<endl;
	}

	
	DistPopular *tarmore=0;

	for(tarmore=target;tarmore;tarmore=tarmore->more) {
		
		if(tarmore->aim==0||tarmore->des==0) {
			tarmore->aim=strdup(a);
			tarmore->des=strdup(b);
			tarmore->resn=n;
			break;
		}
		else if(strcmp(tarmore->aim,a)==0&&strcmp(tarmore->des,b)==0&&tarmore->resn==n) {
			break;
		}
	}

	if(tarmore==0) {
		for(tarmore=target;tarmore->more;tarmore=tarmore->more);
		if(tarmore==0) {cerr<<"tarmore is zero"<<endl;return;}
		tarmore->more=new DistPopular(tr,a,b,n);
		tarmore=tarmore->more;
	}

	 
	if(tarmore==0) {
		cerr<<"tarmore is zero"<<endl;
		return;
	}
	tarmore->addbounds(num0,d0,d1,pop,siz);
}

 
void DistPopular::addbounds(int num0,float d0,float d1,float *pop,int siz) {
 
	dist[0]=d0;
	dist[1]=d1;
	num=num0;

	if(popular) delete [] popular;
	size=siz;
	popular=new float[siz];
	for(int i=0;i<siz;i++) popular[i]=pop[i];
}
 
float DistPopular::findmostpopulardistance() {

	float d=(dist[0]+dist[1])/2;
	float x=-100;
	for(int i=0;i<10;i++) {
		if(popular[i]>x){
			x=popular[i];
			d=dist[0]+(dist[1]-dist[0])/10*(i+0.5);
		}
	}
	return d;
}

float *DistPopular::findmostpopulardp(float cut) {

	if(dist[1]-dist[0]<1) return 0;

	int dp[10];
	float tmp[10];
	int i;

	for(i=0;i<10;i++) dp[i]=-1;
	
	int j=0;
	for(i=0;i<10;i++) {
		if(popular[i]>=20) {
			dp[j]=i;
			tmp[j]=popular[i];
			j++;
		}	
	}

	Qsort cc;

	int order[10];

	cc.sort(tmp,j,order);

	float *de=new float[10];

	for(i=0;i<10;i++) de[i]=-1;

	de[0]=tmp[order[0]];

	float e=(dist[1]-dist[0])/10;
	
	int mm=1;
	int op=dp[order[0]];
	for(i=1;i<j;i++) {
		int n=dp[order[i]];
		if(e*abs(n-op)>cut){
			de[mm++]=tmp[order[i]];
			op=dp[order[n]];
		}
	}
	return de;
}


void DistPopular::setbond() {

        DistPopular *target=0;
	Chn chn;
	char *s=TRES.getcappedsequence('A');
	chn.create(s); 
	chn.header();
	chn.allnearbond(3);
	if(s) delete [] s; 
	
        for(DistPopular *target0=this;target0;target0=target0->next) {        
        	for(target=target0;target;target=target->more) {
			Res *r=0,*r0=0;
			Atm *a=0,*a0=0;
			if(target->tres->name!='A') r=chn.isres(target->tres->name,0);
			else r=chn.isres(target->tres->name,1);

			int n=r->id0+target->resn;
			r0=chn.isres0(n);
			if(r0==0) {
				cerr<<"warning: the previous residue not exist in setbond"<<endl;
				continue;
			}
		 	a=r->isatm(target->aim);
			a0=r0->isatm(target->des);
			
			if(a==0||a0==0) {
				cerr<<"warning: the previous residue not exist in setbond"<<endl;
				continue;
			}
			target->bond=a->isnear(a0,4);			
        	}
        }
}
