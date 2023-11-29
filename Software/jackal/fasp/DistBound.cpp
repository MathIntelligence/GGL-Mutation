#include "source.h"

DistBound::DistBound() {

	atoms=0;
	predt=0;
	natom=0;
	improper=0;
	nimp=0;
	dist=0;
	atmpair=0;
	ndist=0;
	arraysize=0;
	model=0;
	seed=1901;
	flag=0;
	next=0;
	weight=0;
        chiropt=2;
	tag=0;
	distsecond=0;
}

DistBound::~DistBound() {
	if(dist)  delete [] dist;
	if(atmpair) delete [] atmpair;
	if(atoms) delete [] atoms;
	if(improper) delete [] improper;
	if(predt) delete [] predt;
	if(next) delete next;
	if(weight) delete [] weight;
	if(tag) delete [] tag;
	if(model&&(flag&TRES.constant->pdbundeletable)==0)delete model;
	if(distsecond) delete distsecond;
}
	
void DistBound::settag() {

	
	int i,i1;

	cerr<<"set tags ... "<<endl;

	if(natom==0) return; 

	tag=new int[2*natom];

	for(i=0;i<2*natom;i++) tag[i]=-1;

	for(i=0;i<ndist;i++) {
		i1=atmpair[2*i];
		if(tag[2*i1]==-1) {tag[2*i1]=i;tag[2*i1+1]=i;}
		else		  tag[2*i1+1]=i;	
	}
}

void DistBound::setmodel(Pdb *s) {

	model=s;
}
	

void DistBound::ready() {

	Chn *c;
        Res *r;
        Atm *a;

	srandom(seed);

	if(model==0) {
		cerr<<"model is null"<<endl;
		exit(0);
	}

	int n=model->manyatm();
	atoms=new Atm*[n];
	predt=new int[n];
	natom=n;

	int i;

	for(i=0;i<n;i++) {predt[i]=1;atoms[i]=0;}
	
	n=0;
        for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next)
        for(a=r->atm;a;a=a->next) {
		a->allnear(3,0);
		if(a->id0>=natom) {
			cerr<<"size overflow"<<endl;
			cerr<<"the atom id0 is too large:"<<a->name<<a->id0<<endl;
			exit(0);
		}
		if(atoms[a->id0]) {
			cerr<<"two atoms have the same sequential number"<<endl;
			exit(0);
		}
		atoms[a->id0]=a;
		n++;
        }

	//the atom number does not match

	if(n!=natom) {
		cerr<<"atom number does not match"<<endl;
		exit(0);	
	}

	arraysize=10000;
	dist=new float[3*arraysize];
	atmpair=new int[2*arraysize];
        weight=new float[arraysize];
	ndist=0;
}

void DistBound::setsize(int n) {

	arraysize=n;
	dist=(float *)realloc(dist,3*sizeof(float)*n);
	atmpair=(int *)realloc(atmpair,2*sizeof(int)*n);
	weight=(float *)realloc(weight,sizeof(float)*n);
	//predt=(int *)realloc(predt,2*sizeof(int)*n);
}


void DistBound::setstartsize(int n) {

	arraysize=n;

	dist=(float *)realloc(dist,3*sizeof(float)*n);
        atmpair=(int *)realloc(atmpair,2*sizeof(int)*n);
        weight=(float *)realloc(weight,sizeof(float)*n);
        predt=(int *)realloc(predt,1*sizeof(int)*n);

	for(int i=0;i<n;i++) {
		int j=0;
		for(j=0;j<3;j++) dist[i*3+j]=0;
		for(j=0;j<2;j++) atmpair[i*2+j]=0;
		for(j=0;j<1;j++) predt[i*2+j]=1;
	}

	ndist=0;
}

int DistBound::readdist(char *filename) {

	FILE *fp=fopen(filename,"r");

	if(fp==0) {
		cerr<<"file: "<<filename<<" does not exist!"<<endl;
		return ndist;
	}
	
	cerr<<"reading the distance constraint from file:"<<filename<<endl;

	char line[1000];

	int norg=ndist;
	
	char *s;
	Strhandler cc;
	int i1,i2;
	float d1,d2,d3;
	while(fgets(line,256,fp)!=NULL) {
		s=strdup(line);		
		s=cc.clearendchar(s," ");
		if(s==0||s[0]=='#'||s[0]=='!') {
			if(s) delete [] s; s=0;
			continue;
		}
		sscanf(s,"%i %i %f %f %f",&i1,&i2,&d1,&d2,&d3);
		if(s) delete [] s; s=0;
		if(model) {
			Atm *a=model->getatmbyoid(i1);
			if(a) i1=a->id0;
			else  i1=i1;
			a=model->getatmbyoid(i2);
			if(a) i2=a->id0;
			else  i2=i2;
		}
		addbounds(i1,i2,d1,d2,d3,1);
		//if(i1<i2) addbounds(i1-1,i2-1,d1,d2,d3,2);
		//else      addbounds(i2-1,i1-1,d1,d2,d3,2);
	}
	cerr<<"the total "<<ndist<<" distance constraints read from the file "<<filename<<endl;
	return norg;
}

void DistBound::addbounds(Atm *a1,Atm *a2,float d,float w) {

	float d0=TRES.distance(a1,a2);

	float d1=d0-d;
	float d2=d0;
	float d3=d0+d;

	addbounds(a1->id0,a2->id0,d2,d1,d3,w);
	//if(a1->id0<a2->id0) addbounds(a1->id0,a2->id0,d2,d1,d3,w);
	//else addbounds(a2->id0,a1->id0,d2,d1,d3,w);
}

void DistBound::addbounds(int i1,int i2,float d1,float d2,float d3,float w) {

	 int nn=ndist;
	
	 /*
	 if(nn+10>arraysize) {
                setsize(10000+arraysize);
         }
	 */

	 if(i1>i2) {
		int t=i2;
		i1=i2;
		i2=t;
	 }
	 else if(i1==i2) {
		cerr<<"the distance bound constraint should be in the same atoms"<<endl;
		return;
	 }

	 atmpair[nn*2]=i1;
         atmpair[nn*2+1]=i2;
         dist[nn*3]=d1;
         dist[nn*3+1]=d2;
         dist[nn*3+2]=d3;
	 weight[nn]=w;
		
	 
         if(nn+10>arraysize) {
                setsize(10000+arraysize);
         }
         
         nn++;
	 ndist=nn;
	 cerr<<ndist<<" "<<d1<<" "<<d2<<" "<<d3<<endl;
}

void DistBound::getring(Atm **ring,Res *r) {

	int i=0;
	for(Atm *a=r->atm;a;a=a->next) {
		if(a->tatm->ring) ring[i++]=a;
	}
	ring[i]=0;
}

int DistBound::setquasiringdist() {

	Chn *c;
        Res *r;

	int norg=ndist;
	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next){
	
		if(strchr("FYW",r->name)) {
			quasiringdist(r);
		}
	}
	return norg;
}

int DistBound::setmodpdbquasiringdist() {

        Chn *c;
        Res *r;

        int norg=ndist;
        for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next){

                if(strchr("FYW",r->name)) {
                        quasiringdist(r);
                }
        }
        return norg;
}


int DistBound::quasiringdist(Res *r) {

	if(strchr("FYW",r->name)==0) return ndist;
 
	int norg=ndist;
	Atm *ring[100];

	getring(ring,r);

	Atm *ca;

	ca=r->isatm(" CA ");

	for(int i=0;i<100;i++) {
		if(ring[i]==0) break;
		addbounds(ca,ring[i],0.3,1); 
	}	
	return norg;
}

int DistBound::modpdbquasiringdist(Res *r) {

        if(strchr("FYW",r->name)==0) return ndist;

        int norg=ndist;
        Atm *ring[100];

        getring(ring,r);

        Atm *ca;

        ca=r->isatm(" CA ");

        for(int i=0;i<100;i++) {
                if(ring[i]==0) break;
		
                addbounds(ca,ring[i],0.3,1);
        }
        return norg;
}


int DistBound::isnextres(Res *r1,Res *r2){

	if(r2->id-r1->id!=1) return 0;

	Atm *n,*c;

	n=r2->isatmid(0);

	c=r1->isatmid(2);

	if(n==0||c==0) return 0;
	float d=TRES.distance(n->xyz,c->xyz);

	if(d>2.0) return 0;

	return 1;

}


int DistBound::setphipsidist() {

	Chn *c;
        Res *r,*r0;

	int norg=ndist;
	for(c=model->chn;c;c=c->next) {
		r0=0;
        	for(r=c->res;r;r=r->next){
	
		 	if(r0==0) continue;
			
			if(isnextres(r0,r)==0) continue;

			Atm *c0,*c;			
			c0=r0->isatmid(2);
			c=r->isatmid(2);
			addbounds(c0,c,0.3,1); 
			r0=r;
		}		
	}
	return norg;
}


int DistBound::set14dist() {

	Chn *c;
        Res *r;

	int norg=ndist;
	for(c=model->chn;c;c=c->next) {
		 
        	for(r=c->res;r;r=r->next){	
		 	 set14dist(r);
		}		
	}
	return norg;
}

int DistBound::set14dist(Res *r) {

	Atm *a,*a1;

	int norg=ndist;
	for(a=r->atm;a;a=a->next)
	for(a1=a->next;a1;a1=a1->next) {

		if(a->isnear(a1,4)!=4) continue;
		addbounds(a,a1,0.5,1); 
	}
	return norg;
}
	


int DistBound::setn14dist() {

        Chn *c;
        Res *r,*r0;

	int norg=ndist;
        for(c=model->chn;c;c=c->next) {
		r0=0;
                for(r=c->res;r;r=r->next){
			if(r0!=0&&isnextres(r0,r)==1) {
				setn14dist(r0,r);
			}
		 	r0=r;
                }
        }
	return norg;
}

int DistBound::setn14dist(Res *r0,Res *r) {

        Atm *a;
	int norg=ndist;
	Atm *n=r->isatmid(0);
        for(a=r->atm;a;a=a->next){

		if(a->name[1]=='H'&&a->bond[0]->tatm->id>4) {
                	addbounds(a,n,0.5,1);
		}
		else if(a->name[1]!='H'&&a->tatm->id>4) {
			addbounds(a,n,0.5,1);
		}
        }
	return norg;
}


int DistBound::setsidechainrestriction() {

	Chn *c;
        Res *r;

	int norg=ndist;
	for(c=model->chn;c;c=c->next) {
		 
        	for(r=c->res;r;r=r->next){	
		 	 setsidechainrestriction(r);
		}		
	}
	return norg;
}

int DistBound::setsidechainrestriction(Res *r) {

	Atm *a,*o;

	int norg=ndist;
	o=r->isatmid(3);
	if(o==0) return norg;
	for(a=r->atm;a;a=a->next){
	 	if(a->name[1]=='H'&&a->bond[0]->tatm->id>4) {
			addbounds(a,o,0.5,1); 
		}
		else if(a->tatm->id>4) {
			addbounds(a,o,0.5,1); 
		}
	}
	return norg;
}



int DistBound::setinternaldist() {

	if(model==0) return ndist;

	Chn *c;
	Res *r;
	Atm *a;

	int j;
	int norg=ndist;

	Atm *ring[100];

	//1-2 and 1-3 atom pairs

	for(c=model->chn;c;c=c->next)
	for(r=c->res;r;r=r->next)
	for(a=r->atm;a;a=a->next) 
	for(j=0;j<10000&&a->near[j];j++) {
		if(a->id0>=a->near[j]->id0) continue;
		if(a->isbond(a->near[j])) {
			float x=TRES.distance(a,a->near[j]);
			addbounds(a->id0,a->near[j]->id0,x,x-0.03,x+0.03,1);
		}
		else {
			float x=TRES.distance(a,a->near[j]);
                        addbounds(a->id0,a->near[j]->id0,x,x-0.3,x+0.3,1);
		}	
	}

	//ring restrictions
	
	int i,n;

	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		if(strchr("YFWHP",r->name)) {
			for(i=0;i<100;i++) ring[i]=0;
			getring(ring,r);
		}
		else continue;
		
		n=0;
		while(ring[n]) n++;

		for(i=0;i<n;i++) 
		for(j=0;j<n;j++) {
			if(ring[i]->id0>=ring[j]->id0) continue;		 
			if(ring[i]->isnear(ring[j])) continue;
			addbounds(ring[i],ring[j],0.05,1);
		}
	}

	//residue sequential restrictions

	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		if(r->next==0) continue;

		for(i=0;i<100;i++) ring[i]=0;
		getsequential(ring,r);	

		n=0;
                while(ring[n]) n++;

                for(i=0;i<n;i++)
                for(j=0;j<n;j++) {
                        if(ring[i]->id0>=ring[j]->id0) continue;
                        if(ring[i]->isnear(ring[j])) continue;
                        addbounds(ring[i],ring[j],0.05,1);
                }
	}
	return norg;
}



float DistBound::evaluate() {

	float d=0;
	int n=0;
	for(int i=0;i<ndist;i++) {

		int n1,n2;
		n1=atmpair[i*2];
		n2=atmpair[i*2+1];	

		float e=TRES.distance(atoms[n1]->xyz,atoms[n2]->xyz);
		if(e*e>dist[i*3+2]) e=(e-sqrt(dist[i*3+2]));
		else if(e*e<dist[i*3+1]) e=(sqrt(dist[i*3+1])-e);
		else e=0;	
		d+=e*e;		
		n++;
	}
	
	if(n==0) return 0;
	
	d=sqrt(d/n);

	return d;
}

void DistBound::getsequential(Atm **ring,Res *r) {

	int i=0;

	/*
	ring[i++]=r->next->isatm("  C ");
        ring[i++]=r->next->isatm("  O ");
        ring[i++]=r->next->isatm(" CA ");
	*/

	Atm *c,*n;

	c=r->isatm("  C ");
	n=r->next->isatm("  N ");
	if(c==0||n==0) {
		return;
		//cerr<<"warning!..zero distance "
	}
	float d=TRES.distance(c,n);

	if(r->id-r->next->id!=-1||d>2) {
		ring[0]=0;
		ring[1]=0;
		return;
	} 

	ring[i++]=r->isatm("  C ");
        ring[i++]=r->isatm("  O ");
        ring[i++]=r->isatm(" CA ");


	ring[i++]=r->next->isatm("  N ");
	ring[i++]=r->next->isatm(" HN ");
	ring[i++]=r->next->isatm(" CA ");
	i=0; 
	for(int j=0;j<6;j++) {
		if(ring[j]==0) continue;
		ring[i++]=ring[j];	
	}
	ring[i]=0;
}

void DistBound::deletebound(int n) {

	int j=0;
	int k;
	for(int i=0;i<ndist;i++) {
		if(i==n) continue;
		for(k=0;k<3;k++) {
                        dist[j*3+k]=dist[i*3+k];
                }

		for(k=0;k<2;k++) {
                        atmpair[j*2+k]=atmpair[i*2+k];
                }

		j++;
	}
	ndist=j;
}

void DistBound::deletefixedbound(int n) {

        int k;
        int j=0;
        for(int i=0;i<ndist;i++) {
                if(predt[atmpair[2*i]]==0&&predt[atmpair[2*i+1]]==0) continue;
                for(k=0;k<3;k++) {
                        dist[j*3+k]=dist[i*3+k];
                }
                for(k=0;k<2;k++) {
                        atmpair[j*2+k]=atmpair[i*2+k];
                }
                j++;
        }
        ndist=j;
}


void DistBound::reorder() {

//reoder the dist bounds from small to large

	float *temp;
	int   *order;
	int i,j,k,i1,i2;

	cerr<<"sorting the distance constraints...."<<endl;

	if(ndist==0) return; 



	temp=new float[ndist];
	order=new int[ndist];
	
	
	float maxw=-10000;
	for(i=0;i<ndist;i++) {
		if(weight[i]>maxw) maxw=weight[i];
	}
	

	for(i=0;i<ndist;i++) {
		temp[i]=0;
		order[i]=0;
		i1=atmpair[i*2];
		i2=atmpair[i*2+1];
		//temp[i]=i1*natom+i2+i*1.0/ndist;
		temp[i]=i1*natom+i2+weight[i]/maxw/2;
	}
	
	Qsort cc;

	cc.sort(temp,ndist,order);
	
	float *dis=new float[3*ndist];
	int   *atmp=new int[2*ndist];
	//int   *predt0=new int[2*ndist];
	float *weight0=new float[ndist];

	for(i=0;i<ndist;i++) {
		
		weight0[i]=weight[i];
		for(j=0;j<3;j++) {
			dis[i*3+j]=dist[i*3+j];
			dist[i*3+j]=0;
		}

		for(j=0;j<2;j++) {
			atmp[i*2+j]=atmpair[i*2+j];
			atmpair[i*2+j]=0;
			//predt0[i*2+j]=predt[i*2+j];
		}
		//for(j=0;j<1;j++) predt0[i*2+j]=predt[i*2+j];
	}

	for(i=0;i<ndist;i++) {
		j=order[i];
			
		for(k=0;k<3;k++) {
                        dist[i*3+k]=dis[j*3+k];
                }

		for(k=0;k<2;k++) {
                        atmpair[i*2+k]=atmp[j*2+k];
			//predt[i*2+k]=predt0[j*2+k];
                }
		weight[i]=weight0[j];
	}

	delete [] temp;
	delete [] order;
	delete [] dis;
	delete [] atmp;
	delete [] weight0;
	//delete [] predt0;
}

void DistBound::writeoutdistbound(char *s) {

	FILE *fp=fopen(s,"w");

	if(fp==0) {
		cerr<<"could not write file:"<<s<<endl;
		exit(0);
	}

	for(int i=0;i<ndist;i++) {

		fprintf(fp,"%5i %5i %8.3f %8.3f %8.3f\n",atmpair[2*i]+1,atmpair[2*i+1]+1,dist[3*i],dist[3*i+1],dist[3*i+2]);
	} 

	fclose(fp);

}

void DistBound::writeoutdistbound(char *s,int u0,int u1) {

        FILE *fp=fopen(s,"w");

        if(fp==0) {
                cerr<<"could not write file:"<<s<<endl;
                exit(0);
        }

        for(int i=u0;i<min(u1,ndist);i++) {

                fprintf(fp,"%5i %5i %8.3f %8.3f %8.3f\n",atmpair[2*i]+1,atmpair[2*i+1]+1,dist[3*i],dist[3*i+1],dist[3*i+2]);
        }

        fclose(fp);

}


void DistBound::clearredundancy() {

        int i,n;


        cerr<<"clear distance redundancy..."<<endl;

	if(ndist==0) return;

	int *delt=new int[ndist];
	for(i=0;i<ndist;i++) delt[i]=1;

        for(i=0;i<ndist-1;i++) {
	    n=2*i;
            if(atmpair[n]==atmpair[n+2]&&atmpair[n+1]==atmpair[n+3]) {
		  
	          delt[i+1]=0;
            }
            else if(predt[n]==0&&predt[n+1]==0) {
	          delt[i]=0;
            }
        }

	int k;
        int j=0;
        for(i=0;i<ndist;i++) {
                if(delt[i]==0) continue;
                for(k=0;k<3;k++) {
                        dist[j*3+k]=dist[i*3+k];
                }
                for(k=0;k<2;k++) {
                        atmpair[j*2+k]=atmpair[i*2+k];
                }
                j++;
        }

	int m=ndist-j;

        ndist=j;

	delete [] delt;

        cerr<<"the total number of redundancy distance constraints is:"<<m<<endl;
}


void DistBound::setimproper() {
	
	Chn *c;
        Res *r;

	cerr<<"find out the improper dihedral angles..."<<endl;

	nimp=0;
	int n=model->manyatm();
	improper=new Atm*[4*n];
	for(c=model->chn;c;c=c->next) 
	for(r=c->res;r;r=r->next) {

		//make sure we have l- amino acids
		if(r->name!='G') {
			improper[4*nimp]=r->isatm(" CA ");
			improper[4*nimp+1]=r->isatm(" N  ");	
			improper[4*nimp+2]=r->isatm(" C  ");
			improper[4*nimp+3]=r->isatm(" CB ");
			if(improper[4*nimp]==0||improper[4*nimp+1]==0||
			   improper[4*nimp+2]==0||improper[4*nimp+3]==0)
			{
				for(int i=0;i<4;i++) improper[4*nimp+i]=0;
			}
			else {
				nimp++;
			}
		}			
	
		//sidechains

		if(r->name=='L') {

			improper[4*nimp]=r->isatm(" CG ");      
                        improper[4*nimp+1]=r->isatm(" CD2");    
                        improper[4*nimp+2]=r->isatm(" CD1");
                        improper[4*nimp+3]=r->isatm(" CB ");
                        if(improper[4*nimp]==0||improper[4*nimp+1]==0||
                           improper[4*nimp+2]==0||improper[4*nimp+3]==0)
                        {
                                for(int i=0;i<4;i++) improper[4*nimp+i]=0;
                        }
                        else {
                                nimp++;
                        }
		}
		
		if(r->name=='V') {


			improper[4*nimp]=r->isatm(" CB ");      
                        improper[4*nimp+1]=r->isatm(" CG2");    
                        improper[4*nimp+2]=r->isatm(" CG1");
                        improper[4*nimp+3]=r->isatm(" CA ");
                        if(improper[4*nimp]==0||improper[4*nimp+1]==0||
                           improper[4*nimp+2]==0||improper[4*nimp+3]==0)
                        {
                                for(int i=0;i<4;i++) improper[4*nimp+i]=0;
                        }
                        else {
                                nimp++;
                        }

		}

		if(r->name=='I') {


                        improper[4*nimp]=r->isatm(" CB ");      
                        improper[4*nimp+1]=r->isatm(" CG1");    
                        improper[4*nimp+2]=r->isatm(" CG2");
                        improper[4*nimp+3]=r->isatm(" CA ");
                        if(improper[4*nimp]==0||improper[4*nimp+1]==0||
                           improper[4*nimp+2]==0||improper[4*nimp+3]==0)
                        {
                                for(int i=0;i<4;i++) improper[4*nimp+i]=0;
                        }
                        else {
                                nimp++;
                        }

                }

		if(r->name=='T') {


                        improper[4*nimp]=r->isatm(" CB ");      
                        improper[4*nimp+1]=r->isatm(" OG1");    
                        improper[4*nimp+2]=r->isatm(" CG2");
                        improper[4*nimp+3]=r->isatm(" CA ");
                        if(improper[4*nimp]==0||improper[4*nimp+1]==0||
                           improper[4*nimp+2]==0||improper[4*nimp+3]==0)
                        {
                                for(int i=0;i<4;i++) improper[4*nimp+i]=0;
                        }
                        else {
                                nimp++;
                        }
                }
	}

	cerr<<"the total number of improper dihedral is:"<<nimp<<endl;
}

void DistBound::printoutimproper() {

	for(int i=0;i<nimp;i++) {
	     float d=getta(improper[i*4],improper[i*4+1],improper[i*4+2],improper[i*4+3]);
             cerr<<i<<": ";
	     cerr<<improper[i*4]->name<<improper[i*4]->id0<<" ";
	     cerr<<improper[i*4+1]->name<<improper[i*4+1]->id0<<" ";
	     cerr<<improper[i*4+2]->name<<improper[i*4+2]->id0<<" "; 
	     cerr<<improper[i*4+3]->name<<improper[i*4+3]->id0<<" ";
	     cerr<<" "<<d<<endl;
        }
}

int DistBound::optimize() {
      return optimize(500);
}

int DistBound::optimize(int maxnit){

// makes distances correct
 
      float dd[3],coo[3],d,len;
      int nit,k,i,j,i1,i2,nviol,nchk;
      int low;//,oldlow;
      int   *viol,*ip;
      float *imp;

      cerr<<"correcting the distance constraint..."<<endl;
 
      viol=new int[natom];
      ip=new int[ndist];
      imp=new float[nimp];

      //start of iteration

      for(i=0;i<ndist;i++) ip[i]=i;

      low=ndist;

      //oldlow=ndist;

      nchk=ndist;	

      for(nit=0;nit<maxnit;nit++) {
	
         for (i=0;i<natom;i++) viol[i]=0;

         nviol=0;

         //randomise list

         for(i=0;i<nchk;i++) {
	     j=random();
	     j=j%ndist;
	     k=ip[i];
	     ip[i]=ip[j];
	     ip[j]=k; 	
         }	

         // now walk through the list

         for(i=0;i<nchk;i++) {
             k=ip[i];
             i1=atmpair[2*k];
	     i2=atmpair[2*k+1];
             for(j=0;j<3;j++) dd[j]=atoms[i1]->xyz[j]-atoms[i2]->xyz[j];

             len=dd[0]*dd[0]+dd[1]*dd[1]+dd[2]*dd[2];

             //if outside distance bound: make correction

             if (len>dist[k*3+2]||len<dist[k*3+1]) {

		if(len<1) len=1;
	        d=dist[k*3]/len;
		d=sqrt(d)/2;
                nviol=nviol+1;
                if(predt[i1]) viol[i1]=1;
                if(predt[i2]) viol[i2]=1;

                //both to be predicted
	        if(predt[i1]&&predt[i2]) {	

		   for(j=0;j<3;j++)  coo[j]=(atoms[i1]->xyz[j]+atoms[i2]->xyz[j])/2;
		   for(j=0;j<3;j++) { 
			atoms[i1]->xyz[j]=coo[j]+dd[j]*d;
		        atoms[i2]->xyz[j]=coo[j]-dd[j]*d;
                   }	
		   //cerr<<len<<" "<<TRES.distance(atoms[i1],atoms[i2])<<endl;
	        }
                else if(predt[i1]) {  //atom i to be predicted
		    for(j=0;j<3;j++) { 
                        atoms[i1]->xyz[j]=atoms[i2]->xyz[j]+dd[j]*d*2;
                    }

                }
                else if(predt[i2]) { //atom j to be predicted
		    for(j=0;j<3;j++) { 
			atoms[i2]->xyz[j]=atoms[i1]->xyz[j]-dd[j]*d*2;
                    }
                }
            }
         }

         //update list of distances to be checked in next round

         nchk=0;
         for(i=0;i<natom;i++) {
	    if(predt[i]==0) continue;
            if(tag[2*i]==-1) continue;
	    for(j=tag[2*i];j<=tag[2*i+1];j++){
                 i1=atmpair[2*j];
                 i2=atmpair[2*j+1];
                 if(viol[i1]||viol[i2]) {
                      ip[nchk++]=j;
                 }
            }
         }

         low=min(low,nviol);
	 /*
	 if (nit>10&&nviol>ndist/5) goto re200; 
         if (nit>100&&nviol>ndist/10) goto re200;
         if((nit+1)%50==0) {
            if (low>9*oldlow/10&&low>100) goto re200;
            oldlow=low;
         }
	 */

	 //if there are still violations, do another roun

	 if (nchk>0) {
           if (!((nit+1)%10==0&&chiropt==2)) continue;
         }

         //check impropers to get chirality right

         for(i=0;i<nimp;i++) {
	     imp[i]=getta(improper[i*4],improper[i*4+1],improper[i*4+2],improper[i*4+3]);
         }

         j=0;
         for (i=0;i<nimp;i++) {
            if (imp[i]<0.0||imp[i]>70.0) j=j+1;
         }

         //if many impropers are incorrect : try mirror image

         if (j>nimp/2) {
              for(i=0;i<natom;i++) {
	         if(predt[i]==0) continue;
	         atoms[i]->xyz[1]=-atoms[i]->xyz[1];
              }
	      for(i=0;i<nimp;i++) {
                 imp[i]=getta(improper[i*4],improper[i*4+1],improper[i*4+2],improper[i*4+3]);
              }
         }

         //now make the corrections, if necessary

	 if(chiropt==0) goto re200;
         for (i=0;i<nimp;i++) {
             if (imp[i]<0.0||imp[i]>70.0) {
                 swapatm(i,viol);
	     }
         }
	 if(chiropt==1) goto re200;

         //update list of distances to be checked and do another round if
         //required

         if (nchk>0) continue;
     }

     re200:
	
     if(nviol!=0&&chiropt==2) cerr<<"the structure rejected after: "<<nit<<" iterations"<<endl;

     if(viol) delete [] viol;
     if(ip)   delete [] ip;
     if(imp)  delete [] imp;
     
     return nviol;
}

float DistBound::getta(Atm *at1,Atm *at2,Atm *at3,Atm *at4) {

      //calculates torsion angle at1 at2 at3 at4

      float ux,uy,uz,vx,vy,vz,wx,wy,wz,uvx,uvy,uvz,vwx,vwy,vwz,uvl,vwl;
      float cosphi,temp,si,uvlvwl;

      ux = at2->xyz[0]-at1->xyz[0]; 
      uy = at2->xyz[1]-at1->xyz[1];
      uz = at2->xyz[2]-at1->xyz[2];
      vx = at3->xyz[0]-at2->xyz[0];
      vy = at3->xyz[1]-at2->xyz[1];
      vz = at3->xyz[2]-at2->xyz[2];
      wx = at4->xyz[0]-at3->xyz[0];
      wy = at4->xyz[1]-at3->xyz[1];
      wz = at4->xyz[2]-at3->xyz[2];
      
      uvx = uy*vz - uz*vy;
      uvy = uz*vx - ux*vz;
      uvz = ux*vy - uy*vx;

      vwx = vy*wz - vz*wy;
      vwy = vz*wx - vx*wz;
      vwz = vx*wy - vy*wx;
      
      uvl = sqrt(uvx*uvx + uvy*uvy + uvz*uvz);
      vwl = sqrt(vwx*vwx + vwy*vwy + vwz*vwz);
      
      uvlvwl=uvl*vwl;
      if (uvlvwl==0.000000) uvlvwl=1e-12;
      cosphi = (uvx*vwx + uvy*vwy + uvz*vwz) / uvlvwl;
      if (cosphi>1.000000) cosphi=1.000000;
      if (cosphi<-1.000000) cosphi=-1.000000;
      temp=uvx*wx+uvy*wy+uvz*wz;
      if (fabs(temp)==0.000000) {
         si=1.0;
      } else {
         si=temp/fabs(temp);
      }
      float getta = si*(360.0 / (2.0*3.1415927))*acos(cosphi);

      return getta;
}

void DistBound::boundsquare(float n) {

	cerr<<"calculate distance square of distance constraint..."<<endl;

	for(int i=0;i<ndist;i++) {

		dist[i*3]=pow(dist[i*3],n);
		dist[i*3+1]=pow(dist[i*3+1],n);
		dist[i*3+2]=pow(dist[i*3+2],n);
	}
}

void DistBound::searchstructure(int n) {

	int i=0;
	while(1) {
		i++;
		if(i>n) break;
		initialstructure();
		int n=optimize();
		float d=evaluate();
		cerr<<i<<" : "<<n<<" "<<d<<endl;
		//if(n==0) {
			char line[1000];
			sprintf(line,"s%i.pdb",i);
			model->write(line);
			cerr<<line<<endl;
		//}
	} 

}
void DistBound::initialstructure() {

	for(int i=0;i<natom;i++) {	

		if(predt[i]==0) continue;
		
		for(int k=0;k<3;k++) {		
			atoms[i]->xyz[k]=(random())%1000-1.0;
		}
	}
	model->center();
}

void DistBound::mirror(Atm *a1,Atm *a2,Atm *a3,Atm *a4) {

//adapt the coordinates of atom n4 such that it is mirrorred
//with respect to the plane defined by n1,n2,n3 (n1 origin)

      float v1[3],v2[3],v3[3],v4[3],rk,rl;
      float l1,l2,i12,o[3],p[3],pp[3],p3,p4;

      int i;

      for(i=0;i<3;i++) {

		o[i]=a1->xyz[i];
		p[i]=a4->xyz[i];
      }
      	
      //v1 and v2 span the plane

      for(i=0;i<3;i++) {
		v1[i]=a2->xyz[i]-a1->xyz[i];
		v2[i]=a3->xyz[i]-a1->xyz[i];
      }

      l1=TRES.distance(a1,a2);
      l2=TRES.distance(a1,a3);

      i12=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

//we need orthonormal basis set
//l3 is normalized l1
//l4 is orthogonal to it and lies in the plane
//(linear combination of l1 and l2, with a1 and a2 as coefficients)

      rk=0.0-i12/(l1*l1);
      rl=(rk*rk)*(l1*l1)+(l2*l2)+2*rk*i12;

      float aa2=sqrt(rl)/rl;
      float aa1=rk*aa2;

      for(i=0;i<3;i++) {
	v3[i]=v1[i]/l1;
        v4[i]=aa1*v1[i]+aa2*v2[i];
      }

//adapt coordinates

      p3=(p[0]-o[0])*v3[0]+(p[1]-o[1])*v3[1]+(p[2]-o[2])*v3[2];
      p4=(p[0]-o[0])*v4[0]+(p[1]-o[1])*v4[1]+(p[2]-o[2])*v4[2];
      for(i=0;i<3;i++) {
      	pp[i]=o[i]+p3*v3[i]+p4*v4[i];
	a4->xyz[i]=2*pp[i]-p[i];
      }

}

void DistBound::swapatm(Atm *a,Atm *b) {

      //swaps atom a and atom b

      float dum[3];

      for(int i=0;i<3;i++) {
         dum[i]=a->xyz[i];
         a->xyz[i]=b->xyz[i];
         b->xyz[i]=dum[i];
      }
}

void DistBound::swapatm(int n,int *viol) {

      //generic swap routine for atoms involved in incorrect
      //improper dihedrals

      Atm *a1,*a2,*a3,*a4;

      //conversion of d to l amino acid

      if(improper[4*n+2]->tatm->id==2) {
         a1=improper[4*n];
	 a2=improper[4*n+1];
         a3=improper[4*n+2];
         for(a4=a1->res->atm;a4;a4=a4->next) {
            viol[a4->id0]=1;
            if(a4!=a1&&a4!=a2&&a4!=a3) {
		mirror(a1,a2,a3,a4);
            }
         }
      }
      else if (improper[4*n+2]->res->name=='L') {   //swap of CD1 and CD2 in LEU or CG1 and CG2 in VAL
         a1=improper[4*n+2];
	 a2=improper[4*n+1];
         swapatm(a1,a2);
      }
      else if (improper[4*n+2]->res->name=='V') {
	 a1=improper[4*n+2];
         a2=improper[4*n+1];
         swapatm(a1,a2);
      }
      else if (improper[4*n+2]->res->name=='I') {  //mirror ILE or THR side chains
	 a1=improper[4*n+2];
         a2=improper[4*n];
         a3=improper[4*n+3];
         a4=a1->res->isatm(" CD ");
	 if(a4==0) a4=a1->res->isatm(" CD1");
         viol[a4->id0]=1;
         mirror(a1,a2,a3,a4);
         
         a4=improper[4*n+1];
         viol[a4->id0]=1;
         mirror(a1,a2,a3,a4);
      }
      else if (improper[4*n+2]->res->name=='T') {
         a1=improper[4*n+1];
         a2=improper[4*n];
         a3=improper[4*n+3];
         a4=improper[4*n+2];
         mirror(a1,a2,a3,a4);
         viol[a4->id0]=1;
         a4=a1->res->isatm(" HG1");
         if (a4!=0) {
            viol[a4->id0]=1;
            mirror(a1,a2,a3,a4);
         }
      }
      else{
      }
}


void DistBound::adjust() {

      int i1,i2;

      cerr<<"adjusting the distance constraint..."<<endl;

      for(int i=0;i<ndist;i++) {

	i1= atmpair[2*i];
	i2= atmpair[2*i+1];

	if(predt[i1]==0&&predt[i2]==1) {
		
		atmpair[2*i]=i2;
		atmpair[2*i+1]=i1;
	}
      }

}

void DistBound::defaultaction() {

	adjust();
	reorder();
	clearredundancy();
	setimproper();
	settag();
	boundsmooth(5);
	boundsquare(2);
}

int DistBound::findpair(HashTable **hash,int k,int l) {

	long a,i;

	i=2*ndist;
	a=k*l+k+l;
	a=a%i;

	HashTable *c;
	int a1,a2;

	for(c=hash[a];c;c=c->next) {

		if(c->key==-1) return -1;

		a1=atmpair[2*c->key];
		a2=atmpair[2*c->key+1];

		if(a1==k&&a2==l) return c->key;
	}

	return -1;
}

void DistBound::boundsmooth(int tot) {

	int i,j;
	int a,b,c,k;
	int n1,n2,n3;
	//float dum;

	int num=0;

	cerr<<"bound smoothing..."<<endl;
	if(ndist==0) return;
	HashTable **hash=new HashTable*[ndist*2];
	for(i=0;i<ndist*2;i++) hash[i]=new HashTable();	

	long a1,a2,a3;
	a=2*ndist;
	for(i=0;i<ndist;i++) {
		a1=atmpair[2*i];	
		a2=atmpair[2*i+1];
		a3=a1*a2+a1+a2;	
		a3=a3%a;
		hash[a3]->add(i);
	}

	while(1) {
	num++;
	if(num>tot) break;
	    for(i=0;i<ndist;i++) {

         	a=atmpair[2*i];
         	b=atmpair[2*i+1];

         	for(j=tag[a*2];j<i;j++) {
			if(j==-1) break;

            		c=atmpair[2*j+1];
			if(b<c) k=findpair(hash,b,c);
			else    k=findpair(hash,c,b);

			if(k==-1) continue;
			
			n1=i;n2=j;n3=k;

			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
               		dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
			dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);

			n1=j;n2=i;n3=k;
			
			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);

			n1=k;n2=j;n3=k;

			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);
        	} 
	    }
	}

	for(i=0;i<2*ndist;i++) delete hash[i];

	delete [] hash;
}

void DistBound::createnewbound() {

	if(ndist==0) return;

	int i;
	HashTable **hash=new HashTable*[ndist*2];
        for(i=0;i<ndist*2;i++) hash[i]=new HashTable();


	long a1,a2,a3;
        int a=2*ndist;
        for(i=0;i<ndist;i++) {
                a1=atmpair[2*i];
                a2=atmpair[2*i+1];
                a3=a1*a2+a1+a2;
                a3=a3%a;
                hash[a3]->add(i);
        }
}
