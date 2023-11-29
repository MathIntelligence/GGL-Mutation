#include "source.h"

Bound::Bound() {
	npower=1;
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
	modtop=0;
	dssp=0;
	range=0;
	maxnit=500;
	bond=0;
	curve=0;
	atmat=0;
}

Bound::~Bound() {
	if(atmat) delete [] atmat;
	if(bond) delete [] bond;
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
	if(dssp) delete [] dssp;dssp=0;
	if(range) delete [] range;range=0;
	if(curve) delete [] curve;curve=0;
}

void Bound::settag() {


	int i,i1;

	if(TRES.logg>3) cerr<<"set tags ... "<<endl;

	if(natom==0) return;
	if(tag) delete [] tag;tag=0;

	tag=new int[2*natom];

	for(i=0;i<2*natom;i++) tag[i]=-1;

	for(i=0;i<ndist;i++) {
		i1=atmpair[2*i];
		if(tag[2*i1]==-1) {tag[2*i1]=i;tag[2*i1+1]=i;}
		else		  tag[2*i1+1]=i;
	}
}

void Bound::setbuffertag() {

	if(ndist==0) return;

	int *atmpair0;
	float *dist0,*weight0;
	

	dist0=new float[6*ndist];
	atmpair0=new int[4*ndist];
	weight0=new float[2*ndist];

	int i,j;

	for(i=0;i<ndist*2;i++) {
		for(j=0;j<3;j++) dist0[3*i+j]=0;
		for(j=0;j<2;j++) atmpair0[2*i+j]=0;
		weight0[i]=0;
	}

	int k=0;

	for(i=0;i<ndist;i++) {	
		int i1=atmpair[i*2];
		int i2=atmpair[i*2+1];
		if(predt[i1]==1&&predt[i2]==1) {
			atmpair0[k*2]=i1;
			atmpair0[k*2+1]=i2;
			dist0[3*k]=dist[3*i];
			dist0[3*k+1]=dist[3*i+1];
			dist0[3*k+2]=dist[3*i+2];
			weight0[k]=weight[i];
			k++;
			atmpair0[k*2]=i2;
                        atmpair0[k*2+1]=i1;
                        dist0[3*k]=dist[3*i];
                        dist0[3*k+1]=dist[3*i+1];
                        dist0[3*k+2]=dist[3*i+2];
			weight0[k]=weight[i];
			k++;
		}
		else if(predt[i1]==1&&predt[i2]==0) {
			atmpair0[k*2]=i1;
                        atmpair0[k*2+1]=i2;
                        dist0[3*k]=dist[3*i];
                        dist0[3*k+1]=dist[3*i+1];
                        dist0[3*k+2]=dist[3*i+2];
                        weight0[k]=weight[i];
                        k++;
		}
		else if(predt[i1]==0&&predt[i2]==1) {
			atmpair0[k*2]=i2;
                        atmpair0[k*2+1]=i1;
                        dist0[3*k]=dist[3*i];
                        dist0[3*k+1]=dist[3*i+1];
                        dist0[3*k+2]=dist[3*i+2];
                        weight0[k]=weight[i];
                        k++;
                }
	}

	arraysize=2*ndist;
	ndist=k;

	if(dist) delete [] dist;
	if(atmpair)delete [] atmpair;
	if(weight) delete [] weight; 
		
	dist=dist0;
	atmpair=atmpair0;
	weight=weight0;
	
	reorder();
	settag();
	if(TRES.logg>3) writeoutdistbound("out1");
}

void Bound::setmodel(Pdb *s) {

	model=s;
}


void Bound::ready() {

	Chn *c;
        Res *r;
        Atm *a;

	srandom(seed);

	if(model==0) {
		cerr<<"model is null"<<endl;
		exit(0);
	}
	
	model->setallnear();

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

	arraysize=1000;
	dist=new float[3*arraysize];
	atmpair=new int[2*arraysize];
        weight=new float[arraysize];
	ndist=0;
}
void Bound::easyready() {

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
		//a->allnear(3,0);
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

	arraysize=1000;
	dist=new float[3*arraysize];
	atmpair=new int[2*arraysize];
        weight=new float[arraysize];
	ndist=0;
}
void Bound::setsize(int n) {

	arraysize=n;
	dist=(float *)realloc(dist,3*sizeof(float)*n);
	atmpair=(int *)realloc(atmpair,2*sizeof(int)*n);
	weight=(float *)realloc(weight,sizeof(float)*n);
	//predt=(int *)realloc(predt,2*sizeof(int)*n);
}


void Bound::setstartsize(int n) {

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

int Bound::readdist(char *filename) {

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

void Bound::addbounds(Atm *a1,Atm *a2,float d,float w) {

	float d0=TRES.distance(a1,a2);

	float d1=d0-d;
	float d2=d0;
	float d3=d0+d;

	addbounds(a1->id0,a2->id0,d2,d1,d3,w);
	//if(a1->id0<a2->id0) addbounds(a1->id0,a2->id0,d2,d1,d3,w);
	//else addbounds(a2->id0,a1->id0,d2,d1,d3,w);
}

void Bound::addbounds(int i1,int i2,float d1,float d2,float d3,float w) {

	 int nn=ndist;

	 /*
	 if(nn+10>arraysize) {
                setsize(10000+arraysize);
         }
	 */

	 if(i1>i2) {
		int t=i1;
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
	 //cerr<<ndist<<" "<<d1<<" "<<d2<<" "<<d3<<endl;
}

void Bound::getring(Atm **ring,Res *r) {

	int i=0;
	Atm *a;
	for(a=r->atm;a;a=a->next) {
		if(a->tatm->ring) ring[i++]=a;
	}

	if(strchr("YFWH",r->name)) {
		a=r->isatm(" CB ");
		if(a) ring[i++]=a;
	}

	if(strchr("Y",r->name)) {
		a=r->isatm(" OH ");
		if(a) ring[i++]=a;
	}
	
	ring[i]=0;
}

int Bound::setquasiringdist() {

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


int Bound::quasiringdist(Res *r) {

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


int Bound::isnextres(Res *r1,Res *r2){

	if(r2->id-r1->id!=1) return 0;

	Atm *n,*c;

	n=r2->isatmid(0);

	c=r1->isatmid(2);

	if(n==0||c==0) return 0;
	float d=TRES.distance(n->xyz,c->xyz);

	if(d>2.0) return 0;

	return 1;

}


int Bound::setphipsidist() {

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


int Bound::set14dist() {

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

int Bound::set14dist(Res *r) {

	Atm *a,*a1;

	int norg=ndist;
	for(a=r->atm;a;a=a->next)
	for(a1=a->next;a1;a1=a1->next) {

		if(a->isnear(a1,4)!=4) continue;
		addbounds(a,a1,0.5,1);
	}
	return norg;
}



int Bound::setn14dist() {

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

int Bound::setn14dist(Res *r0,Res *r) {

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


int Bound::setsidechainrestriction() {

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

int Bound::setsidechainrestriction(Res *r) {

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

int Bound::setexternalbounds() {

	int norg=ndist;

	StrFmt *mf=modtop->strfmt->getrootStrFmt();
	DistPopular *dst=TRES.popbin->getfaraway();

	Chn *c;
	Res *r,*r0;
	Atm *a,*a0;

	for(c=model->chn;c;c=c->next)
	for(r=c->res;r;r=r->next)
	for(r0=r->next;r0;r0=r0->next) {
		if(abs(r0->id-r->id)<2) continue;
		for(a=r->atm;a;a=a->next)
		for(a0=r0->atm;a0;a0=a0->next) {
		 	DistPopular *p=dst->getDistPopular(r->tres,a->name,a0->name,-1);
			float x=p->findmostpopulardistance();
			addbounds(a->id0,a0->id0,x,p->dist[0],p->dist[1],1);
		}
	}


	for(c=model->chn;c;c=c->next)
	for(r=c->res;r;r=r->next)
	for(r0=r->next;r0;r0=r0->next) {
		if(abs(r0->id-r->id)<2) continue;
		for(a=r->atm;a;a=a->next)
		for(a0=r0->atm;a0;a0=a0->next) {

			float *dp=mf->getfardistbounds(a,a0);
			int ii=0;
                        while(dp&&dp[ii]>0) {
			      addbounds(a->id0,a0->id0,dp[4*ii],dp[4*ii+1],dp[4*ii+2],dp[4*ii+3]);
                              //if(a->id0<a0->id0)  distsecond->addbounds(a->id0,a0->id0,dp[ii],1);
                              //else distsecond->addbounds(a0->id0,a->id0,dp[ii],1);
			      ii++;
                        }
                        if(dp) delete [] dp;
		}
	}

	return norg;
}

int Bound::setnearnextbounds() {


	if(model==0) return ndist;

	int norg=ndist;

	//find internal dist popular;

	DistPopular *dst=TRES.popbin->getnearnext();

	Chn *c;
	Res *r,*r0,*r1;
	Atm *a,*a0;

	for(c=model->chn;c;c=c->next)
	for(r=c->res;r;r=r->next) {

		r0=c->isres(r->id-1);

		r1=c->isres(r->id+1);

		//handle r0
		if(r0) {
			for(a=r->atm;a;a=a->next)
			for(a0=r0->atm;a0;a0=a0->next) {
				DistPopular *p=dst->getDistPopular(r->tres,a->name,a0->name,-1);
				float x=p->findmostpopulardistance();
				addbounds(a->id0,a0->id0,x,p->dist[0],p->dist[1],1);
				float *dp=p->findmostpopulardp(1.0);
				for(int ii=0;dp&&ii<10;ii++) {
					if(dp[ii]<0) continue;
					if(a->id0<a0->id0)  distsecond->addbounds(a->id0,a0->id0,dp[ii],1+ii/10.);
					else distsecond->addbounds(a0->id0,a->id0,dp[ii],1+ii/10.);
				}
				if(dp) delete [] dp;
			}
		}

		//handle r1
		if(r1) {
			for(a=r->atm;a;a=a->next)
                        for(a0=r0->atm;a0;a0=a0->next) {
                                DistPopular *p=dst->getDistPopular(r->tres,a->name,a0->name,1);
                                float x=p->findmostpopulardistance();
                                addbounds(a->id0,a0->id0,x,p->dist[0],p->dist[1],1);
                                float *dp=p->findmostpopulardp(1.0);
                                for(int ii=0;dp&&ii<10;ii++) {
                                        if(dp[ii]<0) continue;
                                        if(a->id0<a0->id0)  distsecond->addbounds(a->id0,a0->id0,dp[ii],1+ii/10.);
                                        else distsecond->addbounds(a0->id0,a->id0,dp[ii],1+ii/10.);
                                }
                                if(dp) delete [] dp;
                        }
		}

		if(r) {

			for(a=r->atm;a;a=a->next)
                        for(a0=r->atm;a0;a0=a0->next) {
				if(a==a0) continue;
                                DistPopular *p=dst->getDistPopular(r->tres,a->name,a0->name,1);
                                float x=p->findmostpopulardistance();
                                addbounds(a->id0,a0->id0,x,p->dist[0],p->dist[1],1);
                                float *dp=p->findmostpopulardp(1.0);
                                for(int ii=0;dp&&ii<10;ii++) {
                                        if(dp[ii]<0) continue;
                                        if(a->id0<a0->id0)  distsecond->addbounds(a->id0,a0->id0,dp[ii],1+ii/10.);
                                        else distsecond->addbounds(a0->id0,a->id0,dp[ii],1+ii/10.);
                                }
                                if(dp) delete [] dp;
                        }
		}
	}


	//set internal cool

	StrFmt *mf=modtop->strfmt->getrootStrFmt();

	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {

                r0=c->isres(r->id-1);

                r1=c->isres(r->id+1);

		//handle r0
                if(r0) {
                        for(a=r->atm;a;a=a->next)
                        for(a0=r0->atm;a0;a0=a0->next) {

				float *dp=mf->getdistbounds(a,a0);
				int ii=0;
                                while(dp&&dp[4*ii]>0) {
					addbounds(a->id0,a0->id0,dp[4*ii],dp[4*ii+1],dp[4*ii+2],dp[4*ii+3]);
                                        //if(a->id0<a0->id0)  distsecond->addbounds(a->id0,a0->id0,dp[ii],1);
                                        //else distsecond->addbounds(a0->id0,a->id0,dp[ii],1);
					ii++;
                                }
                                if(dp) delete [] dp;
                        }
                }

		//handle r1
                if(r1) {
                        for(a=r->atm;a;a=a->next)
                        for(a0=r0->atm;a0;a0=a0->next) {
				float *dp=mf->getdistbounds(a,a0);
                                int ii=0;
                                while(dp&&dp[ii]>0) {
					addbounds(a->id0,a0->id0,dp[4*ii],dp[4*ii+1],dp[4*ii+2],dp[4*ii+3]);
                                        //if(a->id0<a0->id0)  distsecond->addbounds(a->id0,a0->id0,dp[ii],1);
                                        //else distsecond->addbounds(a0->id0,a->id0,dp[ii],1);
                                        ii++;
                                }
                                if(dp) delete [] dp;
                        }
                }


		if(r) {

                        for(a=r->atm;a;a=a->next)
                        for(a0=r->atm;a0;a0=a0->next) {
                                if(a==a0) continue;
				float *dp=mf->getdistbounds(a,a0);
                                int ii=0;
                                while(dp&&dp[ii]>0) {
					addbounds(a->id0,a0->id0,dp[4*ii],dp[4*ii+1],dp[4*ii+2],dp[4*ii+3]);
                                        //if(a->id0<a0->id0)  distsecond->addbounds(a->id0,a0->id0,dp[ii],1);
                                        //else distsecond->addbounds(a0->id0,a->id0,dp[ii],1);
                                        ii++;
                                }
                                if(dp) delete [] dp;
                        }
                }

	}

	return norg;
}

int Bound::setinternaldist() {

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

int Bound::setalldist() {

	int norg=ndist;
	setcloseresidue();  //set constraints in residues
	//setallbackbone();
	//setsecondarybounds();
	//sethbondbounds();
	//setclosedistance(); //set constraint close in distance
	return norg;
}


void Bound::setallbackbone() {

	/*
	//StrFmt *c;

        //StrFmt *sq=fmt->findsequencefmt();
        //if(sq==0||sq->seqngap==0) return;

        //int n=strlen(sq->seqngap);
	if(model==0) return ndist;

	int norg=ndist;

        StrFmt *mf=modtop->strfmt->getrootStrFmt();

        int i,j;
        int n1,n2;
        int m1,m2;
        Atm *a1,*a2;
        Atm *b1,*b2;
	Res *r1,*r2;
        //for(i=0;i<n;i++) {
	for(r1=model->chn->res;r1;r1=r1->next) {
                //for(j=i+1;j<n;j++) {
		for(r2=r1->next;r2;r2=r2->next) {
			if(abs(r2->id-r1->id)<=1) continue;
                        //if(sq->seqngap[i]=='-'||sq->seqngap[j]=='-') continue;
			for(a1=r1;a1;
                        m1=sq->match[i];
                        m2=sq->match[j];
                        //b1=sq->resn[m1]->atm->next;
                        //b2=sq->resn[m2]->atm->next;
                        for(c=fmt;c;c=c->more) {
                                if(c->seqngap[i]=='-'||c->seqngap[j]=='-') continue;
                                if(c->token&&strcmp(c->token,"sequence")==0) continue;
                                n1=c->match[i];
                                n2=c->match[j];
                                for(int ii=0;ii<=4;ii++)
                                for(int jj=0;jj<=4;jj++) {
                                        a1=c->resn[n1]->isatmid(ii);
                                        a2=c->resn[n2]->isatmid(jj);
                                        b1=sq->resn[m1]->isatmid(ii);
                                        b2=sq->resn[m2]->isatmid(jj);
                                        if(a1==0||a2==0) continue;
                                        float d=TRES.distance(a1,a2);
                                        //cerr<<a1->name<<" "<<a2->name<<" "<<c->code<<d<<endl;
                                        if(b1->id0<b2->id0) distbound->addbounds(b1->id0,b2->id0,d,d-5,d+5,2.0);
                                        else                distbound->addbounds(b2->id0,b1->id0,d,d-5,d+5,2.0);
                                }
                        }
                }
        }

	*/
}

int Bound::sethbondbounds() {

	if(model==0) return ndist;

	int norg=ndist;

        StrFmt *mf=modtop->strfmt->getrootStrFmt();

	Chn *c;
    	Res *r,*r0;

	for(c=model->chn;c;c=c->next)
    	for(r=c->res;r;r=r->next) 
	for(r0=r->next;r0;r0=r0->next){
		
		if(r0->id-r->id<=1) continue;

		//if(mf->ishbondexist(r,r0)==0) continue;

		float *dp=mf->gethbondbounds(r,r0);
                int ii=0;
                while(dp&&dp[6*ii]>0) {
                      addbounds((int)dp[6*ii],(int)dp[6*ii+1],dp[6*ii+2],dp[6*ii+3],dp[6*ii+4],dp[6*ii+5]);
                      ii++;
                }
                if(dp) delete [] dp;
	}
	return norg;
}

int Bound::setsecondarybounds() {

        if(model==0) return ndist;

        int norg=ndist;

        StrFmt *mf=modtop->strfmt->getrootStrFmt();

        Chn *c;
        Res *r,*r0;

        for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next)
        for(r0=r->next;r0;r0=r0->next){

                if(r0->id-r->id<=1) continue;

                //if(mf->ishbondexist(r,r0)==0) continue;

                float *dp=mf->gethbondbounds(r,r0);
                int ii=0;
                while(dp&&dp[6*ii]>0) {
                      addbounds((int)dp[6*ii],(int)dp[6*ii+1],dp[6*ii+2],dp[6*ii+3],dp[6*ii+4],dp[6*ii+5]);
                      ii++;
                }
                if(dp) delete [] dp;
        }
        return norg;
}


int Bound::setclosedistance() {

	if(model==0) return ndist;

        int norg=ndist;

        StrFmt *mf=modtop->strfmt->getrootStrFmt();

	StrFmt *c;

        StrFmt *sq=mf->findsequencefmt();

        if(sq==0||sq->seqngap==0) return 0;
	
	DistPopular *dst=TRES.popbin->getfaraway();

        int n=strlen(sq->seqngap);

        int i,j;
        int n1,n2;
        int m1,m2;
        Atm *a1,*a2;
        Atm *b1,*b2;
        for(i=0;i<n;i++) {
                for(j=i+1;j<n;j++) {
                        if(sq->seqngap[i]=='-'||sq->seqngap[j]=='-') continue;
                        m1=sq->match[i];
                        m2=sq->match[j];
                        //b1=sq->resn[m1]->atm->next;
                        //b2=sq->resn[m2]->atm->next;
                        for(c=mf;c;c=c->more) {
                                if(c->seqngap[i]=='-'||c->seqngap[j]=='-') continue;
                                if(c->token&&strcmp(c->token,"sequence")==0) continue;
                                n1=c->match[i];
                                n2=c->match[j];
                                for(int ii=0;ii<=4;ii++)
                                for(int jj=0;jj<=4;jj++) {
                                        a1=c->resn[n1]->isatmid(ii);
                                        a2=c->resn[n2]->isatmid(jj);
                                        b1=sq->resn[m1]->isatmid(ii);
                                        b2=sq->resn[m2]->isatmid(jj);
                                        if(a1==0||a2==0) continue;
					//if(a1->id0>=a2->id0) continue;	
					//if(abs(b1->res->id-b2->res->id)<=1) continue;
					//cerr<<a1->res->id<<" "<<a2->res->id<<" "<<a1->res->sec<<" "<<a2->res->sec<<endl;
					//float dd=3;

					float p=0;
					if(a1->res->sec=='-') p+=2;
					if(a2->res->sec=='-') p+=2;

					if(a1->res->sec=='e') p+=1;
					if(a2->res->sec=='e') p+=1;

					if(a1->res->sec=='h') p+=0.5;
                                        if(a2->res->sec=='h') p+=0.5;


                                        float d=TRES.distance(a1,a2);
					DistPopular *pt=dst->getDistPopular(a1->res->tres,a1->name,a2->name,0);
                                        //cerr<<a1->name<<" "<<a2->name<<" "<<c->code<<d<<endl;
					float x=max(pt->dist[0],d-p);
					x=pt->dist[0];p=pt->dist[1];
                                        if(b1->id0<b2->id0) addbounds(b1->id0,b2->id0,d,x,p,2.0);
                                        else                addbounds(b2->id0,b1->id0,d,x,p,2.0);
                                }
                        }
                }
        }

	return norg;
}

int Bound::setcloseresidue() {

	int norg=ndist;

	if(model==0) return ndist;

	//DistPopular *dst=TRES.popbin->getnearnext();
	DistPopular *dst=TRES.popbin->getDistPopular("internal");
	
	if(dst==0) {
		cerr<<"distance bounds library not exist"<<endl;
		return norg;
	}
	//StrFmt *mf=modtop->strfmt->getrootStrFmt();


	Chn *c;
        Res *r,*r0;
	

	r0=0;
	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {

		Atm *a1,*a2;

		//check if constraint exists.
		int n=0;		
		r0=c->isres(r->id-1);
		for(a1=r->atm;a1;a1=a1->next)n+=predt[a1->id0];
		if(r0) {
			for(a1=r0->atm;a1;a1=a1->next)n+=predt[a1->id0];
		}
		if(n==0) continue;
		
		//find the constraint library
		DistPopular *t=dst->getDistPopular(r->tres);
		DistPopular *d;
		
		if(r0&&r->id-r0->id!=1) r0=0; //not in close contact
		
		for(d=t;d;d=d->more) {
			//if(d->resn==1) continue;
			if(r0==0&&d->resn==-1) continue;
			a1=r->isatm(d->aim);
			if(d->resn==0) {
				a2=r->isatm(d->des);
			}			
			else if(d->resn==-1&&r0) {
				a2=r0->isatm(d->des);
			}
			else {
				continue;
			}
			if(a1==0||a2==0) continue;
			//if(d->resn!=0&&a2->tatm->id>4) continue;
			if(predt[a1->id0]==0&&predt[a2->id0]==0) continue;
			float x=d->findmostpopulardistance();
			addbounds(a1->id0,a2->id0,x,d->dist[0],d->dist[1],1);
			/*
 			if(d->dist[1]-d->dist[0]<1) {
                		addbounds(a1->id0,a2->id0,x,x-0.03,x+0.03,1); 
			}
			else {
                		addbounds(a1->id0,a2->id0,x,d->dist[0],d->dist[1],1); 
			}*/
			if(TRES.logg>3) writeatmbound(stdout,ndist-1);			 			
		}		 
	}
	
	return norg;
}

 

int Bound::setfmtcloseresidue() {

	int norg=ndist;

	if(model==0) return ndist;

	DistPopular *dst=TRES.popbin->getnearnext();
	
	StrFmt *mf=modtop->strfmt->getrootStrFmt();


	Chn *c;
        Res *r,*r0;

	r0=0;
	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		DistPopular *t=dst->getDistPopular(r->tres);
		DistPopular *d;
		r0=c->isres(r->id-1);
		if(r0&&c->isnearnext(r0,r)==0) r0=0;
		Atm *a1,*a2;
		for(d=t;d;d=d->more) {
			if(d->resn==1) break;
			if(r0==0&&d->resn==-1) break;
			a1=r->isatm(d->aim);
			if(d->resn==0) {
				a2=r->isatm(d->des);
			}
			else if(d->resn==1&&r->next) {
				a2=r->next->isatm(d->des);
			}
			else if(d->resn==-1&&r0) {
				a2=r0->isatm(d->des);
			}
			else {
				continue;
			}
			if(a1==0||a2==0) continue;
			if(d->resn!=0&&a2->tatm->id>4) continue;
			float x=d->findmostpopulardistance();
			float z=d->dist[1]-d->dist[0];
			
			if(z<1) {
                		addbounds(a1->id0,a2->id0,x,d->dist[0],d->dist[1],10000);
				writeatmbound(stdout,ndist-1);
			}
			else {
				float *dp=mf->getdistbounds(a1,a2,d->dist[0],d->dist[1]);	
				if(dp==0) addbounds(a1->id0,a2->id0,x,d->dist[0],d->dist[1],1);
				if(dp[0]==0) dp[0]=x;
				int ii=0;
                        	while(dp&&dp[4*ii]>0) {
                                	addbounds(a1->id0,a2->id0,dp[4*ii],dp[4*ii+1],dp[4*ii+2],dp[4*ii+3]);
					writeatmbound(stdout,ndist-1);
                                	ii++;
                        	}
                        	if(dp) delete [] dp;
			}
		}
	}
	return norg;
}




int Bound::setcloseresidueold() {

	if(model==0) return ndist;

	Chn *c;
	Res *r;
	Atm *a,*a0;

	int j;

	int norg=ndist;

	//Atm *ring[100];

	DistPopular *dst=TRES.popbin->getnearnext();
	StrFmt *mf=modtop->strfmt->getrootStrFmt();

	//1-2 and 1-3 and 1-4 atom pairs

	for(c=model->chn;c;c=c->next)
	for(r=c->res;r;r=r->next)
	for(a=r->atm;a;a=a->next)
	for(j=0;j<10000&&a->near[j];j++) { //all atoms within +/- residues

		a0=a->near[j];

		if(a->id0>=a0->id0) continue;  //from low -> high atom id number
		if(abs(a->res->id-a0->res->id)>=2) continue;  //only neighboring residues
		if((a->res!=a0->res)&&(c->isnearnext(a->res,a0->res)==0)) continue; //only neighbors

		int n=a0->res->id-a->res->id;
		DistPopular *p=dst->getDistPopular(r->tres,a->name,a0->name,n);

        	float x=p->findmostpopulardistance();
        	addbounds(a->id0,a0->id0,x,p->dist[0],p->dist[1],1);

		//cerr<<r->name<<" "<<a->name<<" "<<a0->name<<" "<<n<<"  "<<x<<" ";
		//cerr<<p->dist[0]<<" "<<p->dist[1]<<" "<<a->nnear[j]<<endl;

		if(fabs(p->dist[1]-p->dist[0])>=1.0&&a->nnear[j]==3) {  //only 1-4 atom pairs without ring
			float *dp=mf->getdistbounds(a,a0);
			int ii=0;
            		while(dp&&dp[4*ii]>0) {
                 		addbounds(a->id0,a0->id0,dp[4*ii],dp[4*ii+1],dp[4*ii+2],dp[4*ii+3]);
                 		ii++;
            		}
            		if(dp) delete [] dp;
		}
	}


	//set 15 dist

	Atm *a1;
	for(c=model->chn;c;c=c->next)
    	for(r=c->res;r;r=r->next)
	for(a=r->atm;a;a=a->next)
	for(a1=a->next;a1;a1=a1->next) {
		if(a->isnear(a1,4)!=4) continue;
		DistPopular *p=dst->getDistPopular(r->tres,a->name,a1->name,0);
        	float x=p->findmostpopulardistance();
        	addbounds(a->id0,a1->id0,x,p->dist[0],p->dist[1],1);
		if(fabs(p->dist[1]-p->dist[0])>=1) {
			float *dp=mf->getdistbounds(a,a1);
			int ii=0;
            		while(dp&&dp[4*ii]>0) {
                  		addbounds(a->id0,a0->id0,dp[4*ii],dp[4*ii+1],dp[4*ii+2],dp[4*ii+3]);
                  		ii++;
            		}
            		if(dp) delete [] dp;
		}
	}

	return norg;
}


float Bound::evaluate() {

	float d=0;
	int n=0;
	float md=0,mx=0,my=0,me;
	//float rmax=0;
	for(int i=0;i<ndist;i++) {

		int n1,n2;
		n1=atmpair[i*2];
		n2=atmpair[i*2+1];

		float e=TRES.distance(atoms[n1]->xyz,atoms[n2]->xyz);
		cerr<<"Evaluate "<<i<<" "<<weight[i]<<" "<<atoms[n1]->res->id0<<" "<<atoms[n1]->name<<"--";
		cerr<<atoms[n2]->res->id0<<" "<<atoms[n2]->name;
		cerr<<"** "<<e<<" : "<<range[i*3]<<" "<<range[i*3+1]<<" "<<range[i*3+2]<<endl;
		//if(e*e>dist[i*3+2]) e=(e-sqrt(dist[i*3+2]));
		//else if(e*e<dist[i*3+1]) e=(sqrt(dist[i*3+1])-e);
		if(e>range[i*3+2]) e=(e-range[i*3+2]);
		else if(e<range[i*3+1]) e=(range[i*3+1]-e);
		else e=0;
		if(e!=0) cerr<<"Out of Bounds "<<i<<" "<<e<<" : "<<range[i*3]<<" "<<range[i*3+1]<<" "<<range[i*3+2]<<endl;		 
		e=e/(range[i*3+2]-range[i*3+1]);
		if(e>md) {
			md=e;
			me=e*(range[i*3+2]-range[i*3+1]);
			mx=range[i*3+1];
			my=range[i*3+2];
		}
		d+=e*e;
		n++;
	}

	if(n==0) return 0;

	d=sqrt(d/n);
	cerr<<"the maximum distortion is:"<<md<<" "<<me<<" out of:"<<mx<<" "<<my<<endl;
	return d;
}

void Bound::getsequential(Atm **ring,Res *r) {

	int i=0;

	/*
	ring[i++]=r->next->isatm(" C  ");
        ring[i++]=r->next->isatm(" O  ");
        ring[i++]=r->next->isatm(" CA ");
	*/

	Atm *c,*n;

	c=r->isatm(" C  ");
	n=r->next->isatm(" N  ");
	if(c==0||n==0) {
		return;
		//cerr<<"warning!..zero distance "
	}
	float d=TRES.distance(c,n);

	//if(r->id-r->next->id!=-1||d>2) {
	if(r->chn->isnearnext(r,r->next)==0) {
		ring[0]=0;
		ring[1]=0;
		return;
	}

	ring[i++]=r->isatm(" C  ");
        ring[i++]=r->isatm(" O  ");
        ring[i++]=r->isatm(" CA ");


	ring[i++]=r->next->isatm(" N  ");
	ring[i++]=r->next->isatm(" HN ");
	ring[i++]=r->next->isatm(" CA ");
	i=0;
	for(int j=0;j<6;j++) {
		if(ring[j]==0) continue;
		ring[i++]=ring[j];
	}
	ring[i]=0;
}

void Bound::getnewsequential(Atm **ring,Res *r) {

	int i=0;
 	
	if(r->chn->isnearnext(r,r->next)==0) {
		ring[0]=0;
		ring[1]=0;
		return;
	}
	ring[i++]=r->isatm(" C  ");
        ring[i++]=r->isatm(" O  ");
        ring[i++]=r->isatm(" CA ");
	
	ring[i++]=r->next->isatm(" N  ");
	ring[i++]=r->next->isatm(" HN ");
	ring[i++]=r->next->isatm(" CA ");
	i=0;
	for(int j=0;j<6;j++) {
		if(ring[j]==0) continue;
		ring[i++]=ring[j];
	}
	ring[i]=0;
}

void Bound::deletebound(int n) {

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

void Bound::deletefixedbound(int n) {

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


void Bound::predtorder() {
	
	//reoder the dist bounds from small to large

	float *temp;
	int   *order;
	int i,j,k,i1,i2;

	if(TRES.logg>3) cerr<<"sorting the distance constraints...."<<endl;

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
		if(predt[i1]==0||predt[i2]==0) {
			temp[i]=i-2*ndist;
		}	
		else {
			temp[i]=i;
		}		 
	}

	Qsort cc;

	cc.sort(temp,ndist,order);

	float *dis=new float[3*ndist];
	int   *atmp=new int[2*ndist];
	 
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
		}
		 
	}

	for(i=0;i<ndist;i++) {
		j=order[i];

		for(k=0;k<3;k++) {
                        dist[i*3+k]=dis[j*3+k];
                }

		for(k=0;k<2;k++) {
                        atmpair[i*2+k]=atmp[j*2+k];			 
                }
		weight[i]=weight0[j];
	}

	delete [] temp;
	delete [] order;
	delete [] dis;
	delete [] atmp;
	delete [] weight0;

}

void Bound::reorder() {

//reoder the dist bounds from small to large

	float *temp;
	int   *order;
	int i,j,k,i1,i2;

	if(TRES.logg>3) cerr<<"sorting the distance constraints...."<<endl;

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
		temp[i]=i1*natom+i2-weight[i]/maxw/2;		
	}

	Qsort cc;

	cc.sort(temp,ndist,order);

	float *dis=new float[3*ndist];
	int   *atmp=new int[2*ndist];
	 
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
		}
		 
	}

	for(i=0;i<ndist;i++) {
		j=order[i];

		for(k=0;k<3;k++) {
                        dist[i*3+k]=dis[j*3+k];
                }

		for(k=0;k<2;k++) {
                        atmpair[i*2+k]=atmp[j*2+k];			 
                }
		weight[i]=weight0[j];
	}

	delete [] temp;
	delete [] order;
	delete [] dis;
	delete [] atmp;
	delete [] weight0;
	 
}

void Bound::writeoutdistbound(char *s) {

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

void Bound::writeatmbound(char *s) {

        FILE *fp=fopen(s,"w");

        if(fp==0) {
                cerr<<"could not write file:"<<s<<endl;
                exit(0);
        }

        for(int i=0;i<ndist;i++) {
		Atm *a=atoms[atmpair[2*i]];
		Atm *b=atoms[atmpair[2*i+1]];
		Res *ra=a->res;
		Res *rb=b->res;
                fprintf(fp,"%c%i %s %c%i %s %8.3f %8.3f %8.3f %8.3f\n",ra->name,ra->id,a->name,rb->name,rb->id,b->name,weight[i], dist[3*i],dist[3*i+1],dist[3*i+2]);
        }

        fclose(fp);

}

void Bound::writeatmbound(FILE *fp,int nn) {

		int i=nn;
                Atm *a=atoms[atmpair[2*nn]];
		
                Atm *b=atoms[atmpair[2*nn+1]];
                Res *ra=a->res;
                Res *rb=b->res;
                fprintf(fp,"%c%i %s %c%i %s %8.3f %8.3f %8.3f %8.3f\n",ra->name,ra->id,a->name,rb->name,rb->id,b->name,weight[i],dist[3*i],dist[3*i+1],dist[3*i+2]);


}

void Bound::writeoutdistbound(char *s,int u0,int u1) {

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


void Bound::clearredundancy() {

        int i,n;


        cerr<<"clear distance redundancy..."<<endl;

	if(ndist==0) return;

	int *delt=new int[ndist];
	for(i=0;i<ndist;i++) delt[i]=1;

        for(i=0;i<ndist-1;i++) {
	    n=2*i;
            if(atmpair[n]==atmpair[n+2]&&atmpair[n+1]==atmpair[n+3]) {
		  fprintf(stderr,"redundancy:");
		  writeatmbound(stderr,i+1);
	          delt[i+1]=0;
            }
            else if(predt[atmpair[n]]==0&&predt[atmpair[n+1]]==0) {
		  fprintf(stderr,"redundancy:");
		  writeatmbound(stderr,i+1);
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
		weight[j]=weight[i];
                j++;
        }

	int m=ndist-j;

        ndist=j;

	delete [] delt;

        cerr<<"the total number of redundancy distance constraints is:"<<m<<endl;
}

int Bound::calcpredt(Res *r){

	int n=0;
	for(Atm *a=r->atm;a;a=a->next) n+=predt[a->id0];
	return n;
}

void Bound::setimproper() {

	Chn *c;
        Res *r;

	cerr<<"find out the improper dihedral angles..."<<endl;

	nimp=0;
	int n=model->manyatm();
	improper=new Atm*[4*n];
	for(c=model->chn;c;c=c->next)
	for(r=c->res;r;r=r->next) {
		if(calcpredt(r)==0) continue;
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

void Bound::printoutimproper() {

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

int Bound::optimize() {
      if(maxnit==0) maxnit=500;
      int n =optimize(10);
      model->write("ds.pdb");
      n= steepoptimize(500);
      cerr<<"the violated atoms:"<<n<<endl; 
      return n;
}

int Bound::getnumberofpredt() {
     int n=0;
     for(int i=0;i<natom;i++) {
 	if(predt[i]) n++;
     }
     return n;
}
 
int Bound::optimize(int maxnit){

// makes distances correct

      float dd[3],coo[3],d,len;
      int nit,k,i,j,i1,i2,nviol,nchk;
      int low,oldlow;
      int   *viol,*ip;
      float *imp;

      cerr<<"correcting the distance constraint..."<<endl;
      int allp=getnumberofpredt();
      viol=new int[natom];
      ip=new int[ndist];
      imp=new float[nimp];

      //start of iteration

      for(i=0;i<ndist;i++) ip[i]=i;

      low=ndist;

      oldlow=ndist;

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
         int ran=0;
         if (nit>10&&nviol>ndist/5) ran=1;
         if (nit>100&&nviol>ndist/10) ran=1;
         if((nit+1)%50==0) {
            if (low>9*oldlow/10&&low>100) ran=1;
            oldlow=low;
         }
         if(allp==natom&&ran==1) goto re200;

	 /*
         low=min(low,nviol);
	 

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
	 if(allp!=natom) continue;

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

int Bound::bondoptimize(int maxnit,int fff){

// makes distances correct

      float dd[3],coo[3],d,len;
      int nit,k,i,j,i1,i2,nviol,nchk;
      int low,oldlow;
      int   *viol,*ip;
      float *imp;

      cerr<<"correcting the distance constraint..."<<endl;
      int allp=getnumberofpredt();	
      if(allp==natom&&fff==3) return -1;

      viol=new int[natom];
      ip=new int[ndist];
      imp=new float[nimp];

      //start of iteration

      for(i=0;i<ndist;i++) ip[i]=i;

      low=ndist;

      oldlow=ndist;

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
	     if(fff==1&&bond[k]==-1) continue; //1-2,1-3,1-4 bond
	     if(fff==2&&bond[k]!=1&&bond[k]!=2) continue; //1-2,1-3
	     
             i1=atmpair[2*k];
	     i2=atmpair[2*k+1];
	     if(fff==3&&predt[i1]*predt[i2]==1) continue; //only to fixed
	     if(fff==4&&(atoms[i1]->tatm->id>4||atoms[i2]->tatm->id>4)) continue;
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
	 int ran=0;
	 if (nit>10&&nviol>ndist/5) ran=1;//goto re200;
         if (nit>100&&nviol>ndist/10) ran=1;//goto re200;
         if((nit+1)%50==0) {
            if (low>9*oldlow/10&&low>100) ran=1;//goto re200;
            oldlow=low;
         }
	 if(ran==1&&allp==natom) goto re200;
	

	 //if there are still violations, do another roun

	 if (nchk>0) {
           if (!((nit+1)%10==0&&chiropt==2)) continue;
         }
	 if(allp!=natom) continue;

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

int Bound::selectoptimize(int maxnit,int fff){

// makes distances correct

      float dd[3],coo[3],d,len;
      int nit,k,i,j,i1,i2,nviol,nchk;
      int low,oldlow;
      int   *viol,*ip;
      float *imp;

      cerr<<"correcting the distance constraint..."<<endl;
      int allp=getnumberofpredt();	
      if(allp==natom&&fff==3) return -1;

      viol=new int[natom];
      ip=new int[ndist];
      imp=new float[nimp];

      //start of iteration

      for(i=0;i<ndist;i++) ip[i]=i;

      low=ndist;

      oldlow=ndist;

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
	     if(fff==1&&bond[k]==-1) continue; //1-2,1-3,1-4 bond
	     if(fff==2&&bond[k]!=1&&bond[k]!=2) continue; //1-2,1-3
	     
             i1=atmpair[2*k];
	     i2=atmpair[2*k+1];
	     if(fff==3&&predt[i1]*predt[i2]==1) continue; //only to fixed
	     if(fff==4&&(atoms[i1]->tatm->id>4||atoms[i2]->tatm->id>4)) continue;
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
	 int ran=0;
	 if (nit>10&&nviol>ndist/5) ran=1;//goto re200;
         if (nit>100&&nviol>ndist/10) ran=1;//goto re200;
         if((nit+1)%50==0) {
            if (low>9*oldlow/10&&low>100) ran=1;//goto re200;
            oldlow=low;
         }
	 if(ran==1&&allp==natom) goto re200;
	

	 //if there are still violations, do another roun

	 if (nchk>0) {
           if (!((nit+1)%10==0&&chiropt==2)) continue;
         }
	 if(allp!=natom) continue;

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

    
 
int Bound::myoptimize(int maxnit){

// makes distances correct

      float dd[3],coo[3],d,len;
      int nit,k,i,j,i1,i2,nviol,nchk;
      int low,oldlow;
      int   *viol,*ip;
      float *imp;
      int allp;

      allp=getnumberofpredt();
      
      
      cerr<<"correcting the distance constraint..."<<endl;

      viol=new int[natom];
      ip=new int[ndist];
      imp=new float[nimp];

      //start of iteration

      for(i=0;i<ndist;i++) ip[i]=i;

      low=ndist;

      oldlow=ndist;

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
		 
		len=sqrt(len);	 
		float y=range[k*3+2]-range[k*3+1];
		float x=0.3*y;
		if(len>range[k*3+2])   d=(range[k*3+1]+x)/len/2;	
		else if(len<range[k*3+1]) d=(range[k*3+2]-x)/len/2;		
		
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
		   //if(i1==3||i1==4)atoms[i1]->write(stdout); 
		   //if(i2==3||i2==4)atoms[i2]->write(stdout); 
		   //cerr<<len<<" "<<TRES.distance(atoms[i1],atoms[i2])<<endl;
	        }
                else if(predt[i1]) {  //atom i to be predicted
		    for(j=0;j<3;j++) {
                        atoms[i1]->xyz[j]=atoms[i2]->xyz[j]+dd[j]*d*2;
                    }
		    //cout<<len<<" "<<TRES.distance(atoms[i1],atoms[i2])<<" "<<range[k*3]<<endl;
		    //if(i1==3||i1==4)atoms[i1]->write(stdout); 
                }
                else if(predt[i2]) { //atom j to be predicted
		    for(j=0;j<3;j++) {
			atoms[i2]->xyz[j]=atoms[i1]->xyz[j]-dd[j]*d*2;
                    }
		    //cout<<len<<" "<<TRES.distance(atoms[i1],atoms[i2])<<" "<<range[k*3]<<endl;
		    //if(i2==3||i2==4)atoms[i2]->write(stdout); 
                }
            }
         }

    	 //check if quit
         low=min(low,nviol);
	 int ran=0;
	 if(allp==natom) {
	 	if (nit>50&&nviol>ndist/2) ran=1;
         	if (nit>100&&nviol>ndist/5) ran=1;
         	if((nit+1)%50==0) {
            	if (low>oldlow*0.9&&low>100) ran=1;
            		oldlow=low;
         	}
	 }
 	 if(ran) goto re200;

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
 
	 if(nviol>0) continue;
	 //if there are still violations, do another round
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

	 //if(chiropt==0) goto re200;

	 j=0;
         for (i=0;i<nimp;i++) {
             if (imp[i]<0.0||imp[i]>70.0) {
                 swapatm(i,viol);
		 j++;
	     }
         }
	 
	 //if(chiropt==1) goto re200;

         //update list of distances to be checked and do another round if
         //required
	 if(j==0&&nviol==0) goto re200;
         //if (nchk>0) continue;
     }

     re200:

     if(nviol!=0&&chiropt==2) cerr<<"the structure rejected after: "<<nit<<" iterations"<<endl;
     else 	cerr<<"the structure accepted after: "<<nit<<" iterations"<<endl;
     if(viol) delete [] viol;
     if(ip)   delete [] ip;
     if(imp)  delete [] imp;

     return nviol;
}


int Bound::smoothoptimize(int maxnit){

// makes distances correct

      float rms,ent;
      int nit,i,j,i1,i2,nviol,nchk,allp;
      int low,oldlow;
      int   *viol,*ip;
      float *imp;
      float history[1000];
      float *steps; 
      //number of predictable atoms
      allp=getnumberofpredt();
      
      cerr<<"correcting the distance constraint..."<<endl;
	
      //violate atoms space
 
      viol=new int[natom];
      ip=new int[natom];
      imp=new float[nimp];
      steps=new float[natom];
      
      //history=new float[natom*9];

      //start of iteration

      for(i=0;i<natom;i++) ip[i]=i;
      for(i=0;i<natom;i++) steps[i]=0;

      low=ndist;

      oldlow=ndist;

      //maxnit=500;
      nchk=natom;
      for(nit=0;nit<maxnit;nit++) {

	 for(i=0;i<nchk;i++) {	   
	     	j=random();
	     	j=j%natom;
	     	int k=ip[i];
	     	ip[i]=ip[j];
	     	ip[j]=k;	     
         }

         for (i=0;i<natom;i++) viol[i]=0;
	 
         nviol=0; rms=0; ent=0;
  
	 nchk=natom;	
	 for(int kk=0;kk<nchk;kk++) {
		//int jj=ip[kk];
		Atm *a=atoms[kk];
		if(predt[a->id0]==0) continue;
		int k,i1,i2,i3;
		float d,e;		
		Atm *aa;
		i1=tag[2*a->id0];
		i2=tag[2*a->id0+1];
		int tnviol=nviol;
		float trms=rms;
		float tent=ent;
		for(int rr=0;rr<20;rr++) {
		     	float x=0,y=0,z=0;
		     	//float ex=0,ey=0,ez=0;
		     	float dx,dy,dz;
		     	nviol=tnviol;
		     	rms=trms;
			ent=tent;
			float rms0=0;
			float nviol0=0;
			float ent0=0;
			//float enpre=0;
		     	for(k=i1;k<=i2;k++) {
			
				i3=atmpair[k*2+1];
				aa=atoms[i3];
			
				dx=a->xyz[0]-aa->xyz[0];
				dy=a->xyz[1]-aa->xyz[1];
				dz=a->xyz[2]-aa->xyz[2];
				d=sqrt(dx*dx+dy*dy+dz*dz);	
 				if(d<0.1) d=0.1;
				if(d>range[3*k+2]) {
					rms0+=(d-range[3*k+2])*(d-range[3*k+2]);
					if(predt[a->id0]) viol[a->id0]=1;
					if(predt[aa->id0]) viol[aa->id0]=1;
					nviol0++;
				}
				else if(d<range[3*k+1]) {
					rms0+=(d-range[3*k+1])*(d-range[3*k+1]);
					if(predt[a->id0]) viol[a->id0]=1;
					if(predt[aa->id0]) viol[aa->id0]=1;
					nviol0++;
				}
			
				float a,b,c;

				if(d>=range[3*k]) {
					a=curve[3*k+1];
				}
				else {
					a=curve[3*k];	
				}
				b=range[3*k];
				c=curve[3*k+2];
			 
				if(d==0) d=0.1;
				e=2*a*(d-b)/d;

				//if(d>range[3*k+1]&&d<range[3*k+2]) e=0;
				x+=dx*e;		 	
				y+=dy*e;
				z+=dz*e; 
				ent0+=a*(d-b)*(d-b)+c;
			
				//cerr<<x<<" "<<y<<" "<<z<<" "<<dx<<" "<<dy<<" "<<dz<<" "<<d<<endl;
				//e=TRES.engtable->gettablevalue(a->tatm->sim->id,aa->tatm->sim->id,d);
				//en+=TRES.engtable->getenergyvalue(a->tatm->sim->id,aa->tatm->sim->id,d);		      
		   	}
			float gg=0,dgg=0,gam=0;//,h=0;
			float ddx=0,ddy=0,ddz=0;
			float re=0;
			re=x*x+y*y+z*z;
			if(re<0.1) break;	
			if(rr==0){	
				history[0]=-x;
				history[1]=-y;
				history[2]=-z;
				history[3]=-x;
				history[4]=-y;
				history[5]=-z;
				x=history[3];
				y=history[4];
				z=history[5];
				re=sqrt(re)/100;
				if(re>1) re=1;
				ddx=ddy=ddz=0.1;				 
			}
			else {
				//cerr<<"x,y,z:"<<x<<" "<<y<<" "<<z<<endl;
				gg=history[0]*history[0]+history[1]*history[1]+history[2]*history[2];
				dgg=(x+history[0])*x+(y+history[1])*y+(z+history[2])*z;
				gam=dgg/gg;
				if(fabs(history[6])<0.001) ddx=0;
				else ddx=(x+history[0])/history[6];
				if(fabs(history[7])<0.001) ddy=0;
				else ddy=(y+history[1])/history[7];
				if(fabs(history[8])<0.001) ddz=0;
				else ddz=(z+history[2])/history[8];
				if(fabs(ddx)<0.001) ddx=0;
				else ddx=-x/ddx;
				if(fabs(ddy)<0.001) ddy=0;
				else ddy=-y/ddy;
				if(fabs(ddz)<0.001) ddz=0;
				else ddz=-z/ddz;
				ddx=ddy=ddz=0.1;
				//dd=pow(-x/ddx,2)+pow(-y/ddy,2)+pow(-z/ddz,2);
				//dd=sqrt(dd);dd=0.3;
				history[0]=-x;
				history[1]=-y;
				history[2]=-z;				
				history[3]=history[0]+gam*history[3];
				history[4]=history[1]+gam*history[4];
				history[5]=history[2]+gam*history[5];
				x=history[3];
				y=history[4];
				z=history[5];				
			}
			
			//x=-x;y=-y;z=-z;dd=1.0;
			rms+=rms0;
			nviol+=(int)nviol0;
			ent+=ent0;
			if(nviol0) rms0=sqrt(rms0/nviol0);			
		   	
		   	//decides step length
			//stop here.
			e=x*x+y*y+z*z;
			if(e==0) break;
			if(e<0.1) break;
			//if(e>0)  {  
			e=sqrt(e);
			x=fabs(ddx)*x/e;
			y=fabs(ddy)*y/e;
			z=fabs(ddz)*z/e;			
			history[6]=x;			
			history[7]=y;
			history[8]=z;
                	a->xyz[0]+=x;
                	a->xyz[1]+=y;
                	a->xyz[2]+=z;
			if(fabs(ddx)+fabs(ddy)+fabs(ddz)<0.1&&rr>0) break;
			//if(fabs(x)+fabs(y)+fabs(z)<0.1) break;
			//}
                	//float d1=TRES.distance(a->xyz,aa->xyz);
                	cerr<<nit<<" "<<rr<<" "<<a->id0<<" "<<rms0<<" "<<nviol0<<" "<<ent0<<endl;
		}
	 }	
	 
	 nchk=0;
         for(i=0;i<natom;i++) {
            if(predt[i]==0) continue;
	    if(viol[i]==0)  continue;
	    ip[nchk++]=i;          
         }

 
	 cerr<<"the total energy:"<<nit<<" "<<sqrt(rms/nviol)<<" "<<ent<<" "<<nviol<<" "<<nchk<<endl;
	 if(nviol) continue;
         
	 /*
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
		 
		len=sqrt(len);	 
		float y=range[k*3+2]-range[k*3+1];
		float x=0.3*y;
		if(len>range[k*3+2])   d=(range[k*3+1]+x)/len/2;	
		else if(len<range[k*3+1]) d=(range[k*3+2]-x)/len/2;		
		
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
		   //if(i1==3||i1==4)atoms[i1]->write(stdout); 
		   //if(i2==3||i2==4)atoms[i2]->write(stdout); 
		   //cerr<<len<<" "<<TRES.distance(atoms[i1],atoms[i2])<<endl;
	        }
                else if(predt[i1]) {  //atom i to be predicted
		    for(j=0;j<3;j++) {
                        atoms[i1]->xyz[j]=atoms[i2]->xyz[j]+dd[j]*d*2;
                    }
		    //cout<<len<<" "<<TRES.distance(atoms[i1],atoms[i2])<<" "<<range[k*3]<<endl;
		    //if(i1==3||i1==4)atoms[i1]->write(stdout); 
                }
                else if(predt[i2]) { //atom j to be predicted
		    for(j=0;j<3;j++) {
			atoms[i2]->xyz[j]=atoms[i1]->xyz[j]-dd[j]*d*2;
                    }
		    //cout<<len<<" "<<TRES.distance(atoms[i1],atoms[i2])<<" "<<range[k*3]<<endl;
		    //if(i2==3||i2==4)atoms[i2]->write(stdout); 
                }
            }
         }
	 */
	 
    	 //check if quit
         low=min(low,nviol);
	 int ran=0;
	 if(allp==natom) {
	 	if (nit>50&&nviol>ndist/2) ran=1;
         	if (nit>100&&nviol>ndist/5) ran=1;
         	if((nit+1)%50==0) {
            	if (low>oldlow*0.9&&low>100) ran=1;
            		oldlow=low;
         	}
	 }
 	 if(ran) goto re200;

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
 
	 if(nviol>0) continue;
	 //if there are still violations, do another round
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

	 //if(chiropt==0) goto re200;

	 j=0;
         for (i=0;i<nimp;i++) {
             if (imp[i]<0.0||imp[i]>70.0) {
                 swapatm(i,viol);
		 j++;
	     }
         }
	 
	 //if(chiropt==1) goto re200;

         //update list of distances to be checked and do another round if
         //required
	 if(j==0&&nviol==0) goto re200;
         //if (nchk>0) continue;
     }

     re200:

     if(nviol!=0&&chiropt==2) cerr<<"the structure rejected after: "<<nit<<" iterations"<<endl;
     else 	cerr<<"the structure accepted after: "<<nit<<" iterations"<<endl;
     if(viol) delete [] viol;
     if(ip)   delete [] ip;
     if(imp)  delete [] imp;

     return nviol;
}

int Bound::steepoptimize(int maxnit){

// makes distances correct

      float rms,ent;
      int nit,i,j,i1,i2,nviol,nchk,allp;
      int low,oldlow;
      int   *viol,*ip;
      float *imp;
    
      float *history; 
    
      //number of predictable atoms
      allp=getnumberofpredt();
      
      cerr<<"correcting the distance constraint..."<<endl;
	
      //violate atoms space
 
      viol=new int[natom];
      ip=new int[natom];
      imp=new float[nimp];
      history=new float[7*natom];
  
      //start of iteration

      for(i=0;i<natom;i++) ip[i]=i;
      for(i=0;i<7*natom;i++) history[i]=0;

      low=ndist;

      oldlow=ndist;

      //maxnit=500;
      nchk=natom;
      float previous=0;
      for(nit=0;nit<maxnit;nit++) {
	 
	 for(i=0;i<nchk;i++) {	   
	     	j=random();
	     	j=j%natom;
	     	int k=ip[i];
	     	ip[i]=ip[j];
	     	ip[j]=k;	     
         }

         for (i=0;i<natom;i++) viol[i]=0;
	 
         nviol=0; rms=0; ent=0;
  
	 nchk=natom;	
	 float rmax=0;
	 for(int kk=0;kk<nchk;kk++) {
 
		Atm *a0=atoms[kk];
		if(predt[a0->id0]==0) continue;
		int k,i1,i2,i3;
		float d,e;		
		Atm *aa;
		i1=tag[2*a0->id0];
		i2=tag[2*a0->id0+1];
		if(i1==-1||i2==-1) continue;
		int tnviol=nviol;
		float trms=rms;
		float tent=ent;
		int tt=a0->id0*7;
		history[tt+0]=a0->xyz[0];
		history[tt+1]=a0->xyz[1];
		history[tt+2]=a0->xyz[2];
		//int upon=0;
	
		
		for(int rr=0;rr<20;rr++) {
		     	float x=0,y=0,z=0;
		     	float ex=0,ey=0,ez=0;
		     	float dx,dy,dz;
		     	nviol=tnviol;
		     	rms=trms;
			ent=tent;
			float rms0=0;
			float nviol0=0;
			float ent0=0;
		     	for(k=i1;k<=i2;k++) {

				float a,b,c;
				i3=atmpair[k*2+1];
				aa=atoms[i3];
				float len=(range[3*k+2]-range[3*k+1]);
				float mlen=(range[3*k+2]+range[3*k+1])/2;
				dx=a0->xyz[0]-aa->xyz[0];
				dy=a0->xyz[1]-aa->xyz[1];
				dz=a0->xyz[2]-aa->xyz[2];
				d=sqrt(dx*dx+dy*dy+dz*dz);	
 				if(d<0.1) d=0.1;
				if(d>range[3*k+2]) {
					rms0+=(d-range[3*k+2])*(d-range[3*k+2]);
					if(predt[a0->id0]) viol[a0->id0]=1;
					if(predt[aa->id0]) viol[aa->id0]=1;
					nviol0++;

					//energy
					
					a=200/len/mlen;
					
					b=range[3*k+2];
					c=0;
					e=2*a*(d-b)/d;
					x+=dx*e;		 	
					y+=dy*e;
					z+=dz*e; 
					ent0+=a*(d-b)*(d-b)+c;
					if((d-range[3*k+2])/len>rmax) rmax=(d-range[3*k+2])/len;
				}
				else if(d<range[3*k+1]) {
					rms0+=(d-range[3*k+1])*(d-range[3*k+1]);
					if(predt[a0->id0]) viol[a0->id0]=1;
					if(predt[aa->id0]) viol[aa->id0]=1;
					nviol0++;
					
					//energy
					a=100/(range[3*k+2]-range[3*k+1]);
					b=range[3*k+1];
					c=0;
					e=2*a*(d-b)/d;
					
					x+=dx*e;		 	
					y+=dy*e;
					z+=dz*e; 
					ent0+=a*(d-b)*(d-b)+c;
					if((range[3*k+1]-d)/len>rmax) rmax=(range[3*k+1]-d)/len;
				}								 			      
		   	}
 			
			float dd=0;
			float re=0;
			//
			rms+=rms0;
			nviol+=(int)nviol0;
			ent+=ent0;

			//rmsd deviation of this atom
			if(nviol0) rms0=sqrt(rms0/nviol0);	
 			 		
			//vector length	
			re=sqrt(x*x+y*y+z*z);

			if(re!=0) {
				x=x/re; y=y/re; z=z/re;
			}
			
			dd=rms0; 
			if(rr==0) {
				history[tt+0]=a0->xyz[0];
				history[tt+1]=a0->xyz[1];
				history[tt+2]=a0->xyz[2];
				 
				history[tt+3]=x;
				history[tt+4]=y;
				history[tt+5]=z;	
				 
				history[tt+6]=ent0;				 
				 							
				x=-dd*x;	
				y=-dd*y;
				z=-dd*z;
				
			}
			else if(ent0>history[tt+6]) {
				//upon++;
				dx=x*history[tt+3];
                                dy=y*history[tt+4];
                                dz=z*history[tt+5];

				ex=fabs(a0->xyz[0]-history[tt+0]);
                                ey=fabs(a0->xyz[1]-history[tt+1]);
                                ez=fabs(a0->xyz[2]-history[tt+2]);

				a0->xyz[0]=history[tt+0];
				a0->xyz[1]=history[tt+1];
				a0->xyz[2]=history[tt+2];						
				
				x=history[tt+3];
				y=history[tt+4];
				z=history[tt+5];

				if(dx<0) ex=ex/2;
                                if(dy<0) ey=ey/2;
                                if(dy<0) ez=ez/2;
                                
				if(dx>=0&&dy>=0&&dz>=0) {
					ex=ex/2;
					ey=ey/2;
					ez=ez/2;
				}

				x=-ex*x;
                                y=-ey*y;
                                z=-ez*z;
			}		 				
			else{
				dx=x*history[tt+3];
				dy=y*history[tt+4];
				dz=z*history[tt+5];				
				ex=fabs(a0->xyz[0]-history[tt+0]);
				ey=fabs(a0->xyz[1]-history[tt+1]);
				ez=fabs(a0->xyz[2]-history[tt+2]);
				
				
				if(dx<0) {
					ex=ex/2;
				}
				else {
					ex=ex*3;
				}

				if(dy<0) {
                                        ey=ey/2;
                                }
                                else {
                                        ey=ey*3;
                                }

				if(dy<0) {
                                        ez=ez/2;
                                }
                                else {
                                        ez=ez*3;
				}


				history[tt+0]=a0->xyz[0];
				history[tt+1]=a0->xyz[1];
				history[tt+2]=a0->xyz[2];
				history[tt+3]=x;
				history[tt+4]=y;
				history[tt+5]=z;	
				history[tt+6]=ent0;

				x=-ex*x;	
				y=-ey*y;
				z=-ez*z;
			}
      		 	 
			//if(ent0==0) ent0=0.0001;
			if(nviol0==0) break;
			//if(fabs(enpre-ent0)/fabs(ent0)<0.0001&&rr>0) {				 
			//	break;
			//}	
	
			a0->xyz[0]+=x;
                	a0->xyz[1]+=y;
                	a0->xyz[2]+=z;			
		}
	 }	
	 
	 nchk=0; rmax=0;  nviol=0;
	 for(i=0;i<natom;i++) viol[i]=0;
	 for(i=0;i<ndist;i++) {
		int n1,n2;
                n1=atmpair[i*2];
                n2=atmpair[i*2+1];
		float len=range[i*3+2]-range[i*3+1];
		float e=TRES.distance(atoms[n1]->xyz,atoms[n2]->xyz);
		if(e>range[i*3+2]) {
			e=(e-range[i*3+2]);
			viol[n1]=1;
			viol[n2]=1;
			 
		}
		else if(e<range[i*3+1]) {
			e=(range[i*3+1]-e);
			viol[n1]=1;
			viol[n2]=1; 
		}
		else e=0;
		float ee=e/len;
		if(ee>rmax) {
			rmax=ee;
		}
		
	 }
	 
	 nchk=0;
         for(i=0;i<natom;i++) {
            if(predt[i]==0) continue;
            if(viol[i]==0)  continue;
            ip[nchk++]=i;
	    nviol++;
         }

	 
	 cerr<<"the total energy:"<<nit<<" "<<rmax<<" "<<sqrt(rms/nviol)<<" "<<ent<<" "<<nviol<<" "<<nchk<<endl;
	 
	 if(ent==0) ent=0.0001;
	 if(nit>0&&fabs(previous-ent)/fabs(ent)<0.00001) {
		cerr<<previous-ent<<" "<<fabs(previous-ent)/fabs(ent)<<endl;
		break;
	 }
	 previous=ent;
	 if(nviol) continue;
         
	 
    	 //check if quit
         low=min(low,nviol);
	 int ran=0;
	 if(allp==natom) {
	 	if (nit>50&&nviol>ndist/2) ran=1;
         	if (nit>100&&nviol>ndist/5) ran=1;
         	if((nit+1)%50==0) {
            	if (low>oldlow*0.9&&low>100) ran=1;
            		oldlow=low;
         	}
	 }
 	 if(ran) goto re200;

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
 
	 if(nviol>0) continue;
	 //if there are still violations, do another round
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

	 //if(chiropt==0) goto re200;

	 j=0;
         for (i=0;i<nimp;i++) {
             if (imp[i]<0.0||imp[i]>70.0) {
                 swapatm(i,viol);
		 j++;
	     }
         }
	 
	 //if(chiropt==1) goto re200;

         //update list of distances to be checked and do another round if
         //required
	 if(j==0&&nviol==0) goto re200;
         //if (nchk>0) continue;
     }

     re200:

     if(nviol!=0&&chiropt==2) cerr<<"the structure rejected after: "<<nit<<" iterations"<<endl;
     else 	cerr<<"the structure accepted after: "<<nit<<" iterations"<<endl;
     if(viol)   delete [] viol;
     if(ip)     delete [] ip;
     if(imp)    delete [] imp;
     if(history)  delete [] history;
     return nviol;
}

int Bound::steepoptimizeold(int maxnit){

// makes distances correct

      float rms,ent;
      int nit,i,j,i1,i2,nviol,nchk,allp;
      int low,oldlow;
      int   *viol,*ip;
      float *imp;
    
      float *history; 
    
      //number of predictable atoms
      allp=getnumberofpredt();
      
      cerr<<"correcting the distance constraint..."<<endl;
	
      //violate atoms space
 
      viol=new int[natom];
      ip=new int[natom];
      imp=new float[nimp];
      history=new float[7*natom];
  
      //start of iteration

      for(i=0;i<natom;i++) ip[i]=i;
      for(i=0;i<7*natom;i++) history[i]=0;

      low=ndist;

      oldlow=ndist;

      //maxnit=500;
      nchk=natom;
      float previous=0;
      for(nit=0;nit<maxnit;nit++) {

	 for(i=0;i<nchk;i++) {	   
	     	j=random();
	     	j=j%natom;
	     	int k=ip[i];
	     	ip[i]=ip[j];
	     	ip[j]=k;	     
         }

         for (i=0;i<natom;i++) viol[i]=0;
	 
         nviol=0; rms=0; ent=0;
  
	 nchk=natom;	
	 
	 for(int kk=0;kk<nchk;kk++) {
 
		Atm *a0=atoms[kk];
		if(predt[a0->id0]==0) continue;
		int k,i1,i2,i3;
		float d,e;		
		Atm *aa;
		i1=tag[2*a0->id0];
		i2=tag[2*a0->id0+1];
		if(i1==-1||i2==-1) continue;
		int tnviol=nviol;
		float trms=rms;
		float tent=ent;
		int tt=a0->id0*7;
		history[tt+0]=a0->xyz[0];
		history[tt+1]=a0->xyz[1];
		history[tt+2]=a0->xyz[2];
		//int upon=0;
		float enpre=0;
		for(int rr=0;rr<20;rr++) {
		     	float x=0,y=0,z=0;
		     	float ex=0,ey=0,ez=0;
		     	float dx,dy,dz;
		     	nviol=tnviol;
		     	rms=trms;
			ent=tent;
			float rms0=0;
			float nviol0=0;
			float ent0=0;
		     	for(k=i1;k<=i2;k++) {
			
				i3=atmpair[k*2+1];
				aa=atoms[i3];
			
				dx=a0->xyz[0]-aa->xyz[0];
				dy=a0->xyz[1]-aa->xyz[1];
				dz=a0->xyz[2]-aa->xyz[2];
				d=sqrt(dx*dx+dy*dy+dz*dz);	
 				if(d<0.1) d=0.1;
				if(d>range[3*k+2]) {
					rms0+=(d-range[3*k+2])*(d-range[3*k+2]);
					if(predt[a0->id0]) viol[a0->id0]=1;
					if(predt[aa->id0]) viol[aa->id0]=1;
					nviol0++;
				}
				else if(d<range[3*k+1]) {
					rms0+=(d-range[3*k+1])*(d-range[3*k+1]);
					if(predt[a0->id0]) viol[a0->id0]=1;
					if(predt[aa->id0]) viol[aa->id0]=1;
					nviol0++;
				}
			
				float a,b,c;

				if(d>=range[3*k]) {
					a=curve[3*k+1];
				}
				else {
					a=curve[3*k];	
				}
				b=range[3*k];
				c=curve[3*k+2];
			 
				if(d==0) d=0.1;
				e=2*a*(d-b)/d;

				//if(d>range[3*k+1]&&d<range[3*k+2]) e=0; 
				x+=dx*e;		 	
				y+=dy*e;
				z+=dz*e; 
				ent0+=a*(d-b)*(d-b)+c;
								      
		   	}
 
			float dd=0;
			float re=0;
			//
			rms+=rms0;
			nviol+=(int)nviol0;
			ent+=ent0;
			if(nviol0) rms0=sqrt(rms0/nviol0);	
 
			re=sqrt(x*x+y*y+z*z);
			dd=rms0;
			if(rr==0) {

				history[tt+0]=a0->xyz[0];
				history[tt+1]=a0->xyz[1];
				history[tt+2]=a0->xyz[2];
				history[tt+3]=x;
				history[tt+4]=y;
				history[tt+5]=z;	
				history[tt+6]=ent0;

				if(re==0) break;
				if(dd<0.05) dd=0.2;								
				x=-dd*x/re;	
				y=-dd*y/re;
				z=-dd*z/re;
				
			}
			else if(ent0>history[tt+6]) {
				//upon++;
				dx=x*history[tt+3];
                                dy=y*history[tt+4];
                                dz=z*history[tt+5];

				ex=fabs(a0->xyz[0]-history[tt+0]);
                                ey=fabs(a0->xyz[1]-history[tt+1]);
                                ez=fabs(a0->xyz[2]-history[tt+2]);

				a0->xyz[0]=history[tt+0];
				a0->xyz[1]=history[tt+1];
				a0->xyz[2]=history[tt+2];						
				
				x=history[tt+3];
				y=history[tt+4];
				z=history[tt+5];

				if(dx<0) ex=ex/2;
                                if(dy<0) ey=ey/2;
                                if(dy<0) ez=ez/2;
                                
				if(dx>=0&&dy>=0&&dz>=0) {
					ex=ex/2;
					ey=ey/2;
					ez=ez/2;
				}

				x=-ex*x/re;
                                y=-ey*y/re;
                                z=-ez*z/re;
			}		 				
			else{
				dx=x*history[tt+3];
				dy=y*history[tt+4];
				dz=z*history[tt+5];				
				ex=fabs(a0->xyz[0]-history[tt+0]);
				ey=fabs(a0->xyz[1]-history[tt+1]);
				ez=fabs(a0->xyz[2]-history[tt+2]);
				
				
				if(dx<0) {
					ex=ex/2;
				}
				else {
					ex=ex*3;
				}

				if(dy<0) {
                                        ey=ey/2;
                                }
                                else {
                                        ey=ey*3;
                                }

				if(dy<0) {
                                        ez=ez/2;
                                }
                                else {
                                        ez=ez*3;
				}


				history[tt+0]=a0->xyz[0];
				history[tt+1]=a0->xyz[1];
				history[tt+2]=a0->xyz[2];
				history[tt+3]=x;
				history[tt+4]=y;
				history[tt+5]=z;	
				history[tt+6]=ent0;

				x=-ex*x/re;	
				y=-ey*y/re;
				z=-ez*z/re;
			}
      			//if(a->id0%5==0) cerr<<nit<<" "<<rr<<" "<<a->id0<<" "<<rms0<<" "<<nviol0<<" "<<sqrt(x*x+y*y+z*z)<<" "<<ent0<<endl;
			//if(x*x+y*y+z*z<0.001) break; 
			if(ent0==0) ent0=0.0001;
			if(fabs(enpre-ent0)/fabs(ent0)<0.01&&rr>0) {
				//cerr<<nit<<" "<<rr<<" "<<a->id0<<" "<<previous-ent<<" "<<fabs(previous-ent)/fabs(ent)<<endl;
				break;
			}	
			enpre=ent0;
			a0->xyz[0]+=x;
                	a0->xyz[1]+=y;
                	a0->xyz[2]+=z;			
		}
	 }	
	 
	 nchk=0;
         for(i=0;i<natom;i++) {
            if(predt[i]==0) continue;
	    if(viol[i]==0)  continue;
	    ip[nchk++]=i;          
         }
	 
 	
	 cerr<<"the total energy:"<<nit<<" "<<sqrt(rms/nviol)<<" "<<ent<<" "<<nviol<<" "<<nchk<<endl;
	 
	 if(ent==0) ent=0.0001;
	 if(nit>0&&fabs(previous-ent)/fabs(ent)<0.0001) {
		cerr<<previous-ent<<" "<<fabs(previous-ent)/fabs(ent)<<endl;
		break;
	 }
	 previous=ent;
	 if(nviol) continue;
         
	 
    	 //check if quit
         low=min(low,nviol);
	 int ran=0;
	 if(allp==natom) {
	 	if (nit>50&&nviol>ndist/2) ran=1;
         	if (nit>100&&nviol>ndist/5) ran=1;
         	if((nit+1)%50==0) {
            	if (low>oldlow*0.9&&low>100) ran=1;
            		oldlow=low;
         	}
	 }
 	 if(ran) goto re200;

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
 
	 if(nviol>0) continue;
	 //if there are still violations, do another round
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

	 //if(chiropt==0) goto re200;

	 j=0;
         for (i=0;i<nimp;i++) {
             if (imp[i]<0.0||imp[i]>70.0) {
                 swapatm(i,viol);
		 j++;
	     }
         }
	 
	 //if(chiropt==1) goto re200;

         //update list of distances to be checked and do another round if
         //required
	 if(j==0&&nviol==0) goto re200;
         //if (nchk>0) continue;
     }

     re200:

     if(nviol!=0&&chiropt==2) cerr<<"the structure rejected after: "<<nit<<" iterations"<<endl;
     else 	cerr<<"the structure accepted after: "<<nit<<" iterations"<<endl;
     if(viol)   delete [] viol;
     if(ip)     delete [] ip;
     if(imp)    delete [] imp;
     if(history)  delete [] history;
     return nviol;
}

float Bound::getta(Atm *at1,Atm *at2,Atm *at3,Atm *at4) {

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

void Bound::boundsquare(float n) {

	cerr<<"calculate distance square of distance constraint..."<<endl;
	npower=npower*n;
	for(int i=0;i<ndist;i++) {

		dist[i*3]=pow(dist[i*3],n);
		dist[i*3+1]=pow(dist[i*3+1],n);
		dist[i*3+2]=pow(dist[i*3+2],n);
	}
}

void Bound::setrange() {

	if(TRES.logg>3) cerr<<"calculate distance square of distance constraint..."<<endl;
	if(range) delete [] range;
	//if(arraysize<ndist) arraysize=ndist+100;
	range=new float[3*ndist];
	int i=0;
	for(i=0;i<ndist*3;i++) range[i]=0;
	for(i=0;i<ndist;i++) {
		range[i*3]=dist[i*3];
		range[i*3+1]=dist[i*3+1];
		range[i*3+2]=dist[i*3+2];
	}
}


void Bound::searchstructure(int n) {

	//int i=0;
	int ii=0;
	while(1) {
		ii++;
		if(ii>n) break;
		initialstructure();
		char line0[1000];
		sprintf(line0,"x%i.pdb",ii);
		model->write(line0);
		int nv=optimize();
		if(1||nv) {
			sprintf(line0,"y%i.pdb",ii);
			model->write(line0);
			//ii++;	
		}
		float d=evaluate();
		cerr<<ii<<" : "<<nv<<" "<<d<<endl;
		/*
		//if(nv==0) ii++;
		if(nv==0) {
			char line[1000];
			sprintf(line,"s%i.pdb",ii);
			model->write(line);
		}
		*/
	}

}
void Bound::initialstructure() {
	
	int cen=0;
	for(int i=0;i<natom;i++) {

		if(predt[i]==0) continue;
		cen++;
		for(int k=0;k<3;k++) {
			atoms[i]->xyz[k]=(random())%10000*0.01;
		}
	}
	
	if(cen==natom) model->center();
}

void Bound::mirror(Atm *a1,Atm *a2,Atm *a3,Atm *a4) {

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
      if(a4->id0==3||a4->id0==4) {
		a4->res->write("se");
		//a4->write(stdout);
      }
      for(i=0;i<3;i++) {
      	pp[i]=o[i]+p3*v3[i]+p4*v4[i];
	a4->xyz[i]=2*pp[i]-p[i];	
      }

}

void Bound::swapatm(Atm *a,Atm *b) {

      //swaps atom a and atom b

      float dum[3];

      for(int i=0;i<3;i++) {
         dum[i]=a->xyz[i];
         a->xyz[i]=b->xyz[i];
         b->xyz[i]=dum[i];
      }
}

void Bound::swapatm(int n,int *viol) {

      //generic swap routine for atoms involved in incorrect
      //improper dihedrals

      Atm *a1,*a2,*a3,*a4;

      //conversion of d to l amino acid

      if(improper[4*n+2]->tatm->id==2) {
         a1=improper[4*n];
	 a2=improper[4*n+1];
         a3=improper[4*n+2];
         for(a4=a1->res->atm;a4;a4=a4->next) {	    
	    if(a4->tatm->name[1]=='H'&&a4->tatm->bond[0]->id<4) continue;
	    else if(a4->tatm->name[1]!='H'&&a4->tatm->id<4) continue;	
            viol[a4->id0]=1;
            //if(a4!=a1&&a4!=a2&&a4!=a3) {
		mirror(a1,a2,a3,a4);
            //}
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


void Bound::adjust() {

      int i1,i2;

      if(TRES.logg>3) cerr<<"adjusting the distance constraint..."<<endl;

      for(int i=0;i<ndist;i++) {

	i1= atmpair[2*i];
	i2= atmpair[2*i+1];

	if(predt[i1]==0&&predt[i2]==1) {

		atmpair[2*i]=i2;
		atmpair[2*i+1]=i1;
	}
      }
}

void Bound::defaultaction() {

	adjust();
	reorder();
	clearredundancy();
	setimproper();
	settag();
	//boundsmooth(5);
	boundsquare(2);
}

int Bound::findpair(HashTable **hash,int k,int l) {

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

void Bound::boundsmooth(int tot) {

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

			//float pdi=0;

			n1=i;n2=j;n3=k;
			//pdi=dist[3*n1+1];
			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
               		dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
			dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);
			//if(dist[3*n1+1]>dist[3*n1+2]) dist[3*n1+1]=pdi;

			n1=j;n2=k;n3=i;
			//pdi=dist[3*n1+1];
			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);
			//if(dist[3*n1+1]>dist[3*n1+2]) dist[3*n1+1]=pdi;

			n1=k;n2=i;n3=j;
			//pdi=dist[3*n1+1];
			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);
			//if(dist[3*n1+1]>dist[3*n1+2]) dist[3*n1+1]=pdi;
        	}
	    }
	}

	for(i=0;i<2*ndist;i++) delete hash[i];

	delete [] hash;
}

void Bound::createnewbound() {

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


void Bound::detectbadtriangle(){

	int i;
	int n=0;
	for(i=0;i<ndist;i++) {
		
		float a=dist[3*i];
		float b=dist[3*i+1];
		float c=dist[3*i+2];
		if(a>c||a<b||b>c) {
			fprintf(stderr,"warning: bad triangle:");
			writeatmbound(stderr,i);
			dist[3*i]=(dist[3*i+1]+dist[3*i+2])/2;
			n++;
		}
	}
	fprintf(stderr,"the total number of bad triangles:%i\n",n);
}

void Bound::detectbadbound(FILE *fp) {

	int i,j;
	int a,b,c,k;
	//int n1,n2,n3;
	//float dum;
	int tot=1;
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

			//float pdi=0;

			//n1=i;n2=j;n3=k;
			
			/*
			pdi=dist[3*n1+1];
			
			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
               		dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
			dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);
			//if(dist[3*n1+1]>dist[3*n1+2]) dist[3*n1+1]=pdi;

			//n1=j;n2=i;n3=k;
			pdi=dist[3*n1+1];
			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);
			//if(dist[3*n1+1]>dist[3*n1+2]) dist[3*n1+1]=pdi;

			//n1=k;n2=j;n3=k;
			pdi=dist[3*n1+1];
			dist[3*n1+2]=min(dist[3*n1+2],dist[3*n2+2]+dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n2+1]-dist[3*n3+2]);
                        dist[3*n1+1]=max(dist[3*n1+1],dist[3*n3+1]-dist[3*n2+2]);
			//if(dist[3*n1+1]>dist[3*n1+2]) dist[3*n1+1]=pdi;
			*/
        	}
	    }
	}

	for(i=0;i<2*ndist;i++) delete hash[i];

	delete [] hash;
}


void Bound::setdssp(){

	if(dssp) delete [] dssp;

	int n=model->manyres();
	
	dssp=new int[n+1];

	int i;

	for(i=0;i<n+1;i++) dssp[i]=0;

	float *helix=new float[n+1];
	for(i=0;i<n+1;i++) helix[i]=0;

	float *sheet=new float[n+1];
        for(i=0;i<n+1;i++) sheet[i]=0;

	float *coil=new float[n+1];
        for(i=0;i<n+1;i++) coil[i]=0;


	StrFmt *mf=modtop->strfmt->getrootStrFmt();
	
	int ntt=mf->getstructurenumber();

	StrFmt *t,*t0;
	for(t0=mf;t0;t0=t0->next) {

		for(t=t0;t;t=t->more) {

			StrFmt *s=t->getsequenceStrFmt();

			if(strcmp(t->token,"structure")) continue;

			Res *r;

			//working on alphas
			r=t->pdb->chn->res;
			while(r) {
				
				if(r->sec!='h') {r=r->next;continue;}

				Res *r0;
				
				Res *rr,*rr0;

				rr=r;

				for(r0=r;r0;r0=r0->next) {			
				
					if(r0->sec!='h') break;
					rr0=r0;
				}

				int n1=rr->id0;

				int n2=rr0->id0;

				n1=t->compare[n1];

				n2=t->compare[n2];
				
				float x3=(n1+n2)/2.0;
				for(int i=n1;i<=n2;i++) {						
					if(s->seqngap[i]=='-') continue;
					int a=s->match[i];
					if(s->resn[a]==0) continue;
					float x=1-fabs(i-x3)/x3/2;
					x=sqrt(x);
					helix[a]+=t->sitescore[i]*x;
				}
				r=r0;
			}

			//working on sheet
                        r=t->pdb->chn->res;
                        while(r) {

                                if(r->sec!='e') {r=r->next;continue;}

                                Res *r0;

                                Res *rr,*rr0;

                                rr=r;

                                for(r0=r;r0;r0=r0->next) {

                                        if(r0->sec!='e') break;
                                        rr0=r0;
                                }

                                int n1=rr->id0;

                                int n2=rr0->id0;

                                n1=t->compare[n1];

                                n2=t->compare[n2];
				
				float x3=(n1+n2)/2.0;

                                for(int i=n1;i<=n2;i++) {
                                        if(s->seqngap[i]=='-') continue;
                                        int a=s->match[i];
                                        if(s->resn[a]==0) continue;
					float x=1-fabs(i-x3)/x3/2;
                                        x=sqrt(x);
                                        sheet[a]+=t->sitescore[i]*x;
                                }
                                r=r0;
                        }

			r=t->pdb->chn->res;
                        while(r) {

                                if(r->sec!='t') {r=r->next;continue;}

                                Res *r0;

                                Res *rr,*rr0;

                                rr=r;

                                for(r0=r;r0;r0=r0->next) {

                                        if(r0->sec!='t') break;
                                        rr0=r0;
                                }

                                int n1=rr->id0;

                                int n2=rr0->id0;

                                n1=t->compare[n1];

                                n2=t->compare[n2];
				float x3=(n1+n2)/2.0;

                                for(int i=n1;i<=n2;i++) {
                                        if(s->seqngap[i]=='-') continue;
                                        int a=s->match[i];
                                        if(s->resn[a]==0) continue;
					float x=1-fabs(i-x3)/x3/2;
                                        x=sqrt(x);
                                        coil[a]+=t->sitescore[i]*x;
                                }
                                r=r0;
                        }

		}
	}

	//helix
	float *chelix=new float[n+1];
	float *temp=new float[n+1];	
	int *order=new int[n+1];	
	for(i=0;i<n;i++) temp[i]=-helix[i];
	for(i=0;i<n;i++) chelix[i]=0;
	
	Qsort cc;

	cc.sort(temp,n,order);
	
	int n0=0;
	for(i=0;i<n;i++) if(temp[i]==0) break;
	for(i=0;i<n0/2;i++) { 
		int j=order[i];
		float x=temp[i]/temp[0]*100;
		chelix[j]=x;	
	}

	
	//to be done
	

	for(i=0;i<n+1;i++) {
		Res *r=model->isres(i);
		if(r==0) continue;
		if(helix[i]) cerr<<r->name<<r->id0<<" h"<<endl;
		if(sheet[i]) cerr<<r->name<<r->id0<<" e"<<endl;
		if(coil[i]) cerr<<r->name<<r->id0<<" t"<<endl;
		if(helix[i]==0&&sheet[i]==0&&coil[i]==0) cerr<<r->name<<r->id0<<" -"<<endl;
	}

	if(helix) delete [] helix;
	if(sheet) delete [] sheet;
	if(coil) delete [] coil;
}


void Bound::setpredt(int f) {

	for(int i=0;i<natom;i++) predt[i]=f;
}

void Bound::setpredt(int n,int f) {
	predt[n]=f;
}

void Bound::setpredt(Res *from,int end,int f){

	for(Res *r=from;r&&r->id0<end;r=r->next) {
		for(Atm *a=r->atm;a;a=a->next) {
			predt[a->id0]=f;
		}
	}

}

void Bound::setbondold() {

	if(bond) delete [] bond;

	bond=0;

	bond=new int[ndist];

	int i;
	for(i=0;i<ndist;i++) bond[i]=-1;
	
	Chn chn;
	char *s=TRES.getcappedsequence("AAA");
        chn.create(s);
	chn.configure();
        chn.header();
        chn.allnearbond(3);
	if(s) delete [] s;
	
	for(i=0;i<ndist;i++) {

		int n1=atmpair[2*i];
		int n2=atmpair[2*i+1];

		Atm *a1=atoms[n1];
		Atm *a2=atoms[n2];

		if(fabs(a1->res->id-a2->res->id)>1) {
			bond[i]=-1;
			//continue;
		}
		else if(a1->res==a2->res) {
			
			Res *r=chn.isres(a1->res->name,0);
			Atm *b1=r->isatmid(a1->tatm->id);
			Atm *b2=r->isatmid(a2->tatm->id);
			bond[i]=b1->isbondnear(b2);
			if(bond[i]==0) bond[i]=-1;
			//continue;
		}
		else if(a1->res->id-a2->res->id==1) {
			Atm *b1=chn.res->next->isatm(a1->name);
			Atm *b2=chn.res->isatm(a2->name);	
			if(b1==0||b2==0)  continue;
			bond[i]=b1->isbondnear(b2);
			if(bond[i]==0) bond[i]=-1;
			//continue;
		}
		else if(a1->res->id-a2->res->id==-1) {
                        Atm *b1=chn.res->isatm(a1->name);
                        Atm *b2=chn.res->next->isatm(a2->name);       
                        if(b1==0||b2==0)  continue;
                        bond[i]=b1->isbondnear(b2);
			if(bond[i]==0) bond[i]=-1;
                        //continue;
                }
		//cerr<<a1->res->name<<a1->res->id<<" "<<a1->name<<" "<<a2->res->name<<" "<<a2->res->id<<" "<<a2->name<<" "<<bond[i]<<endl;
	}	
}


void Bound::setbond() {

	if(bond) delete [] bond;

	bond=0;

	bond=new int[ndist];

	int i;
	for(i=0;i<ndist;i++) bond[i]=-1;
	 
	for(i=0;i<ndist;i++) {

		int n1=atmpair[2*i];
		int n2=atmpair[2*i+1];

		Atm *a1=atoms[n1];
		Atm *a2=atoms[n2];
		
		if(a1->res->chn!=a2->res->chn) continue;
		
		if(a1==a2) {
			bond[i]=0;
		}
		else if(abs(a1->res->id-a2->res->id)>1) {
			bond[i]=-1;
		}
		else if(a1->res==a2->res) {
			bond[i]=a1->isnewbondnear(a2);									 	 
		}
		else if(a1->res->id-a2->res->id==1) {
			bond[i]=a1->isnewbondnear(a2);	
			
		}
		else if(a1->res->id-a2->res->id==-1) {
                        bond[i]=a1->isnewbondnear(a2);
                }
		else {
			bond[i]=-1;
		}		 
	}	
}

void Bound::setcurve(float e) {
//forgot the meaning

	if(curve) delete [] curve;

	curve=0;

	curve=new float[3*ndist];

	int i;
	for(i=0;i>3*ndist;i++)curve[i]=0;

	for(i=0;i<ndist;i++) {

		float low=dist[3*i+1];
		float high=dist[3*i+2];
		float mid=dist[3*i];
 		
		if(bond[i]==1) {
			curve[3*i]=10;
			curve[3*i+1]=10;
			curve[3*i+2]=e;
		}
		else if(bond[i]==2) {
			curve[3*i]=5;
			curve[3*i+1]=5;
			curve[3*i+2]=e;
		}
		else if(high-low<0.3) {	
			dist[3*i]=(high+low)/2;
			curve[3*i]=0.2;
			curve[3*i+1]=0.2;
			curve[3*i+2]=e;					
		}
		else {		 
			if(mid-low<0.1) {
				mid+=0.1;
				dist[3*i]=mid;	
			}
			if(high-mid<0.1) {
				mid-=0.1;
				dist[3*i]=mid;	
			}
			
			float a=pow(mid-low,2);
			float b=pow(high-mid,2);
			curve[3*i]=-e/a;
			curve[3*i+1]=-e/b;
			curve[3*i+2]=e;					
		}	
		/*
		curve[3*i]=1;
                curve[3*i+1]=1;
                curve[3*i+2]=e;
		*/
		int i1=atmpair[2*i];
		int i2=atmpair[2*i+1];
		if(atoms[i1]==0||atoms[i2]==0) {
			curve[3*i]=10;
                	curve[3*i+1]=10;
                	curve[3*i+2]=10*e;
		}
	}
}


void Bound::setbuddy() {

	int *buddytemp=0;
	if(ndist==0) return;

	if(buddytemp) delete [] buddytemp; buddytemp=0;

	buddytemp=new int[natom];

	int *temp=new int[natom];

	int i;

	for(i=0;i<natom;i++) {
		temp[i]=-1;
		buddytemp[i]=-1;
	}

	for(i=0;i<natom;i++) {
		if(predt[i]==0) buddytemp[i]=0;
	}

	//prepare the seed, the first level to the fixed atoms
	int i1,i2;
	int nn=0;
	for(i=0;i<ndist;i++) {
		if(bond[i]!=1) continue;
		i1=atmpair[2*i];	
		i2=atmpair[2*i+1];
		if(predt[i1]==0) {
			if(buddytemp[i2]==-1) {
				temp[nn++]=i2;
				buddytemp[i2]=1;
			}
		}
		else if(predt[i2]==0) {
			if(buddytemp[i1]==-1) {
                                temp[nn++]=i1;
				buddytemp[i1]=1;
                        }
		}
	}
	
	
	//find all levels starting from the seed
	int j=0, st=0,pre=0;
	int n1,n2,n3;
	while(1) {
		
		pre=nn;
		for(i=st;i<nn;i++) {
			n3=temp[i];
			i1=tag[2*n3];i2=tag[2*n3+1];
			for(j=i1;j<=i2;j++) {		
				if(bond[j]!=1) continue;
				n1=atmpair[2*j];	
				n2=atmpair[2*j+1];
				if(buddytemp[n1]==-1) {
					temp[nn++]=n1;
					buddytemp[n1]=buddytemp[n2]+1;	
				} 
				else if(buddytemp[n2]==-1) {
					temp[nn++]=n2;        
                                        buddytemp[n2]=buddytemp[n1]+1; 
				}
			}
		}
		st=pre;
		if(nn==pre) break; //no new found
	}

	//reorganize the levels with distance constraints
	 
	delete [] temp;
}


void Bound::setbackboneaction() {

	int *action=0;
	if(action) delete [] action;action=0;

	action=new int[ndist];

	int i=0;

	for(i=0;i<ndist;i++) action[i]=0;

	int n=0;
	int i1,i2;
	Atm *a1,*a2;
	for(i=0;i<ndist;i++) {
		i1=atmpair[2*i];
		i2=atmpair[2*i+1];
		a1=atoms[i1];
		a2=atoms[i2];
		if(a1->tatm->id>4&&predt[i1]) continue;
		if(a2->tatm->id>4&&predt[i1]) continue;
		action[i]=1;
		n++;
	}
	cerr<<"the total backbone constraints is:"<<n<<endl;
}

void Bound::delsidechainbound() {

	int *action=0;
        if(action) delete [] action;action=0;

        action=new int[ndist];

        int i=0;

        for(i=0;i<ndist;i++) action[i]=0;

        int n=0;
        int i1,i2;
        Atm *a1,*a2;
        for(i=0;i<ndist;i++) {
                i1=atmpair[2*i];
                i2=atmpair[2*i+1];
                a1=atoms[i1];
                a2=atoms[i2];
                if(a1->tatm->id>4&&predt[i1]) continue;
                if(a2->tatm->id>4&&predt[i1]) continue;
                action[i]=1;
                n++;
        }
	
        cerr<<"the total backbone constraints after deletion is:"<<n<<endl;

	n=0;
	for(i=0;i<ndist;i++) {
		if(action[i]==0) continue;

		int j;
		for(j=0;j<3;j++) dist[3*n+j]=dist[3*i+j];

		for(j=0;j<2;j++) atmpair[2*n+j]=atmpair[2*i+j];
		weight[n]=weight[i];
		n++;
	}
	ndist=n;
}

void Bound::printoutdistance(){

	int i;

	float mid,low,high;
	float x=1/npower;
	for(i=0;i<ndist;i++) {
		if(npower==1) {
			mid=dist[3*i];
			low=dist[3*i+1];
			high=dist[3*i+2];
		}		
		else {
			mid=pow(dist[3*i],x);
                        low=pow(dist[3*i+1],x);
                        high=pow(dist[3*i+2],x);
		}
		int n1=atmpair[2*i];
		int n2=atmpair[2*i+1];
		float d=TRES.distance(atoms[n1],atoms[n2]);
		cerr<<"bound->";
		cerr<<atoms[n1]->res->name<<atoms[n1]->res->id0<<" "<<atoms[n1]->id0<<atoms[n1]->name<<"--";
		cerr<<atoms[n2]->res->name<<atoms[n2]->res->id0<<" "<<atoms[n2]->id0<<atoms[n2]->name<<" : ";				
		cerr<<d<<" "<<mid<<" "<<low<<" "<<high<<" "<<endl;
	}
}


void Bound::setsmallbound() {

	int i;
	for(i=0;i<ndist;i++) {
		int j=3*i;  
		if(bond[i]==1) {
			dist[j+1]=dist[j]-0.05;
			dist[j+2]=dist[j]+0.05;
		}
		else if(bond[i]==2) {
			dist[j+1]=dist[j]-0.1;
			dist[j+2]=dist[j]+0.1;
		}
		 
	}

}

void Bound::setatmat(Chn *chnorg) {

	if(atmat) delete [] atmat;atmat=0;
	atmat=new Atm*[natom];
	int i;
	for(i=0;i<natom;i++) atmat[i]=0;

	Chn *c;
	for(c=model->chn;c;c=c->next) {
		if(c->id!=chnorg->id) continue;
		Res *r;
		for(r=c->res;r;r=r->next) {
			Res *r0=chnorg->isresoid(r->oid);
			if(r0==0) continue;
			Atm *a;
			for(a=r->atm;a;a=a->next) {
				Atm *a0=r0->isatmid(a->tatm->id);
				if(a0==0) continue;
				atmat[a->id0]=a0;
			}
		}
	}
}

int Bound::setpdbfixbound(int flg){
	
	if(model==0) return ndist;

	Chn *c;
	Res *r;
	Atm *a;

	int j;

	int norg=ndist;

	Atm *ring[100];

	DistPopular *dst=TRES.popbin->getDistPopular("internal");
	if(dst==0) {
                cerr<<"distance bounds library not exist"<<endl;
                return norg;
        }
 
	//1-2, 1-3 and 1-4 atom pairs,using existed
	//weight 10, most important
	for(c=model->chn;c;c=c->next)
	for(r=c->res;r;r=r->next)
	for(a=r->atm;a;a=a->next) {
		for(j=0;j<100;j++) {
			if(a->near[j]==0) break;
			if(a->near[j]->res!=r) continue;
			//if(a->id0<=a->near[j]->id0) continue;
			int n=a->isnewbondnear(a->near[j]);
			cerr<<r->name<<r->id0<<" "<<a->name<<" "<<a->near[j]->res->name<<a->near[j]->res->id0<<" "<<a->near[j]->name<<endl;		
			if(n==1&&flg==0) {
				float x=TRES.distance(a,a->near[j]);
				addbounds(a->id0,a->near[j]->id0,x,x-0.03,x+0.03,5);
			}
			else if(n==2&&flg==0) {
				float x=TRES.distance(a,a->near[j]);
                       	 	addbounds(a->id0,a->near[j]->id0,x,x-0.1,x+0.1,5);
			}
			else if(n==3&&flg==0){
				float x=TRES.distance(a,a->near[j]);	
				addbounds(a->id0,a->near[j]->id0,x,x-0.5,x+0.5,5);
			}
			else if(n==1&&flg==1) {
				int m=a->near[j]->res->id-a->res->id;
				DistPopular *d=dst->getDistPopular(r->tres,a->name,a->near[j]->name,m);
				if(d==0) continue;
				float x=d->findmostpopulardistance();							
				addbounds(a->id0,a->near[j]->id0,x,d->dist[0],d->dist[1],5);				 
			}
			else if(n==2&&flg==1) {
				int m=a->near[j]->res->id-a->res->id;
				DistPopular *d=dst->getDistPopular(r->tres,a->name,a->near[j]->name,m);
				if(d==0) continue;
				float x=d->findmostpopulardistance();							
				addbounds(a->id0,a->near[j]->id0,x,d->dist[0],d->dist[1],5);
			}
			else if(n==3&&flg==1) {//database
				int m=a->near[j]->res->id-a->res->id;
				DistPopular *d=dst->getDistPopular(r->tres,a->name,a->near[j]->name,m);
				if(d==0) continue;
				float x=d->findmostpopulardistance();							
				addbounds(a->id0,a->near[j]->id0,x,d->dist[0],d->dist[1],5);
			}
			else if(n==1&&flg==2) {
				float x=TRES.distance(a,a->near[j]);
				addbounds(a->id0,a->near[j]->id0,x,x-0.03,x+0.03,5);
			}
			else if(n==2&&flg==2) {
				float x=TRES.distance(a,a->near[j]);
                       	 	addbounds(a->id0,a->near[j]->id0,x,x-0.1,x+0.1,5);
			}
			else if(n==3&&flg==2){
				float x=TRES.distance(a,a->near[j]);	
				addbounds(a->id0,a->near[j]->id0,x,x-0.5,x+0.5,5);
			}
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
		for(j=i+1;j<n;j++) {	
			if(ring[i]->res!=ring[j]->res) continue;							
			int n=ring[i]->isnewbondnear(ring[j]);
			float x=TRES.distance(ring[i],ring[j]);
			addbounds(ring[i]->id0,ring[j]->id0,x,x-0.05,x+0.05,10);	 				 
		}
	}
	
	
	//residue sequential restrictions
	//omega restrictions

	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		if(r->next==0) continue;
		
		for(i=0;i<100;i++) ring[i]=0;
			
		getnewsequential(ring,r);

		n=0;
                while(ring[n]) n++;

                for(i=0;i<n;i++)
                for(j=i+1;j<n;j++) {	
                        int n=ring[i]->isnewbondnear(ring[j]);
			float x=TRES.distance(ring[i],ring[j]);	
                       	addbounds(ring[i]->id0,ring[j]->id0,x,x-0.0002,x+0.0002,10); 	                       	
                }
	}
	//
	//return;
	//

	//phi-psi 
	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		//
		break;
		//
		if(r->next==0) continue;
		
		if(r->next->id-r->id!=1) continue;
		if(r->chn->isnearnext(r,r->next)==0) continue;
		
		Atm *ca=r->isatm(" CA ");
		Atm *ca0=r->next->isatm(" CA ");

		if(ca&&ca0) {
			int n=ca->isnewbondnear(ca0);				
			if(n==3&&flg==1) {//database								 				 
				DistPopular *d=dst->getDistPopular(r->next->tres,ca0->name,ca->name,-1);
				if(d) {					
				float x=d->findmostpopulardistance();							
				addbounds(ca->id0,ca0->id0,x,d->dist[0],d->dist[1],5);
				}
			}
			else if(n==3&&flg==0){
				float x=TRES.distance(ca,ca0);	
				addbounds(ca->id0,ca0->id0,x,x-0.3,x+0.3,5);
			}
			else if(n==3&&flg==2) {//database								 				 
				DistPopular *d=dst->getDistPopular(r->next->tres,ca0->name,ca->name,-1);
				if(d) {					
				float x=d->findmostpopulardistance();							
				addbounds(ca->id0,ca0->id0,x,d->dist[0],d->dist[1],5);
				}
			}
		}

		ca=r->isatm(" N  ");
		ca0=r->next->isatm(" N  ");

		if(ca&&ca0) {
			int n=ca->isnewbondnear(ca0);				
			if(n==3&&flg==1) {//database								 				 
				DistPopular *d=dst->getDistPopular(r->next->tres,ca0->name,ca->name,-1);
				if(d) {					
				float x=d->findmostpopulardistance();							
				addbounds(ca->id0,ca0->id0,x,d->dist[0],d->dist[1],5);
				}
			}
			else if(n==3&&flg==0){
				float x=TRES.distance(ca,ca0);	
				addbounds(ca->id0,ca0->id0,x,x-0.3,x+0.3,5);
			}
			else if(n==3&&flg==2) {//database								 				 
				DistPopular *d=dst->getDistPopular(r->next->tres,ca0->name,ca->name,-1);
				if(d) {					
				float x=d->findmostpopulardistance();							
				addbounds(ca->id0,ca0->id0,x,d->dist[0],d->dist[1],5);
				}
			}
		}

		ca=r->isatm(" C  ");
		ca0=r->next->isatm(" C  ");

		if(ca&&ca0) {
			int n=ca->isnewbondnear(ca0);				
			if(n==3&&flg==1) {//database								 				 
				DistPopular *d=dst->getDistPopular(r->next->tres,ca0->name,ca->name,-1);
				if(d) {					
				float x=d->findmostpopulardistance();							
				addbounds(ca->id0,ca0->id0,x,d->dist[0],d->dist[1],5);
				}
			}
			
			else if(n==3&&flg==0){
				float x=TRES.distance(ca,ca0);	
				addbounds(ca->id0,ca0->id0,x,x-0.3,x+0.3,5);
			}
			else if(n==3&&flg==2) {//database								 				 
				DistPopular *d=dst->getDistPopular(r->next->tres,ca0->name,ca->name,-1);
				if(d) {					
				float x=d->findmostpopulardistance();							
				addbounds(ca->id0,ca0->id0,x,d->dist[0],d->dist[1],5);
				}
			}
		}
	}

	
	//quasi ring restrictions,ring interacts with ca
	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {
		if(strchr("YFWHP",r->name)) {
			for(i=0;i<100;i++) ring[i]=0;
			getring(ring,r);
		}
		else continue;

		n=0;
		while(ring[n]) n++;
		a=r->isatm(" CA ");
		if(a==0) continue;
		for(i=0;i<n;i++) {
			if(ring[i]->res!=r) continue;
			int n=a->isnewbondnear(ring[i]);
			if(n==3&&flg==1) {//database
				int m=ring[i]->res->id-a->res->id;
				if(m!=-1&&m!=0) continue;
				DistPopular *d=dst->getDistPopular(r->tres,a->name,ring[i]->name,m);	
				if(d==0) continue;
				float x=d->findmostpopulardistance();							
				addbounds(a->id0,ring[i]->id0,x,d->dist[0],d->dist[1],5);
			}
			else if(n==3&&flg==0){
				float x=TRES.distance(a,ring[i]);	
				addbounds(a->id0,ring[i]->id0,x,x-0.3,x+0.3,5);
			}
			else if(n==3&&flg==2) {//database
				int m=ring[i]->res->id-a->res->id;
				if(m!=-1&&m!=0) continue;
				DistPopular *d=dst->getDistPopular(r->tres,a->name,ring[i]->name,m);	
				if(d==0) continue;
				float x=d->findmostpopulardistance();							
				addbounds(a->id0,ring[i]->id0,x,d->dist[0],d->dist[1],5);
			}			
		}
	}
	
	
	//all others	
        Res *r0;
	
	r0=0;
	for(c=model->chn;c;c=c->next)
        for(r=c->res;r;r=r->next) {

		//
		//break;
		//
		Atm *a1,*a2;

		//check if constraint exists.
		 	
		r0=c->isres(r->id-1);
		
		
		//find the constraint library
		DistPopular *t=dst->getDistPopular(r->tres);
		DistPopular *d;
		
		if(r0&&r->id-r0->id!=1) r0=0; //not in close contact
		
		for(d=t;d;d=d->more) {
				
			if(r0==0&&d->resn==-1) continue;
			if(d->resn!=0) continue; //only the same residue
			a1=r->isatm(d->aim);
			if(d->resn==0) {
				a2=r->isatm(d->des);
			}			
			else if(d->resn==-1&&r0) {
				a2=r0->isatm(d->des);
			}
			else {
				continue;
			}
			if(a1==0||a2==0) continue;
			if(flg==1) {
				float x=d->findmostpopulardistance();
				addbounds(a1->id0,a2->id0,x,d->dist[0],d->dist[1],1);			
				if(TRES.logg>3) writeatmbound(stdout,ndist-1);	
			}	
			else if(flg==0) {
				float x=TRES.distance(a1,a2);	
				addbounds(a1->id0,a2->id0,x,x-0.3,x+0.3,1);
			}	
			else if(flg==2) {
				float x=d->findmostpopulardistance();
				addbounds(a1->id0,a2->id0,x,d->dist[0],d->dist[1],1);			
				if(TRES.logg>3) writeatmbound(stdout,ndist-1);	
			}				
		}		 
	}	


	return norg;

}
