#include"source.h"

DistDatabase::DistDatabase()
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
	cutoff=0.05;
}

DistDatabase::DistDatabase(Tres *tr)
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
	cutoff=0.05;
}

DistDatabase::DistDatabase(Tres *tr,char *a,char *b,int n)
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
	cutoff=0.05;
}

DistDatabase::~DistDatabase()
{
	if(popular) delete [] popular;
	if(aim) delete [] aim;
	if(des)  delete [] des;
	if(more) delete more;
	if(next) delete next;
}


void DistDatabase::write(FILE *fp) {

	DistDatabase *target=0;
	float e[10];

	if(tres==0) fprintf(fp,"!start residue internal:all\n");
	else fprintf(fp,"!start residue internal:%c\n",tres->name);

	for(target=this;target;target=target->more) {
	
		int i,n;

		n=0;
		for(i=(int)(target->dist[0]*10);i<=(int)(target->dist[1]*10);i++) n+=(int)(target->popular[i]);

		int nnn=n;

		if(n==0) {
                        cerr<<"n is zero in DistDatabase::write"<<endl;
                        continue;
                }

		//for(i=target->dist[0]*10;i<=target->dist[1]*10;i++) target->popular[i]=target->popular[i]/n*100;
			
		int n1,n2;
		n1=(int)(target->dist[0]*10);n2=(int)(target->dist[1]*10);
		float p1=0;
		for(i=(int)(target->dist[0]*10);i<=(int)(target->dist[1]*10);i++) {
			if(p1+target->popular[i]>cutoff*n) break;
			p1+=target->popular[i];
			n1=i;
		}
		
		float p2=0;

		for(i=(int)(target->dist[1]*10);i>=(int)(target->dist[0]*10);i--) {
                        if(p2+target->popular[i]>cutoff*n) break;
                        p2+=target->popular[i]; 
			n2=i;
                }
		
		target->dist[0]=n1/10.;
		target->dist[1]=n2/10.;
	
		if(target->tres)   
		fprintf(fp,"%c %s %s %2i %4.1f %4.1f %6i",target->tres->name,target->aim,target->des,target->resn,target->dist[0],target->dist[1],nnn);
		else 
		fprintf(fp,"%c %s %s %2i %4.1f %4.1f %6i",'0',target->aim,target->des,target->resn,target->dist[0],target->dist[1],nnn);

			

		int m;
		float d=(target->dist[1]-target->dist[0])/10;
		for(m=0;m<10;m++) {
			p1=target->dist[0]+m*d;
			p2=target->dist[0]+(m+1)*d;
			float ee=0;
			for(i=(int)(p1*10);i<=(int)(p2*10);i++) {
				ee+=target->popular[i];
			}
			e[m]=ee;
			//fprintf(fp," %4.1f",e);
		}

		p1=0;
		for(m=0;m<10;m++) p1+=e[m];

		for(m=0;m<10;m++) {
			fprintf(fp," %4.1f",e[m]/p1*100);
		}
		fprintf(fp,"\n");

		/*
		for(i=target->dist[0]*10;i<=target->dist[1]*10;i++) {
						
			fprintf(fp,"popular: %i %4.1f %4.1f\n",m,i/10.,target->popular[i]*1.0/n*100);
			m++;
		}
		*/
	}	
	if(tres) fprintf(fp,"!end residue internal:%c\n",tres->name);
	else fprintf(fp,"!end residue internal:all\n");
}


void DistDatabase::write(FILE *fp,int ff) {

	DistDatabase *target=0;
	float e[10];

	if(tres==0) fprintf(fp,"!start residue internal:all\n");
	else fprintf(fp,"!start residue internal:%c\n",tres->name);

	for(target=this;target;target=target->more) {
	
		int i,n;

		n=0;
		for(i=(int)(target->dist[0]*10);i<=(int)(target->dist[1]*10);i++) n+=(int)(target->popular[i]);

		int nnn=n;

		if(n==0) {
                        cerr<<"n is zero in DistDatabase::write"<<endl;
                        continue;
                }

		//for(i=target->dist[0]*10;i<=target->dist[1]*10;i++) target->popular[i]=target->popular[i]/n*100;
			
		int n1,n2;
		n1=(int)(target->dist[0]*10);n2=(int)(target->dist[1]*10);
		float p1=0;
		for(i=(int)(target->dist[0]*10);i<=(int)(target->dist[1]*10);i++) {
			if(p1+target->popular[i]>cutoff*n) break;
			p1+=target->popular[i];
			n1=i;
		}
		
		float p2=0;

		for(i=(int)(target->dist[1]*10);i>=(int)(target->dist[0]*10);i--) {
                        if(p2+target->popular[i]>cutoff*n) break;
                        p2+=target->popular[i]; 
			n2=i;
                }
		
		target->dist[0]=n1/10.;
		target->dist[1]=n2/10.;
	
		if(ff) {
		  if(target->tres)   
		    fprintf(fp,"%c %s %s %c %4.1f %4.1f %6i",target->tres->name,target->aim,target->des,char(target->resn),target->dist[0],target->dist[1],nnn);
		  else 
		    fprintf(fp,"%c %s %s %c %4.1f %4.1f %6i",'0',target->aim,target->des,char(target->resn),target->dist[0],target->dist[1],nnn);
		}
		else {
		  if(target->tres)   
		       fprintf(fp,"%c %s %s %c %4.1f %4.1f %6i",target->tres->name,target->aim,target->des,char(target->resn),target->dist[0],target->dist[1],nnn);
	          else 
		      fprintf(fp,"%c %s %s %c %4.1f %4.1f %6i",'0',target->aim,target->des,char(target->resn),target->dist[0],target->dist[1],nnn);
		}
			

		int m;
		float d=(target->dist[1]-target->dist[0])/10;
		for(m=0;m<10;m++) {
			p1=target->dist[0]+m*d;
			p2=target->dist[0]+(m+1)*d;
			float ee=0;
			for(i=(int)(p1*10);i<=(int)(p2*10);i++) {
				ee+=target->popular[i];
			}
			e[m]=ee;
			//fprintf(fp," %4.1f",e);
		}

		p1=0;
		for(m=0;m<10;m++) p1+=e[m];

		for(m=0;m<10;m++) {
			fprintf(fp," %4.1f",e[m]/p1*100);
		}
		fprintf(fp,"\n");

		/*
		for(i=target->dist[0]*10;i<=target->dist[1]*10;i++) {
						
			fprintf(fp,"popular: %i %4.1f %4.1f\n",m,i/10.,target->popular[i]*1.0/n*100);
			m++;
		}
		*/
	}	
	if(tres) fprintf(fp,"!end residue internal:%c\n",tres->name);
	else fprintf(fp,"!end residue internal:all\n");
}

void DistDatabase::addbounds(Tres *tr,char *a, char *b,int n,float d) {

	//find target with the same tr
	DistDatabase *target=0;
	for(target=this;target;target=target->next) {
 		if(target->tres==tr) break;
	}

	//if no such target
	if(target==0) {
		for(target=this;target->next;target=target->next);
		if(target==0) {cerr<<"target is zero"<<endl;return;}
		target->next=new DistDatabase(tr,a,b,n);
		target=target->next;
	}

	if(target==0) {
		cerr<<"target should not be zero..."<<endl;
	}

	
	DistDatabase *tarmore=0;

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
		tarmore->more=new DistDatabase(tr,a,b,n);
		tarmore=tarmore->more;
	}

	 
	if(tarmore==0) {
		cerr<<"tarmore is zero"<<endl;
		return;
	}
	tarmore->addbounds(d);
}

 
void DistDatabase::addbounds(float d) {
 
	dist[0]=min(dist[0],d);
	dist[1]=max(dist[1],d);
 
	int nn=(int)(d*10);
	if(popular==0||size==0) {
		popular=new float[nn+50];
		for(int i=0;i<nn+50;i++) popular[i]=0;
		popular[nn]++;
		size=nn+50;
		return;
	}
	else if(nn<size) {
		popular[nn]++;
		return;
	}
	else {		 
		float *t=new float[nn+50];
		for(int i=0;i<nn+50;i++) t[i]=0;
		for(int ii=0;ii<size;ii++) t[ii]=popular[ii];
		delete [] popular; 
		popular=t;
		popular[nn]++;
		size=nn+50;
		return;
	}
}
 

