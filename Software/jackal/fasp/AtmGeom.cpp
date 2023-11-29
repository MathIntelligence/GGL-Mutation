#include "source.h"

AtmGeom::AtmGeom() {
num=0;
int i;
for(i=0;i<3;i++) coo[i]=0;
for(i=0;i<300;i++) weight[i]=1;
}

AtmGeom::~AtmGeom() {

}

void AtmGeom::addbounds(float *x,float d,float off) {
	addbounds(x,d,off,1);
}

void AtmGeom::addbounds(float *x,float d,float off,int w) {
        int n=num;
        dist[n*3]=d*d;
        dist[n*3+1]=(d-off)*(d-off);
        dist[n*3+2]=(d+off)*(d+off);
        int i;
        for(i=0;i<3;i++) xyz[3*n+i]=x[i];
	weight[n]=w;
        n++;
        num=n;
}

void AtmGeom::findmincoo(float er) {

	float sv[3],svv[3];
	int k,i,j,m;
	
	//find initial
	findcoo();

	//copy
	for(k=0;k<3;k++) sv[k]=coo[k];
	 
	//search the minimal
	float dd=1000;
	int n=0;
	for(i=-5;i<=5;i++)
	for(j=-5;j<=5;j++)
	for(k=-5;k<=5;k++){
		 svv[0]=sv[0]+i*er;
		 svv[1]=sv[1]+j*er;
		 svv[2]=sv[2]+k*er;
		 float d=calcscore(svv);
		 if(n==0) {
			dd=d;
		 }
		 else if(d<dd) {
			for(m=0;m<3;m++) coo[m]=svv[m];
			dd=d;
		 }
		 n++;
	}
}

void AtmGeom::findmincoo(int ne,float er) {

        float sv[3],svv[3];
        int k,i,j,m;

        //find initial
        findcoo();

        //copy
        for(k=0;k<3;k++) sv[k]=coo[k];

        //search the minimal
        float dd=1000;
        int n=0;
        for(i=-ne;i<=ne;i++)
        for(j=-ne;j<=ne;j++)
        for(k=-ne;k<=ne;k++){
                 svv[0]=sv[0]+i*er;
                 svv[1]=sv[1]+j*er;
                 svv[2]=sv[2]+k*er;
                 float d=calcscore(svv);
                 if(n==0) {
                        dd=d;
                 }
                 else if(d<dd) {
                        for(m=0;m<3;m++) coo[m]=svv[m];
                        dd=d;
                 }
                 n++;
        }
}

float AtmGeom::calcscore(float *svv) {

	float dd=0;
	float x,y,z,d;
	int j;
	for(int i=0;i<num;i++) {
		j=3*i;
		x=svv[0]-xyz[j];
                y=svv[1]-xyz[j+1];
                z=svv[2]-xyz[j+2];
                d=x*x+y*y+z*z;
		d=sqrt(d)-sqrt(dist[j]);
		dd+=d*d*weight[i];
	}
	return dd;
}


void AtmGeom::findcoo() {

        int nn,n,i,j,m;

        float d,x,y,z,e;

        for(nn=0;nn<50;nn++) {
           if(nn) {
                coo[0]=random();
                coo[1]=random();
                coo[2]=random();
           }

           for(n=0;n<100;n++) {
		
                m=0;
                for(i=0;i<num;i++) {
                        j=3*i;
                        x=coo[0]-xyz[j];
                        y=coo[1]-xyz[j+1];
                        z=coo[2]-xyz[j+2];
                        d=x*x+y*y+z*z;
                        if(d>dist[j+2]||d<dist[j+1]) {
                                if(d<dist[j]/10) d=dist[j]/10;
                                e=dist[j]/d;
                                e=sqrt(e);
                                coo[0]=xyz[j]+e*x;
                                coo[1]=xyz[j+1]+e*y;
                                coo[2]=xyz[j+2]+e*z;
                                m++;
                        }
                }
                if(m==0) return;
           }
        }
}

 
float AtmGeom::evaluate(){

	float d,x,y,z;
	float dd=0;
	for(int i=0;i<num;i++) {
		int j=3*i;
		x=coo[0]-xyz[j];
                y=coo[1]-xyz[j+1];
                z=coo[2]-xyz[j+2];
                d=x*x+y*y+z*z;
		d=sqrt(d)-sqrt(dist[j]);	
		dd+=d*d;
	}
	return sqrt(dd/num);
}

float AtmGeom::evaluaterange(){

        float d,x,y,z;
        float dd=0;
        for(int i=0;i<num;i++) {
                int j=3*i;
                x=coo[0]-xyz[j];
                y=coo[1]-xyz[j+1];
                z=coo[2]-xyz[j+2];
                d=x*x+y*y+z*z;
		if(d>dist[j+2]||d<dist[j+1]) {
                  d=sqrt(d)-sqrt(dist[j]);
		}
		else d=0;
                dd+=d*d;
        }
        return sqrt(dd/num);
}

float AtmGeom::evaluatemin(){

        float d,x,y,z;
        float dd=0;
        for(int i=0;i<num;i++) {
                int j=3*i;
                x=coo[0]-xyz[j];
                y=coo[1]-xyz[j+1];
                z=coo[2]-xyz[j+2];
                d=x*x+y*y+z*z;
                d=sqrt(d)-sqrt(dist[j]);
                dd+=d*d;
        }
        return sqrt(dd/num);
}

