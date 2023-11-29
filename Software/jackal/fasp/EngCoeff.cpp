#include"source.h"
EngCoeff::EngCoeff()
{
next=0;
flag=0;
dots=0;
coeff=0;
size=0;
}
void EngCoeff::clear(){
	if(dots) delete [] dots;
 	if(coeff) delete [] coeff;
	dots=0;
	coeff=0;
	size=0;
}

void EngCoeff::setdots(float x,float y,float w) {
	
	if(dots==0) {
		dots=new float[300];
		for(int i=0;i<300;i++) dots[i]=0;
		size=0;
	}

	dots[3*size]=x;
	dots[3*size+1]=y;
	dots[3*size+2]=w;
	size++;
}

EngCoeff::~EngCoeff()
{
 if(dots) delete [] dots;
 if(coeff) delete [] coeff;
 if(next) delete next;
 next=0;
}


float EngCoeff::getmatch() {

	int i,j,k;
	if(coeff==0) coeff=new float[100];
	for(i=0;i<100;i++) coeff[i]=0;
	

	float dmin=1000000;
	for(i=0;i<10;i++) 
	for(j=0;j<10;j++)
	for(k=0;k<10;k++) {
		float b[3];
		b[0]=(i+0.5)*0.1;
		b[1]=(j+0.5)*0.1;
		b[2]=(k+0.5)*0.1;

		float a[3];
		float pp[3][3];
		float kk[3];
		
		for(int h=0;h<3;h++) {
		   kk[h]=0;
		   for(int f=0;f<size;f++) {
		     kk[h]+=dots[3*f+2]*exp(-b[h]*dots[3*f]*dots[3*f])*dots[3*f+1];
		   }
		   for(int r=0;r<3;r++) {
		     pp[h][r]=0;
		     for(int f=0;f<size;f++) {
			pp[h][r]+=dots[3*f+2]*exp(-b[h]*dots[3*f]*dots[3*f])*
 					      exp(-b[r]*dots[3*f]*dots[3*f]);
		     }		
		   }	
		}
 		
		float d=getdelta(pp[0][0],pp[1][0],pp[2][0],pp[0][1],pp[1][1],pp[2][1],pp[0][2],pp[1][2],pp[2][2]);

		if(fabs(d)<=0.001) continue;

		float dx=getdelta(kk[0],kk[1],kk[2],pp[0][1],pp[1][1],pp[2][1],pp[0][2],pp[1][2],pp[2][2]);
		float dy=getdelta(pp[0][0],pp[1][0],pp[2][0],kk[0],kk[1],kk[2],pp[0][2],pp[1][2],pp[2][2]);
		float dz=getdelta(pp[0][0],pp[1][0],pp[2][0],pp[0][1],pp[1][1],pp[2][1],kk[0],kk[1],kk[2]);
		a[0]=dx/d;
		a[1]=dy/d;
		a[2]=dz/d;
		float diff=0;
	 	for(int f=0;f<size;f++) {
			float di=0;
			for(int r=0;r<3;r++) {
		  		di+=a[r]*exp(-b[r]*dots[3*f]*dots[3*f]);
			}
		   	di=di-dots[3*f+1];
		   	diff+=dots[3*f+2]*di*di;	
		}
		if(diff<dmin) {
			dmin=diff;
			for(int r=0;r<3;r++) {
				coeff[2*r]=a[r];
				coeff[2*r+1]=b[r];
			}
		}
	}
	cerr<<"dmin.."<<dmin<<endl;
	dmin=1000000;
	float bx=coeff[1];
	float by=coeff[3];
	float bz=coeff[5];
	for(i=-10;i<=10;i++) 
	for(j=-10;j<=10;j++)
	for(k=-10;k<=10;k++) {
		float b[3];
		b[0]=bx+(i+0.5)*0.1/20;
		b[1]=by+(j+0.5)*0.1/20;
		b[2]=bz+(k+0.5)*0.1/20;
		float a[3];
		float pp[3][3];
		float kk[3];
		
		for(int h=0;h<3;h++) {
		   kk[h]=0;
		   for(int f=0;f<size;f++) {
		     kk[h]+=dots[3*f+2]*exp(-b[h]*dots[3*f]*dots[3*f])*dots[3*f+1];
		   }
		   for(int r=0;r<3;r++) {
		     pp[h][r]=0;
		     for(int f=0;f<size;f++) {
			pp[h][r]+=dots[3*f+2]*exp(-b[h]*dots[3*f]*dots[3*f])*
 					      exp(-b[r]*dots[3*f]*dots[3*f]);
		     }		
		   }	
		}
 		
		float d=getdelta(pp[0][0],pp[1][0],pp[2][0],pp[0][1],pp[1][1],pp[2][1],pp[0][2],pp[1][2],pp[2][2]);

		if(fabs(d)<=0.001) continue;

                float dx=getdelta(kk[0],kk[1],kk[2],pp[0][1],pp[1][1],pp[2][1],pp[0][2],pp[1][2],pp[2][2]);
                float dy=getdelta(pp[0][0],pp[1][0],pp[2][0],kk[0],kk[1],kk[2],pp[0][2],pp[1][2],pp[2][2]);
                float dz=getdelta(pp[0][0],pp[1][0],pp[2][0],pp[0][1],pp[1][1],pp[2][1],kk[0],kk[1],kk[2]);

		a[0]=dx/d;
		a[1]=dy/d;
		a[2]=dz/d;
		float diff=0;
	 	for(int f=0;f<size;f++) {
			float di=0;
			for(int r=0;r<3;r++) {
		  		di+=a[r]*exp(-b[r]*dots[3*f]*dots[3*f]);
			}
		   	di=di-dots[3*f+1];
		   	diff+=dots[3*f+2]*di*di;	
		}
		if(diff<dmin) {
			dmin=diff;
			for(int r=0;r<3;r++) {
				coeff[2*r]=a[r];
				coeff[2*r+1]=b[r];
			}
		}
	}

	cerr<<"dmin.."<<dmin<<endl;
	return dmin;
}

float EngCoeff::get2match() {

	int i,j,k;
	if(coeff==0) coeff=new float[100];
	for(i=0;i<100;i++) coeff[i]=0;
	

	float dmin=1000000;
	for(i=0;i<5;i++) 
	for(j=0;j<5;j++)
	{
		float b[2];
		b[0]=(i+0.5)*0.1;
		b[1]=(j+0.5)*0.1;
		 

		float a[2];
		float pp[2][2];
		float kk[2];
		
		for(int h=0;h<2;h++) {
		   kk[h]=0;
		   for(int f=0;f<size;f++) {
		     kk[h]+=dots[3*f+2]*exp(-b[h]*dots[3*f]*dots[3*f])*dots[3*f+1];
		   }
		   for(int r=0;r<2;r++) {
		     pp[h][r]=0;
		     for(int f=0;f<size;f++) {
			pp[h][r]+=dots[3*f+2]*exp(-b[h]*dots[3*f]*dots[3*f])*
 					      exp(-b[r]*dots[3*f]*dots[3*f]);
		     }		
		   }	
		}
 		
		float d=getdelta(pp[0][0],pp[1][0],pp[0][1],pp[1][1]);

		if(fabs(d)<=0.001) continue;

		float dx=getdelta(kk[0],kk[1],pp[0][1],pp[1][1]);
		float dy=getdelta(pp[0][0],pp[1][0],kk[0],kk[1]);
		 
		a[0]=dx/d;
		a[1]=dy/d;
		 
		float diff=0;
	 	for(int f=0;f<size;f++) {
			float di=0;
			for(int r=0;r<2;r++) {
		  		di+=a[r]*exp(-b[r]*dots[3*f]*dots[3*f]);
			}
		   	di=di-dots[3*f+1];
		   	diff+=dots[3*f+2]*di*di;	
		}
		if(diff<dmin) {
			dmin=diff;
			for(int r=0;r<2;r++) {
				coeff[2*r]=a[r];
				coeff[2*r+1]=b[r];
			}
		}
	}
	cerr<<"dmin:"<<dmin<<endl;
	dmin=1000000;
	float bx=coeff[1];
	float by=coeff[3];
	 
	for(i=-5;i<=5;i++) 
	for(j=-5;j<=5;j++)
	 {
		float b[2];
		b[0]=bx+(i+0.5)*0.01;
		b[1]=by+(j+0.5)*0.01;
	 
		float a[2];
		float pp[2][2];
		float kk[2];
		
		for(int h=0;h<2;h++) {
		   kk[h]=0;
		   for(int f=0;f<size;f++) {
		     kk[h]+=dots[3*f+2]*exp(-b[h]*dots[3*f]*dots[3*f])*dots[3*f+1];
		   }
		   for(int r=0;r<2;r++) {
		     pp[h][r]=0;
		     for(int f=0;f<size;f++) {
			pp[h][r]+=dots[3*f+2]*exp(-b[h]*dots[3*f]*dots[3*f])*
 					      exp(-b[r]*dots[3*f]*dots[3*f]);
		     }		
		   }	
		}
 		
		float d=getdelta(pp[0][0],pp[1][0],pp[0][1],pp[1][1]);

		if(fabs(d)<=0.001) continue;

		float dx=getdelta(kk[0],kk[1],pp[0][1],pp[1][1]);
		float dy=getdelta(pp[0][0],pp[1][0],kk[0],kk[1]);
		 
		a[0]=dx/d;
		a[1]=dy/d;
		float diff=0;
	 	for(int f=0;f<size;f++) {
			float di=0;
			for(int r=0;r<2;r++) {
		  		di+=a[r]*exp(-b[r]*dots[3*f]*dots[3*f]);
			}
		   	di=di-dots[3*f+1];
		   	diff+=dots[3*f+2]*di*di;	
		}
		if(diff<dmin) {
			dmin=diff;
			for(int r=0;r<2;r++) {
				coeff[2*r]=a[r];
				coeff[2*r+1]=b[r];
			}
		}
	}
	cerr<<"dmin:"<<dmin<<endl;
	return dmin;
}

float EngCoeff::getdelta(float a1,float a2,float a3,float b1,float b2,float b3,float c1,float c2,float c3) {

	float d=0;
	d=a1*b2*c3+b1*c2*a3+c1*b3*a2-c1*b2*a3-b1*a2*c3-a1*b3*c2;
	return d;
}
float EngCoeff::getdelta(float a1,float a2,float b1,float b2) {

	float d=0;
	d=a1*b2-b1*a2;
	return d;
}

void EngCoeff::setnewcoeff(float t,float n){
// E(vdw) = Aexp(-B*r^2)(1/(r+t)^2n-C/(r+t)^n)
// C=(2*n+B*(1+t))/[B*(1+t)^(n+1)+n*(1+t)^n]
// A=-exp(B)/[1/(1+t)^2n-C/(1+t)^n]

/*
// E(vdw) = Aexp(-B*r^2)(1/(r+t)-C/(r+t)^0.5)
// C=(1+2*B*(1+t))/[2*B*(1+t)^1.5+0.5*(1+t)^0.5]
// A=-exp(B)/[1/(1+t)-C/(1+t)^0.5]
*/
        int i,j,k;
        if(coeff==0) coeff=new float[100];
        for(i=0;i<100;i++) coeff[i]=0;

        float dmin=1000000;
        for(i=0;i<2000;i++) {
                 float b=0.001*(i+0.5);
                 float c=(2*n+b*(1+t))/(b*pow(1+t,1+n)+n*pow(1+t,n));
                 float a=-exp(b)/(1/pow(1+t,2*n)-c/pow(1+t,n));
                 float diff=0;
                 for(int f=0;f<size;f++) {
                        float x=dots[3*f];
                        float y=dots[3*f+1];
                        float w=dots[3*f+2];
                        float di=a*exp(-b*x*x)*(1/pow(x+t,2*n)-c/pow(x+t,n));
                        diff+=w*(di-y)*(di-y);
                }
                if(diff<dmin) {
                        dmin=diff;
                        coeff[0]=a;
                        coeff[1]=b;
                        coeff[2]=c;
                        coeff[3]=t;
                }
        }
        coeff[4]=n;

        if(TRES.logg>3) cerr<<"the minimal difference:"<<sqrt(dmin/size)<<endl;
        if(TRES.logg>3) cerr<<"the coefficient:"<<coeff[0]<<" "<<coeff[1]<<" "<<coeff[2]<<" "<<coeff[3]<<" "<<coeff[4]<<endl;
}


void EngCoeff::setcoeff(float t,float n){
// E(vdw) = Aexp(-B*r^2)(1/(r+t)^2n-C/(r+t)^n)
// C=(2*n+2*B*(1+t))/[2*B*(1+t)^(n+1)+n*(1+t)^n]
// A=-exp(B)/[1/(1+t)^2n-C/(1+t)^n]

/*
// E(vdw) = Aexp(-B*r^2)(1/(r+t)-C/(r+t)^0.5)
// C=(1+2*B*(1+t))/[2*B*(1+t)^1.5+0.5*(1+t)^0.5]
// A=-exp(B)/[1/(1+t)-C/(1+t)^0.5]
*/
        int i,j,k;
        if(coeff==0) coeff=new float[100];
        for(i=0;i<100;i++) coeff[i]=0;

        float dmin=1000000;
        for(i=0;i<2000;i++) {
                 float b=0.001*(i+0.5);
                 float c=(2*n+2*b*(1+t))/(2*b*pow(1+t,1+n)+n*pow(1+t,n));
                 float a=-exp(b)/(1/pow(1+t,2*n)-c/pow(1+t,n));
                 float diff=0;
                 for(int f=0;f<size;f++) {
                        float x=dots[3*f];
                        float y=dots[3*f+1];
                        float w=dots[3*f+2];
                        float di=a*exp(-b*x*x)*(1/pow(x+t,2*n)-c/pow(x+t,n));
                        diff+=w*(di-y)*(di-y);
                }
                if(diff<dmin) {
                        dmin=diff;
                        coeff[0]=a;
                        coeff[1]=b;
                        coeff[2]=c;
                        coeff[3]=t;
                }
        }
	coeff[4]=n;
	
        cerr<<"the minimal difference:"<<sqrt(dmin/size)<<endl;
        cerr<<"the coefficient:"<<coeff[0]<<" "<<coeff[1]<<" "<<coeff[2]<<" "<<coeff[3]<<" "<<coeff[4]<<endl;
}


void EngCoeff::setcoeff(float t){
// E(vdw) = Aexp(-B*r^2)(1/(r+t)-C/(r+t)^0.5)
// C=(1+2*B*(1+t))/[2*B*(1+t)^1.5+0.5*(1+t)^0.5]
// A=-exp(B)/[1/(1+t)-C/(1+t)^0.5]

	int i,j,k;
	if(coeff==0) coeff=new float[100];
	for(i=0;i<100;i++) coeff[i]=0;

	float dmin=1000000;
	for(i=0;i<1000;i++) {
		 float b=0.002*(i+0.5);
		 float c=(1+2*b*(1+t))/(2*b*pow(1+t,1.5)+0.5*pow(1+t,0.5));
		 float a=-exp(b)/(1/(1+t)-c/sqrt(1+t));
		 float diff=0;
	 	 for(int f=0;f<size;f++) {
			float x=dots[3*f];
			float y=dots[3*f+1];
			float w=dots[3*f+2];
			float di=a*exp(-b*x*x)*(1/(x+t)-c/sqrt(x+t));
			diff+=w*(di-y)*(di-y);			
		}
		if(diff<dmin) {
			dmin=diff;
			coeff[0]=a;
			coeff[1]=b;
			coeff[2]=c;
			coeff[3]=t;
		} 		 
	}
	coeff[4]=0.5;
	cerr<<"the minimal difference:"<<sqrt(dmin/size)<<endl;
	cerr<<"the coefficient:"<<coeff[0]<<" "<<coeff[1]<<" "<<coeff[2]<<" "<<coeff[3]<<" "<<coeff[4]<<endl;
}


void EngCoeff::printout() {

	for(float r=0;r<10;r+=0.1) {
		float a=coeff[0];
		float b=coeff[1];
		float c=coeff[2];
		float t=coeff[3];
		float x=r;
		float n=coeff[4];
		float e=a*exp(-b*x*x)*(1/pow(x+t,2*n)-c/pow(x+t,n));
		float d=1/pow(x,12)-2*1/pow(x,6);		
		printf("%f %f %f\n",x,e,d);
	}
}

void EngCoeff::printnewout() {

        for(float r=0;r<10;r+=0.1) {
                float a=coeff[0];
                float b=coeff[1];
                float c=coeff[2];
                float t=coeff[3];
                float x=r;
                float n=coeff[4];
                float e=a*exp(-b*x*x)*(1/pow(x*x+t,2*n)-c/pow(x*x+t,n));
                float d=1/pow(x,12)-2*1/pow(x,6);
		int i=4*TRES.smt;
  		float ta=TRES.smooth[i+0];
  		float tb=TRES.smooth[i+1];
  		float tc=TRES.smooth[i+2];
  		float tn=TRES.smooth[i+3];
		float f=ta*exp(-tb*r*r)*(1/pow(r*r,2*tn)-tc/pow(r*r,tn));
                printf("new old exact:%f %f %f %f\n",x,e,f,d);
        }
}


void EngCoeff::printoutorder(float delt) {

	for(float r=0;r<10;r+=0.1) {
		float a=getvdwvalue(r);
		float a1=getvdwvalue(r+delt);
		float d1=getvdw1storder(r);
		float d2=getvdw2storder(r);
		float a2=a+d1*delt+d2/2*delt*delt;	
		printf("%f %f %f %f %f\n",r,a,a1,a2,a+d1*delt);
	}
}

void EngCoeff::printoutbound(float mid,float low,float high,float k) {
	float e=1;
	for(float r=0;r<20;r+=0.5) {
		//float b=e/2*exp(k*(high-mid));
		//float a=e/2*exp(k*(mid-low)); 
 		//float tt=a*exp(k*(low-r))+b*exp(k*(r-high)); 
		float tt=exp((low-r)*2)+exp(2*(r-high))-0.125*exp(-fabs(r-mid));
		//float tt=exp((low-r)*2)+exp(2*(r-high))-100*(r-mid)*(r-mid);
		/*
		float tt;
		if(r<low) tt=exp((low-r)*2);
		else if(r>high) tt=exp(2*(r-high));
		else tt=-(r-mid)*(r-mid);
		*/
                printf("%f %f\n",r,tt);
        }
}

float EngCoeff::getcrgvalue(float t,float r) {

	float e=1/(1+r);
	return e; 
}

float EngCoeff::getcrg1storder(float t,float r) {

        float e=-1/pow(1+r,2);
        return e; 
}

float EngCoeff::getcrg2storder(float t,float r) {

        float e=2/pow(1+r,3);
        return e; 
}

float EngCoeff::getcrgvalue(float r) {

        float e=1/(1+r);
        return e;
}

float EngCoeff::getcrg1storder(float r) {

        float e=-1/pow(1+r,2);
        return e;
}

float EngCoeff::getcrg2storder(float r) {

        float e=2/pow(1+r,3);
        return e;
}



float EngCoeff::getvdwvalue(float r) {
	float a=coeff[0];
	float b=coeff[1];
	float c=coeff[2];
	float t=coeff[3];
	float n=coeff[4];
	float x=r;
	float e=a*exp(-b*x*x)*(1/pow(x*x+t,2*n)-c/pow(x*x+t,n));
	return e;
}

float EngCoeff::getvdw1storder(float r) {
	float a=coeff[0];
	float b=coeff[1];
	float c=coeff[2];
	float t=coeff[3];
	float n=coeff[4];
	float rr=r*r;
	float re=1/(t+rr);
	float ree=pow(re,n);
	float e=a*exp(-b*rr)*(-b)*(ree*ree-c*ree)+
		a*exp(-b*rr)*(-2*n*ree*ree*re+n*c*ree*re);
	return e;
}

float EngCoeff::getvdw2storder(float r) {
	float a=coeff[0];
	float b=coeff[1];
	float c=coeff[2];
	float t=coeff[3];
	float e=a*exp(-b*r*r)*pow(-2*b*r,2)*(1/(r+t)-c/sqrt(r+t))+
		a*exp(-b*r*r)*(-2*b)*(1/(r+t)-c/pow(r+t,0.5))+
		2*a*exp(-b*r*r)*(-2*b*r)*(-1/pow(r+t,2)+0.5*c/pow(r+t,1.5))+
		a*exp(-b*r*r)*(2/pow(r+t,3)-0.75*c/pow(r+t,2.5));
	return e;
}
 
void EngCoeff::printoutcharge(float t){

	for(float r=0;r<10;r+=0.1) {
		printf("%f %f \n",r,r*r-6*r+1/(r+0.1));
		//printf("%f %f %f\n",r,1/(1+r),exp(-t*r));	
	}
}
