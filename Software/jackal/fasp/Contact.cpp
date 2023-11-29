#include"source.h"

Contact::Contact(int n)
{
atm=0;
next=0;
num=0;
flag=0;
distance=0;
energy=0;
total=0;
setatmarray(n);
}

Contact::Contact()
{
atm=0;
next=0;
num=0;
flag=0;
distance=0;
energy=0;
total=0;
}

Contact::~Contact()
{
 if(next) delete next; next=0;
 if(atm)  delete [] atm; atm=0;
 if(distance) delete [] distance;distance=0;
 if(energy) delete [] energy;energy=0;
}

void Contact::setatmarray(int n){

	atm=new Atm*[n];
	for(int i=0;i<n;i++) atm[i]=0;
}

void Contact::setengarray(int n){

        energy=new float[n];
        for(int i=0;i<n;i++) energy[i]=0;
}

