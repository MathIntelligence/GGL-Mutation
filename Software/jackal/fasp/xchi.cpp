#include<stdio.h>
#include<iostream.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
// calculate the range of chi angle
main(int argc,char *argv[])
{
char line[256];
FILE *fp;
int i,j,k,n,p;
float chi[100000][7];
int range[36];
float average[36];
fp=fopen(argv[1],"r");
p=0;
while(fgets(line,256,fp)!=NULL)
{
 j=strlen(line);
 i=8;
 if(j<10) continue;
 if(argv[2][0]!='0')
 if(line[0]!=argv[2][0]) continue; 
 n=(j+1-8)/8;
 for(k=0;k<n;k++)
 {
   chi[p][k]=atof(line+i);
   if(chi[p][k]>180)chi[p][k]-=360;
   i+=8;
 }
 p++; 
}
float xx;
if(argc>3) n=atoi(argv[3]);
if(argc>4) xx=atof(argv[4]);
cout<<n<<endl;
for(i=0;i<36;i++) {range[i]=0;average[i]=0;}
for(i=0;i<p;i++)
{
//if(!(chi[i][0]<-120&&chi[i][0]>-140))continue;
//if(!(chi[i][1]>120&&chi[i][1]<140))continue;
j=(chi[i][n]+180)/10.;
if(j==36) j=35;
if(argc>4)
{
if(n>0&&fabs(chi[i][n-1]-xx)<10) {range[j]++;average[j]+=chi[i][n];}
}
else {range[j]++;average[j]+=chi[i][n];}
}
for(i=0;i<36;i++)
{
k=range[i];
if(k==0) k=1;
printf("%5i %5i %6.1f\n",i*10-180,range[i],average[i]/k);
}
}
