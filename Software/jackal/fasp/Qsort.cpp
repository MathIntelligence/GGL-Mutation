#include"source.h"

Qsort::Qsort() {}

Qsort::~Qsort(){}

void Qsort::sort(double *area,int n,int *kseq)
{
int k,nx,i,j;
int m,l,ns; 
double x;
if(n==0) return;
for( k=0;k<n;k++) kseq[k]=k;

l=n;
k=n/2;
for(m=k;m>=1;m--)
{
   i=m;
   j=2*i;
   x=area[i-1];
   nx=kseq[i-1];

re10: 
   if(j<=l) 
   {
      if(j<l) 
      {
        if(area[j-1]<area[j]) j=j+1;
      }
      if(x<area[j-1]) 
      {
        area[i-1]=area[j-1];
        kseq[i-1]=kseq[j-1];
        i=j;
        j=2*i;
      }
      else
      {
        goto re20;
      }
      goto re10;
   }
re20: 
   area[i-1]=x;
   kseq[i-1]=nx;
}

for(ns=l;ns>=2;ns--)
{
   x=area[0];
   nx=kseq[0];

   area[0]=area[ns-1];
   kseq[0]=kseq[ns-1];
   area[ns-1]=x;
   kseq[ns-1]=nx;

   i=1;
   j=2;
   x=area[0];
   nx=kseq[0];

   k=ns-1;
   re50:
   if(j<=k) 
   {
     if(j<k) 
     {
       if(area[j-1]<area[j]) j=j+1;
     }
     if(x<area[j-1])
     {
       area[i-1]=area[j-1];
       kseq[i-1]=kseq[j-1];
       i=j;
       j=2*i;
     }
     else
     {
       goto re60;
     }
       goto re50;
   }
   re60:
   area[i-1]=x;
   kseq[i-1]=nx;
} 
}


void Qsort::sort(float *area,int n,int *kseq)
{
int k,nx,i,j;
int m,l,ns; 
float x;
if(n==0) return;
for( k=0;k<n;k++) kseq[k]=k;

l=n;
k=n/2;
for(m=k;m>=1;m--)
{
   i=m;
   j=2*i;
   x=area[i-1];
   nx=kseq[i-1];

re10: 
   if(j<=l) 
   {
      if(j<l) 
      {
        if(area[j-1]<area[j]) j=j+1;
      }
      if(x<area[j-1]) 
      {
        area[i-1]=area[j-1];
        kseq[i-1]=kseq[j-1];
        i=j;
        j=2*i;
      }
      else
      {
        goto re20;
      }
      goto re10;
   }
re20: 
   area[i-1]=x;
   kseq[i-1]=nx;
}

for(ns=l;ns>=2;ns--)
{
   x=area[0];
   nx=kseq[0];

   area[0]=area[ns-1];
   kseq[0]=kseq[ns-1];
   area[ns-1]=x;
   kseq[ns-1]=nx;

   i=1;
   j=2;
   x=area[0];
   nx=kseq[0];

   k=ns-1;
   re50:
   if(j<=k) 
   {
     if(j<k) 
     {
       if(area[j-1]<area[j]) j=j+1;
     }
     if(x<area[j-1])
     {
       area[i-1]=area[j-1];
       kseq[i-1]=kseq[j-1];
       i=j;
       j=2*i;
     }
     else
     {
       goto re60;
     }
       goto re50;
   }
   re60:
   area[i-1]=x;
   kseq[i-1]=nx;
} 
}

void Qsort::back(float *g,int n, int *order)
{
  float *g0;
  int i,j;
  g0=new float[n];
  for(i=0;i<n;i++)
  {
    g0[i]=g[i];
  }

  for(i=0;i<n;i++)
  {
    j=order[i];
    g[j]=g0[i];
  }

  delete [] g0;

}

int Qsort::slnpd(float **g, float eps, int n, int m)
{
int k,j,l,i1,ik,i,k1;
float bmax,t;

for(k=0;k<n;k++) 
{
    bmax=0;
    for(i=k;i<n;i++) 
    {
       if(bmax-fabs(g[i][k])<0) 
       {
         bmax=fabs(g[i][k]);
         l=i;
       }
    }
    if(bmax<1.e-5) return 0;
    if(l==k) goto re301;
    for(j=k;j<m;j++) 
    {
       t=g[l][j]; 
       g[l][j]=g[k][j]; 
       g[k][j]=t;       
    }

re301:

    t=1./g[k][k]; 
    k1=k+1;       
    for(j=k1;j<m;j++) 
    {
      g[k][j]=g[k][j]*t; 
      for(i=k1;i<n;i++) 
       g[i][j]=g[i][j]-g[i][k]*g[k][j];
    }
      
} 

for(ik=1;ik<n;ik++) 
{
    i=m-ik-2;
    i1=i+1;
    for(j=i1;j<n;j++)
    {
      g[i][m-1]=g[i][m-1]-g[i][j]*g[j][m-1];
    }
}
    return 1;
}


int Qsort::slnpd(float **g,int *order,int wide,int length) 
{
//find the maximum coefficent to cancel out. not tested!
int k,j,i,k1,k2,wide0;
float bmax,t;

//find the max coefficient
for(j=0;j<length-1;j++) order[j]=j;
wide0=wide;

for(k=0;k<wide;k++)
{
  bmax=0;
  for(i=k;i<wide;i++)
  {
    for(j=k;j<length-1;j++)
    {
      if(bmax-fabs(g[i][j])<0)
      {
        bmax=fabs(g[i][j]);
        k1=i; k2=j;
      }
    }
  }
  if(bmax<1.e-4) {wide0=k;break;};
  
//put k1 back
  for(j=k;j<length;j++)
  {
    if(k1==k) break;
    t=g[k1][j];
    g[k1][j]=g[k][j];
    g[k][j]=t;
  }
  
//put k2 back
  for(j=0;j<wide;j++)
  {
    if(k2==k) break;
    t=g[j][k2];
    g[j][k2]=g[j][k];
    g[j][k]=t;
  }
  
//unknown order
  j=order[k2];
  order[k2]=order[k];
  order[k]=j; 

  t=1./g[k][k];
  for(j=0;j<length;j++) g[k][j]=g[k][j]*t;
  
  for(i=k+1;i<wide;i++)
  {
    for(j=length-1;j>=k;j--)
    {
     g[i][j]=g[i][j]*g[k][k]-g[i][k]*g[k][j];
    }
  }
}

for(i=wide0-1;i>=0;i--)
{
 for(j=i-1;j>=0;j--)
 { 
   for(k=length-1;k>j;k--)
   {
     g[j][k]=g[j][k]-g[j][i]*g[i][k];
   }
 }
}

for(i=0;i<wide;i++)
{
  for(k=0;k<length;k++)
  {
    if(i>=wide0) g[i][k]=0;
    else if(k<i) g[i][k]=0;
    else if(k==i) g[i][k]=1;
    else g[i][k]=-g[i][k];
  }
}
return wide0;
}

double Qsort::matrix(float **a,int n)
{
  float **g0;
  double d;
  int i,j,k,m,h;
  d=0;
  if(n<1) 
  {
    return 0; 
  }  
  else if(n==1)
  {
    return a[0][0];
  }

  else if(n==2)
  {
    d=a[0][0]*a[1][1]-a[0][1]*a[1][0];
    return d;
  } 
  else if(n==3)
  {
    d=a[0][0]*a[1][1]*a[2][2]+
      a[0][1]*a[1][2]*a[2][0]+
      a[0][2]*a[2][1]*a[1][0]-
      a[0][2]*a[1][1]*a[2][0]-
      a[0][1]*a[1][0]*a[2][2]-
      a[0][0]*a[2][1]*a[1][2];
    return d;
  }
 
  else 
  {
    g0=new float*[n-1];
    for(i=0;i<n-1;i++) g0[i]=new float[n-1];
    d=0;j=-1;m=0;
    for(h=0;h<n;h++)
    {
      m=0;j=-j;
      for(i=0;i<n;i++)
      {
        if(i==h) continue;
        for(k=1;k<n;k++)
        {
          g0[m][k-1]=a[i][k];
        }
        m++;
      }
      d+=j*a[h][0]*matrix(g0,n-1);
    }
    for(i=0;i<n-1;i++) delete [] g0[i]; 
    delete [] g0;
    return d; 
  }
}

