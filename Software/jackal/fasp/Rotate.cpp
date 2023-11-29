#include"source.h"

Rotate::Rotate()
{
 flag=0;
 rota=0;
}

Rotate::~Rotate() {}

void Rotate::link(Res *r1,Res *r2,int n1,int flg)
{
//in this file, CA1 1,C1 2,O1 3,N2 4,CA2 5 
//if flg eq 1 rotate r2; if flg eq 0, rotate r1;   

   Atm *at;
   int i,k;
   float C1[3],CA1[3],O1[3];
   float xo[3],yo[3],zo[3];
   float xx,yx,zx;
   float d45,d46;

   if(r1==0||r2==0){cerr<<"warning: not linked"<<endl; return;} 
   if(r1->id==r2->id) {cerr<<"warning: not linked"<<endl;return;}
   at=(*r1)[" C  "];
   if(at==0) {cerr<<"header not exist! "<<endl;exit(0);}
   for(i=0;i<3;i++) C1[i]=at->xyz[i];
   at=(*r1)[" CA "];
   if(at==0) {cerr<<"header not exist! "<<endl;exit(0);}
   for(i=0;i<3;i++) CA1[i]=at->xyz[i];
   at=(*r1)[" O  "];
   if(at==0) {cerr<<"header not exist! "<<endl;exit(0);}
   for(i=0;i<3;i++) O1[i]=at->xyz[i];
   
   TRES.copy(r2->temp,  xo,3);
   TRES.copy(r2->temp+3,yo,3);
   TRES.copy(r2->temp+6,zo,3);

   Res *rr;
   Atm *a5;
   int j;
   float l[3][3],m[3][3],n[3][3];
   float x[3],y[3],z[3]; 
   float ax[2][3],ay[2][3],az[2][3];  
   float tmp[3]; 

   //calculate the rotation coefficeint for two systems!
   
   for(i=0;i<2;i++)
   {  
     if(i==0)
     {
       TRES.copy(CA1,x,3);
       TRES.copy(C1,y,3); 
       TRES.copy(O1,z,3);
     }
     else if(i==1)
     {
       TRES.copy(xo,x,3);
       TRES.copy(yo,y,3);
       TRES.copy(zo,z,3);
     }
   
     d45=TRES.distance(x,y);
     for(j=0;j<3;j++) ax[i][j]=(x[j]-y[j])/d45;     
     xx=TRES.distance(x,z);
     d46=TRES.distance(z,y);
     yx=(d45*d45+xx*xx-d46*d46)/(2*d45*xx);
     zx=d45/yx/xx; 
     for(j=0;j<3;j++) tmp[j]=x[j]+(z[j]-x[j])*zx;
     zx=TRES.distance(tmp,y);
     for(j=0;j<3;j++) ay[i][j]=(tmp[j]-y[j])/zx;
     az[i][0]=ax[i][1]*ay[i][2]-ax[i][2]*ay[i][1];
     az[i][1]=ax[i][2]*ay[i][0]-ax[i][0]*ay[i][2];
     az[i][2]=ax[i][0]*ay[i][1]-ax[i][1]*ay[i][0];
  
     j=i;
     if(flg==0)
     {
      if(i==0) j=1;
      else if(i==1) j=0;
     }
     
     l[j][0]=ax[i][0];m[j][0]=ax[i][1];n[j][0]=ax[i][2];
     l[j][1]=ay[i][0];m[j][1]=ay[i][1];n[j][1]=ay[i][2];
     l[j][2]=az[i][0];m[j][2]=az[i][1];n[j][2]=az[i][2];
   }
 
//rotate the residues;

   Res *r3;
    
   if(flg) r3=r2;
   else r3=(*r1->chn)[n1];
   
   if(flg==0)
   {
     for(i=0;i<3;i++) 
     {
       tmp[i]=yo[i];
       yo[i]=C1[i];
       C1[i]=tmp[i];
     }
   } 
   
   for(rr=r3;rr;rr=rr->next) 
   {
     if(flg) {if(rr->id0>n1) break;}   
     else    {if(rr->id0>r1->id0) break;}
     for(a5=rr->atm;a5;a5=a5->next)
     {
      for(i=0;i<3;i++) y[i]=a5->xyz[i]-yo[i];
      for(i=0;i<3;i++) //rotate the new axis    
      x[i]=l[1][i]*y[0]+m[1][i]*y[1]+n[1][i]*y[2];
    
      // rotate to back!   
      y[0]=l[0][0]*x[0]+l[0][1]*x[1]+l[0][2]*x[2];
      y[1]=m[0][0]*x[0]+m[0][1]*x[1]+m[0][2]*x[2];
      y[2]=n[0][0]*x[0]+n[0][1]*x[1]+n[0][2]*x[2];
    
      for(i=0;i<3;i++) y[i]+=C1[i];
      a5->transfer(y,1);
     
     } 

     if(rr->nemp==0) continue;

     for(k=0;k<3;k++)
     {
        for(i=0;i<3;i++) y[i]=rr->temp[k*3+i]-yo[i];
        for(i=0;i<3;i++) x[i]=l[1][i]*y[0]+m[1][i]*y[1]+n[1][i]*y[2];
    
	y[0]=l[0][0]*x[0]+l[0][1]*x[1]+l[0][2]*x[2];
        y[1]=m[0][0]*x[0]+m[0][1]*x[1]+m[0][2]*x[2];
        y[2]=n[0][0]*x[0]+n[0][1]*x[1]+n[0][2]*x[2];
         
	for(i=0;i<3;i++) y[i]+=C1[i];
	TRES.copy(y,rr->temp+k*3,3);
     }
   }
}



void Rotate::rotate(Atm *atm0, float angle)
{ 

//rotating within one residue; 
 
  Atm *atm1,*atm2,*atm3;
  int i,j;
  float r[3],r1[3];
  float bod,ca,sa;
  
  if(atm0==0) {cerr<<"no atom for rotating"<<endl;return;}
  if(atm0->tatm->rotate==0) return;
  
  atm1=atm0->bond[0]; 
   
  if(atm1==0) {cerr<<"no atom for rotating"<<endl;return;} 

 // bond length
  bod=TRES.distance(atm0->xyz,atm1->xyz);
  
  if(bod==0) 
  { 
    cerr<<"bond length is 0   "<<atm0->res->name;
    cerr<<atm0->res->id<<" "<<atm0->name<<endl; 
    return;
    //exit(0);
  }
   
  angle=angle/180*3.1415926;
  ca=cos(angle);sa=sin(angle);
   
  //normalization 
  for(i=0;i<3;i++) r[i]=(atm0->xyz[i]-atm1->xyz[i])/bod;
   
  for(atm3=atm0->res->atm;atm3;atm3=atm3->next)
  {
    if(atm3->name[1]=='H') break;

    if(flag==0&&atm3->tatm->id<atm0->tatm->id) continue;
    else if(flag==1&&rota[atm3->id0]==0) continue; //different scheme

    for(i=0;i<4;i++)
    {
       atm2=atm3->bond[i];
       if(atm2==0) continue;
       if(i>0&&atm2->name[1]!='H') continue;
       if(i==0)atm2=atm3;
       for(j=0;j<3;j++) r1[j]=atm2->xyz[j]-atm0->xyz[j];
       rotate(r,r1,ca,sa);
       for(j=0;j<3;j++) atm2->xyz[j]=r1[j]+atm0->xyz[j];
    }
  }

  return;

}

void Rotate::rotate(float *r,float *r1,float ca, float sa)
{
  int j;
  float r2[3],r3[3];
  float d;

  d=r1[0]*r[0]+r1[1]*r[1]+r1[2]*r[2];
  for(j=0;j<3;j++) r2[j]=d*r[j];
  for(j=0;j<3;j++) r3[j]=r1[j]-r2[j];
  r1[0]=r2[0]+ca*r3[0]+sa*(r[1]*r3[2]-r[2]*r3[1]);
  r1[1]=r2[1]+ca*r3[1]+sa*(r[2]*r3[0]-r[0]*r3[2]);
  r1[2]=r2[2]+ca*r3[2]+sa*(r[0]*r3[1]-r[1]*r3[0]);
}

void Rotate::rotate(Atm *atm0, float angle,int n,int flg)
{
  Atm *atm1,*atm2,*atm3;
  Res *rr,*rr1;
  int i,j;
  float r[3],r1[3];
  float bod,ca,sa;

  if(atm0->tatm->id!=1&&atm0->tatm->id!=2) 
  { 
   cerr<<"wrong subroutine to rotate the residue backbones"<<endl;
   return;
  } 

  if(atm0==0) {cerr<<"no atom for rotating"<<endl;return;}
  if(atm0->tatm->rotate==0) return;
  
  atm1=atm0->bond[0];
  if(atm1==0) {cerr<<"no rotating on the backbone performed"<<endl;return;}
 
  //  bod=atm0->bondlen[0];
 
  bod=TRES.distance(atm0->xyz,atm1->xyz);
  if(bod==0)
  {
    cerr<<"bond length is 0   "<<atm0->res->name;
    cerr<<atm0->res->id<<" "<<atm0->name<<endl;
    return;
    //exit(0);
  }
  
  angle=angle/180*3.1415926;
  ca=cos(angle);sa=sin(angle);
  for(i=0;i<3;i++) r[i]=(atm0->xyz[i]-atm1->xyz[i])/bod;

  if(flg==1) rr1=atm0->res;
  else rr1=(*atm0->res->chn)[n];
  if(rr1==0) {cerr<<" no existing of residue of id:"<< n<<endl;return;}

  for(rr=rr1;rr;rr=rr->next)
  {
   if(flg) {if(rr->id0>n) break;}
   else {if(rr->id0>atm0->res->id0) break;}
   //cerr<<rr->id0<<endl;
   for(atm3=rr->atm;atm3;atm3=atm3->next)
   {
    if(atm3->name[1]=='H') break;
    
    if(flag==1&&rota[atm3->id0]==0) continue; //using different scheme!
    else if(flag==0)
    {
      if(rr==atm0->res)
      {  
        //if(atm3->id==atm0->id) continue; 
        //else if(atm3->id==1) continue;
    
        if(flg&&atm0->tatm->id==1) 
        {
           if(atm3->tatm->id==0) continue; 
        }
        else if(flg&&atm0->tatm->id==2)
        {
           if(atm3->tatm->id!=3) continue;
        }
        else if(!flg&&atm0->tatm->id==1)
        {
           if(atm3->tatm->id!=0) continue;
        }   
        else if(!flg&&atm0->tatm->id==2)
        {
           if(atm3->tatm->id==3) continue;
        }
      }
    }
    for(i=0;i<4;i++)
    {
       if(i==0)atm2=atm3;
       else atm2=atm3->bond[i];
       if(atm2==0) continue;
       if(i>0&&atm2->name[1]!='H') continue;
       for(j=0;j<3;j++) r1[j]=atm2->xyz[j]-atm0->xyz[j];
       rotate(r,r1,ca,sa);
       for(j=0;j<3;j++) atm2->xyz[j]=r1[j]+atm0->xyz[j];
    }
   }
   
   for(i=0;i<rr->nemp/3;i++)
   {
       if(atm0->res==rr&&flg==1&&atm0->tatm->id==1) continue;	
       if(atm0->res==rr&&flg==1&&atm0->tatm->id==2) continue;
       if(atm0->res==rr&&atm0->tatm->id>=4) continue;
       for(j=0;j<3;j++) r1[j]=rr->temp[i*3+j]-atm0->xyz[j];
       rotate(r,r1,ca,sa);
       for(j=0;j<3;j++) rr->temp[i*3+j]=r1[j]+atm0->xyz[j];
   }
   
  }
}
void Rotate::rotateomega(Atm *atm0, float angle,int n,int flg)
{
  Atm *atm1,*atm2,*atm3;
  Res *rr,*rr1;
  int i,j;
  float r[3],r1[3];
  float bod,ca,sa;

  if(atm0->tatm->id!=0) 
  { 
   cerr<<"wrong subroutine to rotate the residue backbones"<<endl;
   return;
  } 

  if(atm0==0) {cerr<<"no atom for rotating"<<endl;return;}
  //if(atm0->tatm->rotate==0) return;
  
  atm1=atm0->bond[0];
  //if(atm1==0) {cerr<<"omega:no rotating on the backbone performed"<<endl;return;}
  if(atm1==0) return;
  //  bod=atm0->bondlen[0];
 
  bod=TRES.distance(atm0->xyz,atm1->xyz);
  if(bod==0)
  {
    cerr<<"bond length is 0   "<<atm0->res->name;
    cerr<<atm0->res->id<<" "<<atm0->name<<endl;
    return;
    //exit(0);
  }
  
  angle=angle/180*3.1415926;
  ca=cos(angle);sa=sin(angle);
  for(i=0;i<3;i++) r[i]=(atm0->xyz[i]-atm1->xyz[i])/bod;

  if(flg==1) rr1=atm0->res;
  else rr1=(*atm0->res->chn)[n];
  if(rr1==0) {cerr<<" no existing of residue of id:"<< n<<endl;return;}

  for(rr=rr1;rr;rr=rr->next)
  {
   if(flg) {if(rr->id0>n) break;}
   else {if(rr->id0>atm0->res->id0) break;}
   
   for(atm3=rr->atm;atm3;atm3=atm3->next)
   {
    if(atm3->name[1]=='H') break;
    
    if(flag==1&&rota[atm3->id0]==0) continue; //using different scheme!
    else if(flag==0)
    {
      if(rr==atm0->res)
      {   
        if(flg==0) 
        {
           continue;
        }         
      }
    }
    for(i=0;i<4;i++)
    {
       if(i==0)atm2=atm3;
       else atm2=atm3->bond[i];
       if(atm2==0) continue;
       if(i>0&&atm2->name[1]!='H') continue;
       for(j=0;j<3;j++) r1[j]=atm2->xyz[j]-atm0->xyz[j];
       rotate(r,r1,ca,sa);
       for(j=0;j<3;j++) atm2->xyz[j]=r1[j]+atm0->xyz[j];
    }
   }
   
   for(i=0;i<rr->nemp/3;i++)
   {
       if(atm0->res==rr&&flg==1) continue;
       for(j=0;j<3;j++) r1[j]=rr->temp[i*3+j]-atm0->xyz[j];
       rotate(r,r1,ca,sa);
       for(j=0;j<3;j++) rr->temp[i*3+j]=r1[j]+atm0->xyz[j];
   }
   
  }
}

void Rotate::rotate(Atm *atm0,float angle,int flg)
{
//rotating beyond this residue

  int n;
  Chn *chn;
  chn=atm0->res->chn;

  if(flg) n=100000;
  else       n=chn->res->id0;
  if(atm0->tatm->id>3||atm0->tatm->rotate==0)
  {cerr<<"no rotation!"<<endl;return;}
  if(atm0->tatm->rotate==0){ cerr<<"no rotation!"<<endl;return;}
  rotate(atm0,angle,n,flg);
}

void Rotate::hook(Res *s1,int n1)
{
Res *s2;
s2=new Res(s1->tres);
s2->configure();
hook(s2,s1,n1);
s2->next=0;s2->more=0;
delete s2;
}

void Rotate::hook(Res *s1)
{
Res *s2;
s2=new Res(s1->tres);
s2->configure();
hook(s2,s1);
s2->next=0;s2->more=0;
delete s2;
}


void Rotate::hook(Res *s2,Res *s1,int n1)
{
Atm *a3,*a4,*a5;
Atm *a1,*a2,*b;
float d,dis[2];
float c[2][2],s[2][2],p[2][2];
float l[2][3],m[2][3],n[2][3];
float x[3],y[3]; 
int i,j;
  
//if(s1->regular==0) return;
a2=(*s1)[n1];
if(a2==0) 
{
  cerr<<"warning in hooking:"<<s1->name<<s1->id0<<endl;
  return;
}
a1=a2->bond[0];
if(a1==0) 
{
  cerr<<"warning in hooking:"<<s1->name<<s1->id0<<endl;
  return;
}
a4=(*s2)[n1];
if(a4==0)
{
  cerr<<"warning in hooking:"<<s2->name<<s2->id0<<endl;
  return;
}
a3=a4->bond[0];   
if(a3==0) return;

//calculate the rotation coefficeint for two systems!

for(i=0;i<2;i++)
{  
  if(i==0)
  {
    TRES.copy(a2->xyz,y,3);
    TRES.copy(a1->xyz,x,3); 
  }
  else if(i==1)
  {
    TRES.copy(a4->xyz,y,3);
    TRES.copy(a3->xyz,x,3);
  }
  if(a1->res->id!=a2->res->id) 
  {
    cerr<<"not the same residue no hook performed"<<endl;
    return;
  }

  dis[i]=TRES.distance(y,x);

  if(dis[i]==0) 
  {
     cerr<<"distance is zero between:";
     cerr<<a1->name<<a1->id0<<"**"<<a2->name<<a2->id0<<endl;
     return;
  }
  c[i][0]=(y[2]-x[2])/dis[i];
  c[i][1]=sqrt(1-c[i][0]*c[i][0]);
  d=c[i][1]*dis[i];
  if(d!=0)
  { 
   s[i][0]=-(y[1]-x[1])/d;
   s[i][1]= (y[0]-x[0])/d;
  }
  else
  {
   s[i][0]=1;
   s[i][1]=0;
  }
  p[i][0]=1;
  p[i][1]=0; 
  l[i][0]=s[i][0]*p[i][0]-c[i][0]*s[i][1]*p[i][1];
  l[i][1]=-s[i][0]*p[i][1]-c[i][0]*s[i][1]*p[i][0];
  l[i][2]=c[i][1]*s[i][1];
  m[i][0]=s[i][1]*p[i][0]+c[i][0]*s[i][0]*p[i][1];
  m[i][1]=-s[i][1]*p[i][1]+c[i][0]*s[i][0]*p[i][0];
  m[i][2]=-c[i][1]*s[i][0];
  n[i][0]=c[i][1]*p[i][1];
  n[i][1]=c[i][1]*p[i][0];
  n[i][2]=c[i][0];
}
  d=dis[1]-dis[0];

  for(a5=s2->atm;a5;a5=a5->next)
  {
    if(a5->name[1]=='H') break;
    if(flag==0&&a5->tatm->id<n1) continue;
    else if(flag==1&&rota[a5->id0]==0) continue;//different scheme   

    for(j=0;j<a5->tatm->nbond;j++)
    {
      a4=a5->bond[j];
      if(a4==0) continue;   
      if(j>0&&a4->name[1]!='H') continue;
      if(j==0) a4=a5;
      for(i=0;i<3;i++) x[i]=a4->xyz[i]-a3->xyz[i];
      for(i=0;i<3;i++) //rotate the new axis    
      y[i]=l[1][i]*x[0]+m[1][i]*x[1]+n[1][i]*x[2];
      y[2]=y[2]-d;
      // rotate to back!   
      x[0]=l[0][0]*y[0]+l[0][1]*y[1]+l[0][2]*y[2];
      x[1]=m[0][0]*y[0]+m[0][1]*y[1]+m[0][2]*y[2];
      x[2]=n[0][0]*y[0]+n[0][1]*y[1]+n[0][2]*y[2];
      for(i=0;i<3;i++) x[i]+=a1->xyz[i];
      b=(*s1)[a4->tatm->id];
      if(b==0) continue;
      TRES.copy(x,b->xyz,3);
    }
  }
}
 

void Rotate::hook(Res *from,Res *to)
{
//int this file,CA1,C1,N1 and CA2,C2,N2

   Atm *CA1,*C1,*N1,*CA2,*C2,*N2;
   int i;
   float xo[3],yo[3],zo[3];
   float xx,yx,zx;
   float d45,d46;
   
   N1=from->atm;
   CA1=N1->next;
   C1=CA1->next;

   N2=to->atm;
   CA2=N2->next;
   C2=CA2->next;

   if(N2->tatm->id!=0||CA2->tatm->id!=1||C2->tatm->id!=2)
   {
    cerr<<"not standard residue in hook! "<<to->name<<to->id<<endl;
   }

   TRES.copy(N1->xyz ,xo,3);
   TRES.copy(CA1->xyz,yo,3);
   TRES.copy(C1->xyz ,zo,3);
   
   Atm *a5,*a6;
   int j;
   float l[3][3],m[3][3],n[3][3];
   float x[3],y[3],z[3]; 
   float ax[2][3],ay[2][3],az[2][3];  
   float tmp[3];
   //calculate the rotation coefficeint for two systems!
   for(i=0;i<2;i++)
   {  
     if(i==0)
     {
       TRES.copy(N2->xyz, x,3);
       TRES.copy(CA2->xyz,y,3); 
       TRES.copy(C2->xyz, z,3);
     }
     else if(i==1)
     {
       TRES.copy(xo,x,3);
       TRES.copy(yo,y,3);
       TRES.copy(zo,z,3);
     }
     d45=TRES.distance(x,y);
     for(j=0;j<3;j++) ax[i][j]=(x[j]-y[j])/d45;     
     xx=TRES.distance(x,z);
     d46=TRES.distance(z,y);
     yx=(d45*d45+xx*xx-d46*d46)/(2*d45*xx);
     zx=d45/yx/xx; 
     for(j=0;j<3;j++) tmp[j]=x[j]+(z[j]-x[j])*zx;
     zx=TRES.distance(tmp,y);
     for(j=0;j<3;j++) ay[i][j]=(tmp[j]-y[j])/zx;
     az[i][0]=ax[i][1]*ay[i][2]-ax[i][2]*ay[i][1];
     az[i][1]=ax[i][2]*ay[i][0]-ax[i][0]*ay[i][2];
     az[i][2]=ax[i][0]*ay[i][1]-ax[i][1]*ay[i][0];
  
     l[i][0]=ax[i][0];m[i][0]=ax[i][1];n[i][0]=ax[i][2];
     l[i][1]=ay[i][0];m[i][1]=ay[i][1];n[i][1]=ay[i][2];
     l[i][2]=az[i][0];m[i][2]=az[i][1];n[i][2]=az[i][2];
   }

   
   for(a5=C1->next->next;a5;a5=a5->next)
   {
      if(a5->tatm->name[1]=='H')
      {
        if(a5->tatm->bond[0]->id<4) continue;
      } 
      for(i=0;i<3;i++) y[i]=a5->xyz[i]-yo[i];
      for(i=0;i<3;i++) //rotate the new axis    
      x[i]=l[1][i]*y[0]+m[1][i]*y[1]+n[1][i]*y[2];
    
      // rotate to back!
      y[0]=l[0][0]*x[0]+l[0][1]*x[1]+l[0][2]*x[2];
      y[1]=m[0][0]*x[0]+m[0][1]*x[1]+m[0][2]*x[2];
      y[2]=n[0][0]*x[0]+n[0][1]*x[1]+n[0][2]*x[2];
      for(i=0;i<3;i++) y[i]+=CA2->xyz[i];
      a6=(*to)[a5->tatm->id];
      a6->transfer(y,1);
   }
 }


void Rotate::hookheader(Res *to)
{
//int this file,CA1,C1,N1 and CA2,C2,N2

   Atm *CA1,*C1,*N1,*CA2,*C2,*N2;
   int i;
   float xo[3],yo[3],zo[3];
   float xx,yx,zx;
   float d45,d46;
   Res *from = new Res(to->tres);
   from->nemp=9;
   from->temp=new float[9];
   for(i=0;i<9;i++) from->temp[i]=to->tres->head[i];
    
   N1=from->atm;
   if(N1==0) return;
   CA1=N1->next;
   if(CA1==0) return;
   C1=CA1->next;
   if(C1==0) return;

   N2=to->atm;
   if(N2==0) return;
   CA2=N2->next;
   if(CA2==0) return;
   C2=CA2->next;
   if(C2==0) return; 

   if(N2->tatm->id!=0||CA2->tatm->id!=1||C2->tatm->id!=2)
   {
    cerr<<"not standard residue in hook! "<<to->name<<to->id<<endl;
    return;
    //exit(0);
   }

   TRES.copy(N1->xyz ,xo,3);
   TRES.copy(CA1->xyz,yo,3);
   TRES.copy(C1->xyz ,zo,3);
   
   //Atm *a5,*a6;
   int j;
   float l[3][3],m[3][3],n[3][3];
   float x[3],y[3],z[3]; 
   float ax[2][3],ay[2][3],az[2][3];  
   float tmp[3];
   //calculate the rotation coefficeint for two systems!
   for(i=0;i<2;i++)
   {  
     if(i==0)
     {
       TRES.copy(N2->xyz, x,3);
       TRES.copy(CA2->xyz,y,3); 
       TRES.copy(C2->xyz, z,3);
     }
     else if(i==1)
     {
       TRES.copy(xo,x,3);
       TRES.copy(yo,y,3);
       TRES.copy(zo,z,3);
     }
     d45=TRES.distance(x,y);
     for(j=0;j<3;j++) ax[i][j]=(x[j]-y[j])/d45;     
     xx=TRES.distance(x,z);
     d46=TRES.distance(z,y);
     yx=(d45*d45+xx*xx-d46*d46)/(2*d45*xx);
     zx=d45/yx/xx; 
     for(j=0;j<3;j++) tmp[j]=x[j]+(z[j]-x[j])*zx;
     zx=TRES.distance(tmp,y);
     for(j=0;j<3;j++) ay[i][j]=(tmp[j]-y[j])/zx;
     az[i][0]=ax[i][1]*ay[i][2]-ax[i][2]*ay[i][1];
     az[i][1]=ax[i][2]*ay[i][0]-ax[i][0]*ay[i][2];
     az[i][2]=ax[i][0]*ay[i][1]-ax[i][1]*ay[i][0];
  
     l[i][0]=ax[i][0];m[i][0]=ax[i][1];n[i][0]=ax[i][2];
     l[i][1]=ay[i][0];m[i][1]=ay[i][1];n[i][1]=ay[i][2];
     l[i][2]=az[i][0];m[i][2]=az[i][1];n[i][2]=az[i][2];
   }

   if(to->temp) delete [] to->temp;
   to->nemp=9;
   to->temp=new float[9];


   for(int k=0;k<3;k++) {   
    
	for(i=0;i<3;i++) y[i]=from->temp[k*3+i]-yo[i];
	//rotate the new axis 
        for(i=0;i<3;i++)x[i]=l[1][i]*y[0]+m[1][i]*y[1]+n[1][i]*y[2];
        // rotate to back!
      	y[0]=l[0][0]*x[0]+l[0][1]*x[1]+l[0][2]*x[2];
      	y[1]=m[0][0]*x[0]+m[0][1]*x[1]+m[0][2]*x[2];
      	y[2]=n[0][0]*x[0]+n[0][1]*x[1]+n[0][2]*x[2];
	for(i=0;i<3;i++) y[i]+=CA2->xyz[i];
	for(i=0;i<3;i++) to->temp[k*3+i]=y[i];
   }
   delete from;
 }

