#include"source.h"
Eng::Eng(Tatm *s)
{
tatm=s;
charge=0;
radius=0;
epslon=0;
bond0=0;
kbond=0;
angle0=0;
kangle=0;
chi0=0;
kchi=0;
nchi=0;
mbond=0;
mangle=0;
mchi=0;
chi02=0;
kchi2=0;
nchi2=0;
dipole=0;
charge0=0;
}

void Eng::tank()
{
int i;
float x;
x=999999;
if(TRES.flag/100==1||TRES.flag/100==2||TRES.flag/100==3||TRES.flag/100==4)
{
  mbond=tatm->nbond;
  bond0=new float[mbond];
  kbond=new float[mbond];
  for(i=0;i<mbond;i++)
  {
   bond0[i]=x;
   kbond[i]=x; 
  }
  if(tatm->nbond==1) return;
  if(tatm->nbond==2) mangle=1;
  else if(tatm->nbond==3) mangle=3;
  else if(tatm->nbond==4) mangle=6;

  angle0=new float[mangle];
  kangle=new float[mangle];
  for(i=0;i<mangle;i++) 
  {
    angle0[i]=x;
    kangle[i]=x;
  }

  mchi=(tatm->nbond-1)*(tatm->bond[0]->nbond-1);
  if(mchi>0)
  {
    chi0=new float[mchi];
    kchi=new float[mchi];
    nchi=new float[mchi];
    chi02=new float[mchi];
    kchi2=new float[mchi];
    nchi2=new float[mchi];
    for(i=0;i<mchi;i++)
    {
      chi0[i]=x;
      kchi[i]=x;
      nchi[i]=x;
      chi02[i]=x;
      kchi2[i]=x;
      nchi2[i]=x;

    }
  }
}
//else 
}

Eng::~Eng()
{
if(bond0)   delete [] bond0;
if(kbond)   delete [] kbond;
if(angle0)  delete [] angle0;
if(kangle)  delete [] kangle;
if(kchi)    delete [] kchi;
if(nchi)    delete [] nchi;
if(chi0)    delete [] chi0;
if(kchi2)   delete [] kchi2;
if(nchi2)   delete [] nchi2;
if(chi02)   delete [] chi02;
bond0=0;
kbond=0;
angle0=0;
kangle=0;
chi0=0;
kchi=0;
nchi=0;
chi02=0;
kchi2=0;
nchi2=0;
}

void Eng::charm(char *filnam0,char *filnam1)
{

  FILE *fp;
  int i,j,m,n,neq,bg,ed,mt;
  char *mlib,line[200],line0[200],*buf;
  char equal[5][200];
  Tres *r;
  Tatm *a,*a0,*a1;
  float d;

//read parameter from charm 

  mlib=RCS["library"];
  //if((TRES.flag%100)/10==1) goto remy;
   
  //read charge 
 
  fp=0;i=0;
  sprintf(line,"%s/%s",mlib,filnam0);
  i=strlen(mlib);
  while((fp=fopen(line,"r"))==NULL)
  {
     mlib=mlib+i+1;
     i=strlen(mlib);
     if(*mlib=='\0')
     {cerr<<"warning:could not open: "<<filnam0<<endl;return;}
     sprintf(line,"%s/%s",mlib,filnam0);
  }
  r=0;
  while(fgets(line,200,fp)!=NULL)
  { 
      buf=line;
      while(*buf==' ')buf++;
      if(*buf=='!'||*buf=='\n') continue;
      
      if(strncmp(line,"RESI",4)==0) //new residue
      {
        r=TRES[line+5];
      } 
      if(r==0) continue; //not standard residue
      if(strncmp(line,"ATOM",4)==0) //new atom
      {
        if(line[5]!='H')
        {
          a=r->isatm(line+4); 
          a->eng->charge=atof(line+17);
          strncpy(a->type,line+10,4);
          a->type[4]='\0';
          i=0; //non hydrogen
          continue;
        }
        if(a==0) continue; 
        if(line[5]=='H')
        {
          m=0;
          if(line[6]==' ')line[6]='N';
          for(j=0;j<a->nbond;j++)  
          {           
           if(a->bond[j]->name[1]!='H')continue;   
           if(m==i) 
           {
               a->bond[j]->eng->charge=atof(line+17);
               strncpy(a->bond[j]->type,line+10,4);
               a->bond[j]->type[4]='\0';
               i++;
               break;
           }
           m++;
         }
           if(i==0) a->eng->charge+=atof(line+17);//H not in top file
        } 
      }//if atom 
  }//while
 
 
  fclose(fp);
  //read parameters of charmm
//remy:

  fp=0;i=0;
  sprintf(line,"%s/%s",mlib,filnam1);
  i=strlen(mlib);
  while((fp=fopen(line,"r"))==NULL)
  {
      mlib=mlib+i+1;
      i=strlen(mlib);
      if(*mlib=='\0')
      {cerr<<"warning:could not open: "<<filnam1<<endl;return;}
      sprintf(line,"%s/%s",mlib,filnam1);
   }
   m=0; 
   neq=0;
   while(fgets(line,200,fp)!=NULL)
   {
      buf=line;
      while(*buf==' ')buf++;
      if(*buf=='!'||*buf=='\n') continue;
      strcpy(line,buf);
      if(strncmp(line,"BONDS",5)==0) {m=1;continue;}
      if(strncmp(line,"ANGLES",6)==0){m=2;continue;}
      if(strncmp(line,"DIHEDRALS",9)==0) {m=3;continue;}
      if(strncmp(line,"IMPROPER",8)==0) {m=4;continue;}
      if(strncmp(line,"NONBONDED",9)==0) {m=5;continue;}
      if(strncmp(line,"EQUIVALENCE",11)==0) {m=-1;neq=0;continue;}
      if(strncmp(line,"END",3)==0) {m=0;continue;}
      if(m==-1)
      {
        strcpy(equal[neq++],line);
      }
      if(m==0) continue;
  
      if(m==1)
      {
        //if(strstr(line,"NH1  CT1")) cout<<line<<endl;
        for(r=&TRES;r;r=r->next)
        {
          for(a=r->tatm;a;a=a->next)
          {
            mt=0;
            if(strncmp(line,a->type,4)==0)mt=1;
            if(mt==0) if(strncmp(line+5,a->type,4)==0)mt=2;
            if(mt==0) continue;
            if(strncmp(line,line+5,4)==0) mt=3;           
            for(i=0;i<a->nbond;i++)
            { 
              if(a->eng->bond0[i]<99999)continue;
              if(mt!=2)
              {               
                if(strncmp(line+5,a->bond[i]->type,4)==0)
                {
                  sscanf(line+9,"%f%f",a->eng->kbond+i,a->eng->bond0+i);
                  continue;
                }
               }
               if(mt!=1)
               {                
                if(strncmp(line,a->bond[i]->type,4)==0)
                {
                  sscanf(line+9,"%f%f",a->eng->kbond+i,a->eng->bond0+i);
                  continue;
                }
               }            
            }//i  
          }//a
        }//r
      }//m
      if(m==2)
      {
        //if(strstr(line,"CH-C2-SH")) cerr<<line<<endl;
        for(r=&TRES;r;r=r->next)
        {
          for(a=r->tatm;a;a=a->next)
          {
            if(a->eng->mangle==0) continue;
            if(strncmp(line+5,a->type,4)!=0) continue;
            n=-1;
            for(i=0;i<a->nbond;i++)
            for(j=i+1;j<a->nbond;j++)
            {
                n++;
                if(a->eng->angle0[n]<99999) continue;
                if(n>=a->eng->mangle) continue;
                a0=a->bond[i];
                a1=a->bond[j];
              
                sprintf(line0,"%s %s %s",a0->type,a->type,a1->type);         

                if(strncmp(line,line0,14)==0)
                {
                 sscanf(line+14,"%f%f",a->eng->kangle+n,a->eng->angle0+n);
                 //if(r->name=='A'&&a->id==0) cout<<line0<<endl;
                 continue;
                }
                sprintf(line0,"%s %s %s",a1->type,a->type,a0->type); 
              //if(r->name=='A')
              //if(strstr(line0,"CT1  NH1  C"))
              //cout<<line0;
                if(strncmp(line,line0,14)==0)
                {
                 sscanf(line+14,"%f%f",a->eng->kangle+n,a->eng->angle0+n);
                 //if(r->name=='A'&&a->id==0) cerr<<line0<<endl;
                 continue;
                }
            }//i
          }//a     
        }//r
      }//m=2
      if(m==3)
      {
        //if(strstr(line,"X    CH1E CH1E X")) cerr<<line<<endl;
        if(line[0]=='X'&&line[15]!='X') { bg=5;ed=14;}
        else if(line[0]=='X'&&line[15]=='X') { bg=5;ed=9;}
        else if(line[0]!='X'&&line[15]=='X') { bg=0;ed=14;}
        else { bg=0;ed=19; }
        for(r=&TRES;r;r=r->next)
        {
          for(a0=r->tatm;a0;a0=a0->next)
          {           
            if(a0->eng->mchi==0) continue;
            a1=a0->bond[0];
            mt=0;
            if(strncmp(line+5, a0->type,4)==0&&
               strncmp(line+10,a1->type,4)==0)mt=1;
            if(mt==0)
            {
             if(strncmp(line+5, a1->type,4)==0&&
                strncmp(line+10,a0->type,4)==0)mt=2;
            }
            if(mt==0) continue;
            if(strncmp(a1->type,a0->type,4)==0)mt=3;
            n=-1;
            for(j=0;j<a1->nbond;j++)
            {
              if(a1->bond[j]==a0) continue;
              for(i=1;i<a0->nbond;i++)
              {
                n++;
                if(n>=a0->eng->mchi) continue;
                if(a0->eng->nchi2[n]<9999) continue;
                if(mt!=2)
                {
                   sprintf(line0,"%s %s %s %s",a0->bond[i]->type,a0->type,
                                               a1->type,a1->bond[j]->type);
                   if(strncmp(line+bg,line0+bg,ed)==0)
                   {
                     if(a0->eng->nchi[n]>9999) 
                       sscanf(line+19,"%f%f%f",a0->eng->kchi+n,
                       a0->eng->nchi+n,a0->eng->chi0+n);
                     else 
                       sscanf(line+19,"%f%f%f",a0->eng->kchi2+n,
                       a0->eng->nchi2+n,a0->eng->chi02+n);
                     continue;
                   }
                }
                if(mt!=1)
                { 
                    sprintf(line0,"%s %s %s %s",a1->bond[j]->type,a1->type,
                                                a0->type,a0->bond[i]->type);
                    if(strncmp(line+bg,line0+bg,ed)==0)
                    {
                      if(a0->eng->nchi[n]>9999) 
                       sscanf(line+19,"%f%f%f",a0->eng->kchi+n,
                       a0->eng->nchi+n,a0->eng->chi0+n);
                      else 
                       sscanf(line+19,"%f%f%f",a0->eng->kchi2+n,
                       a0->eng->nchi2+n,a0->eng->chi02+n);
                      continue;
                    }
                }                 
              }//i
            }//j
          }//a0
        }//r
      }//m=3
      if(m==5)
      { 
        for(r=&TRES;r;r=r->next) 
        {
          for(a=r->tatm;a;a=a->next)
          {
            if(a->eng->radius>0) continue;
            if(strncmp(a->type,line,4)==0)
            {
              sscanf(line+16,"%f%f",&a->eng->epslon,&a->eng->radius);
              continue;
            }
            for(i=0;i<neq;i++)
            {
              if(strncmp(line,equal[i],2)==0)
              {
                if(strstr(equal[i],a->type))
                {
                   sscanf(line+6,"%f%f%f",&d,&a->eng->epslon,&a->eng->radius);
                   break;
                }
              }
            } //i          
          }//a
        }//r
      }//m=5
    }//while
}


void Eng::amber(char *filnam0,char *filnam1)
{

FILE *fp;
int i,j,m,n,h,neq,bg,ed,mt;
char *mlib,line[200],line0[200],*buf;
Tres *r;
Tatm *a,*a0,*a1;
char equal[5][200];
//float d;
//read parameter from amber

    //read charge 
 
     
    fp=0;i=0;
    mlib=RCS["library"];
    if((TRES.flag%100)/10==1) goto remy;
    sprintf(line,"%s/%s",mlib,filnam0);
    i=strlen(mlib);
    while((fp=fopen(line,"r"))==NULL)
    {
      mlib=mlib+i+1;
      i=strlen(mlib);
      if(*mlib=='\0')
      {cerr<<"warning:could not open: "<<filnam0<<endl;return;}
      sprintf(line,"%s/%s",mlib,filnam0);
    }
    r=0;
    while(fgets(line,200,fp)!=NULL)
    { 
      buf=line;
      while(*buf==' ')buf++;
      if(*buf=='!'||*buf=='\n') continue;
      //strcpy(line,buf); 
      if(strncmp(line+6,"INT     1",9)==0) //new residue
      {
        r=TRES[line+1];
      } 
      if(strncmp(buf,"DONE",4)==0) r=0;
      if(atoi(buf)<4) continue;
      if(r==0) continue; //not standard residue
      if(line[3]>='0'&&line[4]<='9'&&strlen(line)>70) //new atom
      {
        if(strncmp(line+6,"DUMM",4)==0) continue;
        if(line[6]!='H')
        {
          a=r->isatm(line+5); 
          a->eng->charge=atof(line+63);
          strncpy(a->type,line+12,2);
          a->type[2]='\0';
          i=0; //non hydrogen
          continue;
        }
        if(a==0) continue; 
        if(line[6]=='H')
        { 
          m=0; 
          for(j=0;j<a->nbond;j++)  
          {           
           if(a->bond[j]->name[1]!='H')continue;   
           if(m==i) 
           {
               a->bond[j]->eng->charge=atof(line+63);
               strncpy(a->bond[j]->type,line+12,2);
               a->bond[j]->type[2]='\0';
               i++;
               break;
           }
           m++;
          } 
          if(i==0)a->eng->charge+=atof(line+63);//H not in top file
        } 
      } 
    }
  
    fclose(fp);
  //read parameters of charmm
remy:

    //if(flag%100==0) strcpy(filnam0,"amber_parm96.dat\0");
    //if(flag%100==1||flag%100==2) strcpy(filnam0,"amber_parm91.dat\0");
    fp=0;i=0;
    sprintf(line,"%s/%s",mlib,filnam1);
    i=strlen(mlib);
    while((fp=fopen(line,"r"))==NULL)
    {
      mlib=mlib+i+1;
      i=strlen(mlib);
      if(*mlib=='\0')
      {cerr<<"warning:could not open: "<<filnam1<<endl;return;}
      sprintf(line,"%s/%s",mlib,filnam1);
    }
    m=0; 
    while(fgets(line,200,fp)!=NULL)
    {
      buf=line;
      while(*buf==' ') buf++;
      if(*buf=='!'||*buf=='\n') continue;
      strcpy(line,buf);
      if(strncmp(line,"BONDS",5)==0) {m=1;continue;}
      if(strncmp(line,"ANGLES",6)==0){m=2;continue;}
      if(strncmp(line,"DIHEDRALS",9)==0) {m=3;continue;}
      if(strncmp(line,"IMPROPER",8)==0) {m=4;continue;}
      if(strncmp(line,"NONBONDED",9)==0) {m=5;continue;}
      if(strncmp(line,"EQUIVALENCE",11)==0) {m=-1;neq=0;continue;}
      if(strncmp(line,"END",3)==0) {m=0;continue;}
      if(m==-1)
      {
        strcpy(equal[neq++],line);
      }
      if(m==0) continue;
      if(m==1)
      {
        //if(strstr(line,"CT-H1")) cout<<line<<endl;
        for(r=&TRES;r;r=r->next)
        {
          for(a=r->tatm;a;a=a->next)
          {
            mt=0;
            if(strncmp(line,a->type,2)==0)mt=1;
            if(mt==0) if(strncmp(line+3,a->type,2)==0)mt=2;
            if(mt==0) continue;
            if(strncmp(line,line+3,2)==0) mt=3;           
            for(i=0;i<a->nbond;i++)
            { 
              if(a->eng->bond0[i]<99999)continue;
              if(mt!=2)
              {               
                if(strncmp(line+3,a->bond[i]->type,2)==0)
                {
                  sscanf(line+6,"%f%f",a->eng->kbond+i,a->eng->bond0+i);
                  continue;
                }
               }
               if(mt!=1)
               {                
                if(strncmp(line,a->bond[i]->type,2)==0)
                {
                  sscanf(line+6,"%f%f",a->eng->kbond+i,a->eng->bond0+i);
                  continue;
                }
               }
            }  
          }
        }
      }
      if(m==2)
      {
        //if(strstr(line,"CH-C2-SH")) cout<<line<<endl;
        for(r=&TRES;r;r=r->next)
        {
          for(a=r->tatm;a;a=a->next)
          {
            if(a->eng->mangle==0) continue;
            if(strncmp(line+3,a->type,2)!=0) continue;
            n=-1;
            for(i=0;i<a->nbond;i++)
             for(j=i+1;j<a->nbond;j++)
             { 
              n++;
              if(a->eng->angle0[n]<99999)continue;
              if(n>=a->eng->mangle) continue;
              a0=a->bond[i];
              a1=a->bond[j];
              sprintf(line0,"%s-%s-%s",a0->type,a->type,a1->type);         
              if(strncmp(line,line0,8)==0)
              {
                sscanf(line+8,"%f%f",a->eng->kangle+n,a->eng->angle0+n);
                continue;
              }
              sprintf(line0,"%s-%s-%s",a1->type,a->type,a0->type); 
              if(strncmp(line,line0,8)==0)
              {
                sscanf(line+8,"%f%f",a->eng->kangle+n,a->eng->angle0+n);
                continue;
              }
            }
          }//a     
        }//r
      }//m=2
      if(m==3)
      {
        //sscanf(line+12,"%f%f%f%f",&d,&d,&d,&h);
        //if(h<0) continue;
        //if(strstr(line,"X -C -N -X")) cout<<endl;
        if(line[0]=='X'&&line[9]!='X') { bg=3;ed=8;}
        else if(line[0]=='X'&&line[9]=='X') { bg=3;ed=5;}
        else if(line[0]!='X'&&line[9]=='X') { bg=0;ed=8;}
        else { bg=0;ed=11; }
        for(r=&TRES;r;r=r->next)
        {
          for(a0=r->tatm;a0;a0=a0->next)
          { 
            if(a0->eng->mchi==0) continue;
            a1=a0->bond[0];
            mt=0;
            if(strncmp(line+3, a0->type,2)==0&&
               strncmp(line+6,a1->type,2)==0)mt=1;              
            if(mt==0)
            {               
               if(strncmp(line+3, a1->type,2)==0&&
                  strncmp(line+6,a0->type,2)==0)mt=2;              
            }
            if(mt==0) continue;
            if(strncmp(a1->type,a0->type,2)==0)mt=3;
            n=-1; 
            for(j=0;j<a1->nbond;j++)
            {
              if(a1->bond[j]==a0) continue;
              for(i=1;i<a0->nbond;i++)
              {
                n++;
                if(n>=a0->eng->mchi) continue;
                if(a0->eng->nchi2[n]<99999) continue;
                if(mt!=2)
                {                   
                    sprintf(line0,"%s-%s-%s-%s",a0->bond[i]->type,a0->type,
                            a1->type,a1->bond[j]->type);
                    if(strncmp(line+bg,line0+bg,ed)==0)
                    {
                       if(a0->eng->chi0[n]>9999)
                       {
                           sscanf(line+12,"%i%f%f%f",&h,a0->eng->kchi+n,
                           a0->eng->chi0+n,a0->eng->nchi+n);
                           a0->eng->kchi[n]=a0->eng->kchi[n]/h;
                           continue;
                       }
                       else
                       {
                           sscanf(line+12,"%i%f%f%f",&h,a0->eng->kchi2+n,
                           a0->eng->chi02+n,a0->eng->nchi2+n);
                           a0->eng->kchi2[n]=a0->eng->kchi2[n]/h;
                           if(a0->eng->nchi2[n]>0) 
                           {
                              if(a0->eng->nchi[n]<0) 
                              a0->eng->nchi[n]=a0->eng->nchi2[n];
                           }
                           else
                           {
                              a0->eng->nchi2[n]=999999;
                           }

                           continue;
                       }
                    }
                }
                if(mt!=1)
                {
                    sprintf(line0,"%s-%s-%s-%s",a1->bond[j]->type,a1->type,
                            a0->type,a0->bond[i]->type);
                    if(strncmp(line+bg,line0+bg,ed)==0)
                    {
                       if(a0->eng->chi0[n]>9999)
                       {
                           sscanf(line+12,"%i%f%f%f",&h,a0->eng->kchi+n,
                           a0->eng->chi0+n,a0->eng->nchi+n);
                           a0->eng->kchi[n]=a0->eng->kchi[n]/h;
                           continue;
                       }
                       else
                       {
                           sscanf(line+12,"%i%f%f%f",&h,a0->eng->kchi2+n,
                           a0->eng->chi02+n,a0->eng->nchi2+n);
                           a0->eng->kchi2[n]=a0->eng->kchi2[n]/h;
                           if(a0->eng->nchi2[n]>0)
                           {
                              if(a0->eng->nchi[n]<0)
                              a0->eng->nchi[n]=a0->eng->nchi2[n];
                           }
                           else
                           {
                              a0->eng->nchi2[n]=999999;
                           }

                           continue;
                       }
                    }
                }
              }//i
            }//j
          }//a0
        }//r
      }//m=3
      if(m==5)
      { 
        for(r=&TRES;r;r=r->next) 
        {
          for(a=r->tatm;a;a=a->next)
          {
            if(a->eng->radius>0) continue;
            if(strncmp(a->type,line,2)==0)
            {
              sscanf(line+6,"%f%f",&a->eng->radius,&a->eng->epslon);
              continue;
            }
            for(i=0;i<neq;i++)
            {
              if(strncmp(line,equal[i],2)==0)
              {
                if(strstr(equal[i],a->type))
                {
                  sscanf(line+6,"%f%f",&a->eng->radius,&a->eng->epslon);
                  break;
                }
              }
            }
          }
        }
      }//m=5
    }
}

