#include"source.h"

Tres::Tres()
{
  int i;
  selfvdw=1;
  logf=0;
  name='?';
  strcpy(name3,"???");
  number=0;
  nhydr=0;
  next=0;
  tatm=0;
  id=0;
  area=0;
  flag=100;
  head=new float[9];
  for(i=0;i<9;i++)head[i]=0;
  numrot=0;
  nsmt=0;
  smooth=0;
  smt=1;
  rotamer=0;
  constant=new ModConstant();
  popbin=0;
  simatm=0;
  engtable=0;
  engcoeff=0;
  contable=0;
  logg=0;
  smoothclash=0;
  smoothmin=0;
  hetm=0;
}

Tres::~Tres()
{
if(contable) delete contable;
if(tatm) delete tatm;
if(next) delete next;
if(head) delete [] head;
if(rotamer) delete rotamer;
if(constant) delete constant;
if(popbin) delete popbin;
if(simatm) delete simatm;
if(engtable)delete engtable;
if(engcoeff)delete engcoeff;
if(hetm) delete hetm;
hetm=0;
simatm=0;
tatm=0;
next=0;
head=0;
rotamer=0;
}

void Tres::readmore(char *back,char *side) {
   rotamer=new Pdb;
   rotamer->token=strdup("backbone");
   rotamer->next=new Pdb;
   rotamer->next->token=strdup("sidechain");
   rotamer->readmore(back);
   rotamer->next->readmore(side);
   set3backrotamer(back);
   set3backrotamer(side);
}

Pdb *Tres::findrotamer(char *s) {

    for(Pdb *c=rotamer;c;c=c->next) {
	if(c->token&&strcmp(c->token,s)==0) return c;
    }	

    return 0;
}

Pdb *Tres::findrotamername(char *s) {

    for(Pdb *c=rotamer;c;c=c->next) {
        if(c->name&&strcmp(c->name,s)==0) return c;
    }

    return 0;
}

void Tres::readmore(char *side) {
   Pdb *s;
   for(s=rotamer;s&&s->next;s=s->next);
   if(s==0) {
	rotamer=new Pdb;
	s=rotamer;
   }
   else {
	s->next=new Pdb;
	s=s->next;
   }
   s->readmore(side);
   if(s->chn==0) return;
   for(Res *r=s->chn->res;r;r=r->next){
	for(Res *rr=r;rr;rr=rr->more) {
		rr->chn=s->chn;
	}
   }
}


Tatm *Tres::findfarringatm(Tatm *a) {

	Tatm *t;
	float d0=-1;
	Tatm *t0=0;
	for(t=tatm;t;t=t->next) {
		if(t->name[1]=='H') continue;
		if(t->ring==0) continue;
		if(a->ring==1&&t->ring==2) continue;
		if(a->ring==2&&t->ring==1) continue;
		float d=distance(a->xyz,t->xyz);	
		if(d0<d) {
			d0=d;
			t0=t;
		} 
	}
	return t0;
}

Tatm **Tres::findtwofarringatm(Tatm *a) {
	
	Tatm *t;	 
	Tatm  *t0[100];
	float temp[100];
	int   order[100];
	 
	int m=0;
	for(t=tatm;t;t=t->next) {
		if(t->name[1]=='H') continue;
		if(t->ring==0) continue;
		if(t==a) continue;
		if(name=='W'&&t->ring==1) continue;
		float d=distance(a->xyz,t->xyz);	
		temp[m]=-d;
		t0[m]=t;
		m++;
	}
	t0[m]=0;
	if(m<2) return 0;

	Qsort cc;

	cc.sort(temp,m,order);
	
	int i;
	Tatm **out=new Tatm *[2];
	out[0]=0;
	out[1]=0;	
	for(i=0;i<2;i++) {
		int j=order[i];
		out[i]=t0[j];
	} 

	return out;
}

void Tres::read(char *filnam)
{
  FILE *fp;
  int i,k,j;
  char *mlib,line[200];
  char residue;
  Tres *tres_temp;
  Tatm *tatm_temp,*tatm_temp0,*tatm_temp1,*tatm_temp2;
  float d;

//find the file path and open it!
  mlib=RCS["library"];
  if(mlib==0) {
	cerr<<"library does not exist in jackal.dir"<<endl;
	cerr<<"check the jackal.dir file first in your current directory then in the environment path:JACKALDIR"<<endl;
	cerr<<"something wrong in specifying the library"<<endl;
	cerr<<"make sure the library points to the abosulte path"<<endl;
	exit(0);
  }
  fp=0;i=0;
  //cerr<<filnam<<mlib+2<<endl;
  sprintf(line,"%s/%s",mlib,filnam);
  i=strlen(mlib);
  while((fp=fopen(line,"r"))==NULL)
  {
   mlib=mlib+i+1;
   i=strlen(mlib);
   if(*mlib=='\0')
   {
	cerr<<"warning:could not open: "<<filnam<<endl;
	cerr<<"check the jackal.dir file first in your current directory then in the environment path:JACKALDIR"<<endl;
        cerr<<"something wrong in specifying the library"<<endl;
        cerr<<"make sure the library points to the abosulte path"<<endl;
	exit(0);
   }
   sprintf(line,"%s/%s",mlib,filnam);
  }
//read line by line 
  residue='?';
  while(fgets(line,200,fp)!=NULL)
  {
    if(line[0]=='!')continue;
    if(line[0]==' '||line[0]=='\n') continue;
    if(line[4]!=residue)  //new residue
    {
      if(residue!='?')   // the next residue
      {
        tres_temp->next=new Tres;
        tres_temp=tres_temp->next;    
      }
      else   tres_temp=this; //the first residue        
      residue=line[4];
      tres_temp->name=line[4];
      strncpy(tres_temp->name3,line,3);
      tres_temp->name3[3]='\0';
    } 
    if(tres_temp->tatm==0) // the first atom
    {
      tres_temp->tatm=new Tatm;
      tatm_temp=tres_temp->tatm;
    }
    else //the next atom
    {
      tatm_temp->next=new Tatm; 
      tatm_temp=tatm_temp->next;
    }    
    tatm_temp->tres=tres_temp;
    strncpy(tatm_temp->name,line+6,4);
    tatm_temp->name[4]='\0';
    for(i=0;i<3;i++)
    tatm_temp->xyz[i]=atof(line+11+9*i);
    //tres_temp->number++;
    tatm_temp->rotate=atoi(line+39);
    tatm_temp->balance=atoi(line+41);
    tatm_temp->hbond=atoi(line+43);
    tatm_temp->keep=atoi(line+45);
    tatm_temp->nonp=atoi(line+47);
    if(tatm_temp->nonp==2)tatm_temp->nonp=-1;    
    tatm_temp->ring=atoi(line+49);
    tatm_temp->polarkeep=atoi(line+51);
    tatm_temp->sidecenter=atoi(line+53);
  }
        
  fclose(fp);

  for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
  {
     i=0;
     tatm_temp=tres_temp->tatm;
     while(tatm_temp)
     {
       if(strstr(tatm_temp->name,"HEAD"))
       {
         copy(tatm_temp->xyz,tres_temp->head+i,3);
         i+=3; 
       }
       else 
       {
         tatm_temp0=tatm_temp;
       }
       tatm_temp=tatm_temp->next;
     }
     delete tatm_temp0->next;
     tatm_temp0->next=0;
  }

//kick out not-intended atoms

  if(flag%10!=0) // not all atom model!
  {
    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
    {
      tatm_temp0=tres_temp->tatm;
      tatm_temp=tatm_temp0->next;
      while(tatm_temp)
      {
        if(flag%10==1) //only polar HB
        {
           if(tatm_temp->keep==0&&tatm_temp->name[1]=='H')
           {
             tatm_temp1=tatm_temp->next;
             tatm_temp->next=0;
             delete tatm_temp;
             tatm_temp=tatm_temp1;
             continue;
           }
        }
        
        else if(flag%10==2) //no hydrogen!
        {
           if(tatm_temp->name[1]=='H')
           {
             tatm_temp1=tatm_temp->next;
             tatm_temp->next=0;
             delete tatm_temp;
             tatm_temp=tatm_temp1;
             continue;
           }
        }
        else if(flag%10==3) //only mainchain H
        {
           if(tatm_temp->name[1]=='H') 
           {
             if(tatm_temp->name[2]=='A') goto re200;
             if(tatm_temp->name[2]=='B') goto re200;
             if(tatm_temp->name[2]=='N') goto re200;
             tatm_temp1=tatm_temp->next;
             tatm_temp->next=0;
             delete tatm_temp;
             tatm_temp=tatm_temp1;
             continue;
           }
        }
        else if(flag%10==4) //only side-chain
        {
           if(tatm_temp->name[1]=='H')
           {
             if(tatm_temp->name[2]=='A'||tatm_temp->name[2]=='B'||tatm_temp->name[2]=='N')
             {
               tatm_temp1=tatm_temp->next;
               tatm_temp->next=0;
               delete tatm_temp;
               tatm_temp=tatm_temp1;
               continue;
             }
           }
        }
	else if(flag%10==5) //only polar HB, which can be clearly defined
        {
           if(tatm_temp->polarkeep==0&&tatm_temp->name[1]=='H')
           {
             tatm_temp1=tatm_temp->next;
             tatm_temp->next=0;
             delete tatm_temp;
             tatm_temp=tatm_temp1;
             continue;
           }
        }
	else if(flag%10==6) //only mainchain H
        {
           if(tatm_temp->name[1]=='H')
           {
             //if(tatm_temp->name[2]=='A') goto re200;
             //if(tatm_temp->name[2]=='B') goto re200;
             if(tatm_temp->name[2]=='N') goto re200;
             tatm_temp1=tatm_temp->next;
             tatm_temp->next=0;
             delete tatm_temp;
             tatm_temp=tatm_temp1;
             continue;
           }
        }
	else if(flag%10==7) //only backbone
 	{
    		if(tatm_temp->name[1]=='H'||tatm_temp->id>4)
    		{
     	 		//if(tatm_temp->name[2]=='A') goto re200;
      			//if(tatm_temp->name[2]=='B') goto re200;
      			if(tatm_temp->name[2]=='N') goto re200;
      			tatm_temp1=tatm_temp->next;
      			tatm_temp->next=0;
      			delete tatm_temp;
      			tatm_temp=tatm_temp1;
      			continue;
    		}
 	}
        else  break; //not defined,other case!
        re200:
        tatm_temp0->next=tatm_temp;
        tatm_temp0=tatm_temp;
        tatm_temp=tatm_temp->next;
        continue;
      } 
      tatm_temp0->next=0;
    }
  }

//calculate the number of H and id for atom and residue

  i=0;
  for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
  {
    tres_temp->number=0;
    tres_temp->nhydr=0;
    tres_temp->id=i;
    i++;
    for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
    {
     tatm_temp->id=tres_temp->number;
     tres_temp->number++;
     if(tatm_temp->name[1]=='H')tres_temp->nhydr++;
    }
  }

//connectivity,find covalent bond;

  for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
  {
    for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
    {
      k=0;
      if(tatm_temp->id==0) tatm_temp->bond[k++]=tres_temp->isatm(2);
      for(tatm_temp0=tres_temp->tatm;tatm_temp0;tatm_temp0=tatm_temp0->next)
      {
        if(tatm_temp0==tatm_temp) continue;
        d=distance(tatm_temp->xyz,tatm_temp0->xyz);
        if(d>2.) continue;
        if(tatm_temp->name[1]=='H'||tatm_temp0->name[1]=='H')
        if(d>1.5) continue;
        if(k>4) 
        {  
           if(TRES.logg) cerr<<"warning!"<<tres_temp->name<<":"<<tatm_temp->name;
           if(TRES.logg) cerr<<" has more than four bonds"<<endl;
        }
        tatm_temp->bond[k]=tatm_temp0;
        k++;
      }
      if(tatm_temp->id==2) 
      {
       tatm_temp->bond[k++]=tatm_temp->bond[1];
       tatm_temp->bond[1]=tres_temp->tatm;
      }
      tatm_temp->nbond=k;
    }
    if(tres_temp->name=='P')
    {
      tatm_temp=tres_temp->isatm(" CD ");
      tatm_temp0=tatm_temp->bond[0];
      tatm_temp->bond[0]=tatm_temp->bond[1];
      tatm_temp->bond[1]=tatm_temp0;
    }
  }

//create eng term for each atom;

  for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
  for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  tatm_temp->eng->tank();

  if(flag/100==1)
  {
    tatm->eng->charm("charm_22.top","charm_22.prm");
  }
  else if(flag/100==2)
  {
    tatm->eng->amber("amber_all_amber94.in","amber_parm94.dat");
  }
  else if(flag/100==3)
  {
    tatm->eng->charm("charm_19.top","charm_19.prm");
  }
  else if(flag/100==4)
  {
    tatm->eng->amber("amber_uni_amber.in","amber_parm91.dat");
  }

//reset charge
  

  if((flag%100)/10==1) //turn off all charge
  {
    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
    {
      for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
        tatm_temp->eng->charge=0;
    }
  }

  if((flag%100)/10==2) //turn off charges on mainchain-chain
  {
    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
    {
      for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
      {
        if(tatm_temp->name[1]=='H')
        {
          if(tatm_temp->bond[0]->id>3) continue;
        }
        else
        {
          if(tatm_temp->id>3) continue;
        }
        tatm_temp->eng->charge=0;
      }
    }
  }
  
  if((flag%100)/10==3) //turn off charges on side-chain
  {
    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
    {
      for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
      {
        if(tatm_temp->name[1]=='H')
        {
          if(tatm_temp->bond[0]->id<4) continue;
        }
        else
        {
          if(tatm_temp->id<4) continue;
        }
        tatm_temp->eng->charge=0;
      }
    }
  }
    
  if((flag%100)/10==4)//only consider charges on polar atom
  {

    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next) //turn off charges on sidechain
    {
      for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
      {
        if(tatm_temp->name[1]=='H')
        {
          if(tatm_temp->bond[0]->id<4) continue;
        }
        else
        {
          if(tatm_temp->id<4) continue;
        }
        tatm_temp->eng->charge=0;
      }
    }
/*
    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
    {
      tres_temp->tatm->eng->charge+=tres_temp->tatm->next->eng->charge;
      tres_temp->tatm->next->eng->charge=0;
      tatm_temp=tres_temp->isatm(" HA ");      
      if(tatm_temp)
      { 
       tatm_temp->bond[0]->eng->charge+=tres_temp->tatm->eng->charge;
       tatm_temp->eng->charge=0;
      }      
      
      for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
      {                 
        if(tatm_temp->name[1]=='H')
        {
          if(tatm_temp->bond[0]->id<4) continue;
        }
        else
        {
          if(tatm_temp->id<4) continue;
        }
        tatm_temp->eng->charge=0;
      }
    }
*/

    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
    {
      if(tres_temp->name=='D')
      {
          tatm_temp=tres_temp->isatm(" OD1");  
          if(tatm_temp)tatm_temp->eng->charge=-0.5;
          tatm_temp=tres_temp->isatm(" OD2");  
          if(tatm_temp)tatm_temp->eng->charge=-0.5;
      }

      if(tres_temp->name=='E')
      {
          tatm_temp=tres_temp->isatm(" OE1");  
          if(tatm_temp)tatm_temp->eng->charge=-0.5;
          tatm_temp=tres_temp->isatm(" OE2");  
          if(tatm_temp)tatm_temp->eng->charge=-0.5;
      }
      if(tres_temp->name=='K')
      {
          tatm_temp=tres_temp->isatm(" NZ ");
          tatm_temp0=tres_temp->isatm("1HZ ");
          tatm_temp1=tres_temp->isatm("2HZ ");
          tatm_temp2=tres_temp->isatm("3HZ ");
          if(tatm_temp0&&tatm_temp1&&tatm_temp2)
          {
            tatm_temp0->eng->charge=0.33;
            tatm_temp1->eng->charge=0.33;
            tatm_temp2->eng->charge=0.33;
          }
          else
          {
            if(tatm_temp)tatm_temp->eng->charge=1.0;
          }
      }
      if(tres_temp->name=='R')
      {
          tatm_temp=tres_temp->isatm(" NH2");
          tatm_temp0=tres_temp->isatm("1HH2");
          tatm_temp1=tres_temp->isatm("2HH2");
          if(tatm_temp0&&tatm_temp1&&tatm_temp2)
          {
            tatm_temp0->eng->charge=0.25;
            tatm_temp1->eng->charge=0.25;
          }
          else
          {
           if(tatm_temp)tatm_temp->eng->charge=0.5;
          }
          tatm_temp=tres_temp->isatm(" NH1");
          tatm_temp0=tres_temp->isatm("1HH1");
          tatm_temp1=tres_temp->isatm("2HH1");
          if(tatm_temp0&&tatm_temp1&&tatm_temp2)
          {
            tatm_temp0->eng->charge=0.25;
            tatm_temp1->eng->charge=0.25;
          }
          else
          {
           if(tatm_temp)tatm_temp->eng->charge=0.5;
          }
      }
      if(tres_temp->name=='N')
      {
          tatm_temp=tres_temp->isatm(" ND2");
          if(tatm_temp)tatm_temp->eng->charge=-0.30;
          tatm_temp=tres_temp->isatm(" OD1");
          if(tatm_temp)tatm_temp->eng->charge=-0.30;
          tatm_temp=tres_temp->isatm(" CG ");
          if(tatm_temp)tatm_temp->eng->charge=0.60;
      }
      if(tres_temp->name=='Q')
      {
          tatm_temp=tres_temp->isatm(" NE2");
          if(tatm_temp)tatm_temp->eng->charge=-0.30;
          tatm_temp=tres_temp->isatm(" OE1");
          if(tatm_temp)tatm_temp->eng->charge=-0.30;
          tatm_temp=tres_temp->isatm(" CD ");
          if(tatm_temp)tatm_temp->eng->charge=0.60;
      }
      if(tres_temp->name=='S')
      {
          tatm_temp=tres_temp->isatm(" OG ");
          if(tatm_temp)tatm_temp->eng->charge=-0.25;
          tatm_temp=tres_temp->isatm(" CB ");
          if(tatm_temp)tatm_temp->eng->charge=0.25;
      }
      if(tres_temp->name=='T') 
      {
          tatm_temp=tres_temp->isatm(" OG1");
          if(tatm_temp)tatm_temp->eng->charge=-0.25;
          tatm_temp=tres_temp->isatm(" CB ");
          if(tatm_temp)tatm_temp->eng->charge=0.25;
      }
    }
  }

  if((flag%100)/10==5) //only consider charges on polar atom self defined
  {
    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
    {

      for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
      {
        tatm_temp->eng->charge=0;
        if(strcmp(tatm_temp->name," N  ")==0) tatm_temp->eng->charge= 0.3;
        if(strcmp(tatm_temp->name," O  ")==0) tatm_temp->eng->charge=-0.3;
      }      

      if(tres_temp->name=='D')
      {
          tatm_temp=tres_temp->isatm(" OD1");  
          if(tatm_temp)tatm_temp->eng->charge=-0.5;
          tatm_temp=tres_temp->isatm(" OD2");  
          if(tatm_temp)tatm_temp->eng->charge=-0.5;
      }

      if(tres_temp->name=='E')
      {
          tatm_temp=tres_temp->isatm(" OE1");  
          if(tatm_temp)tatm_temp->eng->charge=-0.5;
          tatm_temp=tres_temp->isatm(" OE2");  
          if(tatm_temp)tatm_temp->eng->charge=-0.5;
      }
      if(tres_temp->name=='K')
      {
          tatm_temp=tres_temp->isatm(" NZ ");         
          if(tatm_temp)tatm_temp->eng->charge=1.0;          
      }
      if(tres_temp->name=='R')
      {
          tatm_temp=tres_temp->isatm(" NH2");        
          if(tatm_temp)tatm_temp->eng->charge=0.5;       
          tatm_temp=tres_temp->isatm(" NH1");        
          if(tatm_temp)tatm_temp->eng->charge=0.5;
          
      }
      if(tres_temp->name=='N')
      {
          tatm_temp=tres_temp->isatm(" ND2");
          if(tatm_temp)tatm_temp->eng->charge= 0.00;
          tatm_temp=tres_temp->isatm(" OD1");
          if(tatm_temp)tatm_temp->eng->charge=-0.00;         
      }
      if(tres_temp->name=='Q')
      {
          tatm_temp=tres_temp->isatm(" NE2");
          if(tatm_temp)tatm_temp->eng->charge= 0.00;
          tatm_temp=tres_temp->isatm(" OE1");
          if(tatm_temp)tatm_temp->eng->charge=-0.00;         
      }   
    }
  }
  
  if((flag%100)/10==6)//switch N charge with CA
  {
    for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
    {
      tatm_temp=tres_temp->isatm(" HN ");
      if(tatm_temp==0)
      {
        if(tres_temp->tatm->eng->charge<=0&&tres_temp->tatm->next->eng->charge>=0)
        {
          d=tres_temp->tatm->eng->charge;  
          tres_temp->tatm->eng->charge=tres_temp->tatm->next->eng->charge;
          tres_temp->tatm->next->eng->charge=d;
        }
      }
    }
  }
  for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
  for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  {
    if(fabs(tatm_temp->eng->charge)<0.0001) tatm_temp->eng->charge=0;
  }

//set for HN
  for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
  for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  {
    if(tatm_temp->name[1]=='H'&&tatm_temp->name[2]=='N')tatm_temp->ishn=1;
    if(tres_temp->name=='P') 
    if(tatm_temp->name[1]=='C'&&tatm_temp->name[2]=='D')tatm_temp->ishn=1;
  }


//set undefined value!
  i=0;j=0;k=0;
  for(tres_temp=this;tres_temp;tres_temp=tres_temp->next)
  {
   tres_temp->charge=0;
  for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  {
   tres_temp->charge+=tatm_temp->eng->charge;
   for(i=0;i<tatm_temp->eng->mbond;i++)
   {
     if(tatm_temp->eng->bond0[i]>9999) {i++;tatm_temp->eng->bond0[i]=0;}
     if(tatm_temp->eng->kbond[i]>9999) tatm_temp->eng->kbond[i]=0;
   }
   for(i=0;i<tatm_temp->eng->mangle;i++)
   {
     if(tatm_temp->eng->angle0[i]>9999) {j++;tatm_temp->eng->angle0[i]=0;}
     if(tatm_temp->eng->kangle[i]>9999) tatm_temp->eng->kangle[i]=0;
   }
   for(i=0;i<tatm_temp->eng->mchi;i++)
   {
     if(tatm_temp->eng->chi0[i]>9999) {k++;tatm_temp->eng->chi0[i]=0;}
     if(tatm_temp->eng->kchi[i]>9999) tatm_temp->eng->kchi[i]=0;
     if(tatm_temp->eng->nchi[i]>9999) tatm_temp->eng->nchi[i]=0;
   }
   if(tatm_temp->eng->radius==0)
   {
     if(tatm_temp->name[1]=='C') {
       tatm_temp->eng->radius=1.95;
       tatm_temp->eng->epslon=-0.11;
     }
     else if(tatm_temp->name[1]=='N'){
       tatm_temp->eng->radius=1.85;
       tatm_temp->eng->epslon=-0.2; 
     }
     else if(tatm_temp->name[1]=='O'){
       tatm_temp->eng->radius=1.75;
       tatm_temp->eng->epslon=-0.12;
     }
     else if(tatm_temp->name[1]=='S'){
       tatm_temp->eng->radius=2.0;
       tatm_temp->eng->epslon=-0.45;
     }
     else if(tatm_temp->name[1]=='H'){
       tatm_temp->eng->radius=0.2;
       tatm_temp->eng->epslon=0.0;
     }
   }
  }
  }
  if(TRES.logg) cerr<<"the number of bond unset:"<<i<<endl;
  if(TRES.logg) cerr<<"the number of angle unset:"<<j<<endl;
  if(TRES.logg) cerr<<"the number of chi unset:"<<k<<endl;

  for(tres_temp=&TRES;tres_temp;tres_temp=tres_temp->next)
  {
    i=0;
    for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
    {
      if(tatm_temp->rotate==1)  
      {
        i++;
        tatm_temp->ntat=i;
      }
    }
    for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
    {
      if(tatm_temp->name[1]=='H') continue;
      if(tatm_temp->id==0) continue;
      tatm_temp0=tatm_temp->bond[0];
      while(tatm_temp0&&tatm_temp0->ntat==0)
      {
        if(tatm_temp0->id!=0)
         tatm_temp0=tatm_temp0->bond[0];
        else break;
      }
      if(tatm_temp0->tres->name!=tres_temp->name) continue;
      tatm_temp->rotm=tatm_temp0->ntat; 
    }
/*
    for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
    {
      if(tatm_temp->name[1]=='H') tatm_temp->rotm=tatm_temp->bond[0]->rotm;
      cerr<<tres_temp->name<<" "<<tatm_temp->name<<" "<<tatm_temp->rotm<<endl;
    }
*/
  }

  for(tres_temp=&TRES;tres_temp;tres_temp=tres_temp->next)
  tres_temp->numrot=tres_temp->rotable(tres_temp->tatm);

  for(tres_temp=&TRES;tres_temp;tres_temp=tres_temp->next)
  for(tatm_temp=tres_temp->tatm;tatm_temp;tatm_temp=tatm_temp->next)
  {
    tatm_temp->bondleng();
  }
  
  for(tres_temp=&TRES;tres_temp;tres_temp=tres_temp->next) tres_temp->ispolar();

  smoothf("smooth");
  if(smt>=nsmt) smt=nsmt-1;
  switchcharge(0);
}

void Tres::smoothf(char *filnam)
{
 char line[200],*mlib,*poi;
 FILE *fp;
 int i,j;
 float smooth0[4000];
 float a,b,c,n;

 mlib=0;poi=0;
 mlib=RCS["library"];
 
 fp=0; i=0;
 sprintf(line,"%s/%s",mlib,filnam);
 i=strlen(mlib);
 while((fp=fopen(line,"r"))==NULL)
 {
   mlib=mlib+i+1;
   i=strlen(mlib);
   if(*mlib=='\0')
   {cerr<<"warning:could not open: "<<filnam<<endl;}
   sprintf(line,"%s/%s",mlib,filnam);
 }

 if(!fp)
 {
   smooth=new float[4];
   nsmt=1;
   smooth[0]=1;smooth[1]=0;smooth[2]=2;smooth[3]=6;
   return;
 }

 i=0;
 while(fgets(line,200,fp)!=NULL)
 {
   poi=line;
   while(*poi==' ')poi++;
   if(*poi=='!') continue;
   if(strlen(poi)<10) continue;
   sscanf(poi,"%f%f%f%f",&a,&b,&c,&n);
   smooth0[i*4]=a;
   smooth0[i*4+1]=b;
   smooth0[i*4+2]=c;
   smooth0[i*4+3]=n;
   i++;
 }
 if(i==0)
 {
   smooth=new float[4];
   nsmt=1;
   smooth[0]=1;smooth[1]=0;smooth[2]=2;smooth[3]=6;
   smt=0;
   return;
 }

 nsmt=i;
 smooth=new float[i*4];
 for(j=0;j<i*4;j++) smooth[j]=smooth0[j];
}

void Tres::ispolar()
{
  Tatm *t,*t1;
  int i,j;

  for(t=tatm;t;t=t->next) 
  {
     if(t->name[1]=='N'||t->name[1]=='O') {t->ispolar=3;continue;}
     j=1;
     for(i=0;i<t->nbond;i++)
     {
       t1=t->bond[i];
       if(t1==0) continue;
       if(t1->name[1]=='N'||t1->name[1]=='O') j=2;
     }
     t->ispolar=j;
  }
}

char Tres::swap(char *s)
{
  Tres *tres_temp;
  for(tres_temp=&TRES;tres_temp;tres_temp=tres_temp->next)
  {
    if(strncmp(tres_temp->name3,s,3)==0) return tres_temp->name;
  }
  return '?';
}

char *Tres::swap(char s)
{
  Tres *tres_temp;
  for(tres_temp=&TRES;tres_temp;tres_temp=tres_temp->next)
  {
    if(tres_temp->name==s) return tres_temp->name3;
  }
  return "???";
}
 
int Tres::getuse(char *s,int m,int n)
{
Tatm *a;
int i;
i=0;
for(a=tatm;a;a=a->next)
if(strncmp(a->name+m,s+m,n)==0) i++;
return i;
}

Tatm *Tres::isatm(int s)
{
Tatm *tatm_temp;
for(tatm_temp=tatm;tatm_temp;tatm_temp=tatm_temp->next)
{
  if(tatm_temp->id==s) return tatm_temp;
}
return 0;
}

 
Tatm *Tres::isatm(char *s)
{
Tatm *a,*a1;
int i;
for(a=tatm;a;a=a->next)
{
  if(strncmp(a->name,s,4)==0)return a;
}
i=0;
for(a=tatm;a;a=a->next)
{
  if(strncmp(a->name+1,s+1,3)==0)
  {
    a1=a;
    i++;
  }
}
if(i==1) return a1;

i=0;
for(a=tatm;a;a=a->next)
{
  if(strncmp(a->name+1,s+1,2)==0)
  {
    a1=a;
    i++;
  }
}
if(i==1) return a1;
return 0;
}


Tatm *Tres::istatm(char *s) {
Tatm *a;
for(a=tatm;a;a=a->next)
{
  if(strncmp(a->name,s,4)==0)return a;
}
return 0;
}

Tatm *Tres::isatm(char *s, int m,int n,int t)
{
Tatm *a;
int i;
i=0;
for(a=tatm;a;a=a->next)
{
  if(strncmp(a->name+m,s+m,n)==0) 
  {
    if(i==t) return a; 
    i++;
  }
}
return 0;
}

Tres *Tres::operator[](char s)
{
Tres *tres_temp;
for(tres_temp=&TRES;tres_temp;tres_temp=tres_temp->next)
{
  if(tres_temp->name==s) return tres_temp;
}
return 0;
}

Tres *Tres::operator[](int c) 
{
Tres *r1;
for(r1=&TRES;r1;r1=r1->next)
if(r1->id==c) return r1;
return 0;
}

Tres *Tres::operator[](char *ch)
{
  Tres *tres_temp;
  for(tres_temp=&TRES;tres_temp;tres_temp=tres_temp->next)
  {
    if(strncmp(tres_temp->name3,ch,3)==0) return tres_temp;
  }
  return 0;
}

float Tres::dihedral(float *x1,float *x2,float *x3,float *x4)
{

float a[3],b[3],c[3],n1[3],n2[3],n[3];
float d1,d2,d3;
float ang;
int i;
// three vecotrs

  for(i=0;i<3;i++)
  {
   a[i]=x2[i]-x1[i]; 
   b[i]=x2[i]-x3[i];
   c[i]=x3[i]-x4[i];
  }  

//n1 = a X b, n2 = c X b 

  n1[0] = a[1]*b[2] - a[2]*b[1];
  n1[1] = a[2]*b[0] - a[0]*b[2];
  n1[2] = a[0]*b[1] - a[1]*b[0];
  n2[0] = c[1]*b[2] - c[2]*b[1];
  n2[1] = c[2]*b[0] - c[0]*b[2];
  n2[2] = c[0]*b[1] - c[1]*b[0];
  n[0] = n2[1]*n1[2] - n2[2]*n1[1];
  n[1] = n2[2]*n1[0] - n2[0]*n1[2];
  n[2] = n2[0]*n1[1] - n2[1]*n1[0];
   

/* angle between n1 and n2 is the dihedral */

  d1 = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
  d2 = sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
  d3 = sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2]);
 
  if(d2==0||d3==0) return 999;
  ang = d1/(d2*d3);
  if(ang>=1) ang=1.;
  else if(ang<=-1) ang=-1.;
  ang = acos(ang);
  ang = ang *180.0/3.1415926;
  d3 = b[0]*n[0] + b[1]*n[1] + b[2]*n[2];
  if (d3 < 0) {
    ang = -ang;
  }
  return ang;
}

void Tres::vectorm(float *x,float *y,float *z)
{
 int i;
 for(i=0;i<3;i++) z[i]=y[i]-x[i];
}

float Tres::vectord(float *x,float *y)
{
float d;
d=x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
return d;
}

void Tres::vectorx(float *x,float *y,float *z)
{
    z[0]=x[1]*y[2]-x[2]*y[1];
    z[1]=x[2]*y[0]-x[0]*y[2];
    z[2]=x[0]*y[1]-x[1]*y[0];
}
float Tres::distance(float *x)
{
float d;
d=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
d=sqrt(d);
return d;
}

float Tres::distance(Atm *x,Atm *y) {
  return distance(x->xyz,y->xyz);
}
float Tres::distance(float *x,float *y)
{
float d;
d=(x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]);
d=sqrt(d);
return d;
}

float Tres::distance(float x0,float x1,float x2,float y0,float y1,float y2)
{
float d;
d=(x0-y0)*(x0-y0)+(x1-y1)*(x1-y1)+(x2-y2)*(x2-y2);
d=sqrt(d);
return d;
}

float Tres::distsqr(float x0,float x1,float x2,float y0,float y1,float y2)
{
float d;
d=(x0-y0)*(x0-y0)+(x1-y1)*(x1-y1)+(x2-y2)*(x2-y2);
return d;
}

float Tres::distsqr(float *x,float *y)
{
float d;
d=(x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]);
return d;
}

float Tres::angle(float *x,float *y,float *z)
{
float a[3],b[3],c,d,e; 
int i;
for(i=0;i<3;i++) {b[i]=z[i]-y[i];a[i]=x[i]-y[i];}
c=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
d=distance(x,y);
e=distance(y,z);
if(d==0||e==0) return 999;
c=acos(c/(d*e))/3.14159*180;
return c;
}

void Tres::copy(float *from, float *to, int n)
{
int i;
for(i=0;i<n;i++)
to[i]=from[i];
}


void Tres::surface(int m,float probe)
{
Tatm *a,*a1;
float xyz[3];
float t1,t2,f1,f2,c1,c2,c3,d;
int i,j,k;
area=0;

//add the probe to atoms radius;

for(a=tatm;a;a=a->next)
a->eng->radius+=probe;

if(m==0) {if(TRES.logg) cerr<<"give integer value in surface calculation\n";exit(0);}

t1=3.1415926/(m+1);
t2=2*3.1415926/(2*m+1);
c3=t1*t2;

//set up the lattice

 area=0;
 for(a=tatm;a;a=a->next)
 {
   a->area=0;
   for(i=0;i<m+1;i++)
   for(j=0;j<2*m+1;j++)
   {   
       f1=i*t1;f2=j*t2;
       c1=a->eng->radius*sin(f1);
       c2=a->eng->radius*a->eng->radius;
       xyz[0]=c1*cos(f2)+a->xyz[0];
       xyz[1]=c1*sin(f2)+a->xyz[1];
       xyz[2]=a->eng->radius*cos(f1)+a->xyz[2];
       k=0;
       for(a1=tatm;a1;a1=a1->next)
       {
         if(a1==a) continue;
         d=TRES.distance(xyz,a1->xyz);
         if(d<a1->eng->radius){k=1;break;} 
       }
       if(k) continue;
       a->area+=c3*c2*sin(f1);
   }
   //cout<<a->name<<" "<<name3<<" "<<a->area<<endl;
   area+=a->area;
 }

//set back the radius

for(a=tatm;a;a=a->next)
a->eng->radius-=probe;

}

void Tres::surface(int m,float probe,int f)
{
Tres *r;
for(r=&TRES;r;r=r->next)
if(f>=0) 
{
 if(r->name==f) r->surface(m,probe);
}
else r->surface(m,probe);
}

float Tres::maxlen(int a,int b)
{
  Tatm *t1,*t2,*t;
  float d,x; 
  t1=isatm(a);
  t2=isatm(b);
  if(t1&&t2)
  {
    d=distance(t1->xyz,t2->xyz);
    mlen=d;
    return d;
  }
  if(t1&&!t2)
  {
    d=-100;
    for(t=tatm;t;t=t->next)
    {
      x=distance(t->xyz,t1->xyz);
      if(d<x)d=x;
    }
    mlen=d;
    return d;
  }
  else if(!t1&&t2)
  {
    d=-100;
    for(t=tatm;t;t=t->next)
    {
      x=distance(t->xyz,t2->xyz);
      if(d<x)d=x;
    }
    mlen=d;
    return d;
  }
  else
  {
    mlen=0;
    return 0;
  }
}

void Tres::maxlen(int a)
{
Tres *t;
for(t=&TRES;t;t=t->next) t->maxlen(a,100);
}


void Tres::transfer(float *co,int f)
{
int i,j;
Tatm *a;
i=0;
for(a=tatm;a;a=a->next)
{
for(j=0;j<3;j++)
if(f) a->xyz[j]=co[i*3+j];
else co[i*3+j]=a->xyz[j];
i++;
}
}

int Tres::rotable(Tatm *a)
{
Tatm *b;
int i;
i=0;
for(b=a;b;b=b->next)
{
  if(b->rotate==1) i++;
}
return i;
}

float Tres::distsqr(float *s)
{

  return s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
} 


float Tres::space(float *a1,float *a2,float *a3)
{
float b1,b2,b3,b;
b1=distance(a1,a2);
b2=distance(a1,a3);
b3=distance(a2,a3);
b=(b1+b2+b3)/2;
return sqrt(b*(b-b1)*(b-b2)*(b-b3));
}

char *Tres::delnonstandard(char *s) {

	if(s==0) return 0;

	int n=strlen(s);

	int i;

	Tres *t;
	
	char *out=new char[n+1];

	int k=0;
		
	for(i=0;i<n;i++) {

		t=TRES[s[i]];

		if(t==0) continue;	

		out[k++]=s[i];		
	}

	out[k]='\0';

	return out;
}


int Tres::ispolaraa() {

 	if(strchr("KRHDENQTS",name)) return 1;
	else return 0;
}

int Tres::isnonpolaraa() {

	if(strchr("AVLIFWCMY",name)) return 1;
	else return 0;
}

int Tres::getaanum() {

	return strlen("ACDEFGHIKLMNPQRSTVWY");
}

int Tres::isaa(char s) {
	char *n=strchr("ACDEFGHIKLMNPQRSTVWY",s);
	if(n) return 1;
	else return 0;
}

int Tres::getorder() {

	int n=0;

	Tres *t;

	for(t=&TRES;t;t=t->next) {
		if(t==this) return n;
		n++;
	}

	return -1;
}

void Tres::setchargezero() {

	for(Tres *t=&TRES;t;t=t->next) 
	for(Tatm *a=t->tatm;a;a=a->next) {
		a->eng->charge=0;
	}
}

void Tres::setdonaldcharge() {

	for(Tres *t=&TRES;t;t=t->next) 
        for(Tatm *a=t->tatm;a;a=a->next) {
		if(t->name=='D'&&strcmp(a->name," CG ")==0) {
                	a->eng->charge=-1;
		}
		else if(t->name=='E'&&strcmp(a->name," CD ")==0) {
			a->eng->charge=-1;
		}
		else if(t->name=='K'&&strcmp(a->name," NZ ")==0) {
			a->eng->charge=1;
		}
		else if(t->name=='R'&&strcmp(a->name," CZ ")==0) {
			a->eng->charge=1;
		}
		else {
			a->eng->charge=0;
		}
        }
}

void Tres::writedelphisize(char *file) {

	FILE *fp = fopen(file,"w");
	
	fprintf(fp,"atom__res_radius_\n");
	Tres *t;
	for(t=&TRES;t;t=t->next)
	for(Tatm *a=t->tatm;a;a=a->next) {
				
		fprintf(fp,"%s  %s %8.3f %8.3f %8.3f\n",a->name,t->name3,a->eng->radius,a->eng->dipole,a->eng->charge);
	}

	for(t=TRES.hetm;t;t=t->next)
        for(Tatm *a=t->tatm;a;a=a->next) {

                fprintf(fp,"%s  %s %8.3f %8.3f %8.3f\n",a->name,t->name3,a->eng->radius,a->eng->dipole,a->eng->charge);
        }
	fclose(fp);
}


void Tres::writedelphicharge(char *file) {

        FILE *fp = fopen(file,"w");

	fprintf(fp,"atom__resnumbc_charge_\n");
        for(Tres *t=&TRES;t;t=t->next)
        for(Tatm *a=t->tatm;a;a=a->next) {
                fprintf(fp,"%s  %s  %8.3f\n",a->name,t->name3,a->eng->charge);
        }

        fclose(fp);
}

void Tres::assigncharge(char *f) {

	FILE *fp;
	int i;
	char *mlib,line[200];
	char *filnam=f;	

	mlib=RCS["library"];
  	fp=0;i=0;
  	sprintf(line,"%s/%s",mlib,filnam);
  	i=strlen(mlib);
  	while((fp=fopen(line,"r"))==NULL)
  	{
   		mlib=mlib+i+1;
   		i=strlen(mlib);
   		if(*mlib=='\0')
   		{cerr<<"warning:could not open: "<<filnam<<endl;return;}
   		sprintf(line,"%s/%s",mlib,filnam);
  	}
			
	if(fp==0) {
		cerr<<"could not open:"<<f<<endl;
		return;
	}

	while(fgets(line,200,fp)!=NULL) {
		if(line[0]=='!')continue;
		if(strlen(line)<10) continue;
		if(strncmp(line,"atom__",6)==0) continue;
		Tres *t=TRES[line+6];
		if(t==0) continue;
		Tatm *a=t->isatm(line);
		if(a==0) continue;
		a->eng->charge=atof(line+10);
	}
	fclose(fp);
}

Tatm *Tres::findanyatmwithname(char *s) {
	
	for(Tres *t=&TRES;t;t=t->next) {
		Tatm *a=t->isatm(s);
		if(a) return a;
	}

	return 0;
}

void Tres::readdistpopular(char *fil,char *ind) {

	if(fil==0||ind==0) return;

	if(popbin==0) popbin=new DistPopularBin();

	DistPopularBin *t=0;

	for(t=popbin;t;t=t->next) {
		if(t->name==0&&t->dist==0) break; 
	}	

	if(t==0) {
		for(t=popbin;t->next;t=t->next);
		t->next=new DistPopularBin();
		t=t->next;
	}

	t->name=strdup(ind);
	t->dist=new DistPopular();
	t->dist->read(fil);
	
}

int Tres::getresid(char c) {

	Algn cc;

	int n=strlen(cc.aaorder);

	for(int i=0;i<n;i++) {
		if(cc.aaorder[i]==c) return i;
	}
	return -1;
}

void Tres::writechargeradius(FILE *fp){

	if(fp==0) return;

	for(Tres *s=&TRES;s;s=s->next) {

		for(Tatm *t=s->tatm;t;t=t->next) {

			fprintf(fp,"%s %s %8.3f %8.3f %8.3f\n",s->name3,t->name,t->eng->charge,t->eng->radius,t->eng->epslon);

		}
	}

}


void Tres::setdefaultsimatm() {
	
	Tres *s;
	Tatm *t;	
	
	if(simatm) delete simatm;
	simatm=0;

	int nid=0;
	SimAtm *ss=0;

	for(s=&TRES;s;s=s->next) {
                for(t=s->tatm;t;t=t->next) {
			t->sim=0;
			if(simatm==0) {
				simatm=new SimAtm(t->name[1],nid,t->eng->charge,t->eng->radius,t->eng->epslon);
				t->sim=simatm;
				ss=simatm;
				nid++;
				continue;
			}
			SimAtm *tt=simatm->ifexist(t->name[1],t->eng->charge,t->eng->radius,t->eng->epslon);
			
			if(tt==0) {
				ss->next=new SimAtm(t->name[1],nid,t->eng->charge,t->eng->radius,t->eng->epslon);
				ss=ss->next;
                                t->sim=ss;
                                nid++;	
				continue;
			}
			else {
				t->sim=tt;
				continue;
			}
                }
        }

	if(TRES.logg) cerr<<"the number of atm types in constraint optimization is:"<<nid<<endl;
}


void Tres::setsimplesimatm(char *ff) {

        Tres *s;
        Tatm *t;

        if(simatm) delete simatm;
        simatm=0;

        int nid=0;
        SimAtm *ss=0;

	float crg=0;
	float rad=0;
	float eps=0;
	char  c; 
        for(s=&TRES;s;s=s->next) {
                for(t=s->tatm;t;t=t->next) {
                        t->sim=0;
			
			//set charge 
			crg=0;rad=0;eps=0;
			if(t->id==0) {
				crg=0.3;
			}
			else if(t->id==3){
				crg=-0.3;
			}
			else if(t->name[1]=='O'&&strchr("DE",s->name)) {
				crg=-0.3;//crg=-0.5;		
			}
			else if(t->name[1]=='O'&&strchr("NQ",s->name)) {
                                crg=-0.3;               
                        }
			else if(t->name[1]=='N'&&strchr("NQ",s->name)) {
                                crg=0.3;               
                        }
			else if(strcmp(" ND1",t->name)==0&&strchr("H",s->name)) {
				crg=0.3; 
			}
			else if(strcmp(" NE2",t->name)==0&&strchr("H",s->name)) {
                                crg=-0.3; 
                        }
			else if(t->name[1]=='N'&&strchr("KWR",s->name)) {
                                crg=0.3;
                        }       
			//end

			c=t->name[1];
			//rad=int(t->eng->radius*10+0.51)/10.;
			//eps=-int(-t->eng->epslon*100+0.5)/100.;
			rad=t->eng->radius;
			eps=t->eng->epslon;
			if(strchr(ff,'c')==0) crg=0;
			if(strchr(ff,'r')==0) rad=0;
			if(strchr(ff,'e')==0) eps=0;
                        if(simatm==0) {
                                simatm=new SimAtm(c,nid,crg,rad,eps);
                                t->sim=simatm;
                                ss=simatm;
                                nid++;
                                continue;
                        }
                        SimAtm *tt=simatm->ifexist(c,crg,rad,eps);

                        if(tt==0) {
                                ss->next=new SimAtm(c,nid,crg,rad,eps);
                                ss=ss->next;
                                t->sim=ss;
                                nid++;
                                continue;
                        }
                        else {
                                t->sim=tt;
                                continue;
                        }
                }
        }

        if(TRES.logg) cerr<<"the number of atm types in constraint optimization is:"<<nid<<endl;
}

void Tres::setsimplecharge() {

        Tres *s;
        Tatm *t;
 
        //int nid=0;
      
	float crg=0;	 
	 
        for(s=&TRES;s;s=s->next) {
                for(t=s->tatm;t;t=t->next) {
                       
			//set charge 
			crg=0; 
			if(t->id==0) {
				crg=0.3;
			}
			else if(t->id==3){
				crg=-0.3;
			}
			else if(t->name[1]=='O'&&strchr("DE",s->name)) {
				crg=-0.3; //crg=-0.5;		
			}
			else if(t->name[1]=='O'&&strchr("NQ",s->name)) {
                                crg=-0.3;               
                        }
			else if(t->name[1]=='N'&&strchr("NQ",s->name)) {
                                crg=0.3;               
                        }
			else if(strcmp(" ND1",t->name)==0&&strchr("H",s->name)) {
				crg=0.3; 
			}
			else if(strcmp(" NE2",t->name)==0&&strchr("H",s->name)) {
                                crg=-0.3; 
                        }
			else if(t->name[1]=='N'&&strchr("KWR",s->name)) {
                                crg=0.3;
                        } 
			t->eng->charge=crg;			
                }
        }
}


void Tres::writesimatm(FILE *fp){
	
	for(Tres *t=&TRES;t;t=t->next) 

	for(Tatm *s=t->tatm;s;s=s->next) {

		fprintf(fp,"%s %s ",t->name3,s->name);
		s->sim->write(fp);
	}

}

void Tres::switchcharge(int f) {
//f=1, save charge0
//f=0, restore;
	for(Tres *t=&TRES;t;t=t->next)
	for(Tatm *s=t->tatm;s;s=s->next) {
		if(f) s->eng->charge=s->eng->charge0;   
		else  s->eng->charge0=s->eng->charge;
	}
}
void Tres::createengtable(int n){
        if(engtable) delete engtable;
        engtable=new EngTable();
	if(n==1) {
		engtable->name=strdup("vdw");
		engtable->setvdwtable();
	}
	/*
	else if(n==2) {
		engtable->name=strdup("vdw");
		engtable->setvdwtable();
		engtable->next=new EngTable();
		engtable->next->name=strdup("vdw6");
		engtable->next->setvdwtable();
	}
	*/
}

char *Tres::getsequence(){

	char *s=new char[100];

	int n=0;

	for(Tres *t=&TRES;t;t=t->next) {
		s[n++]=t->name;
	}
	s[n]='\0'; 
	return s;
}

char *Tres::getcappedsequence(char c){

        char *s=new char[100];

        int n=0;
	s[n++]=c;
        for(Tres *t=&TRES;t;t=t->next) {
                s[n++]=t->name;
        }
	s[n++]=c;
        s[n]='\0';
        return s;
}


char *Tres::getcappedsequence(char *c){

        char *s=new char[100];

	int m=strlen(c);
	int i;

        int n=0;
	for(i=0;i<m;i++) {
        	s[n++]=c[i];
	}

        for(Tres *t=&TRES;t;t=t->next) {
                s[n++]=t->name;
        }

	for(i=0;i<m;i++) {
                s[n++]=c[i];
        }

        s[n]='\0';
        return s;
}

void  Tres::setengcoeff(int n){

	if(engcoeff) delete engcoeff;

	engcoeff=new EngCoeff();

	if(n==1) {
		engcoeff->setdots(0.7,55.247845,1);
        	engcoeff->setdots(0.8,6.922596,1);
        	engcoeff->setdots(0.9,-0.222634,1);
        	engcoeff->setdots(1,-1,1);
        	engcoeff->setdots(1.1,-0.81,1);
        	engcoeff->setdots(1.6,-0.115656,1);
        	engcoeff->setdots(9.,-0.00001,1);
        	engcoeff->setcoeff(0,0.5);
        	engcoeff->printout();
	}
	else if(n==2) {
		engcoeff->setdots(0.7,55.247845,1);
                engcoeff->setdots(0.8,6.922596,1);
                engcoeff->setdots(0.9,-0.222634,1);
                engcoeff->setdots(1,-1,1);
                engcoeff->setdots(1.1,-0.81,1);
                engcoeff->setdots(1.6,-0.115656,1);
                engcoeff->setdots(9.,-0.00001,1);
                engcoeff->setcoeff(0,1);
                engcoeff->printoutorder(0);
	}

}

void  Tres::setengcoeff(float buffer,float npower){

        if(engcoeff) delete engcoeff;

        engcoeff=new EngCoeff();
	engcoeff->setdots(0.7,55.247845,1);
        engcoeff->setdots(0.8,6.922596,1);
        engcoeff->setdots(0.9,-0.222634,1);
        engcoeff->setdots(1,-1,1);
        engcoeff->setdots(1.1,-0.81,1);
        engcoeff->setdots(1.6,-0.115656,1);
        engcoeff->setdots(9.,-0.00001,1);
        engcoeff->setcoeff(buffer,npower);
        engcoeff->printout();
}

void  Tres::setnewengcoeff(float buffer,float npower){

        if(engcoeff) delete engcoeff;
        engcoeff=new EngCoeff();
        engcoeff->setdots(0.7,55.247845,1);
        engcoeff->setdots(0.8,6.922596,1);
        engcoeff->setdots(0.9,-0.222634,1);
        engcoeff->setdots(1,-1,1);
        engcoeff->setdots(1.1,-0.81,1);
        engcoeff->setdots(1.6,-0.115656,1);
        engcoeff->setdots(9.,-0.00001,1);
        engcoeff->setnewcoeff(buffer,npower);
        if(logg>3)engcoeff->printnewout();
}


EngTable *Tres::findengtable(char *s) {

    for(EngTable *c=engtable;c;c=c->next) {
        if(c->name&&strcmp(c->name,s)==0) return c;
    }  

    return 0;
}

ConTable *Tres::findcontable(char *s) {

    for(ConTable *c=contable;c;c=c->next) {
        if(c->name&&strcmp(c->name,s)==0) return c;
    }

    return 0;
}

void Tres::setcontable() {
	if(contable) {
		ConTable *c=contable->setsqrt(400,100);
		contable->add(c);

	}
	else {
		ConTable *c=contable->setsqrt(400,100);
		contable=c;
	}
}

void Tres::prepareengtable(){

	for(EngTable *s=engtable;s;s=s->next) {
		s->preparecall();
	}
}

void Tres::calcbackrotamertorsion(char *s){

        Chn *brt=TRES.findrotamer(s)->chn;
        char nm[4];

        if(brt==0) return;

        Rotate rot;
        Res *r,*t;
        nm[0]='G';
        nm[2]='G';
        nm[3]='\0';
        for(r=brt->res;r;r=r->next) {
                nm[1]=r->name;
                Chn *c=new Chn();
                c->create(nm);
                c->configure();
                //c->header();
                for(t=r->more;t;t=t->more) {
                        c->res->next->transfer(t);
                        c->res->next->transfertemp(t);
                        Res *e=c->res;
                        rot.link(e,e->next,e->id0,0);
                        rot.link(e->next,e->next->next,e->next->next->id0,1);
                        //c->write("out");
                        e->next->setbackbonetorsion();
                        t->atm->next->chi=e->next->atm->next->chi;
                        t->atm->next->next->chi=e->next->atm->next->next->chi;
                        if(TRES.logg) cerr<<r->name<<" "<<t->atm->next->chi<<" "<<t->atm->next->next->chi<<endl;
                }
                delete c;
        }
}

void Tres::set3backrotamer(char *s){

	Chn *brt=TRES.findrotamername(s)->chn;
	char nm[4];

	if(brt==0) return;
	
	Rotate rot;
	Res *r,*t;
	nm[0]='G';
	nm[2]='G'; 
	nm[3]='\0';
	for(r=brt->res;r;r=r->next) {
		nm[1]=r->name;
		Chn *c=new Chn();
		c->create(nm);
		c->configure();
		//c->header();
		for(t=r->more;t;t=t->more) {
			c->res->next->transfer(t);
			c->res->next->transfertemp(t);
			Res *e=c->res;
			rot.link(e,e->next,e->id0,0);
			rot.link(e->next,e->next->next,e->next->next->id0,1);
			//c->write("out");
			/*
			e->next->setbackbonetorsion();
			t->atm->next->chi=e->next->atm->next->chi;
			t->atm->next->next->chi=e->next->atm->next->next->chi;
			cerr<<r->name<<" "<<t->atm->next->chi<<" "<<t->atm->next->next->chi<<endl;
			*/
			e->next->dihedral(0);
			Atm *a;
			for(a=t->atm;a;a=a->next) {
				Atm *b=e->next->isatmid(a->tatm->id);
				if(b==0) continue;
				a->chi=b->chi;
			}
		}
		delete c;
	}	
	
	
	if(1) return;

	int n=0;
	for(r=brt->res;r;r=r->next) {
		if(r->nummore>n) n=r->nummore;
	}

	if(n==0) return;

	Res **tmp=new Res*[2*n];
	int *order=new int[2*n];
	float *value=new float[2*n];
	Res **nemp=new Res*[2*n];

	Qsort cc;

	for(Res *rr=brt->res;rr;rr=rr->next) {
		n=0;
		r=rr;
		for(t=r->more;t;t=t->more) {
			tmp[n++]=t;
		}
		tmp[n]=0;

		int i;

		
		for(i=0;i<n;i++) {
			r=tmp[i];
			
			if(fabs(r->atm->next->chi+60)<25&&fabs(r->atm->next->next->chi+43)<25) {
				value[i]=0;
			}
			else if(fabs(r->atm->next->chi+100)<40&&fabs(r->atm->next->next->chi-130)<40) {
                                value[i]=1;
                        }
			else {
				value[i]=2;
			}
		}
		cc.sort(value,n,order);

		int j;
		for(i=0;i<n;i++) {
			j=order[i];
			nemp[i]=tmp[j];	
		}

		for(i=0;i<n;i++){
			nemp[i]->more=0;
		}

		t=rr;
		for(i=0;i<n;i++) {
			t->more=nemp[i];
			t=t->more;		
		}
	}
	if(nemp) delete [] nemp;
	if(tmp) delete [] tmp;
	if(value) delete [] value;
	if(order) delete [] order;
}

void Tres::setdipolecharge() {

	Tres *s;
	Tatm *a;

	for(s=&TRES;s;s=s->next)  
	for(a=s->tatm;a;a=a->next)
	a->eng->charge=0;

	for(s=&TRES;s;s=s->next) {		
		for(a=s->tatm;a;a=a->next) {
			if(s->name!='P'&&strcmp(a->name," N  ")==0) {
				a->eng->dipole=0.5;//0.3;				 
			}			
			else if(strcmp(a->name," O  ")==0) {
				a->eng->dipole=0.5;				 
			}
			else if(s->name=='D'&&strcmp(a->name," CG ")==0) {
				a->eng->charge=-1;
			}
			else if(s->name=='E'&&strcmp(a->name," CD ")==0) {
				a->eng->charge=-1;
			} 
			else if(s->name=='K'&&strcmp(a->name," NZ ")==0) {
				a->eng->charge=1;
			}
			else if(s->name=='R'&&strcmp(a->name," CZ ")==0) {
				a->eng->charge=1;
			}			
			else if(s->name=='Q'&&strcmp(a->name," NE2")==0) {
				a->eng->dipole=0.5;//0.4;
			}
			else if(s->name=='Q'&&strcmp(a->name," OE1")==0) {
				a->eng->dipole=0.5;
			} 
			else if(s->name=='N'&&strcmp(a->name," OD1")==0) {
				a->eng->dipole=0.5;
			} 
			else if(s->name=='N'&&strcmp(a->name," ND2")==0) {
				a->eng->dipole=0.5;//0.3;
			} 			
			else if(s->name=='R'&&strcmp(a->name," NE ")==0) {
				a->eng->dipole=0.5;//0.4;
			}
			else if(s->name=='R'&&strcmp(a->name," NH1")==0) {
				a->eng->dipole=0.;//0.4;
			}
			else if(s->name=='R'&&strcmp(a->name," NH2")==0) {
				a->eng->dipole=0.;//0.4;
			}
			else if(s->name=='H'&&strcmp(a->name," ND1")==0) {
				a->eng->dipole=0.5;//0.3;
			}
			else if(s->name=='H'&&strcmp(a->name," NE2")==0) {
				a->eng->dipole=0.;
			}
			else if(s->name=='W'&&strcmp(a->name," NE1")==0) {
				a->eng->dipole=0.5;//0.4;
			} 					
			else if(s->name=='T'&&a->name[1]=='O') {
				a->eng->dipole=0;//0.4;
			}
			else if(s->name=='S'&&a->name[1]=='O') {
				a->eng->dipole=0;//0.4;
			}
			else if(s->name=='Y'&&a->name[1]=='O') {
				a->eng->dipole=0;//0.4;
			}					
		}		
	}
}


float *Tres::getdipolevector(Atm *a,int fh) {
	float *xyz;
	int i;
 	Atm *b=0;
	Atm *a1,*a2;
	Res *s=a->res;
	if(s==0) return 0;
	if(a->res->last==0&&a->tatm->id==0) {//it is n terminal charge
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a->xyz[i];
		return xyz;	
	}
	else if(a->res->last==0&&a->tatm->id==3) {//it is c terminal
		b=s->isatmid(2);
		if(b==0) return 0;
		Atm *a1=s->isatmid(1);
		if(a1==0) return 0;
		xyz=new float[6];		
		for(i=0;i<3;i++) xyz[i]=b->xyz[i]+(b->xyz[i]-a1->xyz[i])/2;
		return xyz;
	}
	else if(a->tatm->id==0) {
		b=a->ischildatm(" HN ");
		if(b==0) {			
			Atm *t1=s->isatmid(1); 	
			Atm *t2=s->last->isatmid(2);
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t1->xyz[i]+t2->xyz[i])/2;
			a1=b;a2=a; 							
		}
		else {
			a1=a;a2=b;
		} 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.01) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
		}
		if(b->res==0) delete b;b=0;
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2; 
		return xyz;
        }		
	else if(a->tatm->id==3) {
		b=a->res->isatm(" C  ");
		if(b==0) return 0;
		a1=a;a2=b;
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
		}	
		for(i=3;i<6;i++) xyz[i]=(a1->xyz[i-3]+a2->xyz[i-3])/2; 
		return xyz;	 
	}
	else if(s->name=='K'&&strcmp(a->name," NZ ")==0) {
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a->xyz[i]; 
		return xyz;
	}			
        else if(s->name=='R'&&strcmp(a->name," CZ ")==0) {
		Atm *a1=s->isatm(" NH1");
		Atm *a2=s->isatm(" NH2");
		if(a1==0||a2==0) return 0;
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=(a->xyz[i]+a2->xyz[i]+a1->xyz[i])/3;
		return xyz;		
	}
	else if(s->name=='D'&&strcmp(a->name," CG ")==0) {
		Atm *a1=a->ischildatm(" OD1");
		Atm *a2=a->ischildatm(" OD2");
		if(a1==0||a2==0) return 0;
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=(a->xyz[i]+a2->xyz[i]+a1->xyz[i])/3;
		return xyz;		
	}
	else if(s->name=='E'&&strcmp(a->name," CD ")==0) {
		Atm *a1=a->ischildatm(" OE1");
		Atm *a2=a->ischildatm(" OE2");
		if(a1==0||a2==0) return 0;
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=(a->xyz[i]+a2->xyz[i]+a1->xyz[i])/3;
		return xyz; 
	} 	
	else if(s->name=='W'&&strcmp(a->name," NE1")==0) {
		b=a->res->isatm(" HE1");
		if(b==0) {			
			Atm *t1=a->res->isatm(" CD1"); 	
			Atm *t2=a->res->isatm(" CE2");
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t1->xyz[i]+t2->xyz[i])/2;
			a1=b;a2=a; 							
		}
		else {
			a1=a;a2=b;
		} 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.01) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
		}
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		if(b->res==0) delete b;b=0;
		return xyz;
	} 
	else if(s->name=='Q'&&strcmp(a->name," NE2")==0&&fh==1) {
		b=a->ischildatm("1HE2");
		if(b==0){
			Atm *t1=a->res->isatm(" CG ");
			Atm *t2=a->res->isatm(" CD "); 				
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t2->xyz[i]-t1->xyz[i])+a->xyz[i];
			a1=a;a2=b; 				
		}
		else {
			a1=a;a2=b;
		}		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		if(b->res==0) delete b;b=0;
		return xyz;
	}	
	else if(s->name=='Q'&&strcmp(a->name," NE2")==0&&fh==2) {
		b=a->ischildatm("2HE2");
		if(b==0){
			Atm *t1=a->res->isatm(" OE1");
			Atm *t2=a->res->isatm(" CD "); 	
			
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t2->xyz[i]-t1->xyz[i])+a->xyz[i];
			a1=a;a2=b; 				
		}
		else {
			a1=a;a2=b;
		}		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		if(b->res==0) delete b;b=0;
		return xyz;
	}	
	else if(s->name=='Q'&&strcmp(a->name," NE2")==0) {
		b=a->res->isatm(" CD ");		
		if(b==0) return 0;
		a1=b;a2=a;
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;return 0;
		}
		
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;   
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		return xyz; 
	} 	
	else if(s->name=='Q'&&strcmp(a->name," OE1")==0) {
		b=a->res->isatm(" CD ");		
		if(b==0) return 0;
		a1=a;a2=b;
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;return 0;
		}
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;   
		for(i=3;i<6;i++) xyz[i]=(a1->xyz[i-3]+a2->xyz[i-3])/2; 
		return xyz; 
	} 
	
	else if(s->name=='N'&&strcmp(a->name," ND2")==0&&fh==1) {
		b=a->ischildatm("1HD1");
		if(b==0){
			Atm *t1=a->res->isatm(" CB ");
			Atm *t2=a->res->isatm(" CG "); 	
			
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t2->xyz[i]-t1->xyz[i])+a->xyz[i];
			a1=a;a2=b; 				
		}
		else {
			a1=a;a2=b;
		}		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2; 
		if(b->res==0) delete b;b=0;
		return xyz;
	} 
	else if(s->name=='N'&&strcmp(a->name," ND2")==0&&fh==2) {
		b=a->ischildatm("2HD2");
		if(b==0){
			Atm *t1=a->res->isatm(" OD1");
			Atm *t2=a->res->isatm(" CG "); 				
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t2->xyz[i]-t1->xyz[i])+a->xyz[i];
			a1=a;a2=b; 				
		}
		else {
			a1=a;a2=b;
		}		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		if(b->res==0) delete b;b=0;
		return xyz;
	} 	
	else if(s->name=='N'&&strcmp(a->name," OD1")==0) {
		b=s->isatm(" CG ");
		a1=a;a2=b;
		if(a2==0) return 0;
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;return 0;
		}
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;   
		for(i=3;i<6;i++) xyz[i]=(a1->xyz[i-3]+a2->xyz[i-3])/2;
		return xyz;
	} 
	else if(s->name=='N'&&strcmp(a->name," ND2")==0) {
		b=s->isatm(" CG ");
		a1=b;a2=a;
		if(a2==0) return 0;
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;return 0;
		}
		for(i=0;i<3;i++) xyz[i]=xyz[i]/t;  
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2; 
		return xyz;
	} 
	else if(s->name=='R'&&strcmp(a->name," NE ")==0) {
		b=a->res->isatm(" HE ");
		if(b==0) {			
			Atm *t1=a->res->isatm(" CD "); 	
			Atm *t2=a->res->isatm(" CZ ");
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t1->xyz[i]+t2->xyz[i])/2;
			a1=b;a2=a; 							
		}
		else {
			a1=a;a2=b;
		} 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.01) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
		}
		if(b->res==0) delete b;b=0;
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2; 
		return xyz;
	}	
	
	else if(s->name=='R'&&strcmp(a->name," NH1")==0&&fh==1) {
		b=a->ischildatm("1HH1");
		if(b==0){
			Atm *t1=a->res->isatm(" NH2"); 	
			Atm *t2=a->res->isatm(" CZ ");
						
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t2->xyz[i]-t1->xyz[i])+a->xyz[i];
			a1=a;a2=b; 				
		}
		else {
			a1=a;a2=b;
		}		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		if(b->res==0) delete b;b=0;
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		return xyz;	
	}
	else if(s->name=='R'&&strcmp(a->name," NH1")==0&&fh==2) {
		b=a->ischildatm("1HH1");
		if(b==0){
			Atm *t1=a->res->isatm(" NE "); 	
			Atm *t2=a->res->isatm(" CZ ");
						
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t2->xyz[i]-t1->xyz[i])+a->xyz[i];
			a1=a;a2=b; 				
		}
		else {
			a1=a;a2=b;
		}		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		if(b->res==0) delete b;b=0;
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		return xyz;	
	}
	else if(s->name=='R'&&strcmp(a->name," NH2")==0&&fh==1) {
		b=a->ischildatm("1HH2");
		if(b==0){
			Atm *t1=a->res->isatm(" NH1"); 	
			Atm *t2=a->res->isatm(" CZ ");
						
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t2->xyz[i]-t1->xyz[i])+a->xyz[i];
			a1=a;a2=b; 				
		}
		else {
			a1=a;a2=b;
		}		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		if(b->res==0) delete b;b=0;
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		return xyz;	
	}
	else if(s->name=='R'&&strcmp(a->name," NH2")==0&&fh==2) {
		b=a->ischildatm("2HH2");
		if(b==0){
			Atm *t1=a->res->isatm(" NE "); 	
			Atm *t2=a->res->isatm(" CZ ");
						
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t2->xyz[i]-t1->xyz[i])+a->xyz[i];
			a1=a;a2=b; 				
		}
		else {
			a1=a;a2=b;
		}		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		if(b->res==0) delete b;b=0;
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		return xyz;	
	}
	else if(s->name=='H'&&strcmp(a->name," ND1")==0) {
		b=a->ischildatm(" HD1");
		if(b==0) {			
			Atm *t1=a->res->isatm(" CG "); 	
			Atm *t2=a->res->isatm(" CE1");
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t1->xyz[i]+t2->xyz[i])/2;
			a1=b;a2=a; 							
		}
		else {
			a1=a;a2=b;
		} 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.01) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
		}
		if(b->res==0) delete b;b=0;
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		return xyz;
	}
	else if(s->name=='H'&&strcmp(a->name," NE2")==0) {
		b=a->ischildatm(" HE2");
		if(b==0) {			
			Atm *t1=a->res->isatm(" CD2"); 	
			Atm *t2=a->res->isatm(" CE1");
			if(t1==0||t2==0) return 0;
			b=new Atm;
			for(i=0;i<3;i++)b->xyz[i]=(t1->xyz[i]+t2->xyz[i])/2;
			a1=b;a2=a; 							
		}
		else {
			a1=a;a2=b;
		} 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.01) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t;
		}
		if(b->res==0) delete b;b=0;
		for(i=3;i<6;i++) xyz[i]=a->xyz[i-3]+xyz[i-3]/2;
		return xyz;
	}
	else if(s->name=='T'&&a->name[1]=='O') {		
		b=a->res->isatm(" CB ");
		if(b==0) return 0;
		a1=a;a2=b;		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		for(i=3;i<6;i++) xyz[i]=(a1->xyz[i-3]+a2->xyz[i-3])/2;
		return xyz;
	}
	else if(s->name=='S'&&a->name[1]=='O') {
		b=a->res->isatm(" CB ");
		if(b==0) return 0;
		a1=a;a2=b;		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}
		for(i=3;i<6;i++) xyz[i]=(a1->xyz[i-3]+a2->xyz[i-3])/2;
		return xyz;
	}
	else if(s->name=='Y'&&a->name[1]=='O') {
		b=a->res->isatm(" CZ ");
		if(b==0) return 0;
		a1=a;a2=b;		 
		xyz=new float[6];
		for(i=0;i<3;i++) xyz[i]=a2->xyz[i]-a1->xyz[i];
		float t=distance(xyz);
		if(t<0.1) {
			delete [] xyz;xyz=0;
		}
		else {
			for(i=0;i<3;i++) xyz[i]=xyz[i]/t; 
		}	
		for(i=3;i<6;i++) xyz[i]=(a1->xyz[i-3]+a2->xyz[i-3])/2;
		return xyz;
	}
	return 0;		
}

void Tres::printsurfacehelp(){
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"surface is a program of calculating protein surface\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage:\n");
fprintf(stderr,"surface -prm num -probe float -size num -out str file.pdb\n");
fprintf(stderr,"-prm: set parameter set.\n"); 
fprintf(stderr,"-prm 1: Charmm22 all atom model\n");
fprintf(stderr,"-prm 2: Amber94  all atom model\n");
fprintf(stderr,"-prm 3: Charmm22 heavy atom model\n");
fprintf(stderr,"-prm 4: Amber94  heavy atom model\n");
fprintf(stderr,"-probe: set probe radius. default is 1.4 A\n");
fprintf(stderr,"-size: number of subareas each atom surface subdivided.\n");
fprintf(stderr,"-out output file.default is standard screen output\n");

}

void Tres::printaddhhelp(){

fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"addh is a program of adding protons\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
printf("Usage:\n");
printf("addh -k file.pdb\n");
printf("-k: keep original hydrogen atoms.\n");
printf("-k 1:keep original protons.default is 1\n");
printf("-k 0:delete original protons\n");
}

void Tres::printchihelp(){
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"chi is a program to print out chi angles of a protein\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage:\n");
fprintf(stderr,"chi -a num -o str file.pdb\n");
fprintf(stderr,"-a num: options for torsion angle output.default is 0\n");
fprintf(stderr,"-a   0: print out all torsion angles.\n");
fprintf(stderr,"-a   1: print out only backbone phi,psi and sidechain torsion angles.\n");
fprintf(stderr,"-a   2: print out only backbone phi,psi torsion angles.\n");
fprintf(stderr,"-a   3: print out only sidechain torsion angles.\n");
fprintf(stderr,"-a   4: print out only backbone omega torsion angles.\n");
fprintf(stderr,"-o output file.default is standard screen output\n");
}

void Tres::printalgnhelp(){

printf("Usage:\n");
printf("algn -f num -o output.pdb file.pdb\n");
printf("-f num: options for algnment two structures.default is 0\n");
printf("-f   0: ca only.\n");
printf("-f   1: backbone only.\n");
printf("-f   2: all heavy atoms.\n");
printf("-f   3: all atoms including protons.\n");
printf("-r   0: not printout residue rmsd. default\n");
printf("-r   1: printout residue rmsd\n");
printf("-o output file. \n");
}
void Tres::printnalgnhelp(){
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"nalgn is a program to superimpose proten structure based on the given\n");
fprintf(stderr,"sequence alignment\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
printf("Usage:\n");
printf("nalgn -f num -o output.pdb file.pir\n");
printf("-f num: options for algnment two structures.default is 0\n");
printf("-f   0: ca only.\n");
printf("-f   1: backbone only.\n");
printf("-f   2: all heavy atoms.\n");
printf("-f   3: all atoms including protons.\n");
printf("-o output file. \n");
}
void Tres::printhbondhelp(){
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"hbond is a program to print out hydrogen bond and disulfide bonds of a protein\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage:\n");
fprintf(stderr,"hbond -h num -o str file.pdb\n");
fprintf(stderr,"-h num: options for output hydrogen bond.default is 0\n");
fprintf(stderr,"-h   0: hydrogen bond plus s-s bond.\n");
fprintf(stderr,"-h   1: hydrogen bond only.\n");
fprintf(stderr,"-h   2: s-s bond only.\n");
fprintf(stderr,"-o output file.optional. default is standard screen output\n");
}

void Tres::printseqxyzhelp(){

printf("Usage:\n");
printf("seqxyz -f num file.pir\n");
printf("-f num: options for atom mdoel. default is 0\n");
printf("-f   0: all atom model.\n");
printf("-f   1: heavy atom model.\n");
}

void Tres::printxyzseqhelp(){
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"pdbseq is a program to print out sequnece information of a protein  in pir format\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage:\n");
fprintf(stderr,"pdbseq -a num -c char -o str file.pdb\n");
fprintf(stderr,"-a num: options for reading sequence card. default is 0\n");
fprintf(stderr,"-a   0: read sequence from sequence card in pdb.\n");
fprintf(stderr,"-a   1: read sequence from ATOM record including gap with X.\n");
fprintf(stderr,"-a   2: read sequence from ATOM record including no gap.\n");
fprintf(stderr,"-c char: one character for the chain to be read.default is all chains.\n");
fprintf(stderr,"-o output file.optional. default is standard screen output\n");
}

void Tres::printchnparserhelp(){
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"chnout is a program to print out individual chain coordinates\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage:\n");
fprintf(stderr,"chnout -c char -o str -i num file.pdb\n");
fprintf(stderr,"-c char: one character for the chain to be read.default is all chains.\n");
fprintf(stderr,"-o output file.optional. default is standard screen output\n");
fprintf(stderr,"-i the id of the first residue in the chain\n");
}

void Tres::printchixyzhelp(){
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"rotate is a program to rotate dihedral angles\n");
fprintf(stderr,"questions? please refer to Dr.Jason Z. Xiang at zx11@columbia.edu\n");
fprintf(stderr,"***************************************************************************\n");
fprintf(stderr,"Usage:\n");
fprintf(stderr,"rotate -d num -e num -t name -n num -r num file.pdb\n");
fprintf(stderr,"-t name: chiangle file name. required.\n");
fprintf(stderr,"-d num:  the direction of rotation. default is 1\n");
fprintf(stderr,"-d 1:  the rotation forward. \n");
fprintf(stderr,"-d 0:  the rotation backward. \n");
fprintf(stderr,"-e num:  the end residue where rotation affects. default is to the end\n");
fprintf(stderr,"-n num:  the content of chiangle file. \n");
fprintf(stderr,"-n 0:  the chiangle file includes omega, phi,psi, side-chain torsion angle\n");
fprintf(stderr,"-n 1:  the chiangle file includes  phi,psi, side-chain torsion angle\n");
fprintf(stderr,"-n 2:  the chiangle file includes  phi,psi\n"); 
fprintf(stderr,"-n 3:  the chiangle file includes side-chain torsion angle\n"); 
fprintf(stderr,"-r num: the option to rotate the bond.default 1\n"); 
fprintf(stderr,"-r 0: rotate the bond of degrees specified in chi file.\n");
fprintf(stderr,"-r 1: rotate the bond so that its torsion angle is equal to that specified in chi file.\n");
}

void Tres::printezchixyzhelp(){
printf("Usage:\n");
printf("ezrot -a id -b float -r num file.pdb\n");
printf("-a id: atom id. the bond is defined between the atom and its parent. required.\n");
printf("-b float: rotation angle. required.\n");
printf("-r num: the option to rotate the bond.default 1\n");
printf("-r 0: rotate the bond of degree specified in -b.\n");
printf("-r 1: rotate the bond so that its torsion angle is equal to that specified in -b.\n");
}

Tres *Tres::updatehetm(char *s) {

	Tres *t;

	for(t=this;t;t=t->next) {

		if(strncmp(t->name3,s,3)==0) return t;
	}

	for(t=this;t->next;t=t->next);
	t->next=new Tres;
	t=t->next;
	strncpy(t->name3,s,3);
	t->name3[3]='\0';
	return t;
}
Tres *Tres::addhetm(char *s) {

	Tres *t;
	for(t=this;t->next;t=t->next);
	t->next=new Tres;
	t=t->next;
	strncpy(t->name3,s,3);
	t->name3[3]='\0';
	return t;
}

Tatm *Tres::updatehetmatm(char *s) {

        Tatm *t;

        for(t=tatm;t;t=t->next) {

                if(strncmp(t->name,s,4)==0) return t;
        }

        for(t=tatm;t->next;t=t->next);
        t->next=new Tatm;
        t=t->next;
        strncpy(t->name,s,4);
	t->name[4]='\0';
        return t;
}
Tatm *Tres::addhetmatm(char *s) {

        Tatm *t;
	if(tatm==0) {
		tatm=new Tatm;
		t=tatm;
	}
        else {
		for(t=tatm;t->next;t=t->next);
		t->next=new Tatm;
        	t=t->next;
	}       
        strncpy(t->name,s,4);
	t->name[4]='\0';
        return t;
}
void Tres::setenghetatm(char *filn){

	if(hetm==0) return;
	char *mlib,line[1000];
 	int i;
	FILE *fp;

	mlib=RCS["library"];
	if(mlib==0) {
        	cerr<<"library does not exist in jackal.dir"<<endl;
        	cerr<<"check the jackal.dir file first in your current directory then in the environment path:JACKALDIR"<<endl;
        	cerr<<"something wrong in specifying the library"<<endl;
        	cerr<<"make sure the library points to the abosulte path"<<endl;
        	exit(0);
  	}
	fp=0;i=0;
  	//cerr<<filnam<<mlib+2<<endl;
  	sprintf(line,"%s/%s",mlib,filn);
  	i=strlen(mlib);
  	while((fp=fopen(line,"r"))==NULL)
  	{
   		mlib=mlib+i+1;
   		i=strlen(mlib);
   		if(*mlib=='\0')
   		{
        		cerr<<"warning:could not open: "<<filn<<endl;
        		cerr<<"check the jackal.dir file first in your current directory then in the environment path:JACKALDIR"<<endl;
        		cerr<<"something wrong in specifying the library"<<endl;
        		cerr<<"make sure the library points to the abosulte path"<<endl;
        		exit(0);
   		}
   		sprintf(line,"%s/%s",mlib,filn);
  	}
	char linef[1000];
	strcpy(linef,line);

	cerr<<endl<<"opening file: "<<linef<<endl;
	if(fp==0) {
		cerr<<"warning:could not open: "<<filn<<endl;
        	cerr<<"check the jackal.dir file first in your current directory then in the environment path:JACKALDIR"<<endl;
        	cerr<<"something wrong in specifying the library"<<endl;
        	cerr<<"make sure the library points to the abosulte path"<<endl;
        	exit(0);
	}
	 
	char *lines[10000];
 	Strhandler cc;
	for(i=0;i<10000;i++) lines[i]=0;

	i=0;
	while(fgets(line,200,fp)!=NULL)	 {
		if(i>9900) {
			cerr<<"too many lines!ignore..."<<line<<endl;
			continue;
		}
		if(strlen(line)<16) {
			cerr<<"bad line!ignore..."<<line<<endl;
			continue;
		}
		if(line[0]=='!') continue;
		char *s=strdup(line);
		s=cc.uppercase(s);
		if(s==0)  continue;
		s=cc.clearendchar(s,"\n");
		if(s==0)  continue;
		if(strchr(s,'\t')) {
			cerr<<"line containing tab.ignored: "<<s<<endl;
			cerr<<"the line from in the library file: "<<linef<<endl;
			cerr<<"which is used to assign force field parameters to hetatms"<<endl;
			continue;
		}
		lines[i++]=s;
		
	}
	int nline=i;

	fclose(fp);


	char *limes[10000];
	float  temp[10000];
	int    order[10000];
	 
	
	for(i=0;i<nline;i++) {
		int n=strlen(lines[i]);		 	
		int j;
		temp[i]=0; 
		for(j=12;j<17&&j<n;j++) {
			if(lines[i][j]!=' ') {
				temp[i]+=-100000;break;
			}
		}
			
		for(j=7;j<10&&j<n;j++) {
			if(lines[i][j]!=' ') {
				temp[i]+=-10000;break;
			}
		} 

		for(j=0;j<5&&j<n;j++) {
			if(lines[i][j]!=' ') {
				temp[i]+=-(int)lines[i][j]; 
			}
		} 
	}
	
	Qsort sort;

	sort.sort(temp,nline,order);

	for(i=0;i<nline;i++) {
		int n=order[i];
		limes[i]=lines[n];
	}

	for(i=0;i<nline;i++) lines[i]=limes[i];

	

	Tres *t;
	Tatm *a;
 	

	for(t=hetm;t;t=t->next) 
	for(a=t->tatm;a;a=a->next) {
 
		if(a->eng==0) {
			a->eng=new Eng(a);
		}
		for(i=0;i<nline;i++) {
			int n=strlen(lines[i]);		 	
			int j;
			
			//residue id
			char resid[10];
			int nresid=0;
			for(j=12;j<17&&j<n;j++) {
				if(lines[i][j]!=' '&&lines[i][j]!='\t') {
					resid[nresid++]=lines[i][j]; 					
				}
			}
			resid[nresid++]='\0';
			cc.clearbadascii(resid);
			nresid=strlen(resid);
			if(nresid) {
				if(atoi(resid)!=t->id) continue;
			}
			
			//residue name
			char res[10];
			int nres=0;
			for(j=7;j<10&&j<n;j++) {
				if(lines[i][j]!=' '&&lines[i][j]!='\t') {
					res[nres++]=lines[i][j]; 
				}
			} 
			res[nres++]='\0';
			cc.clearbadascii(res);
			nres=strlen(res);
			if(nres) {
				char *s=strdup(t->name3);
				cc.clearbadascii(s);			
				if(s&&strcmp(s,res)) {
					s=cc.strdel(s);
					continue;
				}	
				s=cc.strdel(s);	
			}

			char atm[10];
			int natm=0;
			for(j=0;j<5&&j<n;j++) {
				if(lines[i][j]!=' '&&lines[i][j]!='\t') {
					atm[natm++]=lines[i][j];	 
				}
			} 
		 	atm[natm++]='\0';
			cc.clearbadascii(atm);
			natm=strlen(atm);
			if(natm) {
				char *s=strdup(a->name);
				cc.clearbadascii(s); 			
				if(s&&strcmp(s,atm)) {
					s=cc.strdel(s);
					continue;
				}
				s=cc.strdel(s);	
			}	
			
			 
			if(n>19)a->eng->charge=atof(lines[i]+19);
			if(n>27)a->eng->radius=atof(lines[i]+27);
			if(n>35)a->eng->epslon=atof(lines[i]+35);
 			// got=1;
			break;				 			 			 
		}
		if(a->eng->radius<0.01) {
			if(a->name[1]=='C') a->eng->radius=1.7;
			else if(a->name[1]=='O') a->eng->radius=1.6;
			else if(a->name[1]=='N') a->eng->radius=1.65;
			else if(a->name[1]=='H') a->eng->radius=1.00;
			else if(a->name[1]=='P') a->eng->radius=1.90;
			else if(a->name[1]=='S') a->eng->radius=1.90;
			else a->eng->radius=1.40;
		} 
		if(fabs(a->eng->epslon)<0.0001) {
			if(a->name[1]=='C') a->eng->epslon=-0.11;
			else if(a->name[1]=='O') a->eng->epslon=-0.12;
			else if(a->name[1]=='N') a->eng->epslon=-0.2;
			else if(a->name[1]=='H') a->eng->epslon=-0.04;
			else if(a->name[1]=='P') a->eng->epslon=-0.45;
			else if(a->name[1]=='S') a->eng->epslon=-0.45;
			else a->eng->radius=-0.04;
		}
	}
	
	for(i=0;i<nline;i++) {
		if(lines[i]) delete [] 	lines[i];
		lines[i]=0;
	}
}
