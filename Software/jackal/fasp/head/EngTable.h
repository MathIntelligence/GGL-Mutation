#ifndef _EngTable
#define _EngTable

class EngTable
{
public:

//function

EngTable();
EngTable(EngTable *);
~EngTable();
void initial();
void preparecall();
void setvdwtable(int a,int b);
void setvdwtable();
int getnumber();
float gettablevalue(int,int,float);
float getenergyvalue(int,int,float);
EngTable **getallengtable();
float calcvdwvalue(float,float,float,float);
float calcvdwtable(float,float,float,float);
//other members
float getvalue(int,int,float);
float getenergy(int,int,float);
float getvalue(float,float);
float getenergy(float,float);
void createdefwalltable();
void createdefwalltable(float);
void createcolonywalltable();
void createcolonywalltable(float);
void createwalltable();
void createwalltable(float);
void setchargezero();
void createsoftvdwtable();
void createsoftvdwtable(int,int);
void createdefsoftvdwtable();
void createdefsoftvdwtable(int,int);
void createcolonysoftvdwtable();
void createcolonysoftvdwtable(int,int);
void createvdwtable();
void createvdwtable(int,int);
void createdeftable();
void createdeftable(int,int);
void createcolonytable(); //table for pair atom interactions of free energy
void createcolonytable(int,int); //table for pair atom interactions
void chargecurve(FILE *fp,float);
void vdwcurve(FILE *fp,float);

//data

float      distance;    //distance
float      resolution;  //the steps
char      *name;        //the name of the table
float     *table;       //actual table value
float     *energy;
float      radius;      //radius of the two atoms
float      charge;      //charge product of the two atoms
float      epslon;      //epslon of the two atoms
int	   numsim;      //total number of the atom type
int        size;        //size of the array for energy and table
int	   id;
EngTable **engtable;    //all compiles of table
EngTable  *more;
EngTable  *next;      //the next table

//other data
float      expand;      //the gauss distribution coeff
float      wall;        //the high of wall between distance constraint
};
#endif
