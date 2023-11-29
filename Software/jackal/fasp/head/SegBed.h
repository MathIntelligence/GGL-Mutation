#ifndef _SegBed
#define _SegBed

class SegBed
{
public:
SegBed();
~SegBed();
void transfer(SegBed *);
void add(float **,int,int);
void clear();
void resize(int);
void add(float *,int,int);
void setsize(int);
int ifexist(int,int);
float **head;
int *start;
int *end;
Res **post;
float **xyzout;
SegBed *next;
int size;
int num;
int first;
int last;
};

#endif
