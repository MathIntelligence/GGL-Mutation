#include"source.h"

SegBed::SegBed()
{
head=0;
post=0;
start=0;
end=0;
xyzout=0;
next=0;
num=0;
size=0;
first=-1;
last=-1;
}

SegBed::~SegBed()
{
 clear();
 if(next) delete next;
}

void SegBed::clear() {
 Strhandler cc;
 xyzout=cc.floatdel(xyzout);
 start=cc.intdel(start);
 end=cc.intdel(end);
 head=cc.floatdel(head);
 if(post) delete [] post;
 post=0;
 num=0;
 next=0;
 size=0;
 first=-1;
 last=-1;
}

void SegBed::transfer(SegBed *s) {
	clear();
	xyzout=s->xyzout;s->xyzout=0;
	start=s->start;s->start=0;
	end=s->end;s->end=0;
	size=s->size;s->size=0;
	num=s->num;s->num=0;
	first=s->first;s->first=-1;
	last=s->last;s->last=-1;
	post=s->post;s->post=0;
	head=s->head;s->head=0;
}
int SegBed::ifexist(int n1,int n2) {

	int i;
	for(i=0;i<num;i++) {
		if(start[i]==n1&&end[i]==n2) return 1;
	}
	return 0;
}

void SegBed::setsize(int n){
 
 start=new int[n];
 end=new int[n];
 xyzout=new float*[n];
 head=new float*[n];
 for(int i=0;i<n;i++) {
   start[i]=-1;end[i]=-1;xyzout[i]=0;
   head[i]=0;
 }
 size=n;
}

void SegBed::add(float *a,int b,int c) {
	
 	if(num<size-2) {
 		start[num]=b;
 		end[num]=c;
 		xyzout[num]=a;
 		num++;
 	}
	else {
		resize(size+1000);
		add(a,b,c);
	}
	a=0;
}

void SegBed::add(float **a,int b,int c) {

	if(a==0) return;
	int kep=0;while(a[kep])kep++;
	int i;
	for(i=0;i<kep;i++) add(a[i],b,c);
	if(a) delete [] a;a=0;
}

void SegBed::resize(int n){
 
	int i;
	int  *start0=start;
   	int  *end0=end;
	float **xyzout0=xyzout;
	float **head0=head;
	start=new int[n];
 	end=new int[n];
 	xyzout=new float*[n];
	head=new float*[n];

	for(i=0;i<n;i++) {
   		head[i]=0; start[i]=-1;end[i]=-1;xyzout[i]=0;
 	}
	
 	for(i=0;i<size;i++) {
   		start[i]=start0[i];
		end[i]=end0[i];
		xyzout[i]=xyzout0[i];
		head[i]=head0[i];
 	}

	size=n;
	Strhandler cc;
 	cc.floatdel(xyzout0);
 	cc.intdel(start0);
 	cc.intdel(end0);
	cc.floatdel(head0);
}
