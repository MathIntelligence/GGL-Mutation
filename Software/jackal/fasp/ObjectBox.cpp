#include "source.h"
void ObjectBox::remove(int i) {

if(element==0&&ne==0) return;
if(ne>0&&element==0) {
ne=0;
return;
}

if(ne==0&&element) {
ne=0;
delete [] element;
return;
}


int t=0;
for(int j=0;j<ne;j++) {
if(element[j]==i) continue;
element[t++]=element[j];
}
ne=t;
if(ne==0) {
if(element) delete [] element;element=0;
}
}

void ObjectBox::add(int i) {

int t=0;
for(int j=0;j<ne;j++) {
if(element[j]==i) {t=1;break;}
}
if(t==1) return;
element[ne++]=i;
}

ObjectBox::ObjectBox() {
ne=0;
max=30;
element=0;
}

ObjectBox::~ObjectBox() {
if(element) delete element;
}

