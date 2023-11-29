#include "source.h"

FastaReader::FastaReader() {

	fileName=0;
	name=0;
	seqn=0;
	num=0; 
	next=0; 
}

FastaReader::~FastaReader() {

	Strhandler cc; 
	fileName=cc.strdel(fileName);
	name=cc.strdel(name);
	seqn=cc.strdel(seqn);
	num=0;
	if(next) delete next;
}

void FastaReader::readFastaAlign(char *f) {
	
		
	Strhandler cc; 

	fileName=f;

	char **lines = cc.opnfilebylinesimple(f);

	if(lines==0) return;

	int n= cc.gettotnum(lines);
	
	name = cc.opnarray(n);
	seqn = cc.opnarray(n);

	//clear lines
	for(int i=0;i<n;i++) { 
		lines[i]=cc.clearfirstchar(lines[i]," ");	
	}

	num=0;
	for(i=0;i<n;i++) {
		if(lines[i]==0) continue;
		if(lines[i][0]=='>') {
			name[num]=cc.opnstr(lines[i]+1);
			name[num]=cc.clearendchar(name[num]," ");
		}
		else {
			seqn[num]=cc.straddup(seqn[num],lines[i]);
		}
		lines[i]=cc.strdel(lines[i]);
	}
	
	 
	lines=cc.strdel(lines);
}

