#include <stdio.h>
#include <stdlib.h>

main(int argc,char *argv[]) {


	char line[1000];
	
	
	sprintf(line,"%s -static -o ./obj/%s.o -c -D_XWINDOWS -c -D_XWINDOWS  ./%s.cpp",argv[1],argv[2],argv[2]);

	printf("%s\n",line);
	system(line);
	sprintf(line,"%s -static -o ./%s -I. ./obj/%s.o ../lib/libm.a -lc -lm",argv[1],argv[2],argv[2]);
	printf("%s\n",line);
	system(line);

}
