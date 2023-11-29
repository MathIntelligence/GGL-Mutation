#include"source.h"

StrFit::StrFit()
{
	next=0;
}

StrFit::~StrFit()
{
	if(next) delete next;
}

