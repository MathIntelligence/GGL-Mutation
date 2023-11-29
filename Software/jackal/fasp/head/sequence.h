#ifndef _Sequence
#define _Sequence

#include"constant.h"

class Sequence{

public:
 
Sequence();
~Sequence();

//member function

char *seq_to_3(char);
char  seq_to_1(char *);

//member data

char  pos_crg_aa[Amino_Acid_Num]; //positive charged aa
char  neg_crg_aa[Amino_Acid_Num];  //negatiev charged aa
char  polar_aa[Amino_Acid_Num]; //polar aa
char  hydro_aa[Amino_Acid_Num];  // hydrophobic aa
char  pro_res_1[Amino_Acid_Num]; //general aa
char  pro_res_3[Amino_Acid_Num][4];  //triple code of protein aa
char  amino_acid_order[Amino_Acid_Num];
char  pro_crg_aa[Amino_Acid_Num];
char  similar_aa[Amino_Acid_Num][Amino_Acid_Num];
char  user_defined[Amino_Acid_Num];
char  medium_aa[Amino_Acid_Num];
char  large_aa[Amino_Acid_Num];
char  small_aa[Amino_Acid_Num];

};

#endif
