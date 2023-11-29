#include"sourceclass.h" 
extern ResourceList resources;

Sequence::Sequence()
{
 int i,flg,j,k;
 FILE *fp;
 char *parameter_directory;
 char line[256],*buffer,*buffer0[Amino_Acid_Num];
 char  pro_1[]="ARNDBCQEZGHILKMFPSTWYV";
 char* pro_3[]={"ALA","ARG","ASN","ASP","ASX","CYS","GLN",
                "GLU","GLX","GLY","HIS","ILE","LEU","LYS",
                "MET","PHE","PRO","SER","THR","TRP","TYR",
                "VAL",0};
 char  amino_acd_ord[] = "ABCDEFGHIKLMNPQRSTVWXYZ";
 
//initialise 

 for(i=0;i<Amino_Acid_Num+1;i++)
 {
   pro_crg_aa[i]='\0';
   pos_crg_aa[i]='\0';
   neg_crg_aa[i]='\0';
   polar_aa[i]='\0';
   hydro_aa[i]='\0';
   amino_acid_order[i]='\0';
   user_defined[i]='\0';
   pro_res_1[i]='\0';
   small_aa[i]='\0';
   large_aa[i]='\0';
   medium_aa[i]='\0';
   for(flg=0;flg<4;flg++) pro_res_3[i][flg]='\0';
   for(flg=0;flg<Amino_Acid_Num+1;flg++)similar_aa[i][flg]='\0';
 }

 parameter_directory=resources["parameter_directory"];
 sprintf(line,"%s/seq.inf",parameter_directory);
 if((fp=fopen(line,"r"))==NULL) {cerr<<"no open:"<<line<<endl;exit(0);}
 
 while(fgets(line,256,fp)!=NULL)
 {
   i=-1;
   while(line[++i]==' ') continue;
   strcpy(line,line+i);
   if(line[0]=='\n') continue;
   if(line[0]=='!')continue;
   if(line[0]=='_')
   { 
     flg=0;j=0;
     if(strstr(line,"positive_charged_residue")) buffer= pos_crg_aa;
     else if(strstr(line,"negative_charged_residue")) buffer=neg_crg_aa;
     else if(strstr(line,"charged_residue")) buffer=pro_crg_aa;
     else if(strstr(line,"polar_residue")) buffer=polar_aa;
     else if(strstr(line,"hydrophobic_residue"))buffer=hydro_aa;
     else if(strstr(line,"large_residue"))buffer=large_aa;
     else if(strstr(line,"small_residue"))buffer=small_aa;
     else if(strstr(line,"medium_residue"))buffer=medium_aa;
     else if(strstr(line,"similar_residue"))
     {for(j=0;j<Amino_Acid_Num;j++)buffer0[j]=similar_aa[j]; flg=1;j=0;}
     else if(strstr(line,"user_defined_group")) buffer=user_defined; 
     continue;
   }
   else
   {
       if(!flg){
         strcpy(buffer,line);
         for(i=0;i<strlen(buffer);i++) 
         if(buffer[i]=='\n') buffer[i]='\0';
       }
       else
       { 
         k=-1;
         for(i=0;i<strlen(line);i++)
         { 
           if(line[i]!=','&&line[i]!='\n') buffer0[j][++k]=line[i];        
           else {k=-1;j++;}
         }
       } 
       
   } 

 }  
 strcpy(amino_acid_order,amino_acd_ord);
 strcpy(pro_res_1,pro_1);
 i=-1;
 while(pro_3[++i]) strcpy(pro_res_3[i],pro_3[i]); 

}



char Sequence::seq_to_1(char *res){
int i=0;
while(strncmp(res,pro_res_3[i++],3)) 
  if(i>=strlen(pro_res_1))  
  { 
   cerr<<res<<":no found\n";
   if(res[2]!=' ') return Nope;
   else return '0';
  }
  return pro_res_1[i-1];
}


char *Sequence::seq_to_3(char res){
 int i=0;
 while(res!=pro_res_1[i++] ) 
 if(i>=strlen(pro_res_1)) 
 {
  cerr<<res<<":no found\n";
  char s[3];
  s[0]=' ';s[1]=' ';s[2]=res; 
  return s; 
 } 
 return pro_res_3[i-1];
}
 Sequence::~Sequence()
{

}
