import java.lang.*;
import java.io.*;

public class SW{

   private static String reps=
      "4,20,0,7,8,17,5,10,14,20,12,15,13,6,20,3,9,11,1,2,20,16,19,20,18,";
   private static String blosums=
     "9,-1,4,-1,1,5,-3,-1,-1,7,0,1,0,-1,4,-3,0,-2,-2,0,6,"+
  "-3,1,0,-2,-2,0,6,-3,0,-1,-1,-2,-1,1,6,-4,0,-1,-1,-1,-2,0,2,5,"+
  "-3,0,-1,-1,-1,-2,0,0,2,5,-3,-1,-2,-2,-2,-2,1,-1,0,0,8,"+
  "-3,-1,-1,-2,-1,-2,0,-2,0,1,0,5,-3,0,-1,-1,-1,-2,0,-1,1,1,-1,2,5,"+
  "-1,-1,-1,-2,-1,-3,-2,-3,-2,0,-2,-1,-1,5,"+
  "-1,-2,-1,-3,-1,-4,-3,-3,-3,-3,-3,-3,-3,1,4,"+
  "-1,-2,-1,-3,-1,-4,-3,-4,-3,-2,-3,-2,-2,2,2,4,"+
  "-1,-2,0,-2,0,-3,-3,-3,-2,-2,-3,-3,-2,1,3,1,4,"+
  "-2,-2,-2,-4,-2,-3,-3,-3,-3,-3,-1,-3,-3,0,0,0,-1,6,"+
  "-2,-2,-2,-3,-2,-3,-2,-3,-2,-1,2,-2,-2,-1,-1,-1,-1,3,7,"+
  "-2,-3,-2,-4,-3,-2,-4,-4,-3,-2,-2,-3,-3,-1,-3,-2,-3,1,2,11,";

   private int n,m;
   private byte[] s1,s2;
   private byte[] rep;
   private int[] table;
   private char[] res;
   private short[][] V,E,F;
   private byte[][] path;
   
   private void init(){
     rep=new byte[25];  int i, end;  String str=reps;
     res=new char[20];
     for(i=0;i<25;i++){  
        end=str.indexOf(',');  
        rep[i]=(byte)(Integer.parseInt(str.substring(0,end)));
        if(rep[i]<20) res[rep[i]]=(char)('A'+i);
        str=str.substring(end+1); 
     }
     str=blosums;
     table=new int[210];
     for(i=0;i<210;i++){
        end=str.indexOf(',');  
        table[i]=Integer.parseInt(str.substring(0,end));
        str=str.substring(end+1); 
     }
   }
   
   private int seqOnly(byte[] s){
     int i=0, current=0;  
     for(;;) if(s[current++]=='\n') break;
     for(;;)
       if(s[current]=='\n' && s[++current]=='\n') return i;
       else s[i++]=rep[s[current++]-'A'];
   }  
       
   
   public SW(String[] files){
     init();
     FileInputStream f=null;
     try{ 
       f=new FileInputStream(files[0]);
     s1=new byte[f.available()];
     f.read(s1); n=seqOnly(s1);
     f=new FileInputStream(files[1]);
     s2=new byte[f.available()];
     f.read(s2); m=seqOnly(s2);
     }catch(FileNotFoundException e){
        System.out.println(e); System.exit(1); }
      catch(IOException e){
        System.out.println(e); System.exit(1); }
   }

   public void traceback(int i,int j){
     byte[] trace=new byte[n+m+2];
     int p=i, q=j, r=0;
     while(path[p][q]>0){
       trace[r++]=path[p][q];
       switch(path[p][q]){
         case 1: p--; q--; break;
         case 2: q--; break;
         case 3: p--; break;
         default: ;
       }
     }
   }       
   
   public void dp(){
     int mmax=0, maxi=0, maxj=0;
     V=new short[n+1][m+1];  
     E=new short[n+1][m+1];  
     F=new short[n+1][m+1];  
     path=new byte[n+1][m+1];  
     for(int i=0;i<=n;i++){ V[i][0]=E[i][0]=F[i][0]=0; path[i][0]=0; }
     for(int j=1;j<=m;j++){ V[0][j]=E[0][j]=F[0][j]=0; path[0][j]=0; }
     for(int i=1;i<=n;i++) for(int j=1;j<=m;j++){
       int max=0; int a=s1[i-1], b=s2[j-1],c,  p=0;
       int sub=(a<=b)?table[b*(b+1)/2+a]:table[a*(a+1)/2+b];
       if((c=V[i-1][j-1]+sub)>max){ max=c; p=1; }
       E[i][j]=(short)((E[i][j-1]>=V[i][j-1]-8)?E[i][j-1]-2:V[i][j-1]-10); 
       if(E[i][j]>max){ max=E[i][j]; p=2; }
       F[i][j]=(short)((F[i-1][j]>=V[i-1][j]-8)?F[i-1][j]-2:V[i-1][j]-10); 
       if(F[i][j]>max){ max=F[i][j]; p=3; }
       V[i][j]=(short)max; path[i][j]=(byte)p;
       if(max>mmax){ mmax=max; maxi=i; maxj=j; }
     }
     System.out.println(mmax);
     System.out.println(maxi);
     System.out.println(maxj);
   }  
   
     
   public static void main(String[] args){
     SW sw=new SW(args);
     sw.dp();
   }
   
}

