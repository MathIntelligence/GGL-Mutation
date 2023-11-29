/*  A GLOBAL ALIGNMENT PROGRAM (GAP):

Copyright (c) 1992 Xiaoqiu Huang
All Rights Reserved.  E-mail: huang@cs.mtu.edu

Permission to use, copy, modify, and distribute this software and its
documentation for educational, research and non-profit purposes, without
fee, and without a written agreement is hereby granted, provided that the
above copyright notice, this paragraph and the following three paragraphs
appear in all copies.

Permission to incorporate this software into commercial products may be
obtained from the Intellectual Property Office, 1400 Townsend Drive, 301
Administration Building, Houghton, MI  49931, phone (906) 487-3429.

IN NO EVENT SHALL THE AUTHOR OR MICHIGAN TECHNOLOGICAL UNIVERSITY BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
DOCUMENTATION, EVEN IF MICHIGAN TECHNOLOGICAL UNIVERSITY HAS BEEN ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

MICHIGAN TECHNOLOGICAL UNIVERSITY SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN
"AS IS" BASIS, AND MICHIGAN TECHNOLOGICAL UNIVERSITY HAS NO OBLIGATIONS TO
PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

	      Xiaoqiu Huang
	      Department of Computer Science
	      Michigan Technological University
	      Houghton, MI 49931
	      E-mail: huang@cs.mtu.edu
	      WWW: http://www.cs.mtu.edu/faculty/huang.html

     Proper attribution of the author as the source of the software would
     be appreciated:
        Huang, X. (1994)
        On Global Sequence Alignment,
        Computer Applications in the Biosciences, 10, 227-235.

    The GAP program computes a global alignment of two sequences
    without penalizing terminal gaps. It delivers the alignment in
    linear space, so long sequences can be aligned. 

    Users supply scoring parameters. In the simplest form, users just
    provide 3 integers: ms, q and r, where ms is the score of a mismatch
    and the score of an i-symbol indel is -(q + r * i). Each match
    automatically receives score 10. This simple scoring scheme may be
    used for DNA sequences. NOTE: all scores are integers.

    In general, users can define an alphabet of characters appearing
    in the sequences and a matrix that gives the substitution score
    for each pair of symbols in the alphabet. The 127 ASCII characters
    are eligible. The alphabet and matrix are given in a file, where
    the first line lists the characters in the alphabet and the lower
    triangle of the matrix comes next. An example file looks as follows:

    ARNDC	       
     13
    -15  19
    -10 -22  11
    -20 -10 -20  18
    -10 -20 -10 -20  12

    Here the -22 at position (3,2) is the score of replacing N by R.
    This general scoring scheme is useful for protein sequences where the
    set of protein characters and Dayhoff matrix are specified in the file.

    The GAP program is written in C and runs under Unix systems on
    Sun workstations and under DOS systems on PCs.
    We think that the program is portable to many machines.

    Sequences to be analyzed are stored in separate files.
    Sequences must be in FASTA format. The first line begins with the symbol '>'
    followed by the name of the sequence. The sequence is on the remaining lines.
    Protein sequences must be in upper case.
    DNA sequences could be in upper or lower case.
    A sample sequence file is shown below.

>DNA sequence
GAATTCTAATCTCCCTCTCAACCCTACAGTCACCCATTTGGTATATTAAA
GATGTGTTGTCTACTGTCTAGTATCCCTCAAGTAGTGTCAGGAATTAGTC
ATTTAAATAGTCTGCAAGCCAGGAGTGGTGGCTCATGTCTGTAATTCCAG
CACTGGAGAGGTAGAAGTG

    To find the best alignment of two sequences in files A and B,
    use a command of form

	   gap  A  B  gs  ms  q  r > result

    where gap is the name of the object code, gs is the minimum length
    of any gap in the short sequence receiving a constant gap penalty,
    ms is a negative integer specifying mismatch weight, q and r are
    non-negative integers specifying gap-open and gap-extend penalties,
    respectively. Output alignment is saved in the file "result".

    For using a scoring matrix defined in file S, use a command of form

	   gap  A  B  gs  S  q  r > result

    Note that ms is replaced by the file S.

    Acknowledgments
    The functions diff2() and display() evolved from those written by Gene Myers.
    We made the following modifications: similarity weights (integer), instead of
    distance weights (float), are used, terminal gaps are not penalized, and
    any gap of length at least gs in the short sequence is given a constant
    penalty.
*/

#include   <stdio.h>
#include <string.h>
#include <stdlib.h>

#define  MATCHSC    10          /* match score */

static int match, mismh;		/* max and min substitution weights */
static char *dhead1, *dhead2;           /* names of sequences */
static int isdna;                       /* 1, DNA; 0, protein */
static int v[128][128];	/* substitution scores */
static int  q, r;       /* gap penalties */
static int  qr;         /* qr = q + r */
static int  gaplen;     /* minimum length for constant-cost insertion */
static int  pay;	/* constant-cost for long insertion */

static int  change;	/* constant-cost for long insertion */

char **cmain(char *a,char *b){

	char *argv[8];
	char **cmain0(int,char**);
	char **out;
	int i;
	argv[0]=strdup("nest");
	argv[1]=strdup(a);
	argv[2]=strdup(b);
	argv[3]=strdup("5");
	argv[4]=strdup("-5"); 
	argv[5]=strdup("5");
	argv[6]=strdup("5"); 
	argv[7]=0;
	
	out= cmain0(7,argv);
	/*
	for(i=0;i<7;i++) {
		free(argv[i]);
	}*/
	return out;
}

char **cmain0(argc, argv) int argc; char *argv[];
{ 
 	char **out;
	int  M, N;				/* Sequence lengths and k     */
  	char  *A,  *B;				/* Storing two sequences      */
  	int  symbol;				/* The next character	      */
  	int  ms;				/* User-supplied weights      */
  	FILE *Bp, *Ap, *Sp, *ckopen();
  	char *ckalloc();			/* space-allocating function  */
  	register int i, j;
  	char  alph[129], *s;			/* alphabet */
  	int  size,nt;				/* size of alphabet */
  	char  *name;				/* temporary pointer */
	char ** GLOBAL(char *,char *,int,int);  
	if ( argc != 7 )
	{ 	
		fprintf(stderr,"Usage: %s Seq1 Seq2 gap_size mismatch gap_open gap_extend\n\n", argv[0]);
  		fprintf(stderr,"Seq1        file of one sequence in FASTA format\n");
  		fprintf(stderr,"Seq2        file of one sequence in FASTA format\n");
  		fprintf(stderr,"gap_size    gap length for a constant gap penalty, a positive integer\n");
  		fprintf(stderr,"mismatch    a negative integer for DNA or PAM250/BLOSUM62 for protein\n");
  		fprintf(stderr,"gap_open    gap open penalty, a non-negative integer \n");
  		fprintf(stderr,"gap_extend  gap extension penalty, a positive integer \n");
  		exit(1);
	}
 
	A=(char *)strdup(argv[1]);
	M=strlen(A);
	
	B=(char *)strdup(argv[2]);
	N=strlen(B);

	if ( M == 0 && N == 0 ) fatal("Both sequences are empty");

	(void) sscanf(argv[3],"%d", &gaplen);

	if ( gaplen < 1 )  fatal("The minimum length for constant-cost insertion is a positive integer");

	(void) sscanf(argv[argc-2],"%d", &q);

	if ( q < 0 )  fatal("The gap-open penalty is a nonnegative integer");

	(void) sscanf(argv[argc-1],"%d", &r);

	if ( r < 0 )  fatal("The gap-extend penalty is a nonnegative integer");

	isdna = 0;
	pay = q + r * gaplen;
	qr = q + r;

	/* check if the argument represents a negative integer */
	
	s = argv[argc-3];

	if ( *s == '-' ) s++;
	for ( ; *s >= '0' && *s <= '9' ; s++ );
	if ( *s == '\0' )
	  { (void) sscanf(argv[argc-3],"%d", &ms);
	    if ( ms >= 0 )
	       fatal("The mismatch weight is a negative integer");
	    match = MATCHSC;
	    mismh = ms;
	    /* set match and mismatch weights */
	    for ( i = 0; i < 128 ; i++ )
	      for ( j = 0; j < 128 ; j++ )
	         if (i == j )
		    v[i][j] = match;
	         else
	            v[i][j] = mismh;
	    v['N']['N'] = mismh;
            v['n']['n'] = mismh;
            v['A']['a'] = v['a']['A'] = match;
            v['C']['c'] = v['c']['C'] = match;
            v['G']['g'] = v['g']['G'] = match;
            v['T']['t'] = v['t']['T'] = match;
            isdna = 1;
	  }
	else
	  { /* read a file containing alphabet and substitution weights */
	    Sp = ckopen(argv[argc-3], "r");
	    (void) fscanf(Sp, "%s", alph);
	    size = strlen(alph);
	    match = mismh = 0;
	    /* Initialize v[][] */
	    for ( i = 0; i < 128 ; i++ )
	      for ( j = 0; j < 128 ; j++ )
	            v[i][j] = 0;
	    for ( i = 0; i < size ; i++ )
	      for ( j = 0; j <= i ; j++ )
		{ (void) fscanf(Sp, "%d", &ms);
		  v[alph[i]][alph[j]] = v[alph[j]][alph[i]] = ms;
		  if ( ms > match ) match = ms;
		  if ( ms < mismh ) mismh = ms;
		}
	  }
	  change = 0;
	  if ( M <= N )
	    out= (char **)GLOBAL(A,B,M,N);
	  else
	   { change = 1;
	    out= (char **)GLOBAL(B,A,N,M);
           }
	   nt=strlen(out[0]);
	   j=strlen(out[1]);
	   i=strlen(out[2]);
	   if(nt<j) nt=j;
	   if(nt<i) nt=i;

	   j=strlen(out[0]);
	   for(i=0;i<j;i++) {
		if(out[0][i]==' ') out[0][i]='-';
	   }
	   for(i=j;i<nt;i++) {
		 out[0][i]='-';
           }

	   j=strlen(out[1]);
           for(i=0;i<j;i++) {
		if(out[1][i]==' ') out[1][i]='-';
	   }
	   for(i=j;i<nt;i++) {
		 out[1][i]='-';
           }
 	   j=strlen(out[2]);
           for(i=0;i<j;i++) {
		if(out[2][i]==' ') out[2][i]='-';
	   }
	   for(i=j;i<nt;i++) {
		 out[2][i]='-';
           }
	   
	   

	   return out;
}

static int *CC, *DD;			/* saving matrix scores */
static int *RR, *SS;		 	/* saving start-points */
static int  *S;				/* saving operations for diff */

/* The following definitions are for function diff() */

int  diff();
static int  zero = 0;				/* int type zero        */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

#define gap2(k)  ((k) <= 0 ? 0 : ((k) <= gaplen ? q+r*(k) : pay))
/* k-symbol insertion score */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

static int no_mat; 				/* number of matches */ 
static int no_mis; 				/* number of mismatches */ 
static int al_len; 				/* length of alignment */
						/* Append "Delete k" op */
#define DEL(k)				\
{ al_len += k;				\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ al_len += k;				\
  if (last > 0)				\
    last = sapp[-1] += (k);		\
  else					\
    last = *sapp++ = (k);		\
}
						/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
  al_len += 1;				\
}

char **GLOBAL(A,B,M,N) char A[],B[]; int M,N;
{ int  score;   		/* the max score in LIST */
  int  i, j;			/* row and column indices */
  char *ckalloc();		/* space-allocating function */
  char **display(char *,char *,int,int,int *,int,int);	
  char **ao;
	    /* allocate space for all vectors */
	    j = (N + 1) * sizeof(int);
	    CC = ( int * ) ckalloc(j);
	    DD = ( int * ) ckalloc(j);
	    RR = ( int * ) ckalloc(j);
	    SS = ( int * ) ckalloc(j);
	    i = (M + 1) * sizeof(int);
	    S = ( int * ) ckalloc(i + j);
	    (void) printf("Max Match   Min Mismatch   Gap-Open Penalty   Gap-Extension Penalty\n");
	    (void) printf("   %d          %d              %d                  %d\n\n", match, mismh, q, r);
	    (void) printf("Upper Sequence: %s\n", dhead1);
	    if ( change )
	      (void) printf("      Length: %d\n", N);
	    else
	      (void) printf("      Length: %d\n", M);
	    (void) printf("Lower Sequence: %s\n", dhead2);
	    if ( change )
	      (void) printf("      Length: %d\n\n", M);
	    else
	      (void) printf("      Length: %d\n\n", N);
            sapp = S;
            last = 0;
            al_len = 0;
            no_mat = 0;
	    no_mis = 0;
	    score = diff2(A,B,M,N,q,q,0,0,0,0);
            /* Output the best alignment */
            (void) printf("      Best Alignment with Constant Cost for Long Insertions\n");
            (void) printf("      Similarity Score : %d\n",score);
            (void) printf("      Match Percentage : %d%%\n", (100*no_mat)/al_len);
            (void) printf("      Number of Matches : %d\n", no_mat);
            (void) printf("      Number of Mismatches : %d\n", no_mis);
            (void) printf("      Total Length of Gaps : %d\n", al_len-no_mat-no_mis);
            (void) printf("      Note that terminal gaps are not penalized\n");
            ao= (char **)display(A,B,M,N,S,1,1);
		
	    /*
	    free(CC);CC=0;
	    free(DD);DD=0;
	    free(RR);RR=0;
	    free(SS);SS=0;
	    free(S);S=0;
	    */
	    return ao;
}

/* diff2(A,B,M,N,tb,te,sc,sr,ec,er) returns the score of an optimum conversion
   between A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script. If sc = 0, then
   the beginning deletion is not penalized; if sr = 0, the beginning insertion is
   not penalized; if ec = 0, the ending deletion is not charged; if er = 0;
   then the ending insertion is not charged. Any insertion of length at least
   gaplen is given a constant cost */

int diff2(A,B,M,N,tb,te,sc,sr,ec,er)
char *A, *B; int M, N; int tb, te, sc, sr, ec, er;

{ int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;
  int  ss,cc;

{ register int   i, j;
  register int c, e, d, s;
           int t, *va;
	   int  g, temp;
  	   char  *ckalloc();

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      if ( !sc || !ec )
	return 0;
      else
        return - gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
          if ( !sr || !er )
    	    return 0;
          else
            return - gap2(N);
        }
      midc = - (sc * (tb + r) + er * gap2(N) );
      midj = -1;
      if ( midc < ( c =  - (ec * (te + r) + sr * gap2(N) ) ) )
	{ midc = c;
	  midj = 0;
	}
      va = v[A[1]];
      for (j = 1; j <= N; j++)
	{ c = va[B[j]] - ( sr * gap2(j-1) + er * gap2(N-j) );
          if (c > midc)
           { midc = c;
             midj = j;
           }
	}
      if (midj == -1)
        { DEL(1) INS(N) }
      else
      if (midj == 0)
        { INS(N) DEL(1) }
      else
        { if (midj > 1) INS(midj-1)
          REP
	  if ( A[1] == B[midj] || isdna && va[B[midj]] == MATCHSC )
	     no_mat += 1;
	  else
	     no_mis += 1;
          if (midj < N) INS(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = - q * sr;
  if ( N <= gaplen )
    for (j = 1; j <= N; j++)
      { CC[j] = t = (t-r) * sr;
        DD[j] = t-q;
      }
  else
   { for (j = 1; j <= gaplen; j++)
      { CC[j] = t = (t-r) * sr;
        DD[j] = t-q;
      }
     for (j = gaplen+1; j <= N; j++)
      { CC[j] = t = -pay * sr;
        DD[j] = t - q;
      }
   }
  if ( !ec ) DD[N] += q;
  t = -tb * sc;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      CC[0] = c = t = (t-r) * sc;
      e = t-q;
      g = t - pay;
      va = v[A[i]];
      for (j = 1; j <= N; j++)
        { if ((c = c - qr) > (e = e - r)) e = c;
	  if ( j == N && !ec )
            { if ((c = CC[j] ) > (d = DD[j] )) d = c;}
	  else
            if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	  c = s+va[B[j]];
          if (c < d) c = d;
          if (c < e) c = e;
	  if ( j - gaplen > 0 )
	    { if ( g < ( temp = CC[j-gaplen-1] - pay ) )
		g = temp;
	      if ( c < g ) c = g;
	    }
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  t = -q * er;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  if ( N <= gaplen )
    for (j = N-1; j >= 0; j--)
      { RR[j] = t = (t-r) * er;
        SS[j] = t-q;
      }
  else
   { temp = N - gaplen;
     for (j = N-1; j >= temp; j--)
      { RR[j] = t = (t-r) * er;
        SS[j] = t-q;
      }
     for (j = temp-1; j >= 0; j--)
      { RR[j] = t = -pay * er;
        SS[j] = t - q;
      }
   }
  if ( !sc ) SS[0] += q;
  t = -te * ec;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      RR[N] = c = t = (t-r) * ec;
      g = t - pay;
      e = t-q;
      va = v[A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c = c - qr) > (e = e - r)) e = c;
	  if ( !j && !sc )
            { if ((c = RR[j] ) > (d = SS[j] )) d = c;}
	  else
            if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	  c =  s+va[B[j+1]];
          if (c < d) c = d;
          if (c < e) c = e;
	  if ( j + gaplen < N )
	    { if ( g < ( temp = RR[j+gaplen+1] - pay ) )
		g = temp;
	      if ( c < g ) c = g;
	    }
          s = RR[j];
          RR[j] = c;
          SS[j] = d;
        }
    }
  SS[N] = RR[N];

  midc = CC[0]+RR[0];		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  for (j = 0; j <= N; j++)
    if ((c = CC[j] + RR[j]) >= midc)
      if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
        { midc = c;
          midj = j;
        }
  for (j = N; j >= 0; j--)
   { if ( j == N )
       d = q * ec;
     else
       if ( j == 0 )
         d = q * sc;
       else
	 d = q;
     if ((c = DD[j] + SS[j] + d) > midc)
       { midc = c;
         midj = j;
         type = 2;
       }
   }
}

/* Conquer: recursively around midpoint */

  cc = midj == N ? ec : 1;
  ss = midj == 0 ? sc : 1;
  if (type == 1)
    { (void) diff2(A,B,midi,midj,tb,q,sc,sr,cc,1);
      (void) diff2(A+midi,B+midj,M-midi,N-midj,q,te,ss,1,ec,er);
    }
  else
    { (void) diff2(A,B,midi-1,midj,tb,zero,sc,sr,cc,1);
      DEL(2);
      (void) diff2(A+midi+1,B+midj,M-midi-1,N-midj,zero,te,ss,1,ec,er);
    }
  return midc;
}
/* Alignment display routine */

static char ALINE[51], BLINE[51], CLINE[51];

char **display(A,B,M,N,S,AP,BP) char A[], B[]; int M, N; int S[], AP, BP;
{ register char *a, *b, *c;
  register int   i,  j, op;
           int   lines, ap, bp;
 
  char **out;
  int nmax; 
  int itp;
  i = j = op = lines = 0;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
  nmax=(M+N+100)*2;
  out=(char **) malloc(4*sizeof(char));
  for(itp=0;itp<3;itp++) {
  	out[itp]=(char *)malloc(nmax*sizeof(char));
	out[itp][0]='\0';
  }
  out[itp]=0;
  
  while (i < M || j < N)
    { if (op == 0 && *S == 0)
        { op = *S++;
          *a = A[++i];
          *b = B[++j];
	  *c++ = (*a == *b || isdna && v[*a][*b] == MATCHSC ) ? '|' : ' ';
	  a++;
	  b++;
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              *b++ = B[++j];
              op--;
            }
          else
            { *a++ = A[++i];
              *b++ = ' ';
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || i >= M && j >= N)
        { *a = *b = *c = '\0';
          (void) printf("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
            (void) printf("    .    :");
          if (b <= a+5)
            (void) printf("    .");
	  if ( change ) {
          	(void) printf("\n%5d %s\n      %s\n%5d %s\n",bp,BLINE,CLINE,ap,ALINE);
		out[0]=strcat(out[0],BLINE);
		out[1]=strcat(out[1],ALINE);
		out[2]=strcat(out[2],CLINE);
		
          }
	  else {
          	(void) printf("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
		out[0]=strcat(out[0],ALINE);
		out[1]=strcat(out[1],BLINE);
		out[2]=strcat(out[2],CLINE);
		
          }
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
	  	
        }
    }
    return out;
}

/* lib.c - library of C procedures. */

/* fatal - print message and die */
fatal(msg)
char *msg;
{
	(void) fprintf(stderr, "%s\n", msg);
	exit(1);
}

/* fatalf - format message, print it, and die */
fatalf(msg, val)
char *msg, *val;
{
	(void) fprintf(stderr, msg, val);
	(void) putc('\n', stderr);
	exit(1);
}
	
/* ckopen - open file; check for success */
FILE *ckopen(name, mode)
char *name, *mode;
{
	FILE *fopen(), *fp;

	if ((fp = fopen(name, mode)) == NULL)
		fatalf("Cannot open %s.", name);
	return(fp);
}

/* ckalloc - allocate space; check for success */
char *ckalloc(amount)
int amount;
{
	char  *p;

	if ((p = (char *)malloc( (unsigned) amount)) == NULL)
		fatal("Ran out of memory.");
	return(p);
}
