/* A PROGRAM FOR LOCAL SIMILARITIES WITH AFFINE WEIGHTS:

    copyright (c) 1990-1997 Xiaoqiu Huang and Webb Miller
    The distribution of the program is granted provided no charge is made
    and the copyright notice is included.
    E-mail: huang@cs.mtu.edu

IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    Please cite the paper
    "A Time-Efficient, Linear-Space Local Similarity Algorithm"
	(Advances in Applied Mathematics, 12: 337-357, 1991)

	      Xiaoqiu Huang
	      Department of Computer Science
	      Michigan Technological University
	      Houghton, MI 49931

	      Webb Miller
	      Department of Computer Science
	      The Pennsylvania State University
	      University Park, PA 16802

     An early version of this program is presented in the paper
     "A Space-Efficient Algorithm for Local Similarities"
	(Computer Applications in the Biosciences, 6(4): 373-381, 1990)

	      Xiaoqiu Huang
	      Department of Computer Science
	      Michigan Technological University
	      Houghton, MI 49931

	      Ross C. Hardison and Webb Miller
	      Department of Molecular and Cell Biology
	      Department of Computer Science
	      The Pennsylvania State University
	      University Park, PA 16802

    The SIM program finds k best non-intersecting alignments between
    two sequences or within one sequence. Using dynamic programming
    techniques, SIM is guaranteed to find optimal alignments. The
    alignments are reported in order of similarity score, with the
    highest scoring alignment first. The k best alignments share no
    aligned pairs. SIM requires space proportional to the sum of the
    input sequence lengths and the output alignment lengths. Thus
    SIM can handle sequences of tens of thousands, even a hundred of
    thousands, of base pairs on a workstation. For example,
    on 73,360-bp and 44,594-bp sequences, SIM took 15 hours to find
    100 best local alignments on a Sun4 workstation.

    Users supply scoring parameters. In the simplest form, users just
    provide 3 integers: ms, q and r, where ms is the score of a mismatch
    and the score of an i-symbol indel is -(q + r * i). Each match
    automatically receives score 10. This simple scoring scheme may be
    used for DNA sequences. NOTE: all scores are integers.

    In general, users can define an alphabet of characters to appear
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

    The SIM program is written in C and runs under Unix systems on
    Sun workstations and under DOS systems on PCs.
    We think that the program is portable to many machines.

    Sequences to be analyzed are stored in separate files.
    Sequences must be in FASTA format. The first line begins with the symbol '>'
    followed by the name of the sequence.  The sequence is on the remaining lines.
    Protein sequences must be in upper case.
    DNA sequences could be in upper or lower case.
    A sample sequence file is shown below.

>DNA sequence
GAATTCTAATCTCCCTCTCAACCCTACAGTCACCCATTTGGTATATTAAA
GATGTGTTGTCTACTGTCTAGTATCCCTCAAGTAGTGTCAGGAATTAGTC
ATTTAAATAGTCTGCAAGCCAGGAGTGGTGGCTCATGTCTGTAATTCCAG
CACTGGAGAGGTAGAAGTG

    The SIM program generates alignments to the standard output.
    Redirection is required to put alignments in a specified file. 
    A sample output file is shown below, where alignment 1 aligns
    the 3240-3313 region of A with the 26300-26374 region of B.

Match   Mismatch   Gap-Open Penalty   Gap-Extension Penalty
 10      -15            30                 5

Upper Sequence: A
        Length: 36741
Lower Sequence: B
        Length: 73360

*********************************************************
      Number 1 Local Alignment
      Similarity Score : 155
      Match Percentage : 70%
      Number of Matches : 54
      Number of Mismatches : 18
      Total Length of Gaps : 5
      Begins at (3240, 26300) and Ends at (3313, 26374)

    0     .    :    .    :    .    :    .    :    .    :
 3240 GAGGGGCATTTGAGGGTGTTTCCAATGTTCCTGTTATTCGGAATAGCGCT
      || || ||||||-|| || ||||||-||   || ||||  |||||| || 
26300 GATGGCCATTTG GGTTGGTTCCAA GTCTTTGGTATTGTGAATAGTGCC

   50     .    :    .    :    .
 3290 GGTGTGAACATTC   TGCACAGGTCT
      |   | ||||| |---|||||| ||||
26348 GCAATAAACATACGTGTGCACATGTCT

*********************************************************
      Number 2 Local Alignment
      Similarity Score : 145
      Match Percentage : 86%
      Number of Matches : 19
      Number of Mismatches : 3
      Total Length of Gaps : 0
      Begins at (4566, 47235) and Ends at (4587, 47256)

    0     .    :    .    :
 4566 GGGGAGGGAGGGGAGCAGAGCA
      |||||||  |||||| ||||||
47235 GGGGAGGTTGGGGAGAAGAGCA


    To find k best non-intersecting alignments of segments from two
    sequences in files A and B, respectively, use a command of form

	   sim  k  A  B  ms  q  r > result

    where sim is the name of the object code, k is a positive integer,
    ms is a negative integer specifying mismatch weight, q and r are
    non-negative integers specifying gap-open and gap-extend penalties,
    respectively. Output alignments are saved in the file "result".

    For using a scoring matrix defined in file S, use a command of form

	   sim  k  A  B  S  q  r > result

    Note that ms is replaced by the file S.

    To find repeats in one sequence in file A, use a command of form

	   sim  k  A  ms  q  r > result
    or
	   sim  k  A  S  q  r > result

    Acknowledgments
    The functions diff() and display() are from Gene Myers. We made
    the following modifications: similarity weights (integer), instead of
    distance weights (float), are used, the aligned pairs already output
    are not used in the subsequent computation, and the positions of
    sequence characters in output alignments are shown by increment of 50.
    T. Mark Reboul found truncation errors in the previous version
    of SIM and also made a few helpful suggestions.
*/

#include   <stdio.h>

#define  MATCHSC    10		/* match score */

static int match, mismh;		/* max and min substitution weights */
static int v[128][128];			/* substitution scores */
static char *dhead1, *dhead2;           /* names of sequences */
static int isdna;          		/* 1, DNA; 0, protein */

main(argc, argv) int argc; char *argv[];
{ int  M, N, K;				/* Sequence lengths and k     */
  char  *A,  *B;			/* Storing two sequences      */
  int  symbol;				/* The next character	      */
  int  nseq;				/* Number of sequences        */
  int  ms, q, r;			/* User-supplied weights      */
  FILE *Bp, *Ap, *Sp, *ckopen();
  char *ckalloc();			/* space-allocating function  */
  register int i, j;
  char  alph[129], *s;			/* alphabet */
  int  size;				/* size of alphabet */

	if ( argc != 6 && argc != 7 )
{ fprintf(stderr,"Usage: %s k Seq1 [Seq2] mismatch gap_open gap_extend\n\n", argv[0]);
  fprintf(stderr,"k           number of local alignments, a positive integer \n");
  fprintf(stderr,"Seq1        file of one sequence in FASTA format\n");
  fprintf(stderr,"Seq2        file of one sequence in FASTA format\n");
  fprintf(stderr,"mismatch    a negative integer for DNA or PAM250/BLOSUM62 for protein\n");
  fprintf(stderr,"gap_open    gap open penalty, a non-negative integer \n");
  fprintf(stderr,"gap_extend  gap extension penalty, a positive integer \n");
  exit(1);
}
        /* read k: the number of local alignments to find */
	(void) sscanf(argv[1],"%d", &K);

	/* determine the sequence lengths */
	Ap = ckopen(argv[2], "r");
	if ( (symbol = getc(Ap) ) != '>' )
	  fatal("Sequence one must be in FASTA format, starting with '>'");
	for (size = 2; ( symbol = getc(Ap)) != EOF ; size++ )
	  if ( symbol == '\n' )
	    break;
	for (M = 0; ( symbol = getc(Ap)) != EOF ; )
	   if ( symbol != '\n' )
	      ++M;
	(void) fclose(Ap);

	/* allocate space for A and dhead1 */
	A = ( char * ) ckalloc( (M + 1) * sizeof(char));
	dhead1 = ( char * ) ckalloc( size * sizeof(char));

	/* read the first sequence into A */
	Ap = ckopen(argv[2], "r");
	symbol = getc(Ap);
	for (i = 0; ( symbol = getc(Ap)) != EOF && symbol != '\n' ; )
	   dhead1[i++] = symbol;
	dhead1[i] = '\0';
	for (M = 0; ( symbol = getc(Ap)) != EOF ; )
	   if ( symbol != '\n' )
		A[++M] = symbol;
	(void) fclose(Ap);

	nseq = 1;

	/* if there is another sequence, read it into B */
	if ( argc == 7 )
	  { nseq = 2;
	    Bp = ckopen(argv[3], "r");
	    if ( (symbol = getc(Bp) ) != '>' )
	     fatal("Sequence two must be in FASTA format, starting with '>'");
	    for (size = 2; ( symbol = getc(Bp)) != EOF ; size++ )
	      if ( symbol == '\n' )
	        break;
	    for (N = 0; ( symbol = getc(Bp)) != EOF ; )
	       if ( symbol != '\n' )
	          ++N;
	    (void) fclose(Bp);
	    B = ( char * ) ckalloc( (N + 1) * sizeof(char));
	    dhead2 = ( char * ) ckalloc( size * sizeof(char));
	    Bp = ckopen(argv[3], "r");
	    symbol = getc(Bp);
	    for (i = 0; ( symbol = getc(Bp)) != EOF && symbol != '\n' ; )
	      dhead2[i++] = symbol;
	    dhead2[i] = '\0';
	    for (N = 0; ( symbol = getc(Bp)) != EOF ; )
	       if ( symbol != '\n' )
	            B[++N] = symbol;
	    (void) fclose(Bp);
	  }

	(void) sscanf(argv[argc-2],"%d", &q);
	if ( q < 0 )
	   fatal("The gap-open penalty is a nonnegative integer");

	(void) sscanf(argv[argc-1],"%d", &r);
	if ( r <= 0 )
	   fatal("The gap-extend penalty is a positive integer");

	isdna = 0;
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
	    for ( i = 0; i < size ; i++ )
	      for ( j = 0; j <= i ; j++ )
		{ (void) fscanf(Sp, "%d", &ms);
		  v[alph[i]][alph[j]] = v[alph[j]][alph[i]] = ms;
		  if ( ms > match ) match = ms;
		  if ( ms < mismh ) mismh = ms;
		}
	  }

  	if ( nseq == 2 )
	  SIM(A,B,M,N,K,q,r,nseq);
	else
	  SIM(A,A,M,M,K,q,r,nseq);
}

static int q, r;			/* gap penalties */
static int qr;				/* qr = q + r */

typedef struct ONE { int COL ;  struct ONE  *NEXT ;} pair, *pairptr;
pairptr *row, z; 			/* for saving used aligned pairs */
static short tt;

typedef struct NODE
	{ int  SCORE;
	  int  STARI;
	  int  STARJ;
	  int  ENDI;
	  int  ENDJ;
	  int  TOP;
	  int  BOT;
	  int  LEFT;
	  int  RIGHT; }  vertex, *vertexptr;
		
vertexptr  *LIST;			/* an array for saving k best scores */
vertexptr  low = 0;			/* lowest score node in LIST */
vertexptr  most = 0;			/* latestly accessed node in LIST */
static int numnode;			/* the number of nodes in LIST */

static int *CC, *DD;			/* saving matrix scores */
static int *RR, *SS, *EE, *FF;	 	/* saving start-points */
static int *HH, *WW;		 	/* saving matrix scores */
static int *II, *JJ, *XX, *YY;	 	/* saving start-points */
static int  m1, mm, n1, nn;		/* boundaries of recomputed area */
static int  rl, cl;			/* left and top boundaries */
static int  min;			/* minimum score in LIST */
static short flag;			/* indicate if recomputation necessary*/

/* DIAG() assigns value to x if (ii,jj) is never used before */
#define DIAG(ii, jj, x, value)				\
{ for ( tt = 1, z = row[(ii)]; z != 0; z = z->NEXT )	\
    if ( z->COL == (jj) )				\
      { tt = 0; break; }				\
  if ( tt )						\
    x = ( value );					\
}

/* replace (ss1, xx1, yy1) by (ss2, xx2, yy2) if the latter is large */
#define ORDER(ss1, xx1, yy1, ss2, xx2, yy2)		\
{ if ( ss1 < ss2 )					\
    { ss1 = ss2; xx1 = xx2; yy1 = yy2; }		\
  else							\
    if ( ss1 == ss2 )					\
      { if ( xx1 < xx2 )				\
	  { xx1 = xx2; yy1 = yy2; }			\
	else						\
	  if ( xx1 == xx2 && yy1 < yy2 )		\
	    yy1 = yy2;					\
      }							\
}

/* The following definitions are for function diff() */

int diff();
static int  zero = 0;				/* int type zero        */

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

static int I, J;				/* current positions of A ,B */
static int no_mat; 				/* number of matches */ 
static int no_mis; 				/* number of mismatches */ 
static int al_len; 				/* length of alignment */
						/* Append "Delete k" op */
#define DEL(k)				\
{ I += k;				\
  al_len += k;				\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ J += k;				\
  al_len += k;				\
  if (last < 0)				\
    { sapp[-1] = (k); *sapp++ = last; }	\
  else					\
    last = *sapp++ = (k);		\
}

						/* Append "Replace" op */
#define REP 				\
{ last = *sapp++ = 0; 			\
  al_len += 1;				\
}

/* SIM(A,B,M,N,K,Q,R,nseq) reports K best non-intersecting alignments of
   the segments of A and B in order of similarity scores, where
   v[a][b] is the score of aligning a and b, and -(Q+R*i) is the score
   of an i-symbol indel.  						*/

SIM(A,B,M,N,K,Q,R,nseq) char A[],B[]; int M,N,K,nseq; int Q,R;
{ int endi, endj, stari, starj;	/* endpoint and startpoint */ 
  int  score;   			/* the max score in LIST */
  int count;				/* maximum size of list */	
  register  int  i, j;			/* row and column indices */
  char *ckalloc();			/* space-allocating function */
  int  *S;				/* saving operations for diff */
  vertexptr cur; 			/* temporary pointer */
  vertexptr findmax();	 		/* return the largest score node */
	
	/* allocate space for all vectors */
	j = (N + 1) * sizeof(int);
	CC = ( int * ) ckalloc(j);
	DD = ( int * ) ckalloc(j);
	RR = ( int * ) ckalloc(j);
	SS = ( int * ) ckalloc(j);
	EE = ( int * ) ckalloc(j);
	FF = ( int * ) ckalloc(j);
	i = (M + 1) * sizeof(int);
	HH = ( int * ) ckalloc(i);
	WW = ( int * ) ckalloc(i);
	II = ( int * ) ckalloc(i);
	JJ = ( int * ) ckalloc(i);
	XX = ( int * ) ckalloc(i);
	YY = ( int * ) ckalloc(i);
	S = ( int * ) ckalloc(i + j);
	row = ( pairptr * ) ckalloc( (M + 1) * sizeof(pairptr));

	/* set up list for each row */
	for ( i = 1; i <= M; i++ )
	  if ( nseq == 2 )
	     row[i] = 0;
	  else
	    { row[i] = z = ( pairptr ) ckalloc( (int) sizeof(pair));
              z->COL = i;			
              z->NEXT = 0;
	    }

	q = Q;
	r = R;
	qr = q + r;

	LIST = ( vertexptr * ) ckalloc( K * sizeof(vertexptr));
	for ( i = 0; i < K ; i++ )
	   LIST[i] = ( vertexptr ) ckalloc( (int) sizeof(vertex));

 (void) printf("Match   Mismatch   Gap-Open Penalty   Gap-Extension Penalty\n");
 (void) printf(" %d      %d            %d                 %d\n\n",
			match, mismh, q, r);
	if ( nseq == 2 )
	  { (void) printf("Upper Sequence: %s\n", dhead1);
	    (void) printf("        Length: %d\n", M);
	    (void) printf("Lower Sequence: %s\n", dhead2);
	    (void) printf("        Length: %d\n", N);
	  }
	else
	  { (void) printf("Single Sequence: %s\n", dhead1);
	    (void) printf("         Length: %d\n", M);
	  }

	numnode = min = 0;
	big_pass(A,B,M,N,K,nseq);

        /* Report the K best alignments one by one. After each alignment is
           output, recompute part of the matrix. First determine the size
	   of the area to be recomputed, then do the recomputation         */

	for ( count = K - 1; count >= 0 ; count-- )
	  { if ( numnode == 0 )
	      fatalf("There are no %d alignments between two sequences", K);
            cur = findmax();	/* Return a pointer to a node with max score*/
            score = cur->SCORE;
      	    stari = ++cur->STARI;
            starj = ++cur->STARJ;
            endi = cur->ENDI;
            endj = cur->ENDJ;
            m1 = cur->TOP;
            mm = cur->BOT;
            n1 = cur->LEFT;
            nn = cur->RIGHT;
            rl = endi - stari + 1;
            cl = endj - starj + 1;
            I = stari - 1;
            J = starj - 1;
            sapp = S;
            last = 0;
            al_len = 0;
            no_mat = 0;
	    no_mis = 0;
            (void) diff(&A[stari]-1, &B[starj]-1,rl,cl,q,q);
            /* Output the best alignment */
 (void) printf("\n*********************************************************\n");
            (void) printf("      Number %d Local Alignment\n", K - count);
            (void) printf("      Similarity Score : %d\n",score);
          (void) printf("      Match Percentage : %d%%\n", (100*no_mat)/al_len);
            (void) printf("      Number of Matches : %d\n", no_mat);
            (void) printf("      Number of Mismatches : %d\n", no_mis);
      (void) printf("      Total Length of Gaps : %d\n", al_len-no_mat-no_mis);
            (void) printf("      Begins at (%d, %d) and Ends at (%d, %d)\n",
	                stari,starj, endi,endj);
            display(&A[stari]-1,&B[starj]-1,rl,cl,S,stari,starj);
	    (void) fflush(stdout);

            if ( count )
	      { flag = 0;
                locate(A,B,nseq);
                if ( flag )
		   small_pass(A,B,count,nseq);
              }
	  }
}

/* A big pass to compute K best classes */

big_pass(A,B,M,N,K,nseq) char A[],B[]; int M,N,K,nseq;
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  *va;				/* pointer to v(A[i], B[j]) */
  int   addnode();			/* function for inserting a node */

	
	/* Compute the matrix and save the top K best scores in LIST
	   CC : the scores of the current row
	   RR and EE : the starting point that leads to score CC
	   DD : the scores of the current row, ending with deletion
	   SS and FF : the starting point that leads to score DD        */
 	/* Initialize the 0 th row */
	for ( j = 1; j <= N ; j++ )
	  {  CC[j] = 0;
	     RR[j] = 0;
	     EE[j] = j;
	     DD[j] = - (q);
	     SS[j] = 0;
	     FF[j] = j;
	  }
	for ( i = 1; i <= M; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     ci = fi = i;
	     va = v[A[i]];
	     if ( nseq == 2 )
	       { p = 0;
	         pi = i - 1;
	         cj = fj = pj = 0;
	       }
	     else
	       { p = CC[i];
		 pi = RR[i];
		 pj = EE[i];
	         cj = fj = i;
	       }
	     for ( j = (nseq == 2 ? 1 : (i+1)) ; j <= N ; j++ )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  DIAG(i, j, c, p+va[B[j]])		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > min )	/* add the score into list */
		    min = addnode(c, ci, cj, i, j, K, min);
	        }
	  }
}

/* Determine the left and top boundaries of the recomputed area */

locate(A,B,nseq) char A[],B[]; int nseq;
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  short  cflag, rflag;			/* for recomputation */
  int  *va;				/* pointer to v(A[i], B[j]) */
  int   addnode();			/* function for inserting a node */
  int  limit;				/* the bound on j */

	/* Reverse pass
	   rows
	   CC : the scores on the current row
	   RR and EE : the endpoints that lead to CC
	   DD : the deletion scores 
	   SS and FF : the endpoints that lead to DD

	   columns
	   HH : the scores on the current columns
	   II and JJ : the endpoints that lead to HH
	   WW : the deletion scores
	   XX and YY : the endpoints that lead to WW
	*/
	for ( j = nn; j >= n1 ; j-- )
          {  CC[j] = 0;
	     EE[j] = j;
	     DD[j] = - (q);
	     FF[j] = j;
	     if ( nseq == 2 || j > mm )
                RR[j] = SS[j] = mm + 1;
	     else
                RR[j] = SS[j] = j;
	  }

        for ( i = mm; i >= m1; i-- )
	  {  c = p = 0;
	     f = - (q);
	     ci = fi = i;
	     pi = i + 1;
	     cj = fj = pj = nn + 1;
	     va = v[A[i]];
	     if ( nseq == 2 || n1 > i )
		limit = n1;
	     else
		limit = i + 1;
	     for ( j = nn; j >= limit ; j-- )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  DIAG(i, j, c, p+va[B[j]])		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > min )
		    flag = 1;
	        }
	     if ( nseq == 2 || i < n1 )
	       { HH[i] = CC[n1];
	         II[i] = RR[n1];
	         JJ[i] = EE[n1];
	         WW[i] = f;
	         XX[i] = fi;
	         YY[i] = fj;
	       }
	  }
      
  for ( rl = m1, cl = n1; ; )
    { for ( rflag = cflag = 1; ( rflag && m1 > 1 ) || ( cflag && n1 > 1 ) ;  )
        { if ( rflag && m1 > 1 )	/* Compute one row */
            { rflag = 0;
	      m1--;
      	      c = p = 0;
	      f = - (q);
	      ci = fi = m1;
	      pi = m1 + 1;
	      cj = fj = pj = nn + 1;
	      va = v[A[m1]];
	      for ( j = nn; j >= n1 ; j-- )  
	        { f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  DIAG(m1, j, c, p+va[B[j]])		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = m1; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > min )
		     flag = 1;
		  if ( ! rflag && ( ci > rl && cj > cl || di > rl && dj > cl
	 		                            || fi > rl && fj > cl ) )
		      rflag = 1;
	        }
	      HH[m1] = CC[n1];
	      II[m1] = RR[n1];
	      JJ[m1] = EE[n1];
	      WW[m1] = f;
	      XX[m1] = fi;
	      YY[m1] = fj;
	      if ( ! cflag && ( ci > rl && cj > cl || di > rl && dj > cl
			     || fi > rl && fj > cl ) )
	         cflag = 1;
	    }

	  if ( nseq == 1 && n1 == (m1 + 1) && ! rflag )
	     cflag = 0;
	  if ( cflag && n1 > 1 )	/* Compute one column */
	    { cflag = 0;
	      n1--;
	      c = 0;
	      f = - (q);
	      cj = fj = n1;
	      va = v[B[n1]];
	      if ( nseq == 2 || mm < n1 )
		{ p = 0;
	          ci = fi = pi = mm + 1;
	          pj = n1 + 1;
		  limit = mm;
		}
	      else
		{ p = HH[n1];
		  pi = II[n1];
		  pj = JJ[n1];
	          ci = fi = n1;
		  limit = n1 - 1;
		}
	      for ( i = limit; i >= m1 ; i-- )  
	        { f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = HH[i] - qr; 
		  ci = II[i];
		  cj = JJ[i];
		  d = WW[i] - r;
		  di = XX[i];
		  dj = YY[i];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
	          DIAG(i, n1, c, p+va[A[i]])
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = n1; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = HH[i];
		  HH[i] = c;
		  pi = II[i];
		  pj = JJ[i];
		  II[i] = ci;
		  JJ[i] = cj;
		  WW[i] = d;
		  XX[i] = di;
		  YY[i] = dj;
		  if ( c > min )
		     flag = 1;
	          if ( ! cflag && ( ci > rl && cj > cl || di > rl && dj > cl
		               || fi > rl && fj > cl ) )
		     cflag = 1;
	        }
	      CC[n1] = HH[m1];
	      RR[n1] = II[m1];
	      EE[n1] = JJ[m1];
	      DD[n1] = f;
	      SS[n1] = fi;
	      FF[n1] = fj;
	      if ( ! rflag && ( ci > rl && cj > cl || di > rl && dj > cl
		                                 || fi > rl && fj > cl ) )
	         rflag = 1;
	    }
	}
      if ( m1 == 1 && n1 == 1 || no_cross() )
	 break;
   }
  m1--;
  n1--;
}

/* recompute the area on forward pass */
small_pass(A,B,count,nseq) char A[], B[]; int count, nseq;
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  *va;				/* pointer to v(A[i], B[j]) */
  int   addnode();			/* function for inserting a node */
  int  limit;				/* lower bound on j */

	for ( j = n1 + 1; j <= nn ; j++ )
	  {  CC[j] = 0;
	     RR[j] = m1;
	     EE[j] = j;
	     DD[j] = - (q);
	     SS[j] = m1;
	     FF[j] = j;
	  }
	for ( i = m1 + 1; i <= mm; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     ci = fi = i;
	     va = v[A[i]];
	     if ( nseq == 2 || i <= n1 )
	       { p = 0;
	         pi = i - 1;
	         cj = fj = pj = n1;
		 limit = n1 + 1;
	       }
	     else
	       { p = CC[i];
		 pi = RR[i];
		 pj = EE[i];
	         cj = fj = i;
		 limit = i + 1;
	       }
	     for ( j = limit ; j <= nn ; j++ )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  DIAG(i, j, c, p+va[B[j]])		/* diagonal */
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > min )	/* add the score into list */
		    min = addnode(c, ci, cj, i, j, count, min);
	        }
	  }
}

/* Add a new node into list.  */

int addnode(c, ci, cj, i, j, K, cost)  int c, ci, cj, i, j, K, cost;
{ short found;				/* 1 if the node is in LIST */
  register int d;

  found = 0;
  if ( most != 0 && most->STARI == ci && most->STARJ == cj )
    found = 1;
  else
     for ( d = 0; d < numnode ; d++ )
	{ most = LIST[d];
	  if ( most->STARI == ci && most->STARJ == cj )
	    { found = 1;
	      break;
	    }
        }
  if ( found )
    { if ( most->SCORE < c )
        { most->SCORE = c;
          most->ENDI = i;
          most->ENDJ = j;
        }
      if ( most->TOP > i ) most->TOP = i;
      if ( most->BOT < i ) most->BOT = i;
      if ( most->LEFT > j ) most->LEFT = j;
      if ( most->RIGHT < j ) most->RIGHT = j;
    }
  else
    { if ( numnode == K )	/* list full */
	 most = low;
      else
         most = LIST[numnode++];
      most->SCORE = c;
      most->STARI = ci;
      most->STARJ = cj;
      most->ENDI = i;
      most->ENDJ = j;
      most->TOP = most->BOT = i;
      most->LEFT = most->RIGHT = j;
    }
  if ( numnode == K )
    { if ( low == most || ! low ) 
        { for ( low = LIST[0], d = 1; d < numnode ; d++ )
            if ( LIST[d]->SCORE < low->SCORE )
              low = LIST[d];
	}
      return ( low->SCORE ) ;
    }
  else
    return cost;
}

/* Find and remove the largest score in list */

vertexptr findmax()
{ vertexptr  cur;
  register int i, j;

  for ( j = 0, i = 1; i < numnode ; i++ )
    if ( LIST[i]->SCORE > LIST[j]->SCORE )
       j = i;
  cur = LIST[j];
  if ( j != --numnode )
    { LIST[j] = LIST[numnode];
      LIST[numnode] =  cur;
    }
  most = LIST[0];
  if ( low == cur ) low = LIST[0];
  return ( cur );
}

/* return 1 if no node in LIST share vertices with the area */

no_cross()
{ vertexptr  cur;
  register int i;

      for ( i = 0; i < numnode; i++ )
	{ cur = LIST[i];
	  if ( cur->STARI <= mm && cur->STARJ <= nn && cur->BOT >= m1-1 && 
	       cur->RIGHT >= n1-1 && ( cur->STARI < rl || cur->STARJ < cl ))
	     { if ( cur->STARI < rl ) rl = cur->STARI;
	       if ( cur->STARJ < cl ) cl = cur->STARJ;
	       flag = 1;
	       break;
	     }
	}
      if ( i == numnode )
	return 1;
      else
	return 0;
}

/* diff(A,B,M,N,tb,te) returns the score of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

int diff(A,B,M,N,tb,te) char *A, *B; int M, N; int tb, te;

{ int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;

{ register int   i, j;
  register int c, e, d, s;
           int t, *va;
  	   char  *ckalloc();

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      return - gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
          return - gap(N);
        }
      if (tb > te) tb = te;
      midc = - (tb + r + gap(N) );
      midj = 0;
      va = v[A[1]];
      for (j = 1; j <= N; j++)
        {  for ( tt = 1, z = row[I+1]; z != 0; z = z->NEXT )	
              if ( z->COL == j+J )			
	         { tt = 0; break; }		
           if ( tt )			
            { c = va[B[j]] - ( gap(j-1) + gap(N-j) );
              if (c > midc)
               { midc = c;
                 midj = j;
               }
	    }
	}
      if (midj == 0)
        { INS(N) DEL(1) }
      else
        { if (midj > 1) INS(midj-1)
          REP
	  if ( A[1] == B[midj] || isdna && va[B[midj]] == MATCHSC )
	     no_mat += 1;
	  else
	     no_mis += 1;
	  /* mark (A[I],B[J]) as used: put J into list row[I] */	
          I++; J++;
	  z = ( pairptr ) ckalloc( (int) sizeof(pair));
          z->COL = J;			
          z->NEXT = row[I];				
	  row[I] = z;
          if (midj < N) INS(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  t = -q;
  for (j = 1; j <= N; j++)
    { CC[j] = t = t-r;
      DD[j] = t-q;
    }
  t = -tb;
  for (i = 1; i <= midi; i++)
    { s = CC[0];
      CC[0] = c = t = t-r;
      e = t-q;
      va = v[A[i]];
      for (j = 1; j <= N; j++)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	  DIAG(i+I, j+J, c, s+va[B[j]])
          if (c < d) c = d;
          if (c < e) c = e;
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  t = -q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  for (j = N-1; j >= 0; j--)
    { RR[j] = t = t-r;
      SS[j] = t-q;
    }
  t = -te;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      RR[N] = c = t = t-r;
      e = t-q;
      va = v[A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c = c - qr) > (e = e - r)) e = c;
          if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	  DIAG(i+1+I, j+1+J, c, s+va[B[j+1]])
          if (c < d) c = d;
          if (c < e) c = e;
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
    if ((c = DD[j] + SS[j] + q) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
}

/* Conquer: recursively around midpoint */

  if (type == 1)
    { (void) diff(A,B,midi,midj,tb,q);
      (void) diff(A+midi,B+midj,M-midi,N-midj,q,te);
    }
  else
    { (void) diff(A,B,midi-1,midj,tb,zero);
      DEL(2);
      (void) diff(A+midi+1,B+midj,M-midi-1,N-midj,zero,te);
    }
  return midc;
}

/* Alignment display routine */

static char ALINE[51], BLINE[51], CLINE[51];

display(A,B,M,N,S,AP,BP) char A[], B[]; int M, N; int S[], AP, BP;
{ register char *a, *b, *c;
  register int   i,  j, op;
           int   lines, ap, bp;

  i = j = op = lines = 0;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
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
          (void) printf("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
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
	char *malloc(), *p;

	if ((p = malloc( (unsigned) amount)) == NULL)
		fatal("Ran out of memory.");
	return(p);
}
