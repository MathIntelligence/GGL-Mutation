#ifndef _Gap
#define _Gap
#define  MATCHSC    10
#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */
#define gap2(k)  ((k) <= 0 ? 0 : ((k) <= gaplen ? q+r*(k) : pay))

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


class Gap
{
public:
Gap();
~Gap();
 
int match, mismh;		/* max and min substitution weights */
char *dhead1, *dhead2;           /* names of sequences */
int isdna;                       /* 1, DNA; 0, protein */
int v[128][128];	/* substitution scores */
int  q, r;       /* gap penalties */
int  qr;         /* qr = q + r */
int  gaplen;     /* minimum length for constant-cost insertion */
int  pay;	/* constant-cost for long insertion */

int  change;	/* constant-cost for long insertion */
int *CC, *DD;			/* saving matrix scores */
int *RR, *SS;		 	/* saving start-points */
int  *S;				/* saving operations for diff */
int  zero = 0;				/* int type zero        */

int *sapp;				/* Current script append ptr */
int  last;				/* Last script op appended */

int no_mat; 				/* number of matches */ 
int no_mis; 				/* number of mismatches */ 
int al_len; 				/* length of alignment */
						/* Append "Delete k" op */
};

#endif
