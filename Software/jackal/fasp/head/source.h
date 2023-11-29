#ifndef _Source
#define _Source
#include<stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include<iostream>
#include <assert.h>
#include<unistd.h>
#include<sys/time.h>
#include <sys/resource.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<new.h>
#include<signal.h>
#include"Res.h"
#include"Chn.h"
#include"Atm.h"
#include"Rcs.h"
#include"SurfSide.h"
#include "ConTable.h"
#include"Pdb.h"
#include"Tres.h"
#include"Tatm.h"
#include"Lattice.h"
#include"Cell.h"
#include"Scap.h"
#include"Eng.h"
#include"Qsort.h" 
#include"Chiangle.h"
#include"Segen.h"
#include"Loopy.h"
#include"Map.h"
#include"Rotate.h"
#include"Decoy.h"
#include"Fmps.h"
#include"DistMutateBound.h"
#include"DistBound.h"
#include"Contact.h"
#include"ModPdb.h"
#include"Mutate.h"
#include"Stralg.h"
#include"ModAlgnFmt.h"
#include"NearDistCut.h"
#include"EngTable.h"
#include"Strhandler.h"
#include"AtmGeom.h"
#include "SimAtm.h"
#include"Algn.h"
#include"HashTable.h"
#include"HBondList.h"
#include"ModConstant.h"
#include"DistSecond.h"
#include "Vector.h"
#include"Icosahedron.h"
#include "SurfvSurface.h"
#include "LatSurfv.h"
#include "Atom.h"
#include "MapNew.h"
#include "Grid.h"
#include "GridNew.h"
#include "DistDatabase.h"
#include "DistProtein.h"
#include "DistPopular.h"
#include "EngCoeff.h"
#include "DistPopularBin.h"
#include "SegBed.h"
#include "StrFmt.h"
#include "Bound.h"
#include "ModTop.h"
#include "Dipole.h"
#include "Charge.h"
#include "Disc.h"
#include"Mfmt.h"
#include "PdbFix.h"
#include "ProFix.h"
#include "Rotamer.h"
#include "Dace.h"
extern "C" {
	char **cmain(char *,char *); 	
}

#ifndef _min
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#endif

#ifndef _max
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#endif

//#define max(x,y) ( (x) > (y) ? (x) : (y) )
//#define min(x,y) ( (x) < (y) ? (x) : (y) )
extern Tres TRES;
extern Rcs RCS;
#endif

