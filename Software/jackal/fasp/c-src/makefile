LIB=../lib/libm.a

LIBFILES= ../obj/Rotamer.o ../obj/Mfmt.o ../obj/ProFix.o ../obj/Loopy.o ../obj/PdbFix.o ../obj/Dace.o ../obj/SegBed.o ../obj/Lattice.o ../obj/Res.o ../obj/Atm.o  ../obj/Chn.o ../obj/EngTable.o ../obj/ConTable.o\
         ../obj/Cell.o ../obj/Scap.o ../obj/Eng.o ../obj/Algn.o ../obj/NearDistCut.o\
         ../obj/Rcs.o ../obj/Pdb.o ../obj/Tres.o ../obj/HashTable.o ../obj/SimAtm.o\
         ../obj/Tatm.o  ../obj/Strhandler.o ../obj/DistBound.o ../obj/ModPdb.o ../obj/AtmGeom.o\
         ../obj/Qsort.o ../obj/Stralg.o ../obj/Fmps.o ../obj/HBondList.o ../obj/MapNew.o\
         ../obj/Segen.o ../obj/Rotate.o ../obj/ModAlgnFmt.o ../obj/DistSecond.o ../obj/EngCoeff.o\
         ../obj/Map.o  ../obj/Decoy.o ../obj/Chiangle.o ../obj/ModConstant.o ../obj/Vector.o ../obj/Icosahedron.o\
	 ../obj/Atom.o ../obj/Grid.o ../obj/SurfvSurface.o ../obj/LatSurfv.o  ../obj/GridNew.o ../obj/DistDatabase.o\
	 ../obj/DistProtein.o ../obj/DistPopular.o ../obj/DistMutateBound.o ../obj/DistPopularBin.o\
	 ../obj/StrFmt.o ../obj/Bound.o ../obj/ModTop.o ../obj/Mutate.o ../obj/SurfSide.o ../obj/Contact.o\
	 ../obj/Dipole.o ../obj/Charge.o ../obj/Disc.o
         
CINCLUDE=-I.
COPT=-c -g


$(LIB): $(LIBFILES) 
	ar rv $(LIB) $?

$(LIBFILES):  ./$$(@F:.o=.cpp) ./$$(@F:.o=.h) *.h
	CC -g -o  $@ $(COPT) $(CINCLUDE) ../fasp/$(@F:.o=.cpp)

