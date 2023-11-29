#ifndef TROLLCONSTANTS
#define TROLLCONSTANTS


// PDB file parsing.

#define ATOMNAME_COLUMN 12
#define RESIDUENAME_COLUMN 17
#define CHAIN_COLUMN 21
#define RESIDUE_COLUMN 22
#define X_COORDINATE_COLUMN 30
#define Y_COORDINATE_COLUMN 38
#define Z_COORDINATE_COLUMN 46
#define BFACTOR_COLUMN 61
#define ICODE_COLUMN 26


// Selection oprations

#define SEL_UNION 0
#define SEL_INTERSECTION 1


//  Flags

#define YES 1
#define NO 0

// Surfaces

#define Surf_Save_Vertices 0



// Secondary structure.

#define ALPHA 0
#define THREE_TEN 1




// Area types

enum { 
A_Current,
A_Buried,
A_PercentBuried,
A_Extended,
A_Hydrophobic,
A_Polar,
A_BuriedHydrophobic,
A_BuriedPolar
} ;


#define TC_SP2 0
#define TC_SP3 1



// Operate on selection or not.

#define UseSelection 0
#define IgnoreSelection 1
#define TC_And 0
#define TC_Or 1

// Object types

#define TC_DefaultSubset 257
#define TC_Subset 258
#define TC_Transformation 259
#define TC_AtomList 300
#define TC_HighlightSequence 301
#define TC_BackgroundColor 302
#define TC_ViewState 303
#define TC_UpdateModelView 304
#define TC_CPK 0
#define TC_ByPotential 1
#define TC_Solid 4
#define TC_ByAccessibility 5
#define TC_ByGaussianCurvature 6
#define TC_ByDistance 7


// Graphics

#define TC_SphereDL 1

#define TC_CATrace ((1L)<<0)
#define TC_Wires ((1L)<<1)
#define TC_Spheres ((1L)<<2)
#define TC_Worm ((1L)<<3)
#define TC_BallAndStick ((1L)<<4)
#define TC_Show ((1L)<<5)
#define TC_CPKSurface ((1L)<<6)
#define TC_SheetBar ((1L)<<7)
#define TC_HelixTube ((1L)<<8)
#define TC_HBModel ((1L)<<9)
#define TC_MolecularSurface ((1L)<<10)
#define TC_Phimap ((1L)<<12)
#define TC_MolecularModel ((1L)<<13)
#define TC_Ribbon ((1L)<<14)
#define TC_Secondary ((1L)<<15)
#define TC_All (~0L)

#define TC_Helix 0
#define TC_Sheet 1

#define TC_Add 0
#define TC_Replace 1
#define TC_Remove 2

#define TC_CPK 0
#define TC_Single 1
#define TC_BFactor 2


// Sequence display

#define TC_Aligned 0
#define TC_Listed 1

#define TC_StructureView 0
#define TC_AlignmentView 1


// Math

#define PI (float)3.14159265358979393846


// Memory allocation

#define MAX_IDENTIFIERS 256


// File types;

#define PDB 0


#ifdef WIN32
#define separator '\\'
#else
#define separator '/'
#endif
#endif

