#ifndef _Matrix
#define _Matrix

//#ifndef _MMatrix
//#define _MMatrix

char *amino_acid_order = "ABCDEFGHIKLMNPQRSTVWXYZ";

short blosum30mt[]={
  4,
  0,  5,
 -3, -2, 17,
  0,  5, -3,  9,
  0,  0,  1,  1,  6,
 -2, -3, -3, -5, -4, 10,
  0,  0, -4, -1, -2, -3,  8,
 -2, -2, -5, -2,  0, -3, -3, 14,
  0, -2, -2, -4, -3,  0, -1, -2,  6,
  0,  0, -3,  0,  2, -1, -1, -2, -2,  4,
 -1, -1,  0, -1, -1,  2, -2, -1,  2, -2,  4,
  1, -2, -2, -3, -1, -2, -2,  2,  1,  2,  2,  6,
  0,  4, -1,  1, -1, -1,  0, -1,  0,  0, -2,  0,  8,
 -1, -2, -3, -1,  1, -4, -1,  1, -3,  1, -3, -4, -3, 11,
  1, -1, -2, -1,  2, -3, -2,  0, -2,  0, -2, -1, -1,  0,  8,
 -1, -2, -2, -1, -1, -1, -2, -1, -3,  1, -2,  0, -2, -1,  3,  8,
  1,  0, -2,  0,  0, -1,  0, -1, -1,  0, -2, -2,  0, -1, -1, -1,  4,
  1,  0, -2, -1, -2, -2, -2, -2,  0, -1,  0,  0,  1,  0,  0, -3,  2,  5,
  1, -2, -2, -2, -3,  1, -3, -3,  4, -2,  1,  0, -2, -4, -3, -1, -1,  1,  5,
 -5, -5, -2, -4, -1,  1,  1, -5, -3, -2, -2, -3, -7, -3, -1,  0, -3, -5, -3, 20,
  0, -1, -2, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0, -1,  0, -1,  0,  0,  0, -2, -1,
 -4, -3, -6, -1, -2,  3, -3,  0, -1, -1,  3, -1, -4, -2, -1,  0, -2, -1,  1,  5, -1,  9,
  0,  0,  0,  0,  5, -4, -2,  0, -3,  1, -1, -1, -1,  0,  4,  0, -1, -1, -3, -1,  0, -2,  4};


short blosum50mt[]={
  5,
 -2,  5,
 -1, -3, 13,
 -2,  5, -4,  8,
 -1,  1, -3,  2,  6,
 -3, -4, -2, -5, -3,  8,
  0, -1, -3, -1, -3, -4,  8,
 -2,  0, -3, -1,  0, -1, -2, 10,
 -1, -4, -2, -4, -4,  0, -4, -4,  5,
 -1,  0, -3, -1,  1, -4, -2,  0, -3,  6,
 -2, -4, -2, -4, -3,  1, -4, -3,  2, -3,  5,
 -1, -3, -2, -4, -2,  0, -3, -1,  2, -2,  3,  7,
 -1,  4, -2,  2,  0, -4,  0,  1, -3,  0, -4, -2,  7,
 -1, -2, -4, -1, -1, -4, -2, -2, -3, -1, -4, -3, -2, 10,
 -1,  0, -3,  0,  2, -4, -2,  1, -3,  2, -2,  0,  0, -1,  7,
 -2, -1, -4, -2,  0, -3, -3,  0, -4,  3, -3, -2, -1, -3,  1,  7,
  1,  0, -1,  0, -1, -3,  0, -1, -3,  0, -3, -2,  1, -1,  0, -1,  5,
  0,  0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  2,  5,
  0, -4, -1, -4, -3, -1, -4, -4,  4, -3,  1,  1, -3, -3, -3, -3, -2,  0,  5,
 -3, -5, -5, -5, -3,  1, -3, -3, -3, -3, -2, -1, -4, -4, -1, -3, -4, -3, -3, 15,
 -1, -1, -2, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, -2, -1, -1, -1,  0, -1, -3, -1,
 -2, -3, -3, -3, -2,  4, -3,  2, -1, -2, -1,  0, -2, -3, -1, -1, -2, -2, -1,  2, -1,  8,
 -1,  2, -3,  1,  5, -4, -2,  0, -3,  1, -3, -1,  0, -1,  4,  0,  0, -1, -3, -2, -1, -2,  5};

short blosum65mt[]={
  4,
 -2,  4,
  0, -3,  9,
 -2,  4, -4,  6,
 -1,  1, -4,  2,  5,
 -2, -3, -2, -4, -3,  6,
  0, -1, -3, -1, -2, -3,  6,
 -2,  0, -3, -1,  0, -1, -2,  8,
 -1, -3, -1, -3, -3,  0, -4, -3,  4,
 -1,  0, -3, -1,  1, -3, -2, -1, -3,  5,
 -2, -4, -1, -4, -3,  0, -4, -3,  2, -3,  4,
 -1, -3, -2, -3, -2,  0, -3, -2,  1, -2,  2,  6,
 -2,  3, -3,  1,  0, -3, -1,  1, -3,  0, -4, -2,  6,
 -1, -2, -3, -2, -1, -4, -2, -2, -3, -1, -3, -3, -2,  8,
 -1,  0, -3,  0,  2, -3, -2,  1, -3,  1, -2,  0,  0, -1,  6,
 -1, -1, -4, -2,  0, -3, -2,  0, -3,  2, -2, -2,  0, -2,  1,  6,
  1,  0, -1,  0,  0, -2,  0, -1, -2,  0, -3, -2,  1, -1,  0, -1,  4,
  0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,
  0, -3, -1, -3, -3, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4,
 -3, -4, -2, -5, -3,  1, -3, -2, -2, -3, -2, -2, -4, -4, -2, -3, -3, -3, -3, 10,
 -1, -1, -2, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, -2, -1, -1, -1, -1, -1, -2, -1,
 -2, -3, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -2, -2, -2, -2, -1,  2, -1,  7,
 -1,  1, -4,  1,  4, -3, -2,  0, -3,  1, -3, -2,  0, -1,  3,  0,  0, -1, -2, -3, -1, -2,  4};


short blosum80mt[]={
  7,
 -3,  6,
 -1, -6, 13,
 -3,  6, -7, 10,
 -2,  1, -7,  2,  8,
 -4, -6, -4, -6, -6, 10,
  0, -2, -6, -3, -4, -6,  9,
 -3, -1, -7, -2,  0, -2, -4, 12,
 -3, -6, -2, -7, -6, -1, -7, -6,  7,
 -1, -1, -6, -2,  1, -5, -3, -1, -5,  8,
 -3, -7, -3, -7, -6,  0, -7, -5,  2, -4,  6,
 -2, -5, -3, -6, -4,  0, -5, -4,  2, -3,  3,  9,
 -3,  5, -5,  2, -1, -6, -1,  1, -6,  0, -6, -4,  9,
 -1, -4, -6, -3, -2, -6, -5, -4, -5, -2, -5, -4, -4, 12,
 -2, -1, -5, -1,  3, -5, -4,  1, -5,  2, -4, -1,  0, -3,  9,
 -3, -2, -6, -3, -1, -5, -4,  0, -5,  3, -4, -3, -1, -3,  1,  9,
  2,  0, -2, -1, -1, -4, -1, -2, -4, -1, -4, -3,  1, -2, -1, -2,  7,
  0, -1, -2, -2, -2, -4, -3, -3, -2, -1, -3, -1,  0, -3, -1, -2,  2,  8,
 -1, -6, -2, -6, -4, -2, -6, -5,  4, -4,  1,  1, -5, -4, -4, -4, -3,  0,  7,
 -5, -8, -5, -8, -6,  0, -6, -4, -5, -6, -4, -3, -7, -7, -4, -5, -6, -5, -5, 16,
 -1, -3, -4, -3, -2, -3, -3, -2, -2, -2, -2, -2, -2, -3, -2, -2, -1, -1, -2, -5, -2,
 -4, -5, -5, -6, -5,  4, -6,  3, -3, -4, -2, -3, -4, -6, -3, -4, -3, -3, -3,  3, -3, 11,
 -2,  0, -7,  1,  6, -6, -4,  0, -6,  1, -5, -3, -1, -2,  5,  0, -1, -2, -4, -5, -1, -4,  6};

short pam60mt[]={
  5,
 -2,  5,
 -5, -9,  9,
 -2,  5,-10,  7,
 -1,  2,-10,  3,  7,
 -6, -8, -9,-11,-10,  8,
  0, -2, -7, -2, -2, -7,  6,
 -5,  0, -6, -2, -3, -4, -6,  8,
 -3, -4, -4, -5, -4, -1, -7, -6,  7,
 -5, -1,-10, -2, -3,-10, -5, -4, -4,  6,
 -4, -7,-11, -9, -7, -1, -8, -4,  0, -6,  6,
 -3, -6,-10, -7, -5, -2, -6, -7,  1,  0,  2, 10,
 -2,  5, -7,  2,  0, -6, -1,  1, -4,  0, -5, -6,  6,
  0, -4, -6, -5, -3, -7, -4, -2, -6, -4, -5, -6, -4,  7,
 -3, -1,-10, -1,  2, -9, -5,  2, -5, -1, -3, -2, -2, -1,  7,
 -5, -5, -6, -6, -6, -7, -7,  0, -4,  2, -6, -2, -3, -2,  0,  8,
  1,  0, -1, -2, -2, -5,  0, -4, -4, -2, -6, -4,  1,  0, -3, -2,  5,
  1, -2, -5, -3, -4, -6, -3, -5, -1, -2, -5, -2, -1, -2, -4, -4,  1,  6,
 -1, -5, -4, -6, -4, -5, -4, -5,  3, -6, -1,  0, -5, -4, -5, -5, -4, -1,  6,
-10, -8,-12,-11,-12, -3,-11, -5,-10, -8, -4, -9, -6,-10, -9,  0, -4, -9,-11, 13,
 -2, -3, -6, -3, -3, -5, -3, -3, -3, -3, -4, -3, -2, -3, -3, -4, -2, -2, -3, -8, -3,
 -6, -5, -2, -8, -7,  3,-10, -2, -4, -7, -5, -7, -3,-10, -8, -8, -5, -5, -5, -3, -5,  9,
 -2,  1,-10,  2,  5,-10, -3,  0, -4, -2, -5, -4, -1, -2,  6, -2, -3, -4, -5,-11, -3, -7,  5};

short pam120mt[]={
  3,
  0,  4,
 -3, -6,  9,
  0,  4, -7,  5,
  0,  3, -7,  3,  5,
 -4, -5, -6, -7, -7,  8,
  1,  0, -4,  0, -1, -5,  5,
 -3,  1, -4,  0, -1, -3, -4,  7,
 -1, -3, -3, -3, -3,  0, -4, -4,  6,
 -2,  0, -7, -1, -1, -7, -3, -2, -3,  5,
 -3, -4, -7, -5, -4,  0, -5, -3,  1, -4,  5,
 -2, -4, -6, -4, -3, -1, -4, -4,  1,  0,  3,  8,
 -1,  3, -5,  2,  1, -4,  0,  2, -2,  1, -4, -3,  4,
  1, -2, -4, -3, -2, -5, -2, -1, -3, -2, -3, -3, -2,  6,
 -1,  0, -7,  1,  2, -6, -3,  3, -3,  0, -2, -1,  0,  0,  6,
 -3, -2, -4, -3, -3, -5, -4,  1, -2,  2, -4, -1, -1, -1,  1,  6,
  1,  0,  0,  0, -1, -3,  1, -2, -2, -1, -4, -2,  1,  1, -2, -1,  3,
  1,  0, -3, -1, -2, -4, -1, -3,  0, -1, -3, -1,  0, -1, -2, -2,  2,  4,
  0, -3, -3, -3, -3, -3, -2, -3,  3, -4,  1,  1, -3, -2, -3, -3, -2,  0,  5,
 -7, -6, -8, -8, -8, -1, -8, -3, -6, -5, -3, -6, -4, -7, -6,  1, -2, -6, -8, 12,
 -1, -1, -4, -2, -1, -3, -2, -2, -1, -2, -2, -2, -1, -2, -1, -2, -1, -1, -1, -5, -2,
 -4, -3, -1, -5, -5,  4, -6, -1, -2, -5, -2, -4, -2, -6, -5, -5, -3, -3, -3, -2, -3,  8,
 -1,  2, -7,  3,  4, -6, -2,  1, -3, -1, -3, -2,  0, -1,  4, -1, -1, -2, -3, -7, -1, -5,  4};

short pam160mt[]={
  2,
  0,  3,
 -2, -4,  9,
  0,  3, -5,  4,
  0,  2, -5,  3,  4,
 -3, -4, -5, -6, -5,  7,
  1,  0, -3,  0,  0, -4,  4,
 -2,  1, -3,  0,  0, -2, -3,  6,
 -1, -2, -2, -3, -2,  0, -3, -3,  5,
 -2,  0, -5,  0, -1, -5, -2, -1, -2,  4,
 -2, -4, -6, -4, -3,  1, -4, -2,  2, -3,  5,
 -1, -3, -5, -3, -2,  0, -3, -3,  2,  0,  3,  7,
  0,  2, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  3,
  1, -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -2, -1,  5,
 -1,  1, -5,  1,  2, -5, -2,  2, -2,  0, -2, -1,  0,  0,  5,
 -2, -1, -3, -2, -2, -4, -3,  1, -2,  3, -3, -1, -1, -1,  1,  6,
  1,  0,  0,  0,  0, -3,  1, -1, -2, -1, -3, -2,  1,  1, -1, -1,  2,
  1,  0, -2, -1, -1, -3, -1, -2,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,
  0, -2, -2, -3, -2, -2, -2, -2,  3, -3,  1,  1, -2, -2, -2, -3, -1,  0,  4,
 -5, -5, -7, -6, -7, -1, -7, -3, -5, -4, -2, -4, -4, -5, -5,  1, -2, -5, -6, 12,
  0, -1, -3, -1, -1, -3, -1, -1, -1, -1, -2, -1,  0, -1, -1, -1,  0,  0, -1, -4, -1,
 -3, -3,  0, -4, -4,  5, -5,  0, -2, -4, -2, -3, -2, -5, -4, -4, -3, -3, -3, -1, -3,  8,
  0,  2, -5,  2,  3, -5, -1,  1, -2,  0, -3, -2,  1, -1,  3,  0, -1, -1, -2, -6, -1, -4,  3};

short pam250mt[]={
  2,
  0,  3,
 -2, -4, 12,
  0,  3, -5,  4,
  0,  3, -5,  3,  4,
 -3, -4, -4, -6, -5,  9,
  1,  0, -3,  1,  0, -5,  5,
 -1,  1, -3,  1,  1, -2, -2,  6,
 -1, -2, -2, -2, -2,  1, -3, -2,  5,
 -1,  1, -5,  0,  0, -5, -2,  0, -2,  5,
 -2, -3, -6, -4, -3,  2, -4, -2,  2, -3,  6,
 -1, -2, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6,
  0,  2, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  2,
  1, -1, -3, -1, -1, -5,  0,  0, -2, -1, -3, -2,  0,  6,
  0,  1, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4,
 -2, -1, -4, -1, -1, -4, -3,  2, -2,  3, -3,  0,  0,  0,  1,  6,
  1,  0,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1, -1,  0,  2,
  1,  0, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,
  0, -2, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4,
 -6, -5, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17,
  0, -1, -3, -1, -1, -2, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1,  0,  0, -1, -4, -1,
 -3, -3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, -2, 10,
  0,  2, -5,  3,  3, -5,  0,  2, -2,  0, -3, -2,  1,  0,  3,  0,  0, -1, -2, -6, -1, -4,  3};

short pam350mt[]={
  2,
  1,  3,
 -2, -5, 18,
  1,  3, -6,  4,
  1,  3, -6,  4,  4,
 -4, -5, -5, -6, -6, 13,
  2,  1, -4,  1,  1, -6,  5,
 -1,  1, -4,  1,  1, -2, -2,  7,
  0, -2, -3, -2, -2,  2, -2, -2,  5,
 -1,  1, -6,  1,  0, -6, -1,  1, -2,  5,
 -2, -4, -7, -4, -4,  3, -4, -2,  4, -3,  8,
 -1, -2, -6, -3, -2,  1, -3, -2,  3,  0,  5,  6,
  0,  2, -4,  2,  2, -4,  1,  2, -2,  1, -3, -2,  2,
  1,  0, -3,  0,  0, -5,  0,  0, -2, -1, -3, -2,  0,  6,
  0,  2, -6,  2,  3, -5, -1,  3, -2,  1, -2, -1,  1,  1,  4,
 -1,  0, -4, -1,  0, -5, -2,  2, -2,  4, -3,  0,  1,  0,  2,  7,
  1,  1,  0,  1,  0, -4,  1, -1, -1,  0, -3, -2,  1,  1,  0,  0,  1,
  1,  0, -2,  0,  0, -3,  1, -1,  0,  0, -2, -1,  1,  1,  0, -1,  1,  2,
  0, -2, -2, -2, -2, -1, -1, -2,  4, -2,  3,  2, -2, -1, -2, -3, -1,  0,  5,
 -7, -6,-10, -8, -8,  1, -8, -3, -6, -4, -2, -5, -5, -7, -5,  4, -3, -6, -7, 27,
  0,  0, -3, -1,  0, -2, -1,  0,  0, -1, -1,  0,  0,  0,  0, -1,  0,  0,  0, -5, -1,
 -4, -4,  1, -5, -5, 11, -6,  0,  0, -5,  0, -2, -3, -6, -5, -5, -3, -3, -2,  1, -2, 14,
  0,  2, -6,  3,  3, -6,  0,  2, -2,  1, -3, -2,  2,  0,  3,  1,  0,  0, -2, -7,  0, -5,  3};

short md_40mt[]={
  9,
  0,  0,
 -7,  0, 16,
 -6,  0,-13, 11,
 -5,  0,-15,  3, 11,
-11,  0, -5,-15,-16, 13,
 -3,  0, -7, -4, -4,-15, 10,
 -9,  0, -6, -4, -8, -7,-10, 14,
 -6,  0,-11,-12,-12, -5,-13,-11, 11,
 -8,  0,-12, -8, -3,-16, -9, -6,-11, 11,
 -9,  0,-10,-14,-13, -1,-14, -7, -1,-12,  9,
 -6,  0, -9,-12,-11, -7,-12, -9,  1, -7,  1, 14,
 -6,  0, -8,  1, -5,-12, -5,  0, -8, -1,-12, -9, 12,
 -2,  0,-11,-11,-11,-11, -9, -4,-11,-10, -5,-10, -9, 12,
 -7,  0,-12, -6,  0,-14, -9,  2,-12, -1, -6, -8, -5, -3, 12,
 -7,  0, -5,-10, -8,-15, -4,  0,-10,  3, -9, -8, -6, -6,  0, 11,
  0,  0, -2, -6, -8, -6, -2, -6, -8, -7, -7, -8,  1, -1, -7, -5,  9,
  1,  0, -7, -8, -8,-11, -7, -7, -2, -5, -9, -2, -2, -4, -7, -6,  1, 10,
 -1,  0, -7, -9, -8, -6, -8,-12,  4,-12, -2,  0,-10, -9,-11,-11, -7, -4, 10,
-14,  0, -4,-15,-15, -7, -7,-13,-13,-13, -8,-11,-14,-14,-11, -4, -9,-12,-10, 18,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
-13,  0, -2, -8,-14,  2,-13,  2, -9,-13, -9,-11, -6,-13, -9,-10, -7,-10,-11, -6,  0, 14,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

short md_120mt[]={
  6,
  0,  0,
 -3,  0, 14,
 -2,  0, -7,  8,
 -2,  0, -8,  5,  8,
 -6,  0, -2, -9,-10, 11,
  0,  0, -3,  0, -1, -9,  8,
 -4,  0, -2, -1, -3, -2, -4, 11,
 -1,  0, -5, -7, -7, -1, -6, -6,  7,
 -4,  0, -6, -2,  0, -9, -4, -1, -6,  8,
 -4,  0, -5, -8, -8,  2, -8, -4,  2, -6,  7,
 -2,  0, -5, -7, -6, -2, -6, -5,  3, -4,  3, 10,
 -1,  0, -3,  3, -1, -6, -1,  2, -4,  1, -6, -5,  8,
  0,  0, -5, -5, -5, -5, -4, -1, -5, -4, -2, -5, -3,  9,
 -3,  0, -6, -1,  2, -7, -4,  4, -6,  2, -3, -4, -1,  0,  9,
 -3,  0, -2, -4, -3, -8, -1,  2, -6,  4, -5, -4, -2, -2,  2,  8,
  2,  0,  0, -2, -3, -3,  0, -2, -3, -3, -3, -3,  2,  1, -3, -2,  5,
  2,  0, -3, -3, -4, -6, -2, -3,  0, -2, -4,  0,  1,  0, -3, -3,  2,  6,
  1,  0, -3, -5, -5, -2, -4, -6,  5, -6,  1,  2, -5, -4, -6, -6, -3,  0,  7,
 -8,  0,  0, -9, -9, -3, -3, -6, -7, -6, -4, -6, -8, -8, -6, -1, -5, -7, -6, 17,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 -7,  0,  2, -4, -7,  5, -8,  4, -5, -7, -4, -6, -2, -7, -4, -5, -3, -6, -6, -2,  0, 12,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

short md_250mt[]={
  2,
  0,  0,
 -1,  0, 11,
 -1,  0, -3,  5,
 -1,  0, -4,  4,  5,
 -3,  0,  0, -5, -5,  8,
  1,  0, -1,  1,  1, -5,  5,
 -2,  0,  0,  0,  0,  0, -2,  6,
  0,  0, -2, -3, -3,  0, -3, -3,  4,
 -1,  0, -3,  0,  1, -5, -1,  1, -3,  5,
 -1,  0, -2, -4, -4,  2, -4, -2,  2, -3,  5,
  0,  0, -2, -3, -3,  0, -3, -2,  3, -2,  3,  6,
  0,  0, -1,  2,  1, -3,  0,  1, -2,  1, -3, -2,  3,
  1,  0, -2, -2, -2, -2, -1,  0, -2, -1,  0, -2, -1,  6,
 -1,  0, -3,  0,  2, -4, -1,  3, -3,  2, -2, -2,  0,  0,  5,
 -1,  0, -1, -1,  0, -4,  0,  2, -3,  4, -3, -2,  0, -1,  2,  5,
  1,  0,  1,  0, -1, -2,  1, -1, -1, -1, -2, -1,  1,  1, -1, -1,  2,
  2,  0, -1, -1, -1, -2,  0, -1,  1, -1, -1,  0,  1,  1, -1, -1,  1,  2,
  1,  0, -2, -3, -2,  0, -2, -3,  4, -3,  2,  2, -2, -1, -3, -3, -1,  0,  4,
 -4,  0,  1, -5, -5, -1, -1, -3, -4, -3, -2, -3, -4, -4, -3,  0, -3, -4, -3, 15,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 -3,  0,  2, -2, -4,  5, -4,  4, -2, -3, -1, -3, -1, -3, -2, -2, -1, -3, -3,  0,  0,  9,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

short md_350mt[]={
  1,
  0,  0,
  0,  0,  9,
  0,  0, -2,  3,
  0,  0, -2,  3,  3,
 -2,  0,  1, -3, -4,  6,
  1,  0,  0,  1,  1, -3,  4,
 -1,  0,  0,  0,  0,  0, -1,  3,
  0,  0, -1, -2, -2,  1, -2, -2,  3,
 -1,  0, -1,  0,  1, -3,  0,  1, -2,  3,
 -1,  0, -1, -3, -3,  2, -2, -1,  2, -2,  3,
  0,  0, -1, -2, -2,  1, -2, -1,  2, -2,  2,  3,
  0,  0, -1,  1,  1, -2,  0,  1, -1,  1, -2, -1,  2,
  1,  0, -1, -1, -1, -2, -1,  0, -1, -1,  0, -1,  0,  4,
 -1,  0, -2,  1,  1, -2,  0,  2, -2,  2, -1, -1,  0,  0,  3,
 -1,  0,  0,  0,  0, -3,  0,  1, -2,  3, -2, -1,  0,  0,  2,  3,
  1,  0,  0,  0,  0, -1,  1,  0, -1,  0, -1, -1,  1,  1,  0,  0,  1,
  1,  0,  0,  0, -1, -1,  0, -1,  0,  0, -1,  0,  0,  1, -1,  0,  1,  1,
  0,  0, -1, -2, -2,  0, -1, -2,  2, -2,  1,  2, -1, -1, -2, -2,  0,  0,  2,
 -3,  0,  1, -4, -3,  0, -1, -2, -3, -2, -1, -2, -3, -3, -2,  0, -2, -3, -2, 14,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 -2,  0,  2, -2, -2,  5, -3,  3, -1, -2,  0, -1, -1, -2, -1, -1, -1, -2, -2,  0,  0,  7,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

short idmat[]={
10,
 0, 10,
 0, 0, 10,
 0, 0, 0, 10,
 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,10,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,10};

short gon40mt[]={
  92,
   0,   0,
 -31,   0, 163,
 -56,   0,-135, 111,
 -37,   0,-140,  16, 105,
 -92,   0, -64,-152,-143, 126,
 -32,   0, -91, -51, -76,-152, 105,
 -65,   0, -67, -41, -40, -50, -81, 145,
 -76,   0, -87,-150,-106, -39,-158, -94, 104,
 -54,   0,-132, -47, -13,-127, -79, -34, -86, 103,
 -68,   0, -85,-155,-108, -13,-141, -85,   5, -85,  89,
 -45,   0, -63,-130, -80, -16,-114, -60,  10, -57,  16, 140,
 -62,   0, -83,   6, -38,-104, -40,  -7, -99, -20,-112, -91, 115,
 -37,   0,-137, -69, -60,-128, -87, -71,-108, -62, -83,-119, -78, 124,
 -43,   0,-113, -32,  10,-100, -71,   0, -91,   2, -60, -35, -25, -46, 118,
 -61,   0, -86, -77, -50,-130, -69, -31,-103,  19, -84, -81, -47, -73,  -6, 112,
   0,   0, -35, -36, -41,-111, -37, -48, -95, -43, -95, -64, -11, -35, -35, -51,  99,
 -25,   0, -59, -47, -52, -90, -85, -46, -51, -34, -78, -44, -27, -42, -39, -52,  13, 100,
 -22,   0, -43,-133, -74, -58,-122, -98,  28, -82, -18, -22,-103, -86, -79, -88, -74, -25,  97,
-120,   0, -68,-171,-131,  -6,-108, -70, -93,-127, -71, -72,-119,-149, -87, -63, -98,-120,-115, 181,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -95,   0, -56, -98,-107,  31,-129,   5, -76, -88, -64, -66, -62,-106, -81, -75, -69, -87, -73,   1,   0, 135,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

short gon80mt[]={
  75,
   0,   0,
 -10,   0, 154,
 -31,   0, -93,  96,
 -17,   0, -94,  31,  88,
 -64,   0, -39,-111,-102, 114,
 -11,   0, -61, -26, -47,-115,  97,
 -39,   0, -43, -17, -17, -26, -53, 127,
 -43,   0, -54,-106, -73, -15,-114, -64,  86,
 -30,   0, -88, -21,   4, -89, -50, -12, -59,  85,
 -43,   0, -55,-109, -75,   7,-104, -57,  22, -58,  77,
 -26,   0, -39, -88, -53,   3, -83, -38,  25, -37,  31, 117,
 -34,   0, -55,  21, -13, -75, -18,   9, -71,  -2, -79, -62,  97,
 -16,   0, -93, -42, -35, -93, -58, -45, -75, -37, -58, -78, -48, 114,
 -22,   0, -76,  -9,  23, -70, -44,  14, -60,  17, -39, -19,  -6, -24,  95,
 -36,   0, -60, -44, -23, -90, -43, -10, -71,  33, -58, -53, -22, -45,  11,  97,
  14,   0, -15, -14, -19, -77, -16, -25, -62, -20, -64, -41,   5, -14, -15, -27,  78,
  -5,   0, -34, -24, -27, -62, -52, -24, -28, -15, -49, -25,  -7, -20, -18, -27,  25,  81,
  -6,   0, -21, -89, -51, -31, -86, -65,  41, -54,   3,   1, -69, -57, -51, -60, -43,  -9,  80,
 -87,   0, -43,-124, -98,  16, -81, -43, -63, -89, -44, -45, -86,-112, -62, -41, -72, -87, -80, 173,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -65,   0, -32, -69, -74,  49, -94,  21, -47, -60, -35, -37, -39, -76, -53, -50, -46, -58, -47,  23,   0, 123,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

short gon120mt[]={
  59,
   0,   0,
  -1,   0, 144,
 -18,   0, -69,  82,
  -9,   0, -68,  35,  72,
 -48,   0, -26, -87, -78, 102,
  -3,   0, -45, -14, -31, -92,  90,
 -26,   0, -31,  -7,  -6, -14, -37, 110,
 -27,   0, -36, -80, -55,  -3, -87, -48,  72,
 -19,   0, -64,  -8,  11, -67, -34,  -2, -44,  69,
 -30,   0, -39, -82, -57,  15, -82, -42,  28, -44,  66,
 -17,   0, -26, -64, -40,  11, -65, -28,  29, -27,  34,  95,
 -20,   0, -41,  26,  -1, -58,  -7,  14, -55,   5, -61, -46,  80,
  -6,   0, -68, -28, -22, -72, -41, -31, -56, -24, -44, -56, -32, 105,
 -12,   0, -56,   1,  25, -53, -30,  17, -43,  20, -30, -14,   1, -14,  74,
 -23,   0, -45, -27, -10, -68, -30,  -1, -53,  36, -44, -38, -10, -30,  16,  83,
  16,   0,  -7,  -5,  -9, -58,  -6, -14, -44, -10, -47, -29,  10,  -5,  -7, -15,  60,
   2,   0, -21, -13, -15, -47, -35, -14, -17,  -6, -34, -16,   0, -10,  -9, -16,  26,  64,
   0,   0, -11, -65, -38, -17, -65, -47,  42, -39,  13,  10, -50, -42, -36, -44, -28,  -3,  65,
 -68,   0, -29, -96, -78,  27, -66, -28, -46, -68, -29, -31, -68, -89, -49, -30, -57, -67, -59, 166,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -48,   0, -20, -53, -56,  55, -74,  26, -31, -44, -20, -22, -28, -59, -38, -37, -35, -42, -33,  33,   0, 111,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

short gon160mt[]={
  46,
   0,   0,
   3,   0, 135,
 -11,   0, -53,  70,
  -4,   0, -52,  34,  59,
 -38,   0, -18, -70, -62,  91,
   2,   0, -34,  -7, -21, -76,  82,
 -18,   0, -23,  -1,  -1,  -7, -27,  93,
 -18,   0, -25, -62, -43,   3, -70, -37,  59,
 -12,   0, -48,  -1,  13, -53, -24,   2, -35,  55,
 -22,   0, -29, -65, -45,  19, -67, -32,  30, -34,  57,
 -12,   0, -19, -50, -31,  14, -52, -21,  29, -21,  34,  76,
 -12,   0, -31,  26,   5, -47,  -2,  15, -44,   8, -48, -36,  65,
  -1,   0, -52, -19, -14, -58, -30, -22, -43, -16, -35, -42, -22,  96,
  -7,   0, -42,   6,  23, -41, -21,  17, -32,  20, -24, -12,   5,  -8,  56,
 -16,   0, -35, -16,  -3, -53, -21,   3, -41,  35, -35, -29,  -4, -21,  17,  71,
  16,   0,  -2,   0,  -3, -45,  -1,  -8, -33,  -4, -36, -23,  11,   0,  -2,  -9,  44,
   5,   0, -14,  -6,  -8, -36, -24,  -8, -12,  -2, -24, -11,   3,  -4,  -4,  -9,  23,  50,
   1,   0,  -6, -49, -30,  -8, -52, -35,  40, -30,  17,  14, -38, -32, -27, -34, -20,   0,  53,
 -55,   0, -21, -78, -64,  32, -55, -19, -34, -54, -20, -22, -55, -74, -40, -24, -47, -54, -45, 158,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -37,   0, -13, -42, -44,  56, -60,  27, -20, -35, -11, -13, -22, -48, -29, -29, -28, -32, -24,  38,   0, 100,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

short gon250mt[]={
  24,
   0,   0,
   5,   0, 115,
  -3,   0, -32,  47,
   0,   0, -30,  27,  36,
 -23,   0,  -8, -45, -39,  70,
   5,   0, -20,   1,  -8, -52,  66,
  -8,   0, -13,   4,   4,  -1, -14,  60,
  -8,   0, -11, -38, -27,  10, -45, -22,  40,
  -4,   0, -28,   5,  12, -33, -11,   6, -21,  32,
 -12,   0, -15, -40, -28,  20, -44, -19,  28, -21,  40,
  -7,   0,  -9, -30, -20,  16, -35, -13,  25, -14,  28,  43,
  -3,   0, -18,  22,   9, -31,   4,  12, -28,   8, -30, -22,  38,
   3,   0, -31,  -7,  -5, -38, -16, -11, -26,  -6, -23, -24,  -9,  76,
  -2,   0, -24,   9,  17, -26, -10,  12, -19,  15, -16, -10,   7,  -2,  27,
  -6,   0, -22,  -3,   4, -32, -10,   6, -24,  27, -22, -17,   3,  -9,  15,  47,
  11,   0,   1,   5,   2, -28,   4,  -2, -18,   1, -21, -14,   9,   4,   2,  -2,  22,
   6,   0,  -5,   0,  -1, -22, -11,  -3,  -6,   1, -13,  -6,   5,   1,   0,  -2,  15,  25,
   1,   0,   0, -29, -19,   1, -33, -20,  31, -17,  18,  16, -22, -18, -15, -20, -10,   0,  34,
 -36,   0, -10, -52, -43,  36, -40,  -8, -18, -35,  -7, -10, -36, -50, -27, -16, -33, -35, -26, 142,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -22,   0,  -5, -28, -27,  51, -40,  22,  -7, -21,   0,  -2, -14, -31, -17, -18, -19, -19, -11,  41,   0,  78,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

short gon300mt[]={
  16,
   0,   0,
   5,   0, 104,
  -1,   0, -24,  37,
   1,   0, -23,  23,  27,
 -18,   0,  -5, -37, -31,  60,
   5,   0, -15,   3,  -4, -42,  58,
  -6,   0, -10,   5,   4,   0, -10,  45,
  -6,   0,  -7, -30, -21,  11, -36, -16,  33,
  -2,   0, -21,   6,  11, -26,  -7,   5, -17,  24,
  -9,   0, -10, -32, -22,  19, -36, -14,  25, -17,  33,
  -5,   0,  -6, -24, -16,  15, -28, -10,  22, -11,  24,  31,
  -1,   0, -14,  18,   9, -25,   5,  10, -22,   8, -24, -17,  27,
   3,   0, -23,  -4,  -2, -30, -11,  -8, -20,  -3, -18, -19,  -6,  66,
  -1,   0, -18,   9,  14, -20,  -6,   9, -15,  13, -13,  -8,   7,  -1,  18,
  -4,   0, -17,   0,   5, -25,  -6,   6, -19,  22, -18, -13,   4,  -6,  13,  37,
   8,   0,   1,   5,   3, -22,   4,  -1, -14,   2, -17, -11,   7,   4,   2,   0,  15,
   5,   0,  -3,   1,   1, -17,  -7,  -1,  -4,   2,  -9,  -5,   4,   2,   1,  -1,  11,  17,
   0,   0,   1, -23, -15,   4, -26, -15,  26, -13,  17,  15, -17, -14, -12, -15,  -8,   0,  26,
 -29,   0,  -7, -42, -36,  36, -34,  -5, -13, -28,  -4,  -6, -30, -41, -23, -14, -27, -28, -19, 132,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -17,   0,  -3, -22, -22,  46, -33,  18,  -3, -17,   3,   1, -12, -25, -14, -14, -15, -15,  -7,  40,   0,  67,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

#endif
