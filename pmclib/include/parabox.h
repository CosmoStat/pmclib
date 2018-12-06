/* Parameter boxes */

#ifdef __PLANCK__
#include "HL2_likely/tools/errorlist.h"
#include "HL2_likely/tools/io.h"
#else
#include "errorlist.h"
#include "io.h"
#endif

#ifndef __PARABOX__
#define __PARABOX__

/* errors */
#define pb_base       -12000
#define pb_allocate   -1 + pb_base
#define pb_infnan     -2 + pb_base


// 1) A parameter box is the intersection of half spaces (faces)
// 2) A parameter box is initialized as empty (no constraint)
// 3) Constraints (i.e. half spaces aka faces) are added one by one with function "add_face"
// 4) For convenience, functions add_lobound, add_hibound and add_slab 
//    can be used to add faces parallel to the coordinate axes.


typedef struct {
  int ndim ;      // The dimension of the parameter space
  int nfaces ;    // The number of constraints 
  double *faces;  // Pointer to an nfaces * (ndim+1) array of doubles
} parabox;


// Creates a parameter box for vectors of parameters of dimension "parasize"
parabox* init_parabox(int parasize,error **err) ;

// Returns 1 if the vector pointed by pos is in the box, 0 otherwise
int isinBox (const parabox* bx, double* pos,error **err);


// Adds a lower bound for a given coordinate (the face is parallel to one coord. axis)
void add_lobound(parabox* bx, int coord, double vmin ,error **err); 

// Adds an upper bound for a given coordinate (the face is parallel to one coord. axis)
void add_hibound(parabox* bx, int coord, double vmax ,error **err); 

// Adds a slab to the box: a pair of lower and upper bounds for one coordinate
void add_slab(parabox* bx, int coord, double vmin, double vmax ,error **err); 

// Adds an aribtrary (non axis-parallel) face to the box.
// A vector theta is in the box if
//   \sum_{i=1}^ndim  theta_i dirvec_i <= threshold 
void add_face(parabox* bx, double* dirvec, double threshold ,error **err); 

void add_halfspace(parabox* bx, int coord, double sign, double val ,error **err);

// Releases the box
void free_parabox(parabox **bx);

#endif
