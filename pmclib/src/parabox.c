/* Boxes of parameters */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __PLANCK__
#include "HL2_likely/pmclib/parabox.h"
#else
#include "parabox.h"
#endif



// Creates a parameter box for vectors of parameters of dimension "parasize"
parabox* init_parabox(int parasize,error **err) 
{
  parabox* newbox ;

  newbox = malloc_err(sizeof(parabox),err) ;
  forwardError(*err,__LINE__,NULL);

  newbox->ndim   = parasize ;
  newbox->nfaces = (int) 0  ;
  newbox->faces  = NULL ;
  return newbox ;
}


void add_face(parabox* bx, double* dirvec, double threshold ,error **err)
{
  int nfaces  = bx->nfaces ;
  int dimpara = bx->ndim   ;
  double* oldfaces = bx->faces ;
  double* newfaces = malloc_err( sizeof(double)*(nfaces+1)*(dimpara+1) ,err) ;
  forwardError(*err,__LINE__,);
  
  int i, offset  ;
  // copy old to new
  for (i=0; i< nfaces*(dimpara+1); i++)
    newfaces[i] = oldfaces[i] ;
  offset = nfaces*(dimpara+1) ;

  // add directional vector
  for (i= 0 ; i < dimpara ; i++) 
    newfaces[offset+i] = dirvec[i] ;

  // add thrshold
  newfaces[offset+dimpara] = threshold ;

  // Substitute new to old
  free(oldfaces);
  bx->faces  = newfaces ;
  bx->nfaces = 1+nfaces ;

  return  ; // unmitigated happiness
}


// Adds a lower bound for a given coordinate (the face is parallel to one coord. axis)
void add_lobound(parabox* bx, int coord, double vmin ,error **err)
{
  add_halfspace(bx, coord, -1.0, -vmin ,err);
	forwardError(*err,__LINE__,);
}

// Adds an upper bound for a given coordinate (the face is parallel to one coord. axis)
void add_hibound(parabox* bx, int coord, double vmax ,error **err)
{
  add_halfspace(bx, coord, (double) +1.0,  vmax ,err) ;
	forwardError(*err,__LINE__,);
}

// Adds a slab to the box: a pair of lower and upper bounds for one coordinate
void add_slab(parabox* bx, int coord, double vmin, double vmax ,error **err)
{

  // add min
	add_halfspace(bx, coord, (double) -1.0, -vmin ,err) ;
	forwardError(*err,__LINE__,);

  // add max
  add_halfspace(bx, coord, (double) +1.0,  vmax ,err) ;
	forwardError(*err,__LINE__,);

  return ;
}


void add_halfspace(parabox* bx, int coord, double sign, double val ,error **err)
{
  int dimpara = bx->ndim   ;
  double* dirvec ;

  dirvec = calloc_err( dimpara, sizeof(double) ,err) ;  // creates a zero vector
  forwardError(*err,__LINE__,);

  dirvec[coord] = sign ;  

  add_face(bx, dirvec, val,err) ;
  if (isError(*err)) {
     free(dirvec) ;
     forwardError(*err,__LINE__,);		
  }
  free(dirvec) ;

  return  ;
}


// Returns 1 if the vector pointed by pos is in the box, 0 otherwise
int isinBox(const parabox *bx, double* pos, error **err)
{
  int i, j, offset ;
  int dim = bx->ndim ;
  double* coeffs = bx->faces ;
  double threshold ;

  for (j=0; j<dim; j++) {
    testErrorRetVA(!finite(pos[j]), pb_infnan, "Parameter #%d is not finite", *err, __LINE__, 0, j);
  }

  offset = 0 ;
  for (i=0; i< bx->nfaces; i++) { // we go thru all constraints

    threshold = coeffs[offset+dim] ; // The weighted sum should not exceed this
    for (j=0; j<dim; j++)
      threshold -= coeffs[offset+j] * pos[j] ; // removing bit by bit 
    if (threshold < 0) {
       //printf("not in Box %d/%d\n", i, bx->nfaces);
       return 0 ;  // Any violated constraint causes immediate outrage
    }
    offset += dim+1 ;  // ready for next constraint
  }
  return 1 ;  // We got there, so all constraints are satisfied.
}


void free_parabox(parabox **bx)
{
  free((*bx)->faces);
  free(*bx) ;
	bx=NULL;
}


void print_parabox(parabox *bx,FILE* mstream)
{
  int i, j, offset ;
  int dim = bx->ndim ;
  double* coeffs = bx->faces ;
	FILE* out;
	
	out=mstream;
	if (mstream==NULL) {
		out=stdout;
	}
	
  //  printf("Dimension of parameter space..: %i\n", dim) ;
  // printf("Number of faces...............: %i\n", bx->nfaces) ;

  for (j=0; j<dim; j++) {
		fprintf(out,"----------"); 
	} 
	fprintf(out,"------------------\n");   		
  fprintf(out,"Face | ") ;
  for (j=0; j<dim; j++)
    fprintf(out,"w_%i       ", j ) ;
  fprintf(out,"| threshold \n");
  for (j=0; j<dim; j++) {
		fprintf(out,"----------");	
	} 
	fprintf(out,"------------------\n");   
  offset = 0 ;
  for (i=0; i< bx->nfaces; i++) {
    fprintf(out,"#%2i  |", i) ;
    for (j=0; j<dim; j++)
      fprintf(out," %+8.2e", coeffs[offset+j] ) ;
    fprintf(out," | %+8.2e\n", coeffs[offset+dim] ) ;
    offset += dim+1 ;
  }
  for (j=0; j<dim; j++) 
		fprintf(out,"----------"); 
	fprintf(out,"------------------\n");   
  fprintf(out,"theta is in box if \\sum_i theta_i w_i <= threshold for all faces\n");
}
