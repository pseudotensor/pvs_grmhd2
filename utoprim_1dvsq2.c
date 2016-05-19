/*
-------------------------------------------------------------------------------
    Copyright 2005 Scott C. Noble, Charles F. Gammie, 
                   Jonathan C. McKinney, and Luca Del Zanna


    This file is part of PVS-GRMHD.

    PVS-GRMHD is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    PVS-GRMHD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PVS-GRMHD; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

-------------------------------------------------------------------------------
*/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_1dvsq2.c: 
---------------

    Uses the 1D^*_{v^2} method: 
       -- solves for one independent variable (v^2, or vsq) via a 1D 
          Newton-Raphson method 
       -- like the 1D_{v^2} method, except can be used (in principle)
          with a general equation of state. The main difference is how
          it calculates the intermediary value of W: by performing a 
          nested set of Newton-Rapshon iterations to get the best W value
          for the current v^2 value.  

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

#include "u2p_util.h"

#define NEWT_DIM 1

/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D ;

FTYPE vsq_for_gnr2, W_for_gnr;

// Declarations: 
static FTYPE vsq_calc(FTYPE W);
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE vsq);

/**********************************************************************/
/******************************************************************

  Utoprim_1dvsq2():
  
  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may 
     wish to alter the translation as they see fit.  

     It assumes that on input/output:

              /  rho u^t           \
         U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))  
              |  T^t_\mu           |  
              \   B^i              /


             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

     ala HARM. 

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/

int Utoprim_1dvsq2(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
		   FTYPE gdet, FTYPE prim[NPR])
{

  static int Utoprim_new_body(FTYPE U[], FTYPE gcov[NDIM][NDIM], 
			      FTYPE gcon[NDIM][NDIM], FTYPE gdet,  FTYPE prim[]);

  FTYPE U_tmp[NPR], prim_tmp[NPR];
  int i, j, ret; 
  FTYPE alpha;


  /* First update the primitive B-fields */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] / gdet ;


  /* Set the geometry variables: */
  alpha = 1.0/sqrt(-gcon[0][0]);
  
  /* Transform the CONSERVED variables into the new system */
  U_tmp[RHO] = alpha * U[RHO] / gdet;
  U_tmp[UU]  = alpha * (U[UU] - U[RHO])  / gdet ;
  for( i = UTCON1; i <= UTCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }

  /* Transform the PRIMITIVE variables into the new system */
  for( i = 0; i < BCON1; i++ ) {
    prim_tmp[i] = prim[i];
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    prim_tmp[i] = alpha*prim[i];
  }

  ret = Utoprim_new_body(U_tmp, gcov, gcon, gdet, prim_tmp);

  /* Transform new primitive variables back if there was no problem : */ 
  if( ret == 0 ) {
    for( i = 0; i < BCON1; i++ ) {
      prim[i] = prim_tmp[i];
    }
  }

  return( ret ) ;

}


/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the 
        Newton-Raphson routine. 

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu    |
             \  alpha B^i        /



               /    rho        \
	prim = |    uu         |
               | \tilde{u}^i   |
               \  alpha B^i   /


return:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;
	     
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
	     3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], 
			    FTYPE gcon[NDIM][NDIM], FTYPE gdet,  FTYPE prim[NPR])
{
  static void func_1d_gnr(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			  FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static void func_1d_gnr2(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			   FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);

  static int general_newton_raphson( FTYPE x[], int n, 
				     void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
						    FTYPE [][NEWT_DIM], FTYPE *, 
						    FTYPE *, int) );

  static int gnr2( FTYPE x[], int n, 
		   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
				  FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int) );
  
  FTYPE x_1d[1];
  FTYPE QdotB,Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qsq,Qtcon[NDIM];
  FTYPE rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,utsq,vsq,tmpdiff ;
  int    i,j, retval, retval2, i_increase ;


  // Assume ok initially:
  retval = 0 ;

  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  lower_g(Bcon,gcov,Bcov) ;

  for(i=0;i<4;i++) Qcov[i] = U[QCOV0+i] ;
  raise_g(Qcov,gcon,Qcon) ;


  Bsq = 0. ;
  for(i=1;i<4;i++) Bsq += Bcon[i]*Bcov[i] ;

  QdotB = 0. ;
  for(i=0;i<4;i++) QdotB += Qcov[i]*Bcon[i] ;
  QdotBsq = QdotB*QdotB ;

  ncov_calc(gcon,ncov) ;
  raise_g(ncov,gcon,ncon);

  Qdotn = Qcon[0]*ncov[0] ;

  Qsq = 0. ;
  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;

  Qtsq = Qsq + Qdotn*Qdotn ;

  D = U[RHO] ;

  /* calculate W from last timestep and use  for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  gammasq = 1. + utsq ;
  gamma  = sqrt(gammasq);
	
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = D / gamma ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0,u) ;
  w = rho0 + u + p ;

  W_last = w*gammasq ;


  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*Bsq ) 
	    - QdotBsq*(2.*W_last + Bsq) ) <= W_last*W_last*(Qtsq-Bsq*Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }

  // Initialize independent variables for Newton-Raphson:
  W_for_gnr = W_last;
  x_1d[0] = vsq_for_gnr2 = 1. - 1. / gammasq ; 


  // Find vsq via Newton-Raphson:
  retval = general_newton_raphson( x_1d, 1, func_1d_gnr) ; 

  // Find W from this vsq:
  vsq_for_gnr2 = x_1d[0] ; 
  x_1d[0] = W_for_gnr;
  retval2 = gnr2( x_1d, 1, func_1d_gnr2 ) ; 
	
  W = x_1d[0];

  retval += retval2 ;

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  }
  else{

    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
      return(retval) ;
    }
  }

  // Calculate v^2 :
  vsq = vsq_calc(W) ;
  if( vsq >= 1. ) {
    retval = 4;
    return(retval) ;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = D * gtmp;

  w = W * (1. - vsq) ;
  p = pressure_rho0_w(rho0,w) ;
  u = w - (rho0 + p) ;

  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( (rho0 <= 0.) || (u <= 0.) ) { 
    retval = 5;
    return(retval) ;
  }

  prim[RHO] = rho0 ;
  prim[UU] = u ;


  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn;
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;
	
  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  /* done! */
  return(retval) ;
    
}
  

/**********************************************************************/
/****************************************************************************
   vsq_calc(): 
    
      -- evaluate v^2 (spatial, normalized velocity) from 
            W = \gamma^2 w 

****************************************************************************/
static FTYPE vsq_calc(FTYPE W)
{
	FTYPE Wsq,Xsq;
	
	Wsq = W*W ;
	Xsq = (Bsq + W) * (Bsq + W);

	return(  ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq) );
}

/**********************************************************************/
/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0] is physical, based upon its definition;
    -- "corrects" it if it is not physical;

*********************************************************************/

static void validate_x(FTYPE x[1], FTYPE x0[1] ) 
{
  
  FTYPE small = 1.e-10;

  x[0] = (x[0] > 1.0)    ?  ( 0.5*(x0[0] + 1.) )    : x[0];
  x[0] = (x[0] < -small) ?  ( 0.5*x0[0] )           : x[0];
  x[0] = fabs(x[0]);

  return;

}


/**********************************************************************/
/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, 
				   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
						  FTYPE [][NEWT_DIM], FTYPE *, 
						  FTYPE *, int) )
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df =  f = 1.;
  i_extra = doing_extra = 0;

  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  vsq_old = vsq = W = W_old = 0.;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */


    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific
    W_old = W;
    W = W_for_gnr;
    errx  = (W==0.) ?  fabs(W-W_old) : fabs((W-W_old)/W);

    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    // METHOD specific
    validate_x( x, x_old );


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*   before stopping                                                         */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    // See if we've done the extra iterations, or have done too many iterations:
    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) 
	|| (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0) ) {
    return(2);
  }

  // Return in different ways depending on whether a solution was found:
  if( fabs(errx) > MIN_NEWT_TOL){
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}



/************************************************************/
/********************************************************************** 

  gnr2()

    -- used for the new 1D^*_vsq2 method, much like general_newton_raphson()
       but with a few differences (like error-checking
       is different since errant W values do not equal errant vsq values);
       A separate NR routine is most likely not necessary (differences 
       can be handled by adding a "type" flag to the arg. list), but 
       one was made to prevent confusion during the debugging process.

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt(), and should 
       be nearly identical to general_newton_raphson() (above);

*****************************************************************/
static int gnr2( FTYPE x[], int n, 
			    void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
					 FTYPE [][NEWT_DIM],FTYPE *,FTYPE *,int) )
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id,jd, i_extra, doing_extra;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id] ;
    }

    /* Calculate the convergence criterion */
    for( id = 0; id < n ; id++) {
      errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]);
    }
    errx /= 1.*n;

    /* Make sure that the new x[] is physical : */
    // METHOD specific:
    x[0] = fabs(x[0]);


    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*   before stopping                                                         */
    if( (fabs(errx) <= NEWT_TOL2) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 0;
    }

    if( doing_extra == 1 ) i_extra++ ;

    // See if we've done the extra iterations, or have done too many iterations:
    if( ((fabs(errx) <= NEWT_TOL2)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }  


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0)  ) {
    return(2);
  }

  // Return in different ways depending on whether a solution was found:
  if( fabs(errx) > MIN_NEWT_TOL){
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}



/********************************************************************************/
/********************************************************************** 
   func_1d_gnr(): 

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/

static void func_1d_gnr(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE vsq,W,W0,dpdW,Wsq,W3,p_tmp,dWdvsq , dpdvsq, fact_tmp ;
  int retval, iters; 

  static void func_1d_gnr2(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			   FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
  static int gnr2( FTYPE x[], int n, 
		   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
				  FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int) );

  static FTYPE dpdvsq_calc(FTYPE W, FTYPE vsq);


  vsq = x[0];

  // Calculate best value for W given current guess for vsq: 
  vsq_for_gnr2 = vsq;
  W0           = W_for_gnr; 
  retval = 1 ; 
  iters = 0;
  while ( (retval != 0) && (++iters < 6) )  { 
    x[0]         = W0;
    retval = gnr2( x, 1, func_1d_gnr2 ) ;
    W0 *= 10.;
  }
    
  W = W_for_gnr = x[0]; 
  Wsq = W*W;
  W3 = W*Wsq;

  x[0] = vsq;


  dpdW   = dpdW_calc_vsq( W, vsq );
  dpdvsq = dpdvsq_calc(   W, vsq );

  dWdvsq = ( dpdvsq - 0.5*Bsq )  /  ( 1. - dpdW + QdotBsq/W3 ) ;

  fact_tmp = (Bsq + W) ;

  resid[0] = Qtsq  -  vsq * fact_tmp * fact_tmp  +  QdotBsq * ( Bsq + 2.*W ) / Wsq ; 
  jac[0][0] =  -fact_tmp * ( fact_tmp +   2. * dWdvsq * ( vsq + QdotBsq/W3 ) ) ; 

  dx[0] = -resid[0]/jac[0][0];

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

}

/**********************************************************************/
/********************************************************************** 
   func_1d_gnr2(): 

        -- calculates the residuals, and Newton step for gnr2();
        -- for this method, x = W here;

     Arguments:
          x  = current value of independent var's (on input & output);
         dx  = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac = Jacobian matrix based on x (on output);
         f  =  resid.resid/2  (on output)
        df  = -2*f;  (on output)
************************************************************************/

static void func_1d_gnr2(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			 FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE vsq, W,drdW,dpdW,Wsq,p_tmp;
  FTYPE  wtmp;
  int id;
  static FTYPE pressure_W_vsq(FTYPE W, FTYPE vsq) ;


  W = x[0];
  vsq = vsq_for_gnr2;
  
  p_tmp = pressure_W_vsq( W,  vsq) ;
  
  Wsq = W*W;
  dpdW = dpdW_calc_vsq( W, vsq );

    
  resid[0] = 
    +W 
    + 0.5 * Bsq * ( 1. + vsq )
    - 0.5*QdotBsq/Wsq
    + Qdotn
    - p_tmp;

  jac[0][0] = drdW = 1. - dpdW + QdotBsq / (Wsq*W) ;

  dx[0] = -resid[0]/jac[0][0];

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

}


/********************************************************************** 
 ********************************************************************** 
   
 The following routines specify the equation of state.  All routines 
  above here should be indpendent of EOS.  If the user wishes 
  to use another equation of state, the below functions must be replaced 
  by equivalent routines based upon the new EOS. 

 **********************************************************************
**********************************************************************/


/**********************************************************************/
/********************************************************************** 
  pressure_W_vsq():  
  
        -- Gamma-law equation of state;
        -- pressure as a function of W, vsq, and D:
 **********************************************************************/
static FTYPE pressure_W_vsq(FTYPE W, FTYPE vsq) 
{
  FTYPE gtmp;
  
  gtmp = 1. - vsq;
  
  return(  (GAMMA - 1.) * ( W * gtmp  -  D * sqrt(gtmp) ) / GAMMA  );
}



/**********************************************************************/
/********************************************************************** 
  dpdW_calc_vsq(): 
 
      -- partial derivative of pressure with respect to W;
 **********************************************************************/
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE vsq)
{
  return( (GAMMA - 1.) * (1. - vsq) /  GAMMA ) ;
}


/**********************************************************************/
/********************************************************************** 
  dpdvsq_calc(): 
 
      -- partial derivative of pressure with respect to vsq
 **********************************************************************/
static FTYPE dpdvsq_calc(FTYPE W, FTYPE vsq)
{
 return(  (GAMMA - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / GAMMA  ) ;
}


/****************************************************************************** 
             END   OF   UTOPRIM_1DVSQ2.C
 ******************************************************************************/
