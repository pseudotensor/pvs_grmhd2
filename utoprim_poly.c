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

utoprim_poly.c: 
---------------

    Uses the Polynomial method: 
       -- solves for one independent variable (W) via Laguerre's method 
          on the 8th-order polynomial which is derived assuming a Gamma-law
          equation of state;

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/


#include "u2p_util.h"

#define NEWPOLY (0) 



#define POLY_ORDER 8               /* order of the polynomial that we need to solve */


#define MAX_ORDER 8                /* Maximum allowed Order of polynomial for which  
                                        we are root-finding */

#define REDUCED_ROOT_TOL NEWT_TOL  /* Thresh. below which the reduced poly. must be 
                                     for a guess to be chosen as a root */

#define ROOT_TOL NEWT_TOL          /* Threshold below which the poly. must be 
                                      for a guess to be chosen as a root*/

#define MIN_ROOT_TOL MIN_NEWT_TOL  /* Threshold below which the poly. must be for 
                                      a guess to be chosen as a root*/

#define IMAG_FLOOR 1.0e-10         /* Rel. thresh. below which the imag. part of 
                                      a root is considered to be zero */      

#define IMAG_FLR_PHYS  1.0e-6      /* in determining physical roots, rel. factor 
                                      to Re. part that Im. part must 
                                      be g.t. to not be considered negligible */

#define N_CYCLE_BREAKS 20          /* Total number of cycle breaks to attempt */

#define ROOT_MAX_ITER 30           /* Total number of Laguerre iterations */

#define N_EXTRA 2 


/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D ;

/* Used to divide "fake" roots from residual */
FTYPE roots[POLY_ORDER] = {0., 0., 0., 0., 0., 0., 0., 0.};  


/* Used to divide "fake" roots from residual */
FTYPE complex roots_z[POLY_ORDER] = {0., 0., 0., 0., 0., 0., 0., 0.};  

/* Coefficients for polynomial of residual1 */
FTYPE coefs[POLY_ORDER+1] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};  

/* Coefficients for polynomial of residual1 */
FTYPE complex coefs_z[POLY_ORDER+1] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};  


// Declarations: 
static FTYPE vsq_calc(FTYPE W);
static FTYPE utsq_calc(FTYPE W);

/**********************************************************************/
/******************************************************************

  Utoprim_poly():
  
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

int Utoprim_poly(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
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
  
  /* Calculate the transform the CONSERVATIVE variables into the new system */
  U_tmp[RHO] = alpha * U[RHO] / gdet;
  U_tmp[UU]  = alpha * (U[UU] - U[RHO])  / gdet ;
  for( i = UTCON1; i <= UTCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }

  /* Calculate the transform the PRIMITIVE variables into the new system */
  for( i = 0; i < BCON1; i++ ) {
    prim_tmp[i] = prim[i];
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    prim_tmp[i] = alpha*prim[i];
  }

  ret = Utoprim_new_body(U_tmp, gcov, gcon, gdet, prim_tmp);

  /* Transform new primitive variables back : */ 
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
  FTYPE x_1d[1];
  FTYPE QdotB,Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qsq,Qtcon[NDIM];
  FTYPE rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,utsq,vsq,tmpdiff ;
  int i,j, retval, i_increase ;

  static FTYPE find_root_1D(FTYPE guess) ;


  // Assume ok initially:
  retval = 0;

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

  /* calculate W from last timestep and use for guess */
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

  // METHOD specific:
  W = find_root_1D(0.);


  /* Problem with solver, so return denoting error before doing anything further */
  if( W == FAIL_VAL ) {
    retval = 100+1;
    return(retval);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
      return(retval) ;
    }
  }


  // Calculate utsq

  vsq = vsq_calc(W) ;
  if( vsq >= 1. ) {
    retval = 4;
    return(retval) ;
  }

  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = D * gtmp;

  w = W * (1. - vsq) ;
  p = pressure_rho0_w(rho0,w) ;
  u = w - (rho0 + p) ;

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
  return(0) ;

}

/*************************************************************************/
/************************************************************************* 

  residual0():

     -- a residual function derived from the definitions
        of the conserved variables, that should be satisfied by the 
        right conserved/primitive variable set;
   
     -- used to verify roots;

*************************************************************************/

static FTYPE residual0(FTYPE W) 
{
	FTYPE utsq,gamma,rho0,w,resid ;

	utsq = fabs(utsq_calc(W)) ;
	gamma = sqrt(1. + utsq) ;

	rho0 = D/gamma ;
	w = W/(gamma*gamma) ;

	resid = 
		+ 1.
		- QdotBsq/(2.*W*W*W)
		+ Bsq/W
		- Bsq/(2.*gamma*gamma*W)
		+ Qdotn/W
		- pressure_rho0_w(rho0,w)/W ;

	return(resid) ;

}

/***********************************************************************************/
/***********************************************************************************

   calc_poly_coefs(): 

       -- calculates the coefficients of the polynomial to be root-ed; (surprise!)

***********************************************************************************/

static void calc_poly_coefs( ) {

  FTYPE G, max_coef, cftmp; 
  int i;

  G = ( GAMMA - 1. ) / GAMMA ;

  coefs[0] = -4*Bsq*Bsq*G*G*QdotBsq*(Bsq*D*D + QdotBsq);

  coefs[1] = QdotBsq*(Bsq*(4 - 16*G)*G*QdotBsq + 
		      Bsq*Bsq*(Bsq*(-4*Bsq*G - 8*G*Qdotn) + 
			       G*(-16*D*D*G - 4*Qtsq)));

  coefs[2] = -((1 - 4*G)*(1 - 4*G)*QdotBsq*QdotBsq) + 
    Bsq*(QdotBsq*(2*Qtsq - 4*G*(5*D*D*G + 2*(1 + G)*Qtsq)) + 
	 Bsq*(-4*(-1 + 8*G)*QdotBsq*Qdotn - 
	      Qtsq*(4*D*D*G*G + Qtsq) + 
	      Bsq*((2 + 8*(-3 + G)*G)*QdotBsq - 4*Qdotn*Qtsq - 
		   Bsq*(-4*D*D*G*G + (Bsq + 2*Qdotn)*(Bsq + 2*Qdotn) 
			+ 2*Qtsq))));

  coefs[3] = -4*G*QdotBsq*(2*D*D*G + (-1 + 4*G)*Qtsq) + 
    Bsq*(-8*(-1 + 5*G)*QdotBsq*Qdotn - 
	 4*G*Qtsq*(2*D*D*G + Qtsq) + 
	 Bsq*(8*(1 + G*(-7 + 4*G))*QdotBsq - 8*(1 + G)*Qdotn*Qtsq + 
	      4*Bsq*(Bsq*Bsq*(-2 + G) + 4*D*D*G*G + 
		     2*Bsq*(-3 + G)*Qdotn - 2*(2*Qdotn*Qdotn + Qtsq))));

  coefs[4] = -4*(-1 + 4*G)*QdotBsq*Qdotn - 
    4*G*G*Qtsq*(D*D + Qtsq) + 
    Bsq*((10 + 8*G*(-7 + 5*G))*QdotBsq - 4*(1 + 4*G)*Qdotn*Qtsq + 
	 Bsq*(24*D*D*G*G - 24*Qdotn*Qdotn + 
	      Bsq*(Bsq*(-26 - 4*(-6 + G)*G) - 8*(7 - 4*G)*Qdotn) + 
	      (-10 + 8*(-1 + G)*G)*Qtsq));

  coefs[5] = 4*(-1 + G)*(-1 + 4*G)*QdotBsq - 8*G*Qdotn*Qtsq + 
    Bsq*(16*D*D*G*G - 16*Qdotn*Qdotn + 
	 Bsq*(Bsq*(-44 + 8*(7 - 2*G)*G) - (64 - 48*G)*Qdotn) + 
	 4*(-1 + 4*(-1 + G)*G)*Qtsq);

  coefs[6] = 4*D*D*G*G - 4*Qdotn*Qdotn + 
    Bsq*(Bsq*(-41 + 8*(8 - 3*G)*G) - 4*(9 - 8*G)*Qdotn) + 8*(-1 + G)*G*Qtsq;

  coefs[7] = -4*(-1 + G)*(Bsq*(-5 + 4*G) - 2*Qdotn);

  coefs[8] = -4*(-1 + G)*(-1 + G);


  /* Normalize coefficients to largest one: */

  max_coef = -1.e30;
  for( i = 0; i <= POLY_ORDER; i++ ) {
    cftmp = fabs(coefs[i]);
    if( cftmp > max_coef )  max_coef = cftmp;
  }
  for( i = 0; i <= POLY_ORDER; i++ )    coefs[i] /= max_coef;


  return; 
}


/***********************************************************************************/
/***********************************************************************************

   utsq_calc(): 

       -- calculates \tilde{u}^2;

***********************************************************************************/
static FTYPE utsq_calc(FTYPE W)
{
  FTYPE Wsq,W4,utsq ;
	
  Wsq = W*W ;
  W4 = Wsq*Wsq ;

  utsq = (Bsq*QdotBsq + W*(2.*QdotBsq + Qtsq*W))/
    (W4 + 2.*Bsq*Wsq*W + Wsq*(Bsq*Bsq - Qtsq) 
     - QdotBsq*(2.*W + Bsq)) ; 

  return(utsq) ;	
}

/***********************************************************************************/
/***********************************************************************************

   vsq_calc(): 

       -- calculates v^2, where v^i is the physical velocity w.r.t. normal observer;

***********************************************************************************/
static FTYPE vsq_calc(FTYPE W)
{
  FTYPE Wsq,Xsq;
	
  Wsq = W*W ;
  Xsq = (Bsq + W) * (Bsq + W);

  return(  ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq) );

}


/*********************************************************************/
/*********************************************************************

  find_root_1D(void) : 

       -- returns with the physical root;
       -- uses the Laguerre's method for finding all complex roots
           of a given polynomial;
*************************************************************************/
static FTYPE find_root_1D(FTYPE guess) 
{
	int i, j, k, n_real_roots , ret, n;
	FTYPE final_root, min_resid;
	FTYPE residtmp;
	static int find_all_roots( FTYPE complex roots_z[], 
				   FTYPE complex coefs_z[], int n );


	n = POLY_ORDER;

	/* Reset roots[]: */
	for( k = 0; k < n ; k++ ) {
	  roots_z[k] = guess;
	  roots[k]   = 0.;
	}

	/* Calculate coefficients of the polynomial : */
	calc_poly_coefs();

	for( i = 0; i <= n; i++ ) coefs_z[i] = coefs[i]  ; 

	/* Find all the roots : */
	if( (ret = find_all_roots( roots_z, coefs_z, n ) ) < 0 ) {
	  return( FAIL_VAL );
	}

	/* Find all the positive, real roots  : */
	n_real_roots = 0;
	for( i = 0; i < n; i++ ) {
	  if( (fabs(cimag(roots_z[i])) < fabs(IMAG_FLR_PHYS * creal(roots_z[i]) ) ) 
	      && ( creal(roots_z[i]) >= 0. ) )  {
	    roots[n_real_roots] = creal(roots_z[i]);
	    ++n_real_roots;
	  }
	}

	if( n_real_roots == 0 ) {
	  return( FAIL_VAL );
	}

	/***********************************************************************
	    Return the real root that leads to the smallest value of the original 
           residual function:
	************************************************************************/

	min_resid  =  1.e30;
	final_root = -1.e30;
	   
	for( i = 0 ; i < n_real_roots; i++ ) {

	  residtmp = fabs( residual0( roots[i] ) );

	  if( residtmp < min_resid ) {
	    min_resid  = residtmp;
	    final_root = roots[i] ;
	  }
	}

	return( final_root );
}


/*******************************************************

find_all_roots(): 
    
  -- inspired a great deal by the "Numerical Recipes in C"
     routine called "zroots()"

  -- uses find_root_laguerre() to find all the roots 
     (real and complex) of the polynomial of coefficients
     coefs[0..n];

  -- if failure, returns negative integer whose abs. value
     is the number of roots it successively found before
     the error occurred.

  -- roots[] must be set to intended guess prior

*******************************************************/

static int find_all_roots( FTYPE complex roots[], FTYPE complex coefs[], int n )
{

  int nroots,i, j,ret;
  FTYPE tol, tol_out;
  FTYPE complex x;
  FTYPE complex red_coefs[MAX_ORDER+1];  /* the coef's of the reduced poly. */
  int ntries = 1; /* Number of times to attempt Laguerre solves while 
		     increasing the tolerance by 10 each time */
  int itries;

  static void div_poly( FTYPE complex coefs[], FTYPE complex x0, int n ) ;
  static int find_root_laguerre( FTYPE complex *x, FTYPE complex coefs[], 
				 FTYPE tol, FTYPE *tol_out, int n );


  if( n > MAX_ORDER ) {
    fprintf(stderr, "find_all_roots(): n > MAX_ORDER,   n = %d ... \n", n);
    return(-1); 
  }


  /* start with the full poly. :  */
  for( i = 0; i <= n; i++ ) red_coefs[i] = coefs[i] ;
  
  
  /*******************************************************************
    Find roots successively, starting from the smallest, and reduce
     the polynomial as you go in order to avoid finding a root 
     multiple times:
  ******************************************************************/
  for( i = n; i >= 1; i-- ) {
    x = roots[i-1] ;

    tol = REDUCED_ROOT_TOL;
    itries = 0;
    /* Find the root */
    while( (itries < ntries)
	   && ((ret = find_root_laguerre( &x, red_coefs, tol, &tol_out, n )) < 0)){
      tol *= 10;
      ++itries;
    }

    if( (ret < 0) && ( tol_out > MIN_ROOT_TOL ) ) { 
      return( i - n - 1 );
    }

    /*Floor the imaginary part of root */
    if( fabs(cimag(x)) < IMAG_FLOOR*fabs(creal(x)) ) x = creal(x);   
    
    /* set root estimate */
    roots[i-1] = x;    
    
  /* Reduce the polynomial to remove root */
    div_poly( red_coefs, x, n ); 

  }

  
  /* Polish the roots by using the current estimates as close guesses for 
     the full, non-reduced polynomial : */
  for( i = 0; i < n; i++ ) {
    tol = ROOT_TOL;
    itries = 0;
    /* Find the root */
    while( (itries < ntries)
	   && ( (ret = find_root_laguerre( &roots[i], coefs, tol,&tol_out,n)) < 0)) {
      tol *= 10;
      ++itries;
    }
    if( (ret < 0) && ( tol_out > MIN_ROOT_TOL ) ) { 
      return( -i - 1 );
    }
  }


  /* Sort the roots in order of increasing real parts : */
  for( i = 1; i < n ; i++ ) {       /* Start from the beginning */
    x = roots[i];

   /* Make sure that all the interior points are ordered */
    for( j = (i-1); j >= 0; j-- ) { 
      if( creal(roots[j]) <= creal(x) ) break;
      roots[j+1] = roots[j];
    }
    roots[j+1] = x;
  }


  return(0);

}

/*******************************************************

find_root_laguerre(): 
    
  -- uses Laguerre's method for the root of a complex polynomial.
      nearest the initial guess "x". 

  -- coefs[0..n] = array of poly. coefficients;

  -- if failure, returns negative integer whose abs. value
     is the number of iterations it tried before 
     an error occurred.

  -- inspired a great deal by "Num. Recipes in C"

*******************************************************/

static int find_root_laguerre( FTYPE complex *x, FTYPE complex coefs[], 
			       FTYPE tol, FTYPE *tol_out, int n )
{

  int istep, i;
  FTYPE complex fpoly, fpoly_x, fpoly_xx;  /*Value of poly. and first 2 derivatives*/
  FTYPE complex x_init, x_abs, g, h, gsq, sqrt_term, den_p, den_m, x_new, dx;
  FTYPE resid, den_p_abs, den_m_abs ;
  int doing_extra, i_extra, keep_iterating;
  static FTYPE rand_vect[N_CYCLE_BREAKS] = {0.51, 0.24, 0.76, 0.23, 0.63, 
					     0.8, 0.1, 0.65, 0.34, 0.4,
					     0.04, 0.96, 0.83, 0.32, 0.78, 
					     0.19, 0.45, 0.61, 0.3, 0.2 };

  x_init = *x;
  i_extra = doing_extra = 0 ;

  istep = 0;
  keep_iterating = 1;
  while (keep_iterating) { 
    
    fpoly = coefs[n];
    fpoly_x = fpoly_xx = 0. ;

    x_abs = cabs(*x);
    resid = cabs(fpoly);

    for( i = (n-1); i >= 0; i-- ) {
      fpoly_xx = *x * fpoly_xx + fpoly_x ;   /* really f'' / 2 */
      fpoly_x  = *x * fpoly_x  + fpoly   ;
      fpoly    = *x * fpoly    + coefs[i];

      resid =  cabs(fpoly) + x_abs*resid; 
    }

    g   = fpoly_x / fpoly;
    gsq = g*g;
    h   = gsq  -  (2. * fpoly_xx) / fpoly;

    sqrt_term = csqrt( (n-1)*( n*h - gsq ) );
    den_p = g + sqrt_term;
    den_m = g - sqrt_term;

    den_p_abs = cabs(den_p);
    den_m_abs = cabs(den_m);
    
    /* Choose largest denominator :  */
    if( den_p_abs < den_m_abs ) {
      den_p = den_m;
      den_p_abs = den_m_abs;
    }
    
    dx = (   (den_p_abs > 0.0) 
	     ?  ( (1.*n) / den_p ) 
	     : ( exp(log(1+x_abs)) * (cos((FTYPE) istep) + I*sin((FTYPE) istep))));

    x_new = *x - dx;


    /* Check to see if dx is negligibly small and return if it is: */
    if( *x == x_new ) {
      return(istep) ; 
    }

    *x = x_new;

    *tol_out = cabs(dx)/cabs(*x);

    /* Check to see if root has been found */
    if( (*tol_out < tol) && (doing_extra == 0) && (N_EXTRA > 0)  ) {
      doing_extra = 1;
    }
    
    if( doing_extra == 1 ) i_extra++ ;
    

    if( ((*tol_out < tol) && (doing_extra == 0)) 
	|| (i_extra > N_EXTRA) || ( istep >= (ROOT_MAX_ITER+1) )  ) {
      return( istep );
    }


    istep++;

  }
  
  return(-istep);

}


/**********************************************************************************
  Divides the polynomial of coefficients coefs[9]  by (x-x0) where x0 is most 
  likely a root.
***********************************************************************************/

static void div_poly( FTYPE complex coefs[], FTYPE complex x0, int n ) 
{

  FTYPE complex coeftmp, remainder;
  int k;

  remainder = coefs[n];
  coefs[n] = 0.;

  for( k = (n-1); k >= 0; k-- ) {
    coeftmp = coefs[k];
    coefs[k] = remainder;
    remainder = remainder * x0 + coeftmp ;
  }

  return;
}
