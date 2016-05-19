
/**********************************************************************

   main.c:
 
   Test routine for the primitive variable inversion methods in order
   to demonstrate how to use/call the inversion routines.  

   This test routine sets the metric, inverse metric, metric's 
   determinant, and conserved variables.  It takes a given set of 
   primitive variables P_i, calculates U(P_i), perturbs them with 
   random factor to get P_g (P_i*R), and then uses P_g as the 
   guess for the inversion routines.  The calculated primitive 
   variables (P_f) are then outputted and compared with the original set.

   It prompts the user to specify which inversion method to use at runtime.


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


*********************************************************************/

#include "u2p_util.h"


/* The inversion methods: */
int Utoprim_1d(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
               FTYPE gdet, FTYPE prim[NPR]);

int Utoprim_2d(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
               FTYPE gdet, FTYPE prim[NPR]);

int Utoprim_1dvsq1(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
               FTYPE gdet, FTYPE prim[NPR]);

int Utoprim_1dvsq2(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
               FTYPE gdet, FTYPE prim[NPR]);

int Utoprim_5d(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
               FTYPE gdet, FTYPE prim[NPR]);

int Utoprim_poly(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
               FTYPE gdet, FTYPE prim[NPR]);


/**********************************************************************
/**********************************************************************
/**********************************************************************
 The program:
*********************************************************************/

int main()
{

  FTYPE U[NPR], P_i[NPR], P_g[NPR], P_f[NPR];
  FTYPE gdet, gcon[NDIM][NDIM], gcov[NDIM][NDIM];
  int   i,j,k,answer, bad_answer, ret;



  /*********************************************************************
    First determine which inversion method to use, let the user specify:
  ***********************************************************************/
  
  bad_answer = 1;
  while( bad_answer ) { 
    bad_answer = 0;
    printf("Please specify which inversion routine to use: \n");
    printf("  1) 1D_W   \n  2) 2D  \n  3) 1D_{v^2}  \n  4) 1D^*_{v^2} \n");
    printf("  5) 5D   \n  6) Polynomial \n Answer: ");
    scanf("%d",&answer);

    if( (answer < 1) || (answer > 6) ) { 
      bad_answer = 1;
      printf("Invalid answer, please try again \n\n");
    }
  }

  /*********************************************************************
    Set the metric to Minkowski (aka flat spacetime, special relativistic): 
  ***********************************************************************/

  gdet = -1.;
  for( i = 0; i < NDIM; i++) { 
    for( j = 0; j < NDIM; j++) { 
      gcov[i][j] = gcon[i][j] = 0.;
    }
  }

  gcov[0][0] = gcon[0][0] = -1.;

  for( i = 1; i < NDIM; i++) { 
    gcov[i][i] = gcon[i][i] = 1.;
  }

  
  /*********************************************************************
    Set the primitive variables: 
  ***********************************************************************/
  P_i[RHO]    = 1.0e-1;
  P_i[UU]     = 1.0e-2;
  P_i[UTCON1] = 1.0e-1;
  P_i[UTCON2] = -2.0e-1;
  P_i[UTCON3] = 3.0e-1;
  P_i[BCON1]  = 3.0e-1;
  P_i[BCON2]  = 2.0e-1;
  P_i[BCON3]  = -1.0e-1;


  /*********************************************************************
    Calculate the conserved variables associated w/ these primitive var's:
  ***********************************************************************/

  primtoU_g( P_i, gcov, gcon, gdet, U );


  /*********************************************************************
    Perturb the primitive variables by random factors to get the "guess":
  ***********************************************************************/

  for( i = 0; i < NPR; i++ ) { 
    P_g[i] = P_i[i] * 2.*( (1.*rand())/(1.*RAND_MAX) )  ;
    P_f[i] = P_g[i];
  }

  /*********************************************************************
    Feed the guess into the desired inversion routine:
  ***********************************************************************/

  switch( answer ) { 
    
  case 1:  
    printf("Using  Utoprim_1d...\n");      
    ret = Utoprim_1d(    U, gcov, gcon, gdet, P_f); 
    break;

  case 2:  
    printf("Using  Utoprim_2d...\n");
    ret = Utoprim_2d(    U, gcov, gcon, gdet, P_f); 
    break;

  case 3:  
    printf("Using  Utoprim_1dvsq1...\n");  
    ret = Utoprim_1dvsq1(U, gcov, gcon, gdet, P_f); 
    break;

  case 4:  
    printf("Using  Utoprim_1dvsq2...\n");  
    ret = Utoprim_1dvsq2(U, gcov, gcon, gdet, P_f); 
    break;

  case 5:  
    printf("Using  Utoprim_5d...\n");      
    ret = Utoprim_5d(    U, gcov, gcon, gdet, P_f); 
    break;

  case 6:  
    printf("Using  Utoprim_poly...\n");    
    ret = Utoprim_poly(  U, gcov, gcon, gdet, P_f); 
    break;

  default: fprintf(stderr,"Should not be here...\n");exit(1);
  }


  /*********************************************************************
    Output the resultant primtiive variables from the specified method:
  ***********************************************************************/

  printf("\nHere are the results: \n\n");

  printf("The inversion routine returned with status = %d \n", ret );
  printf("U[0-7]   = ");
  for( i = 0 ; i < NPR; i++ ) { 
    printf(" %25.16e", U[i] );
  }
  printf("\n");

  printf("P_i[0-7] = ");
  for( i = 0 ; i < NPR; i++ ) { 
    printf(" %25.16e", P_i[i] );
  }
  printf("\n");

  printf("P_g[0-7] = ");
  for( i = 0 ; i < NPR; i++ ) { 
    printf(" %25.16e", P_g[i] );
  }
  printf("\n");

  printf("P_f[0-7] = ");
  for( i = 0 ; i < NPR; i++ ) { 
    printf(" %25.16e", P_f[i] );
  }
  printf("\n");


  return(0);
  
}

