
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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#define FTYPE double     /* your choice of floating-point data type */
#define NPR 8
#define NDIM 4
#define GAMMA	(4./3.)  /* Adiabatic index used for the state equation */


#define MAX_NEWT_ITER 30     /* Max. # of Newton-Raphson iterations for find_root_2D(); */
#define NEWT_TOL   1.0e-10    /* Min. of tolerance allowed for Newton-Raphson iterations */
#define MIN_NEWT_TOL  1.0e-10    /* Max. of tolerance allowed for Newton-Raphson iterations */
#define EXTRA_NEWT_ITER 2

#define NEWT_TOL2     1.0e-15      /* TOL of new 1D^*_{v^2} gnr2 method */
#define MIN_NEWT_TOL2 1.0e-10  /* TOL of new 1D^*_{v^2} gnr2 method */

#define W_TOO_BIG	1.e9	/* \gamma^2 (\rho_0 + u + p) is assumed
                                  to always be smaller than this.  This
				  is used to detect solver failures */
#define UTSQ_TOO_BIG	1.e9    /* \tilde{u}^2 is assumed to be smaller
                                  than this.  Used to detect solver
				  failures */

#define FAIL_VAL  1.e30    /* Generic value to which we set variables when a problem arises */

#define NUMEPSILON (2.2204460492503131e-16)


/* some mnemonics */
/* for primitive variables */
#ifndef RHO
#define RHO 	0 
#endif

#ifndef UU
#define UU 	1 
#endif

#define UTCON1 	2
#define UTCON2 	3
#define UTCON3 	4
#define BCON1	5
#define BCON2	6
#define BCON3	7

/* for conserved variables */
#define QCOV0	1
#define QCOV1	2
#define QCOV2	3
#define QCOV3	4


#define MYMAX(a,b) ( ((a) > (b)) ? (a) : (b) )

#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

#define delta(i,j) (((i) == (j)) ? 1. : 0.)

