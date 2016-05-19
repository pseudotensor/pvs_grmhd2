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


#include "u2p_defs.h"

extern void primtoU_g( FTYPE prim[], FTYPE gcov[][4], FTYPE gcon[][4], FTYPE gdet,  FTYPE U[] );
extern void ucon_calc_g(FTYPE prim[],FTYPE gcov[][4],FTYPE gcon[][4],FTYPE ucon[]);
extern void raise_g(FTYPE vcov[], FTYPE gcon[][4], FTYPE vcon[]);
extern void lower_g(FTYPE vcon[], FTYPE gcov[][4], FTYPE vcov[]);
extern void ncov_calc(FTYPE gcon[][4],FTYPE ncov[]) ;
extern void bcon_calc_g(FTYPE prim[],FTYPE ucon[],FTYPE ucov[],FTYPE ncov[],FTYPE bcon[]); 
extern FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
extern FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);
extern FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);
extern int gamma_calc_g(FTYPE *pr, FTYPE gcov[NDIM][NDIM], FTYPE *gamma);

