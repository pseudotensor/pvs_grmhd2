.IGNORE:


#-------------------------------------------------------------------------------
#    Copyright 2005 Scott C. Noble, Charles F. Gammie, 
#                   Jonathan C. McKinney, and Luca Del Zanna
#
#
#    This file is part of PVS-GRMHD.
#
#    PVS-GRMHD is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    PVS-GRMHD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PVS-GRMHD; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

#Intel C++ Compiler v.8.0 (these parameters were used in the parameter space survey):
#CC = icc
#CCFLAGS = -O3 -axW -tpp7
#CCCFLAGS = -c
#CCLFLAGS = -lm

#GNU/Default compiler options:
CC = cc
CCFLAGS = -g
CCCFLAGS = -c
CCLFLAGS = -lm

CC_COMPILE  =  $(CC) $(CCFLAGS) $(CCCFLAGS) 
CC_LOAD     = $(CC) $(CCFLAGS) $(CCLFLAGS)

.c.o:
	$(CC_COMPILE) $*.c

INC = u2p_defs.h u2p_util.h

SRC = main.c u2p_util.c utoprim_1d.c utoprim_2d.c utoprim_poly.c \
utoprim_5d.c utoprim_1dvsq1.c utoprim_1dvsq2.c 

OBJS = main.o u2p_util.o utoprim_1d.o utoprim_2d.o utoprim_poly.o \
utoprim_5d.o utoprim_1dvsq1.o utoprim_1dvsq2.o 

EXECUTABLES = pvstest

all: $(EXECUTABLES)

$(OBJS) : $(INC)

pvstest: $(OBJS) 
	$(CC_LOAD) $(OBJS) -o pvstest

clean:
	/bin/rm *.o
	/bin/rm pvstest
