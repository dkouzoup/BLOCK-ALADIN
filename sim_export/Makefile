#
#    This file was auto-generated by ACADO Code Generation Tool.
#    
#    ACADO Code Generation tool is a sub-package of ACADO toolkit --
#    A Toolkit for Automatic Control and Dynamic Optimization.
#    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
#    Milan Vukov, Rien Quirynen, KU Leuven.
#    Developed within the Optimization in Engineering Center (OPTEC)
#    under supervision of Moritz Diehl. All rights reserved.
#    
#    ACADO Toolkit is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 3 of the License, or (at your option) any later version.
#    
#    ACADO Toolkit is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#    
#    You should have received a copy of the GNU Lesser General Public
#    License along with ACADO Toolkit; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#    

UNAME := $(shell uname)

LDLIBS = -lm -lstdc++
ifeq ($(UNAME), Linux)
	LDLIBS += -lrt
endif

CCACHE_APP := $(shell which ccache 2>/dev/null)

CFLAGS = -O3 -I.
CXXFLAGS = -O3 -finline-functions -I. 
CC     = $(CCACHE_APP) gcc
CXX    = $(CCACHE_APP) g++

OBJECTS = 	\
			acado_auxiliary_sim_functions.o \
			acado_integrator.o 

.PHONY: all
all: acado_test acado_compare 

acado_test: ${OBJECTS} acado_test.o

acado_compare: ${OBJECTS} acado_compare.o

acado_auxiliary_sim_functions.o          : 	acado_common.h \
											acado_auxiliary_sim_functions.h
acado_integrator.o          : 	acado_common.h \
								acado_auxiliary_sim_functions.h
acado_test.o                : 	acado_common.h \
								acado_auxiliary_sim_functions.h
acado_compare.o             : 	acado_common.h \
								acado_auxiliary_sim_functions.h

${OBJECTS} : acado_auxiliary_sim_functions.h

.PHONY : clean
clean :
	-rm -f *.o *.a acado_test acado_compare

