##
##	This file is part of qp42.
##
##	qp42 -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al. 
##	All rights reserved.
##
##	qp42 is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qp42 is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qp42; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

################################################################################
#
# Description:
#	qpDUNES configuration file
#
# Authors:
#	Milan Vukov, milan.vukov@esat.kuleuven.be
#
# Year:
#	2013.
#
# NOTE:
#	- Linux/Unix only.
#
# Usage:
#	- Linux - Ubuntu:
#		* Users are supposed to source this file into ~/.bashrc file.
#
################################################################################

################################################################################
#
# Definitions for both users and developers.
#
################################################################################

# 
# Tell the user project where to find our headers, libraries and external
# packages, etc.
#
export qpDUNES_ENV_INCLUDE_DIRS="/Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev;/Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/include;/Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/include/qp;/Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/interfaces;/Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/interfaces/mpc"
export qpDUNES_ENV_LIBRARY_DIRS="/Users/elusiv/Documents/Repositories/GIT/qpDUNES-dev/build/lib"

export qpDUNES_ENV_BUILD_DIR=""

#
# List of qpDUNES libraries
#
export qpDUNES_ENV_STATIC_LIBRARIES="qpdunes"

#
# Sources and headers
#
export qpDUNES_ENV_SOURCES=""
export qpDUNES_ENV_HEADERS=""
