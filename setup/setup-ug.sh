#!/bin/sh

if test "$1" = "help"
then
	echo "usage: setup.sh <block size> <max blocks> <solver> <objdir> <site dir>"
	exit 1
fi

# Excerpt from FLASH 4.3 User Manual
#   > The Uniform Grid has the same resolution in all the blocks throughout the
#   > domain, and each processor has *exactly one block*. The uniform grid can operate
#   > in either of two modes: fixed block size (FIXEDBLOCKSIZE) mode, and non-fixed
#   > block size (NONFIXEDBLOCKSIZE) mode. The default fixed block size grid is
#   > statically defined at compile time and can therefore take advantage of
#   > compile-time optimizations. The non-fixed block size version uses dynamic
#   > memory allocation of grid variables.
#
# We will only use 'fixed block size' mode.
#
#   > In this mode, the block size is specified at compile time as NXB x NYB x NZB. 
# 
# Distribution of blocks on processors is determined via {i,j,k}procs:
# 
#    iprocs * jprocs * kprocs == nprocs
#
# Global domain size:
#
#    nxb*iprocs * nyb*jprocs * nzb*kprocs == gridsize
#
# Eg., for gridsize = 32^3:
#
#      (16 * 2)^3  =  32^3  --> 2^3 =   8 ... processes
#      ( 8 * 4)^3  =  32^3  --> 4^3 =  64 ... processes
#
# Eg., for gridsize = 64^3:
#
#      (32 * 2)^3  =  64^3  --> 2^3 =   8 ... processes
#      (16 * 4)^3  =  64^3  --> 4^3 =  64 ... processes
#
# Eg., for gridsize = 128^3:
#
#      (32 * 4)^3  = 128^3  --> 4^3 =  64 ... processes
#
# Eg., for gridsize = 256^3:
#
#      (64 * 4)^3  = 256^3  --> 4^3 =  64 ... processes
#      (32 * 8)^3  = 256^3  --> 8^3 = 512 ... processes


BLOCKSIZE="${1:?No block size given!}"
MAXBLOCKS="${2:?No maxblocks given}"
   SOLVER="${3:?No solver given!}"
   OBJDIR="${4:?No source directory given!}"
  SITEDIR="${5:?No site directory given!}"

declare -A SOLVERUNIT
SOLVERUNIT['pm']='--with-unit=physics/Hydro/HydroMain/split/PPM'
SOLVERUNIT['8w']='--with-unit=physics/Hydro/HydroMain/split/MHD_8Wave'
SOLVERUNIT['b3']='--with-unit=physics/Hydro/HydroMain/split/Bouchut3'
SOLVERUNIT['b5']='--with-unit=physics/Hydro/HydroMain/split/Bouchut5'
SOLVERUNIT['es']='--with-unit=physics/Hydro/HydroMain/split/ES'

./setup Girichidis-StirTurb -3d -auto -portable -opt +ug \
	-nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=$MAXBLOCKS \
	-site="$SITEDIR" -objdir="$OBJDIR" ${SOLVERUNIT[$SOLVER]}
