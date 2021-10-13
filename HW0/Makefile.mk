## Simple makefile that you can use for all assignments for AMS213A (WINTER 2019)
## You only need to change the variables in the beginning (except for ALLOBJ, that's set automatically)

EXECBIN    = pi_approx                      							 # Put all executables here
FSOURCE    = pi_approx.f90     											 # Put all source files (.f90) here
ALLOBJ     = ${FSOURCE:.f90=.o}             							 # Put all object files here (automatically set with substitution)
#SHAREDOBJ  = Lib.o                         							 # Put all shared object files (fortran modules) here
GFORTRAN   =  gfortran -fdefault-real-8 -fdefault-double-8 -Wall -Wextra # Put compilation command here

# Builds all executables
all : ${EXECBIN}

# Links object files into executables
.SECONDEXPANSION:
${EXECBIN} : ${SHAREDOBJ} $$@.o
	${GFORTRAN} -o $@ ${SHAREDOBJ} $@.o

# Builds all source files into object files
%.o : %.f90
	${GFORTRAN} -c $<

# Removes all executables and object files
clean :
	rm -f ${ALLOBJ} ${EXECBIN} *.mod
