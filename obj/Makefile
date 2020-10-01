# Remember module load netcdf and module load dfftpack first.
#
# Makefile
#
# The compiler options
FC = gfortran
CC = 
#
FFLAGS = -g -fbounds-check
#
FCFLAGS += $(NETCDF_FFLAGS)
#
LDFLAGS = $(NETCDF_LDFLAGS)
#
LIBPATH = -L$$DFFTPACK_HOME/lib
#
INCLUDEPATH =
#
LIB =-ldfftpack
#
PROG = run
# The object files
OBJECTS = structure.o main.o
#
$(PROG) : $(OBJECTS)
	$(FC) $(FFLAGS) $(FCFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBPATH) $(LIB)
	./run 

structure.o : ../src/structure.f90
	$(FC) $(FFLAGS) $(FCFLAGS)  -c ../src/structure.f90

main.o : ../src/main.f90
	$(FC) $(FFLAGS) $(FCFLAGS)  -c ../src/main.f90

clean :
	rm -f *.o
	rm -r *.mod
	rm -f $(PROG)
