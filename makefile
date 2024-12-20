#############################
# Intel Compiler
#############################


F90 = mpifort

F90FLAGS = -cpp -w -O3 -ffixed-line-length-none


MPIDIR=/cm/shared/apps/openmpi/openmpi-4.0.1
MPILIB=$(MPIDIR)/lib
MPIINC=$(MPIDIR)/include

SRC = src
DST = build

LIBPATH = -L$(MPILIB) -lmpi
INCPATH = -I$(MPIINC)


OBJ= $(patsubst %,$(DST)/%, EvolvingYieldBoundary.o)

$(DST)/%.o: $(SRC)/%.f90 
		$(F90) $(F90FLAGS) $(INCPATH) -c $(filter-out $(SRC)/,$^) -o $(DST)/$*.o 

all: EvolvingYieldBoundary


EvolvingYieldBoundary: $(OBJ)
		$(F90) $(F90FLAGS) -o $@ $^ $(LIBPATH)

clean:
		rm -f ./EvolvingYieldBoundary  $(DST)/*.o $(DST)/*genmod.f90 $(DST)/*genmod.mod $(DST)/*.mod


