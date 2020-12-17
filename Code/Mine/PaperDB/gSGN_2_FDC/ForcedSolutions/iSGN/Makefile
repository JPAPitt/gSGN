FF = f77
FFLAGS = -c -Wall -freal-4-real-8 -std=legacy

objects = MultipleExperiments.out

all : $(objects)
 
%.out : FD2C.o InitialConditionProblems.o MiscFuncs.o %.a
	$(FF) $^ -o $@

%.a : %.f FD2C.o InitialConditionProblems.o MiscFuncs.o
	$(FF) $(FFLAGS) $< -o $@
	
%.o : %.f
	$(FF) $(FFLAGS) $< -o $@

clean:
	rm -f ./*.o *~ ./*.out ./*.mod

fresh : clean all
