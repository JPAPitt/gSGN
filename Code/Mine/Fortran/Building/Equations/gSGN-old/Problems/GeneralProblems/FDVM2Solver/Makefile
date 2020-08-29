FF = gfortran
FFLAGS = -c -Wall -fdefault-real-8 -fdefault-integer-8 -fdefault-double-8

objects = IndividualExperiments.out MultipleExperiments.out

all : $(objects)
 
%.out : FDVM2.o InitialConditionProblems.o %.intout 
	$(FF) $^ -o $@

%.intout : %.for FDVM2.o InitialConditionProblems.o
	$(FF) $(FFLAGS) $< -o $@
	
%.o : %.for
	$(FF) $(FFLAGS) $< -o $@

clean:
	rm -f ./*.o *~ ./*.out ./*.mod

fresh : clean all
