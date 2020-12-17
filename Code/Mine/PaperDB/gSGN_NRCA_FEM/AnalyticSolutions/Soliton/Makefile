FF = f77
FFLAGS = -c -Wall -freal-4-real-8 -std=legacy 
objects =   MultipleExperiments.out TestNewMethodsRecons.out
LINKFLAGS = -L/home/jp/libf77 -llapack -lblas

all : $(objects)
 
%.out : FDVM2.o InitialConditionProblems.o MiscFuncs.o %.a
	$(FF) $^ $(LINKFLAGS) -o $@

%.a : %.f FDVM2.o InitialConditionProblems.o MiscFuncs.o
	$(FF) $(FFLAGS) $< -o $@
	
%.o : %.f
	$(FF) $(FFLAGS) $< -o $@

clean:
	rm -f ./*.o *~ ./*.out ./*.mod

fresh : clean all
