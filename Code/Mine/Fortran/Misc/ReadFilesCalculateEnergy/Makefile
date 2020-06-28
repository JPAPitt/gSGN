FF = f77
FFLAGS = -c -Wall -freal-4-real-8 -std=legacy

objects = ReadAndCalc.out ReadAndCalcNew.out ReadAndCalc-DB.out ReadAndCalc-DBBeta.out

all : $(objects)
 
%.out :  %.a
	$(FF) $^ -o $@

%.a : %.f
	$(FF) $(FFLAGS) $< -o $@
	
%.o : %.f
	$(FF) $(FFLAGS) $< -o $@

clean:
	rm -f ./*.o *~ ./*.out ./*.mod

fresh : clean all
