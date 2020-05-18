FF = gfortran
FFLAGS = -c -Wall -fdefault-real-8 -fdefault-integer-8 -fdefault-double-8
LINKFLAGS = -L/home/jp/libf77 -llapack -lblas
SRCS := $(wildcard *.for)
BINS := $(SRCS:%.for=%.out)

all: ${BINS}

%.out : %.o
	$(FF) $<  $(LINKFLAGS) -o $@

%.o : %.for
	$(FF) $(FFLAGS) $< -o $@

fresh : clean all
clean:
	rm -f ./*.o *~ ./*.out

