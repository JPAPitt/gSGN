FF = gfortran
FFLAGS = -Wall

SRCS := $(wildcard *.for)
BINS := $(SRCS:%.for=%.o)

all: ${BINS}

%.o : %.for
	$(FF) $(FFLAGS) $< -o $@

clean:
	rm -f ./*.o *~ 

