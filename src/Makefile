CC=g++
CCFLAGS=-Wall -O3 #-pedantic
LIBS=-lm -lumfpack

#INCLUDES = -I./headers
# INCLUDES=-I./packages/muparser/include
INCLUDES=

#SRC=gnuplot_i.c ini.c expression_parser.c input.c alloc.c math.c struct.c solvers.c fdm.c modules.c ccmSolve.c tests.c main.c

BASE=estimator.cpp field.cpp log.cpp model.cpp modules.cpp PDE.cpp solver.cpp sparse.cpp stencil.cpp

MUP_SRC=./packages/muparser/*.cpp
INIH_SRC=./packages/inih/ini.c ./packages/inih/INIReader.cpp

EXTRA_SRC= $(MUP_SRC) $(INIH_SRC)

SRC=$(BASE) $(TEST_BASE) $(EXTRA_SRC) main.cpp

OBJ=$(SRC:.cpp=.o)

EXE=run

all: $(EXE)

$(EXE): $(SRC)
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ $^ $(LIBS)

.PHONY: test
test:	$(BASE) test_main.c
	$(CC) $(CCFLAGS) $(INCLUDES)  -o test.x $^ $(LIBS)

.PHONY: rebuild
rebuild: all
	./$(EXE) solve inputs.ini

%.o: %.c
	$(CC) $(CCFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) run
	$(RM) outputs/* logs/*
