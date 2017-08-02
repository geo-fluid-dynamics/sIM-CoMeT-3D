CC=g++
CCFLAGS=-O3 #-pedantic
LIBS=-lm -lumfpack
INCLUDES=

BASE=estimator.cpp field.cpp log.cpp model.cpp modules.cpp PDE.cpp plot.cpp solver.cpp sparse.cpp stencil.cpp

MUP_SRC=./packages/muparser/*.cpp
INIH_SRC=./packages/inih/ini.c ./packages/inih/INIReader.cpp
GNUP_SRC=
EXTRA_SRC= $(MUP_SRC) $(INIH_SRC) $(GNUP_SRC)

SRC=$(BASE) $(TEST_BASE) $(EXTRA_SRC) main.cpp
OBJ=$(SRC:.cpp=.o)

EXE=simcomet

all: $(EXE)

$(EXE): $(SRC)
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ $^ $(LIBS)

.PHONY: run
run: all
	./$(EXE) 

%.o: %.c
	$(CC) $(CCFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(EXE)
	$(RM) outputs/* logs/*

# .PHONY: test
# test:	$(BASE) test_main.c
# 	$(CC) $(CCFLAGS) $(INCLUDES)  -o test.x $^ $(LIBS)
