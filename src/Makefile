CC=g++
CCFLAGS=-Wall -O3 #-pedantic
LIBS=-lm -lumfpack -lnlopt
INCLUDES=

BASE=estimator.cpp field.cpp functions.cpp log.cpp model.cpp modules.cpp optimization.cpp PDE.cpp plot.cpp solver.cpp sparse.cpp stencil.cpp operators.cpp

# INIH_SRC=./packages/inih/ini.cpp ./packages/inih/INIReader.cpp
# EXTRA_SRC= $(INIH_SRC) 

PACKAGES   := packages/inih packages/lepton
PACKAGES_SRC := $(foreach mdir,$(PACKAGES),$(wildcard $(mdir)/*.cpp) )

SRC=$(BASE) $(TEST_BASE) $(PACKAGES_SRC) main.cpp
OBJ := $(addsuffix .o, $(basename $(SRC)))

TEST=$(BASE) $(PACKAGES_SRC) tests.cpp
OBJTEST := $(addsuffix .o, $(basename $(TEST)))

EXE=simcomet
TESTEXE = test.x

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ $^ $(LIBS)

.PHONY: run
run: all
	./$(EXE) 

.PHONY:test
test: $(OBJTEST)
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $(TESTEXE) $^ $(LIBS)

%.o: %.cpp
	$(CC) $(CCFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(OBJ)
	$(RM) $(EXE)
	$(RM) -r outputs/* logs/*

# .PHONY: test
# test:	$(BASE) test_main.c
# 	$(CC) $(CCFLAGS) $(INCLUDES)  -o test.x $^ $(LIBS)
