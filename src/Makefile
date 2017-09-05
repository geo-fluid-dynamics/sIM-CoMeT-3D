CC=g++
CCFLAGS=-Wall -O3 #-pedantic
LIBS=-lm -lumfpack
INCLUDES=

BASE=estimator.cpp field.cpp log.cpp model.cpp modules.cpp PDE.cpp plot.cpp solver.cpp sparse.cpp stencil.cpp operators.cpp

# INIH_SRC=./packages/inih/ini.cpp ./packages/inih/INIReader.cpp
# EXTRA_SRC= $(INIH_SRC) 

PACKAGES   := packages/inih packages/lepton
PACKAGES_SRC := $(foreach mdir,$(PACKAGES),$(wildcard $(mdir)/*.cpp) )

SRC=$(BASE) $(TEST_BASE) $(PACKAGES_SRC) main.cpp
OBJ := $(addsuffix .o, $(basename $(SRC)))

EXE=simcomet

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ $^ $(LIBS)

.PHONY: run
run: all
	./$(EXE) 

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
