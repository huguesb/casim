
EXECUTABLE := casim
CU_FILES   := casim.cu
CU_DEPS    :=
CC_FILES   := main.cpp pgm.cpp

###########################################################

ARCH=$(shell uname | sed -e 's/-.*//g')
OBJDIR=objs
CXX=g++ -m64
CXXFLAGS=-O3 -Wall -g
LDFLAGS=-L/opt/cuda-toolkit/lib64

LIBS       :=
FRAMEWORKS := 

LIBS += cudart

LDLIBS  := $(addprefix -l, $(LIBS))
LDFRAMEWORKS := $(addprefix -framework , $(FRAMEWORKS))

NVCC=nvcc
NVCCFLAGS=-O3 -m64 -arch compute_11

OBJS=$(OBJDIR)/main.o $(OBJDIR)/pgm.o $(OBJDIR)/casim.o

.PHONY: dirs clean

default: casim pgmtrim

dirs:
	mkdir -p $(OBJDIR)/

clean:
	rm -rf $(OBJDIR) *~ $(EXECUTABLE)

casim: dirs $(OBJDIR)/main.o $(OBJDIR)/pgm.o $(OBJDIR)/casim.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $(LDFRAMEWORKS) -o casim $(OBJDIR)/main.o $(OBJDIR)/pgm.o $(OBJDIR)/casim.o

pgmtrim: dirs $(OBJDIR)/pgmtrim.o $(OBJDIR)/pgm.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $(LDFRAMEWORKS) -o pgmtrim $(OBJDIR)/pgmtrim.o $(OBJDIR)/pgm.o

$(OBJDIR)/main.o: pgm.h casim.h main.cpp
	$(CXX) main.cpp $(CXXFLAGS) -c -o $(OBJDIR)/main.o

$(OBJDIR)/pgm.o: pgm.h pgm.cpp
	$(CXX) pgm.cpp $(CXXFLAGS) -c -o $(OBJDIR)/pgm.o

$(OBJDIR)/pgmtrim.o: pgm.h pgmtrim.cpp
	$(CXX) pgmtrim.cpp $(CXXFLAGS) -c -o $(OBJDIR)/pgmtrim.o

$(OBJDIR)/casim.o: casim.h casim.cu
	$(NVCC) casim.cu $(NVCCFLAGS) -c -o $(OBJDIR)/casim.o
