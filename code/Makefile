
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

default: $(EXECUTABLE)

dirs:
	mkdir -p $(OBJDIR)/

clean:
	rm -rf $(OBJDIR) *~ $(EXECUTABLE)

$(EXECUTABLE): dirs $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(LDFRAMEWORKS)

$(OBJDIR)/%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

$(OBJDIR)/%.o: %.cu
	$(NVCC) $< $(NVCCFLAGS) -c -o $@
