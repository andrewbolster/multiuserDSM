

# Defaults
N_THREADS=4
ISB_THREADS=2
COMP=INTEL
ROOT=..
MATH_LIB=EIGEN
SCEN=ADSL_DOWNSTREAM

ifdef PROFILE
DEBUG = -pg
COMP = GNU
OPT_LEVEL= -O0
endif

ifdef XC
N_THREADS=16
MKL_INCLUDE=/contrib2/cmkl/10.0.2.018/include/
MKL_LIBDIR=/contrib2/cmkl/10.0.2.018/lib/64
EIGEN_INCLUDE=/genusers4/alastair/eigen-2.0.0
BOOST_LIBDIR=/genusers4/alastair/boost_1_39_0_install/lib/
BOOST_DIR=/genusers4/alastair/boost_1_39_0_install/include/boost-1_39/
DEFINES = -D$(MATH_LIB) -DXC -DN_THREADS=$(N_THREADS) -D$(SCEN) -DISB_THREADS=$(ISB_THREADS) -DEIGEN_NO_DEBUG
else
MKL_INCLUDE=/opt/intel/mkl/10.0.1.014/include
MKL_LIBDIR=/opt/intel/mkl/10.0.1.014/lib/32
BOOST_DIR=/home/alastair/boost_1_39_0/
BOOST_LIBDIR=/home/alastair/boost_1_39_0/stage/lib/
EIGEN_INCLUDE=/home/alastair/eigen-2.0.0
DEFINES = -D$(MATH_LIB) -DN_THREADS=$(N_THREADS) -DISB_THREADS=$(ISB_THREADS) -D$(SCEN) -DEIGEN_NO_DEBUG

endif

INCLUDE = -I$(ROOT)/include -I$(MKL_INCLUDE) -I$(EIGEN_INCLUDE) -I$(BOOST_DIR)
srcfiles = $(shell ls *.c)
objects = $(shell ls *.c | perl -pe 's/\.c/\.o/')
headers = $(shell ls $(ROOT)/include/*.h)
class_headers =  $(shell ls *.h)
COMPILE_DATE = $(shell date +%d%m_%H%M)

ifeq ($(COMP),GNU)

ifdef DBG
DEBUG = -g
OPT_LEVEL = -O0
else
DEBUG=
OPT_LEVEL = -O99
endif


outfile = multiuser-gcc
CXX = g++

LIBDIRS=-L$(MKL_LIBDIR)

CXXFLAGS = $(DEBUG) $(OPT_LEVEL) -Wfatal-errors -Wall $(INCLUDE) $(LIBDIRS) $(DEFINES) -funroll-all-loops

BOOST_THREADLIB=boost_thread-gcc42-mt

endif

ifeq ($(COMP),INTEL)

ifdef DBG
DEBUG = -g
OPT_LEVEL = -O0
else
DEBUG=
OPT_LEVEL = -O2
endif

outfile = multiuser-icc
ifdef XC
MTUNE=itanium2
INTEL_LIBDIR=/contrib2/intel/cc/10.1.011/lib/
CXX = /contrib2/intel/cc/10.1.011/bin/icpc
outfile = multiuser-icc-$(COMPILE_DATE)
else
MTUNE=core2
INTEL_LIBDIR=/opt/intel/Compiler/11.0/081/lib/
CXX=/opt/intel/Compiler/11.0/081/bin/ia32/icpc
endif

LIBDIRS=-L$(INTEL_LIBDIR) -L$(MKL_LIBDIR)

CXXFLAGS = -x c++ $(INCLUDE) $(LIBDIRS) $(DEFINES) $(DEBUG) $(OPT_LEVEL) -mtune=$(MTUNE) -unroll-agressive -no-prec-div -fp-model precise -Wno-deprecated -std=c99

BOOST_THREADLIB=boost_thread-il-mt

endif

ifdef PROFILE
DEBUG = -pg
COMP = GNU
OPT_LEVEL= -O0
endif

LIBDIRS+=-L$(BOOST_LIBDIR)

CXXLIBS = -lm -lgsl -lgslcblas -lpthread -l$(BOOST_THREADLIB)

multiuser: $(srcfiles) $(headers) $(class_headers)
	$(CXX) $(CXXFLAGS) -o $(outfile) $(srcfiles) $(CXXLIBS)
