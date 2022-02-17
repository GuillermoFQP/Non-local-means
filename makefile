COMPILER    = ifort
PREFIX      = /usr/local/src/Healpix_3.80
SUFFIX      = ifort
EXEC_PREFIX = ${PREFIX}/bin${SUFFIX}
LIBDIR      = ${PREFIX}/lib${SUFFIX}
SHARPDIR    = ${PREFIX}/lib
INCLUDEDIR  = ${PREFIX}/include${SUFFIX}

# Fortran compiler flags
FCFLAGS1    = -O3 -r8 -parallel -heap-arrays 256 -mcmodel=small
FCFLAGS2    = -I$(INCLUDEDIR) -cm -w -sox -qopt-report=0 -qopenmp -fPIC
FCFLAGS3    = -L$(LIBDIR) -L/usr/lib -L$(SHARPDIR) -lhealpix -lhpxgif -lsharp -lcfitsio -Wl,-R/usr/lib -Wl,-R$(SHARPDIR) -Wl,-R$(LIBDIR) -lcurl

# Program name
PROGRAM     = nlmeans

# Source file
SRCS        = nlmeans.f90

$(PROGRAM): $(SRCS)
	$(COMPILER) $(FCFLAGS1) $^ $(FCFLAGS2) -o $@ $(FCFLAGS3)
