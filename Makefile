CC = g++
CFLAGS = -O3 -march=znver2 -funroll-loops -ffast-math -faggressive-loop-optimizations -ftree-vectorize
LDFLAGS = -lgsl -lgslcblas -lm

export LD_LIBRARY_PATH = /opt/spack/0.16.1/opt/spack/linux-centos7-x86_64/gcc-10.2.0/gsl-2.6-iw2ra42l6umymgghwu74o5di5t42gx6c/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH = /opt/spack/0.16.1/opt/spack/linux-centos7-x86_64/gcc-10.2.0/gsl-2.6-iw2ra42l6umymgghwu74o5di5t42gx6c/lib:$LIBRARY_PATH

SOURCES := halo_energy.cpp recentering.cpp get_halo_props_PIC_v2.cpp
OBJECTS := $(SOURCES:.cpp=.o)
TARGET = get_halo_props_PIC_v2

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	$(RM) $(TARGET) $(OBJECTS)
