CC = g++
CFLAGS = -O3 -fopenmp -std=c++11
LDFLAGS = -lgsl -lgslcblas -lm

HOST    := $(shell hostname)
ifeq ("$(HOST)" , "zen2")  # zen2 x86 family Rome
        CFLAGS += -march=znver2 # For alumnos node
else
ifeq ("$(HOST)" ,"clemente") # clemente x86 family broadwell
        CFLAGS +=  -march=broadwell # For clemente
endif
endif

SOURCES := halo_energy.cpp recentering.cpp compute_profile.cpp get_halo_props_PIC_v2.cpp moment_of_inertia.cpp calculate_shapes.cpp project_particles.cpp transform_coordinates.cpp io.cpp
OBJECTS := $(SOURCES:.cpp=.o)
TARGET = get_halo_props_PIC_v2

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(LDFLAGS) $< -o $@

clean:
	$(RM) $(TARGET) $(OBJECTS)
