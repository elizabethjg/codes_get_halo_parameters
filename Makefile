CC = g++
CFLAGS += -fopenmp -std=c++20
LDFLAGS = -lgsl -lgslcblas -lm -lhdf5

HOST    := $(shell hostname)
ifeq ("$(HOST)" , "rome01")  # zen2 x86 family Rome
   CC = clang++
   CFLAGS += -mfma -mavx2 -funroll-loops
else
ifeq  ("$(HOST)" , "rome02")
   CC = clang++
   CFLAGS += -mfma -mavx2 -funroll-loops
else
ifeq ("$(HOST)" ,"clemente")  # clemente x86 family broadwell
    CFLAGS +=  -march=broadwell  # For clemente
endif
endif
endif

SOURCES := halo_energy.cpp recentering.cpp compute_profile.cpp \
moment_of_inertia.cpp calculate_shapes.cpp project_particles.cpp transform_coordinates.cpp \
io.cpp pos_to_z.cpp make_z_table.cpp read_sidm_simu.cpp
OBJECTS := $(SOURCES:.cpp=.o)
TARGET = get_halo_props_PIC_v3

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(TARGET).cpp $(CFLAGS) $^ $(LDFLAGS) -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	$(RM) $(TARGET) $(OBJECTS)
