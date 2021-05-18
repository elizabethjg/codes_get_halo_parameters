CC = g++
CFLAGS = -O3 -fopenmp -ffast-math
LDFLAGS = -lm -lgsl -lgslcblas
SOURCES := halo_energy.cpp recentering.cpp get_halo_props_PIC_v2.cpp
OBJECTS := $(SOURCES:.cpp=.o)
TARGET = get_halo_props_PIC_v2

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS)-o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) $(TARGET) $(OBJECTS)
